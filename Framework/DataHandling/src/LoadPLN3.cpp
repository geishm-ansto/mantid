// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "MantidDataHandling/LoadPLN3.h"
#include "MantidAPI/Axis.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/LogManager.h"
#include "MantidAPI/RegisterFileLoader.h"
#include "MantidAPI/Run.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"
#include "MantidKernel/MandatoryValidator.h"
#include "MantidKernel/OptionalBool.h"
#include "MantidKernel/PhysicalConstants.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidNexus/NexusClasses.h"

#include <boost/filesystem.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <regex>
#include <tuple>

namespace Mantid {
namespace DataHandling {

using namespace Kernel;

namespace {

// number of physical detectors
constexpr size_t MONITORS = 8;
constexpr size_t DETECTOR_TUBES = 200;
constexpr size_t HISTO_BINS_X = DETECTOR_TUBES + MONITORS;
constexpr size_t HISTO_BINS_Y_DENUMERATOR = 16;
constexpr size_t TUBE_RESOLUTION = 1024;
constexpr size_t PIXELS_PER_TUBE = TUBE_RESOLUTION / HISTO_BINS_Y_DENUMERATOR;
constexpr size_t DETECTOR_SPECTRA = DETECTOR_TUBES * PIXELS_PER_TUBE;
constexpr size_t HISTOGRAMS = DETECTOR_SPECTRA + MONITORS;

// File loading progress boundaries
constexpr size_t Progress_LoadBinFile = 48;
constexpr size_t Progress_ReserveMemory = 4;
constexpr size_t Progress_Total = 2 * Progress_LoadBinFile + Progress_ReserveMemory;

// Algorithm parameter names
constexpr char FilenameStr[] = "Filename";
constexpr char MaskStr[] = "Mask";
constexpr char SelectDetectorTubesStr[] = "SelectDetectorTubes";
constexpr char SelectDatasetStr[] = "SelectDataset";
constexpr char FilterByTimeStartStr[] = "FilterByTimeStart";
constexpr char FilterByTimeStopStr[] = "FilterByTimeStop";
constexpr char TOFBiasStr[] = "TimeOfFlightBias";
constexpr char CalibrateTOFStr[] = "CalibrateTOFBias";
constexpr char LambdaOnTwoStr[] = "LambdaOnTwoMode";

// Common pairing of limits
using TimeLimits = std::pair<double, double>;

template <typename Type>
void AddSinglePointTimeSeriesProperty(API::LogManager &logManager, const std::string &time, const std::string &name,
                                      const Type value) {
  // create time series property and add single value
  auto p = new Kernel::TimeSeriesProperty<Type>(name);
  p->addValue(time, value);

  // add to log manager
  logManager.addProperty(p);
}

// Utility functions for loading values with defaults
// Single value properties only support int, double, string and bool
template <typename Type>
Type GetNeXusValue(NeXus::NXEntry &entry, const std::string &path, const Type &defval, int32_t index) {
  try {
    NeXus::NXDataSetTyped<Type> dataSet = entry.openNXDataSet<Type>(path);
    dataSet.load();

    return dataSet()[index];
  } catch (std::runtime_error &) {
    return defval;
  }
}

// string is a special cases
template <>
std::string GetNeXusValue<std::string>(NeXus::NXEntry &entry, const std::string &path, const std::string &defval,
                                       int32_t /*unused*/) {

  try {
    NeXus::NXChar dataSet = entry.openNXChar(path);
    dataSet.load();

    return std::string(dataSet(), dataSet.dim0());
  } catch (std::runtime_error &) {
    return defval;
  }
}

template <typename T>
void MapNeXusToProperty(NeXus::NXEntry &entry, const std::string &path, const T &defval, API::LogManager &logManager,
                        const std::string &name, const T &factor, int32_t index) {

  T value = GetNeXusValue<T>(entry, path, defval, index);
  logManager.addProperty<T>(name, value * factor);
}

// sting and float are special cases
template <>
void MapNeXusToProperty<std::string>(NeXus::NXEntry &entry, const std::string &path, const std::string &defval,
                                     API::LogManager &logManager, const std::string &name,
                                     const std::string & /*unused*/, int32_t index) {

  std::string value = GetNeXusValue<std::string>(entry, path, defval, index);
  logManager.addProperty<std::string>(name, value);
}

template <>
void MapNeXusToProperty<float>(NeXus::NXEntry &entry, const std::string &path, const float &defval,
                               API::LogManager &logManager, const std::string &name, const float &factor,
                               int32_t index) {

  float value = GetNeXusValue<float>(entry, path, defval, index);
  logManager.addProperty<double>(name, value);
}

template <typename T>
void MapNeXusToSeries(NeXus::NXEntry &entry, const std::string &path, const T &defval, API::LogManager &logManager,
                      const std::string &time, const std::string &name, const T &factor, int32_t index) {

  auto value = GetNeXusValue<T>(entry, path, defval, index);
  AddSinglePointTimeSeriesProperty<T>(logManager, time, name, value * factor);
}

template <>
void MapNeXusToSeries(NeXus::NXEntry &entry, const std::string &path, const float &defval, API::LogManager &logManager,
                      const std::string &time, const std::string &name, const float &factor, int32_t index) {

  auto value = GetNeXusValue<float>(entry, path, defval, index);
  AddSinglePointTimeSeriesProperty<double>(logManager, time, name, value * factor);
}

// Extract the relevant timestamp data. Get the timestamp first to
// determine the index limits (it assumed the timestamp is ordered).
// Then extract the value from that range. If extn is not zero it will
// n elements before and after the range.
template <typename T>
int extractEvents(NeXus::NXEntry &entry, const std::string &path, uint64_t startTime, uint64_t endTime,
                  std::vector<uint64_t> &times, std::vector<T> &events, int extn = 0) {

  auto timeStamp = entry.openNXDataSet<uint64_t>(path + "/time");
  timeStamp.load();
  int maxn = timeStamp.size();
  int startIx{0}, endIx{0};
  auto itt = timeStamp();
  for (int i = 0; i < maxn; i++) {
    auto v = itt[i];
    if (v <= startTime)
      startIx = i;
    if (v < endTime)
      endIx = i + 2;
  }
  startIx = std::max(startIx - extn, 0);
  endIx = std::min(endIx + extn, maxn);
  times.assign(itt + startIx, itt + endIx);

  auto values = entry.openNXDataSet<T>(path + "/value");
  values.load();
  auto itv = values();
  events.assign(itv + startIx, itv + endIx);

  return endIx - startIx;
}

int64_t beamMonitorCounts(NeXus::NXEntry &entry, const std::string &path, uint64_t startTime, uint64_t endTime,
                          bool sumValues) {
  std::vector<uint64_t> times;
  std::vector<int64_t> values;
  auto n = extractEvents<int64_t>(entry, path, startTime, endTime, times, values, 0);

  if (sumValues) {
    // take a fraction of the first and last based on the timestamp
    // note that values[k] is events between times[k-1] and times[k]
    // so the lead in value needs to be subtracted and the tail added
    int64_t count{0};
    for (int i = 1; i < n - 1; i++)
      count += values[i];
    if (n > 1) {
      if (startTime > times[0]) {
        int64_t delta = startTime - times[0];
        int64_t frac = values[1] * delta;
        frac /= (times[1] - times[0]);
        count -= frac;
      }
      if (endTime > times[n - 2]) {
        int64_t delta = endTime - times[n - 2];
        int64_t frac = values[n - 1] * delta;
        frac /= (times[n - 1] - times[n - 2]);
        count += frac;
      }
    }

    return count;
  } else {
    return values[n - 1];
  }
}

template <class IEventHandler>
void ReadEventData(NeXus::NXEntry &entry, IEventHandler &handler, uint64_t start_nsec, uint64_t end_nsec,
                   const std::string &neutron_path) {

  // the detector event time zero is the actual chopper time and all the events are
  // relative to this base

  // Get the event index, base values and values and
  auto eventID = entry.openNXDataSet<uint32_t>(neutron_path + "/event_id");
  eventID.load();
  auto eventIndex = entry.openNXDataSet<uint32_t>(neutron_path + "/event_index");
  eventIndex.load();
  auto zeroOffset = entry.openNXDataSet<uint64_t>(neutron_path + "/event_time_zero");
  zeroOffset.load();
  auto offsetValues = entry.openNXDataSet<uint32_t>(neutron_path + "/event_time_offset");
  offsetValues.load();
  uint32_t numPulses = eventIndex.size();
  uint32_t totalEvents = offsetValues.size();

  // The chopper times are monotonically increasing but there may be duplicate
  // pulse times when a 'efu' buffer is full mid pulse. In this case the buffer
  // is sent but the next buffer uses the same pulse time. Only send the frame
  // event when it changes.

  uint64_t lastFrameTS{0};

  // Run through the data and forward the event data if the pulse time occurs between
  // start and end time. To be clear the start and end time relates to the pulse times.
  for (uint32_t ix = 0; ix < numPulses; ix++) {
    auto pulseTime = zeroOffset[ix];
    if (start_nsec <= pulseTime && pulseTime < end_nsec) {
      if (pulseTime > lastFrameTS) {
        handler.newFrame();
        lastFrameTS = pulseTime;
      }
      auto baseIndex = eventIndex[ix];
      auto lastIndex = (ix + 1 < numPulses ? eventIndex[ix + 1] : totalEvents);
      for (uint32_t j = baseIndex; j < lastIndex; j++) {
        // convert tof to microseconds as double and pixel as (x,y)
        // and send to handler
        double tof = offsetValues[j] * 1.0e-3;
        auto pixel = eventID[j];
        size_t y = pixel % TUBE_RESOLUTION;
        size_t x = (pixel - y) / TUBE_RESOLUTION;
        handler.addEvent(x, y, tof, 0.0);
      }
    }
  }
}

bool addTimeZone(const std::string &ins, std::string &outs) {
  const std::regex fmt("^(\\d{4}-\\d{1,2}-\\d{1,2}) (\\d{2}:\\d{2}:\\d{2})");
  std::smatch match;
  bool ok = false;
  if (std::regex_match(ins, match, fmt)) {
    auto date = match[1].str();
    auto time = match[2].str();

    // add the ansto time +10:00
    outs = date + "T" + time + "+10:00";
    ok = true;
  }
  return ok;
}

// Extract the start and end time in nsecs from the nexus file
// based on entry/scan_dataset/[time, value]
//
// The time is the start time value is duration in nsec for each dataset
//
std::pair<uint64_t, uint64_t> getTimeScanLimits(NeXus::NXEntry &entry, int datasetIx) {

  auto timestamp = entry.openNXDataSet<uint64_t>("scan_dataset/time");
  timestamp.load();
  auto offset = entry.openNXDataSet<int64_t>("scan_dataset/value");
  offset.load();
  try {
    auto start = timestamp[datasetIx];
    auto end = start + offset[datasetIx];
    return {start, end};
  } catch (std::runtime_error &) {
    return {0, 0};
  }
}

// map the comma separated range of indexes to the vector via a lambda function
// throws an exception if it is outside the vector range
//
template <typename T, typename F> void mapRangeToIndex(const std::string &line, std::vector<T> &result, const F &fn) {

  std::stringstream ss(line);
  std::string item;
  size_t index = 0;
  while (std::getline(ss, item, ',')) {
    auto const k = item.find('-');

    size_t p0, p1;
    if (k != std::string::npos) {
      p0 = boost::lexical_cast<size_t>(item.substr(0, k));
      p1 = boost::lexical_cast<size_t>(item.substr(k + 1, item.size() - k - 1));
    } else {
      p0 = boost::lexical_cast<size_t>(item);
      p1 = p0;
    }

    if (p1 < result.size() && p0 <= p1) {
      while (p0 <= p1) {
        result[p0++] = fn(index);
        index++;
      }
    } else if (p0 < result.size() && p1 < p0) {
      do {
        result[p0] = fn(index);
        index++;
      } while (p1 < p0--);
    } else
      throw std::invalid_argument("invalid range specification");
  }
}

// Simple reader that is compatible with the ASNTO event file loader
class FileLoader {
  std::ifstream _ifs;
  size_t _size;

public:
  explicit FileLoader(const char *filename) : _ifs(filename, std::ios::binary | std::ios::in) {
    if (!_ifs.is_open() || _ifs.fail())
      throw std::runtime_error("unable to open file");

    _ifs.seekg(0, _ifs.end);
    _size = _ifs.tellg();
    _ifs.seekg(0, _ifs.beg);
  }

  bool read(char *s, std::streamsize n) { return static_cast<bool>(_ifs.read(s, n)); }

  size_t size() { return _size; }

  size_t position() { return _ifs.tellg(); }

  size_t selected_position() { return _ifs.tellg(); }
};

} // anonymous namespace

namespace PLN3 {

//
// In the future the ANSTO helper and event file loader will be generalized to
// handle the instruments consistently.

// Simple 1D histogram class
class SimpleHist {
  std::vector<size_t> m_hist;
  double m_M;
  double m_B;
  size_t m_peak;
  size_t m_count;

public:
  SimpleHist(size_t N, double minVal, double maxVal) : m_hist(N, 0) {
    m_M = (static_cast<double>(N) / (maxVal - minVal));
    m_B = -m_M * minVal;
    m_peak = 0;
    m_count = 0;
  }

  inline double ival(double val) const { return m_M * val + m_B; }

  inline double xval(double ix) const { return (ix - m_B) / m_M; }

  inline void add(double val) {
    auto ix = static_cast<size_t>(std::floor(ival(val)));
    if (ix < m_hist.size()) {
      m_hist[ix]++;
      m_count++;
      if (m_hist[ix] > m_peak)
        m_peak = m_hist[ix];
    }
  }

  const std::vector<size_t> &histogram() const { return m_hist; }

  inline size_t peak() const { return m_peak; }
  inline size_t count() const { return m_count; }
};

class EventProcessor {
protected:
  // fields
  const std::vector<bool> &m_roi;
  const std::vector<size_t> &m_mapIndex;
  const double m_framePeriod;
  const double m_gatePeriod;

  // number of frames
  size_t m_frames;
  size_t m_framesValid;

  // time boundaries
  const TimeLimits m_timeBoundary; // seconds

  virtual void addEventImpl(size_t id, size_t x, size_t y, double tof) = 0;

public:
  EventProcessor(const std::vector<bool> &roi, const std::vector<size_t> &mapIndex, const double framePeriod,
                 const double gatePeriod, const TimeLimits &timeBoundary)
      : m_roi(roi), m_mapIndex(mapIndex), m_framePeriod(framePeriod), m_gatePeriod(gatePeriod), m_frames(0),
        m_framesValid(0), m_timeBoundary(timeBoundary) {}

  void newFrame() {
    m_frames++;
    if (validFrame())
      m_framesValid++;
  }

  inline bool validFrame() const {
    double frameTime = m_framePeriod * static_cast<double>(m_frames) * 1.0e-6;
    return (frameTime >= m_timeBoundary.first && frameTime <= m_timeBoundary.second);
  }

  double duration() const {
    // length test in seconds
    return m_framePeriod * static_cast<double>(m_frames) * 1.0e-6;
  }

  inline int64_t frameStart() const {
    // returns time in nanoseconds from start of test
    auto start = m_framePeriod * static_cast<double>(m_frames);
    return static_cast<int64_t>(start * 1.0e3);
  }

  void addEvent(size_t x, size_t p, double tof, double /*taux*/) {

    // check if in time boundaries
    if (!validFrame())
      return;

    // group pixels
    auto y = static_cast<size_t>(p / HISTO_BINS_Y_DENUMERATOR);

    // determine detector id and check limits
    if (x >= HISTO_BINS_X || y >= PIXELS_PER_TUBE)
      return;

    // map the raw detector index to the physical model
    size_t xid = m_mapIndex[x];

    // take the modules of the tof time to account for the
    // longer background chopper rate
    double mtof = tof < 0.0 ? fmod(tof + m_gatePeriod, m_gatePeriod) : fmod(tof, m_gatePeriod);

    size_t id = xid < DETECTOR_TUBES ? PIXELS_PER_TUBE * xid + y : DETECTOR_SPECTRA + xid;
    if (id >= m_roi.size())
      return;

    // check if neutron is in region of interest
    if (!m_roi[id])
      return;

    // finally pass to specific handler
    addEventImpl(id, xid, y, mtof);
  }
};

// The class determines the number of counts linked to the detectors and the
// tof correction.
class EventCounter : public EventProcessor {
protected:
  // fields
  std::vector<size_t> &m_eventCounts;
  double m_L1;
  double m_V0;
  const std::vector<double> &m_L2;
  SimpleHist m_histogram;

  void addEventImpl(size_t id, size_t /*x*/, size_t /*y*/, double tobs) override {
    m_eventCounts[id]++;
    // the maximum occurs at the elastic peak
    double deltaT = 1.0e6 * (m_L1 + m_L2[id]) / m_V0 - tobs;
    m_histogram.add(deltaT);
  }

public:
  // construction
  EventCounter(const std::vector<bool> &roi, const std::vector<size_t> &mapIndex, const double framePeriod,
               const double gatePeriod, const TimeLimits &timeBoundary, std::vector<size_t> &eventCounts,
               const double L1, const double V0, const std::vector<double> &vecL2)
      : EventProcessor(roi, mapIndex, framePeriod, gatePeriod, timeBoundary), m_eventCounts(eventCounts), m_L1(L1),
        m_V0(V0), m_L2(vecL2), m_histogram(5000, -2500.0, 2500.0) {}

  size_t numFrames() const { return m_framesValid; }

  size_t numEvents() const {
    size_t sum(0);
    return std::accumulate(m_eventCounts.begin(), m_eventCounts.end(), sum);
  }

  // clips the histogram above 25% and takes the mean of the values
  double tofCorrection() {

    // determine the points above the 25% threshold
    auto minLevel = static_cast<size_t>(m_histogram.peak() / 4);
    auto hvec = m_histogram.histogram();
    double sum = 0.0;
    size_t count = 0;
    for (size_t i = 0; i < hvec.size(); i++) {
      if (hvec[i] >= minLevel) {
        auto ix = static_cast<double>(i);
        sum += static_cast<double>(hvec[i]) * m_histogram.xval(ix + 0.5);
        count += hvec[i];
      }
    }

    return (count > 0 ? sum / static_cast<double>(count) : 0.0);
  }
};

class EventAssigner : public EventProcessor {
protected:
  // fields
  std::vector<EventVector_pt> &m_eventVectors;
  double m_tofMin;
  double m_tofMax;
  int64_t m_startTime;
  double m_tofCorrection;
  double m_sampleTime;

  void addEventImpl(size_t id, size_t /*x*/, size_t /*y*/, double tobs) override {

    // get the absolute time for the start of the frame
    auto const offset = m_startTime + frameStart();

    // adjust the the tof to account for the correction and allocate events
    // that occur before the sample time as slow events from the previous pulse
    double tof = tobs + m_tofCorrection - m_sampleTime;
    if (tof < 0.0)
      tof = fmod(tof + m_gatePeriod, m_gatePeriod);
    tof += m_sampleTime;
    if (m_tofMin > tof)
      m_tofMin = tof;
    if (m_tofMax < tof)
      m_tofMax = tof;

    auto ev = Types::Event::TofEvent(tof, Types::Core::DateAndTime(offset));
    m_eventVectors[id]->emplace_back(ev);
  }

public:
  EventAssigner(const std::vector<bool> &roi, const std::vector<size_t> &mapIndex, const double framePeriod,
                const double gatePeriod, const TimeLimits &timeBoundary, std::vector<EventVector_pt> &eventVectors,
                int64_t startTime, double tofCorrection, double sampleTime)
      : EventProcessor(roi, mapIndex, framePeriod, gatePeriod, timeBoundary), m_eventVectors(eventVectors),
        m_tofMin(std::numeric_limits<double>::max()), m_tofMax(std::numeric_limits<double>::min()),
        m_startTime(startTime), m_tofCorrection(tofCorrection), m_sampleTime(sampleTime) {}

  double tofMin() const { return m_tofMin <= m_tofMax ? m_tofMin : 0.0; }
  double tofMax() const { return m_tofMin <= m_tofMax ? m_tofMax : 0.0; }
};

template <typename EP>
void loadEvents(API::Progress &prog, const char *progMsg, EP &eventProcessor, NeXus::NXEntry &entry,
                uint64_t start_nsec, uint64_t end_nsec) {

  using namespace ANSTO;

  prog.doReport(progMsg);

  // FileLoader loader(eventFile.c_str());

  // for progress notifications
  // ANSTO::ProgressTracker progTracker(prog, progMsg, 1000U,
  //                                   Progress_LoadBinFile);

  // ReadEventFile(loader, eventProcessor, progTracker, 100, false);
  const std::string neutronPath{"instrument/detector_events"};
  ReadEventData(entry, eventProcessor, start_nsec, end_nsec, neutronPath);
}
} // namespace PLN3

/// Initialise the algorithm and declare the properties for the
/// nexus descriptor.
void LoadPLN3::init() {

  // Specify file extensions which can be associated with a specific file.
  std::vector<std::string> exts;

  // Declare the Filename algorithm property. Mandatory. Sets the path to the
  // file to load.
  exts.clear();
  exts.emplace_back(".anxs");
  declareProperty(std::make_unique<API::FileProperty>(FilenameStr, "", API::FileProperty::Load, exts),
                  "The input filename of the stored data");

  // mask
  exts.clear();
  exts.emplace_back(".xml");
  declareProperty(std::make_unique<API::FileProperty>(MaskStr, "", API::FileProperty::OptionalLoad, exts),
                  "The input filename of the mask data");

  declareProperty(SelectDetectorTubesStr, "",
                  "Comma separated range of detectors tubes to be loaded,\n"
                  "  eg 16,19-45,47");

  declareProperty(
      std::make_unique<API::WorkspaceProperty<API::IEventWorkspace>>("OutputWorkspace", "", Kernel::Direction::Output));

  declareProperty(SelectDatasetStr, 0, "Select the index for the dataset to be loaded.");

  declareProperty(TOFBiasStr, 0.0, "Time of flight correction in micro-sec.");

  declareProperty(CalibrateTOFStr, false, "Calibrate the TOF correction from the elastic pulse.");

  declareProperty(LambdaOnTwoStr, false, "Instrument is operating in Lambda on Two mode.");

  declareProperty(FilterByTimeStartStr, 0.0,
                  "Only include events after the provided start time, in "
                  "seconds (relative to the start of the run).");

  declareProperty(FilterByTimeStopStr, EMPTY_DBL(),
                  "Only include events before the provided stop time, in "
                  "seconds (relative to the start of the run).");

  std::string grpOptional = "Filters";
  setPropertyGroup(FilterByTimeStartStr, grpOptional);
  setPropertyGroup(FilterByTimeStopStr, grpOptional);
}

/// Creates an event workspace and sets the \p title.
void LoadPLN3::createWorkspace(const std::string &title) {

  // Create the workspace
  m_localWorkspace = std::make_shared<DataObjects::EventWorkspace>();
  m_localWorkspace->initialize(HISTOGRAMS, 2, 1);

  // set the units
  m_localWorkspace->getAxis(0)->unit() = Kernel::UnitFactory::Instance().create("TOF");
  m_localWorkspace->setYUnit("Counts");

  // set title
  m_localWorkspace->setTitle(title);
}

/// <summary>
/// Execute the algorithm through the following:
///   Create the workspace
///   Get the instrument properties and load options
///   Load the instrument from the IDF
///   Load the data values and adjust TOF
///   Set up the masks
/// </summary>
void LoadPLN3::exec() {

  namespace fs = boost::filesystem;

  // Create workspace
  // ----------------
  std::string nxsFile = getPropertyValue(FilenameStr);
  fs::path p = nxsFile;
  for (; !p.extension().empty();)
    p = p.stem();
  std::string title = p.generic_string();
  createWorkspace(title);
  API::LogManager &logManager = m_localWorkspace->mutableRun();
  API::Progress prog(this, 0.0, 1.0, Progress_Total);

  // Get the root entry for the nexs file and get the required dataset
  // TBD dataset index in future nxs?
  NeXus::NXRoot root(nxsFile);
  NeXus::NXEntry nxsEntry = root.openFirstEntry();
  m_datasetIndex = getProperty(SelectDatasetStr);
  logManager.addProperty(SelectDatasetStr, m_datasetIndex);
  std::tie(m_startTime, m_endTime) = getTimeScanLimits(nxsEntry, m_datasetIndex);

  loadParameters(nxsEntry, logManager);
  prog.doReport("creating instrument");
  loadInstrument();

  // Get the region of interest and filters and save to log
  std::string const maskfile = getPropertyValue(MaskStr);
  std::string const seltubes = getPropertyValue(SelectDetectorTubesStr);
  logManager.addProperty(SelectDetectorTubesStr, seltubes);
  logManager.addProperty(MaskStr, maskfile);

  std::vector<bool> roi = createRoiVector(seltubes, maskfile);
  double timeMaxBoundary = getProperty(FilterByTimeStopStr);
  if (isEmpty(timeMaxBoundary))
    timeMaxBoundary = std::numeric_limits<double>::infinity();
  TimeLimits timeBoundary(getProperty(FilterByTimeStartStr), timeMaxBoundary);

  // get the detector map from raw input to a physical detector
  auto instr = m_localWorkspace->getInstrument();
  std::string dmapStr = instr->getParameterAsString("DetectorMap");
  std::vector<size_t> detMapIndex = std::vector<size_t>(HISTO_BINS_X, 0);
  mapRangeToIndex(dmapStr, detMapIndex, [](size_t n) { return n; });

  // Load the events file. First count the number of events to reserve
  // memory and then assign the events to the detectors
  size_t numberHistograms = m_localWorkspace->getNumberHistograms();
  std::vector<EventVector_pt> eventVectors(numberHistograms, nullptr);
  std::vector<size_t> eventCounts(numberHistograms, 0);

  double masterRpm = fabs(logManager.getTimeSeriesProperty<double>("FermiChopperFreq")->firstValue());
  double slaveRpm = fabs(logManager.getTimeSeriesProperty<double>("OverlapChopperFreq")->firstValue());
  double framePeriod = 1.0e6 / masterRpm;

  // if fermi chopper freq equals the overlap freq then the gate period is
  // half the frame period
  double gatePeriod = (std::round(masterRpm / slaveRpm) == 1.0 ? 0.5 * framePeriod : framePeriod);
  AddSinglePointTimeSeriesProperty<double>(logManager, m_startRun, "GatePeriod", gatePeriod);

  // count total events per pixel and reserve necessary memory
  loadDetectorL2Values();
  double sourceSample = fabs(instr->getSource()->getPos().Z());
  double wavelength = logManager.getTimeSeriesProperty<double>("Wavelength")->firstValue();
  double velocity = PhysicalConstants::h / (PhysicalConstants::NeutronMass * wavelength * 1e-10);
  double sampleTime = 1.0e6 * sourceSample / velocity;
  PLN3::EventCounter eventCounter(roi, detMapIndex, framePeriod, gatePeriod, timeBoundary, eventCounts, sourceSample,
                                  velocity, m_detectorL2);
  PLN3::loadEvents(prog, "loading neutron counts", eventCounter, nxsEntry, m_startTime, m_endTime);
  ANSTO::ProgressTracker progTracker(prog, "creating neutron event lists", numberHistograms, Progress_ReserveMemory);
  prepareEventStorage(progTracker, eventCounts, eventVectors);

  AddSinglePointTimeSeriesProperty<int32_t>(logManager, m_startRun, "TotalCounts", (int32_t)eventCounter.numEvents());

  // now perform the actual event collection and TOF convert if necessary
  // if a phase calibration is required then load it as raw doppler time
  // perform the calibration and then convert to TOF
  Types::Core::DateAndTime startTime(m_startRun);
  auto const start_nanosec = startTime.totalNanoseconds();
  bool const calibrateTOF = getProperty(CalibrateTOFStr);
  double tofCorrection = getProperty(TOFBiasStr);
  if (calibrateTOF) {
    tofCorrection = eventCounter.tofCorrection();
  }
  logManager.addProperty("CalibrateTOF", (calibrateTOF ? 1 : 0));
  AddSinglePointTimeSeriesProperty<double>(logManager, m_startRun, "TOFCorrection", tofCorrection);
  PLN3::EventAssigner eventAssigner(roi, detMapIndex, framePeriod, gatePeriod, timeBoundary, eventVectors,
                                    start_nanosec, tofCorrection, sampleTime);
  PLN3::loadEvents(prog, "loading neutron events (TOF)", eventAssigner, nxsEntry, m_startTime, m_endTime);

  // perform a calibration and then TOF conversion if necessary
  // and update the tof limits
  auto minTOF = eventAssigner.tofMin();
  auto maxTOF = eventAssigner.tofMax();

  // just to make sure the bins hold it all and setup the detector masks
  m_localWorkspace->setAllX(HistogramData::BinEdges{std::max(0.0, floor(minTOF)), maxTOF + 1});
  setupDetectorMasks(roi);

  // set log values
  auto frame_count = static_cast<int>(eventCounter.numFrames());
  AddSinglePointTimeSeriesProperty<int>(logManager, m_startRun, "frame_count", frame_count);

  std::string filename = getPropertyValue(FilenameStr);
  logManager.addProperty("filename", filename);

  // Types::Core::time_duration duration = boost::posix_time::microseconds(
  //    static_cast<boost::int64_t>(eventCounter.duration() * 1.0e6));
  Types::Core::time_duration duration =
      boost::posix_time::microseconds(static_cast<boost::int64_t>((m_endTime - m_startTime) * 1.0e-3));
  Types::Core::DateAndTime endTime(startTime + duration);
  logManager.addProperty("start_time", startTime.toISO8601String());
  logManager.addProperty("end_time", endTime.toISO8601String());
  logManager.addProperty<double>("dur", duration.total_milliseconds() * 1.0e-3);

  // Finally add the time-series evironment parameters explicitly
  loadEnvironParameters(nxsEntry, logManager);

  setProperty("OutputWorkspace", m_localWorkspace);
}

/// Recovers the L2 neutronic distance for each detector.
void LoadPLN3::loadDetectorL2Values() {

  m_detectorL2 = std::vector<double>(HISTOGRAMS);
  const auto &detectorInfo = m_localWorkspace->detectorInfo();
  auto detIDs = detectorInfo.detectorIDs();
  for (const auto detID : detIDs) {
    auto ix = detectorInfo.indexOf(detID);
    double l2 = detectorInfo.l2(ix);
    m_detectorL2[detID] = l2;
  }
}

/// Set up the detector masks to the region of interest \p roi.
void LoadPLN3::setupDetectorMasks(std::vector<bool> &roi) {

  // count total number of masked bins
  size_t maskedBins = 0;
  for (size_t i = 0; i != roi.size(); i++)
    if (!roi[i])
      maskedBins++;

  if (maskedBins > 0) {
    // create list of masked bins
    std::vector<size_t> maskIndexList(maskedBins);
    size_t maskIndex = 0;

    for (size_t i = 0; i != roi.size(); i++)
      if (!roi[i])
        maskIndexList[maskIndex++] = i;

    API::IAlgorithm_sptr maskingAlg = createChildAlgorithm("MaskDetectors");
    maskingAlg->setProperty("Workspace", m_localWorkspace);
    maskingAlg->setProperty("WorkspaceIndexList", maskIndexList);
    maskingAlg->executeAsChildAlg();
  }
}

/// Allocate space for the event storage in \p eventVectors after the
/// \p eventCounts have been determined.
void LoadPLN3::prepareEventStorage(ANSTO::ProgressTracker &progTracker, std::vector<size_t> &eventCounts,
                                   std::vector<EventVector_pt> &eventVectors) {

  size_t numberHistograms = eventCounts.size();
  for (size_t i = 0; i != numberHistograms; ++i) {
    DataObjects::EventList &eventList = m_localWorkspace->getSpectrum(i);

    eventList.setSortOrder(DataObjects::PULSETIME_SORT);
    eventList.reserve(eventCounts[i]);

    eventList.setDetectorID(static_cast<detid_t>(i));
    eventList.setSpectrumNo(static_cast<detid_t>(i));

    DataObjects::getEventsFrom(eventList, eventVectors[i]);

    progTracker.update(i);
  }
  progTracker.complete();
}

/// Region of interest is defined by the \p selected detectors and the
/// \p maskfile.
std::vector<bool> LoadPLN3::createRoiVector(const std::string &selected, const std::string &maskfile) {

  std::vector<bool> result(HISTOGRAMS, true);

  // turn off pixels linked to missing tubes
  if (!selected.empty()) {
    std::vector<bool> tubes(HISTO_BINS_X, false);
    mapRangeToIndex(selected, tubes, [](size_t) { return true; });
    for (size_t i = 0; i < DETECTOR_TUBES; i++) {
      if (tubes[i] == false) {
        for (size_t j = 0; j < PIXELS_PER_TUBE; j++) {
          result[i * PIXELS_PER_TUBE + j] = false;
        }
      }
    }
    for (size_t i = 0; i < MONITORS; i++) {
      result[DETECTOR_SPECTRA + i] = tubes[DETECTOR_TUBES + i];
    }
  }

  if (maskfile.length() == 0)
    return result;

  std::ifstream input(maskfile.c_str());
  if (!input.good())
    throw std::invalid_argument("invalid mask file");

  std::string line;
  while (std::getline(input, line)) {
    auto i0 = line.find("<detids>");
    auto iN = line.find("</detids>");

    if ((i0 != std::string::npos) && (iN != std::string::npos) && (i0 < iN)) {
      line = line.substr(i0 + 8, iN - i0 - 8); // 8 = len("<detids>")
      mapRangeToIndex(line, result, [](size_t) { return false; });
    }
  }

  return result;
}

/// Load parameters from input \p nxsFile and save to the log manager, \p logm.
void LoadPLN3::loadParameters(NeXus::NXEntry &entry, API::LogManager &logm) {

  MapNeXusToProperty<std::string>(entry, "sample/name", "unknown", logm, "SampleName", "", 0);
  MapNeXusToProperty<std::string>(entry, "sample/description", "unknown", logm, "SampleDescription", "", 0);

  // if dataset index > 0 need to add an offset to the start time
  Types::Core::DateAndTime startTime(GetNeXusValue<std::string>(entry, "start_time", "2000-01-01T00:00:00", 0));
  if (m_datasetIndex > 0) {
    /*
     * Need to extract the timestamp
     */
    auto baseTime = GetNeXusValue<int32_t>(entry, "instrument/detector/start_time", 0, 0);
    auto nthTime = GetNeXusValue<int32_t>(entry, "instrument/detector/start_time", 0, m_datasetIndex);

    Types::Core::time_duration duration =
        boost::posix_time::microseconds(static_cast<boost::int64_t>((nthTime - baseTime) * 1.0e6));
    Types::Core::DateAndTime startDataset(startTime + duration);
    m_startRun = startDataset.toISO8601String();
  } else {
    m_startRun = startTime.toISO8601String();
  }

  // Add support for instrument running in lambda on two mode.
  // Added as UI option as the available instrument parameters
  // cannot be reliably interpreted to predict the mode (as
  // advised by the instrument scientist).
  bool const lambdaOnTwoMode = getProperty(LambdaOnTwoStr);
  float lambdaFactor = (lambdaOnTwoMode ? 0.5f : 1.0f);
  logm.addProperty("LambdaOnTwoMode", (lambdaOnTwoMode ? 1 : 0));

  MapNeXusToSeries<double>(entry, "instrument/fermi_chopper/mchs", 6000.0, logm, m_startRun, "FermiChopperFreq",
                           1.0 / 60, 0);
  MapNeXusToSeries<double>(entry, "instrument/fermi_chopper/schs", 3000.0, logm, m_startRun, "OverlapChopperFreq",
                           1.0 / 60, 0);
  MapNeXusToSeries<float>(entry, "instrument/crystal/wavelength/value", 0.0, logm, m_startRun, "Wavelength",
                          lambdaFactor, 0);
  MapNeXusToSeries<float>(entry, "instrument/detector/stth/value", 0.0, logm, m_startRun, "DetectorTankAngle", 1.0, 0);
  MapNeXusToSeries<float>(entry, "instrument/detector/tofw/value", 5.0, logm, m_startRun, "ChannelWidth", 1, 0);
  MapNeXusToSeries<float>(entry, "sample/mscor/value", 0.0, logm, m_startRun, "SampleRotation", 1, 0);

  auto bmc = beamMonitorCounts(entry, "monitor/bm1_counts", m_startTime, m_endTime, false);
  AddSinglePointTimeSeriesProperty<int32_t>(logm, m_startRun, "MonitorCounts", (int32_t)bmc);
}

/// Load the environment variables from the \p nxsFile and save as
/// time series to the log manager, \p logm.
void LoadPLN3::loadEnvironParameters(NeXus::NXEntry &entry, API::LogManager &logm) {

  auto time_str = logm.getPropertyValueAsType<std::string>("end_time");

  // load the environment variables for the dataset loaded
  std::vector<std::string> tags = {"P01PS05", "P01PSP05", "T01S00", "T01S06", "T01S07", "T01S08", "T01SP00", "T01SP06",
                                   "T01SP07", "T01SP08",  "T2S1",   "T2S2",   "T2S3",   "T2S4",   "T2SP1",   "T2SP2"};

  for (const auto &tag : tags) {
    MapNeXusToSeries<float>(entry, "control/" + tag + "/value", 0.0, logm, time_str, "env_" + tag, 1.0, m_datasetIndex);
  }
}

/// Load the instrument definition.
void LoadPLN3::loadInstrument() {

  // loads the IDF and parameter file
  API::IAlgorithm_sptr loadInstrumentAlg = createChildAlgorithm("LoadInstrument");
  loadInstrumentAlg->setProperty("Workspace", m_localWorkspace);
  loadInstrumentAlg->setPropertyValue("InstrumentName", "PELICAN");
  loadInstrumentAlg->setProperty("RewriteSpectraMap", Mantid::Kernel::OptionalBool(false));
  loadInstrumentAlg->executeAsChildAlg();
}

/// Algorithm's version for identification. @see Algorithm::version
int LoadPLN3::version() const { return 1; }

/// Similar algorithms. @see Algorithm::seeAlso
const std::vector<std::string> LoadPLN3::seeAlso() const { return {"Load", "LoadEMU"}; }
/// Algorithm's category for identification. @see Algorithm::category
const std::string LoadPLN3::category() const { return "DataHandling\\Nexus;DataHandling\\ANSTO "; }

/// Algorithms name for identification. @see Algorithm::name
const std::string LoadPLN3::name() const { return "LoadPLN3"; }

/// Algorithm's summary for use in the GUI and help. @see Algorithm::summary
const std::string LoadPLN3::summary() const { return "Loads an ANSTO nexus like file into a workspace."; }

/// Return the confidence as an integer value that this algorithm can
/// load the file \p descriptor.
int LoadPLN3::confidence(Kernel::NexusDescriptor &descriptor) const {
  if (descriptor.extension() != ".nxs")
    return 0;

  if (descriptor.pathExists("/entry1/program_name") && descriptor.pathExists("/entry1/experiment/gumtree_version") &&
      descriptor.pathExists("/entry1/instrument/detector_events") &&
      descriptor.pathExists("/entry1/instrument/aperture/sh1") &&
      descriptor.pathExists("/entry1/instrument/ag1010/MEAS/Temperature/value") &&
      descriptor.pathExists("/entry1/instrument/detector/daq_dirname") &&
      descriptor.pathExists("/entry1/instrument/detector/dataset_number/value") &&
      descriptor.pathExists("/entry1/scan_dataset/value") && descriptor.pathExists("/entry1/monitor/bm1_counts") &&
      descriptor.pathExists("/entry1/monitor/time")) {
    return 95;
  } else {
    return 0;
  }
}

// register the algorithms into the AlgorithmFactory
DECLARE_NEXUS_FILELOADER_ALGORITHM(LoadPLN3)

} // namespace DataHandling
} // namespace Mantid
