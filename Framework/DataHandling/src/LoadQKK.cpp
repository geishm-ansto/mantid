// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataHandling/LoadQKK.h"

#include "MantidAPI/Axis.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/RegisterFileLoader.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidDataObjects/WorkspaceCreation.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/RectangularDetector.h"
#include "MantidGeometry/Objects/ShapeFactory.h"
#include "MantidIndexing/IndexInfo.h"
#include "MantidKernel/OptionalBool.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidNexus/NexusClasses.h"

#include <Poco/File.h>

#include <fstream>
#include <regex>

using namespace Mantid::DataHandling;
using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::NeXus;

namespace {

template <typename TYPE>
void AddSinglePointTimeSeriesProperty(LogManager &logManager, const std::string &time, const std::string &name,
                                      const TYPE value) {
  // create time series property and add single value
  auto p = new TimeSeriesProperty<TYPE>(name);
  p->addValue(time, value);

  // add to log manager
  logManager.addProperty(p);
}

int GetRunNumber(NXEntry &entry) {
  std::string name = entry.name();
  const std::regex pattern("^QKK0+(?!$)");
  auto str = std::regex_replace(name, pattern, "");
  return stoi(str);
}

} // namespace

namespace Mantid {
namespace DataHandling {

// Register the algorithm into the AlgorithmFactory
DECLARE_NEXUS_FILELOADER_ALGORITHM(LoadQKK)

/**
 * Return the confidence with with this algorithm can load the file
 * @param descriptor A descriptor for the file
 * @returns An integer specifying the confidence level. 0 indicates it will not
 * be used
 */
int LoadQKK::confidence(Kernel::NexusDescriptor &descriptor) const {
  const auto &firstEntryName = descriptor.firstEntryNameType().first;
  std::regex qkkRE("QKK([0-9]*)");
  if (std::regex_match(firstEntryName, qkkRE) && descriptor.pathExists("/" + firstEntryName + "/data/hmm_xy"))
    return 90;
  else
    return 0;
}

/**
 * Initialise the algorithm. Declare properties which can be set before
 * execution (input) or
 * read from after the execution (output).
 */
void LoadQKK::init() {
  // Declare the Filename algorithm property. Mandatory. Sets the path to the
  // file to load.
  declareProperty(std::make_unique<API::FileProperty>("Filename", "", API::FileProperty::Load, ".nx.hdf"),
                  "The input filename of the stored data");
  // Declare the OutputWorkspace property. This sets the name of the workspace
  // to be filled with the data
  // from the file.
  declareProperty(std::make_unique<API::WorkspaceProperty<>>("OutputWorkspace", "", Kernel::Direction::Output));
}

/**
 * Execute the algorithm.
 */
void LoadQKK::exec() {
  using namespace Mantid::API;
  // Get the name of the file.
  std::string filename = getPropertyValue("Filename");

  // Open the root.
  NeXus::NXRoot root(filename);
  // Open the first NXentry found in the file.
  NeXus::NXEntry entry = root.openFirstEntry();

  createWorkspace(entry);
  setFinalProperties(filename);
  loadMetaData(entry);
  loadInstrument();

  // add the ILL data from here for now, ideally need a transformation package
  loadILLMetaData(entry);

  setProperty("OutputWorkspace", m_localWorkspace);
}

void LoadQKK::createWorkspace(NeXus::NXEntry &entry) {
  // Naturally, the relevant data is in the group "data"
  NeXus::NXData data = entry.openNXData("data");
  double wavelength = static_cast<double>(data.getFloat("wavelength"));
  // open and load the data set with the counts
  NeXus::NXInt hmm = data.openIntData();
  hmm.load();

  // Get the wavelength spread
  double wavelength_spread = static_cast<double>(entry.getFloat("instrument/velocity_selector/wavelength_spread"));
  double wavelength0 = wavelength - wavelength_spread / 2;
  double wavelength1 = wavelength + wavelength_spread / 2;

  // hmm is a 3d array with axes: sample_x, y_pixel_offset, x_pixel_offset
  size_t ny = hmm.dim1(); // second dimension
  size_t nx = hmm.dim2(); // third dimension
  size_t nHist = ny * nx; // number of spectra in the dataset
  if (nHist == 0) {
    throw std::runtime_error("Error in data dimensions: " + std::to_string(ny) + " X " + std::to_string(nx));
  }

  // Create a workspace with nHist spectra and a single y bin.
  m_localWorkspace =
      DataObjects::create<DataObjects::Workspace2D>(Indexing::IndexInfo(nHist), HistogramData::BinEdges(2));
  // Set the units of the x axis as Wavelength
  m_localWorkspace->getAxis(0)->unit() = UnitFactory::Instance().create("Wavelength");
  // Set the units of the data as Counts
  m_localWorkspace->setYUnitLabel("Counts");

  using namespace HistogramData;
  const BinEdges binEdges = {wavelength0, wavelength1};
  for (size_t index = 0; index < nHist; ++index) {
    auto x = static_cast<int>(index % nx);
    auto y = static_cast<int>(index / nx);
    auto c = hmm(0, x, y);

    Counts yValue = {static_cast<double>(c)};
    m_localWorkspace->setHistogram(index, binEdges, yValue);
  }
}

/// Load the instrument definition.
void LoadQKK::loadInstrument() {

  // loads the IDF and parameter file
  auto loadInstrumentAlg = createChildAlgorithm("LoadInstrument");
  loadInstrumentAlg->setProperty("Workspace", m_localWorkspace);
  loadInstrumentAlg->setPropertyValue("InstrumentName", "QUOKKA");
  loadInstrumentAlg->setProperty("RewriteSpectraMap", Mantid::Kernel::OptionalBool(true));
  loadInstrumentAlg->executeAsChildAlg();
}

void LoadQKK::setFinalProperties(const std::string &filename) {
  API::Run &runDetails = m_localWorkspace->mutableRun();
  NXhandle nxHandle;
  NXstatus nxStat = NXopen(filename.c_str(), NXACC_READ, &nxHandle);

  if (nxStat != NX_ERROR) {
    m_loadHelper.addNexusFieldsToWsRun(nxHandle, runDetails);
    NXclose(&nxHandle);
  }
}

void LoadQKK::loadMetaData(NeXus::NXEntry &entry) {
  // the IDF logfile parameter are required to be time series entries
  NeXus::NXData params = entry.openNXData("instrument/parameters");
  double offsetX = params.getFloat("BeamCenterX");
  double offsetY = params.getFloat("BeamCenterZ");
  double offsetL1 = params.getFloat("L1");
  double offsetL2 = params.getFloat("L2");
  API::LogManager &logManager = m_localWorkspace->mutableRun();

  std::string startDate = entry.getString("start_time");
  // std::string isoStart = m_loadHelper.dateTimeInIsoFormat(startDate);

  AddSinglePointTimeSeriesProperty(logManager, startDate, "L1", offsetL1);
  AddSinglePointTimeSeriesProperty(logManager, startDate, "L2", offsetL2);
  AddSinglePointTimeSeriesProperty(logManager, startDate, "BeamCenterX", offsetX);
  AddSinglePointTimeSeriesProperty(logManager, startDate, "BeamCenterY", offsetY);

  // end_time is not set, so add duration to the start time
  API::Run &runDetails = m_localWorkspace->mutableRun();
  double duration = entry.openNXData("instrument/detector").getFloat("time");
  Types::Core::DateAndTime startTime(entry.getString("start_time"));
  auto endTime = startTime + duration;
  runDetails.addProperty<std::string>("start_time", startTime.toISO8601String(), true);
  runDetails.addProperty<double>("duration", duration, true);
  runDetails.addProperty<std::string>("end_time", endTime.toISO8601String(), true);

  // set the facility
  runDetails.addProperty<std::string>("Facility", std::string("ANSTO"));
}

void LoadQKK::loadILLMetaData(NeXus::NXEntry &entry) {
  g_log.debug("Loading ILL metadata...");
  API::Run &runDetails = m_localWorkspace->mutableRun();
  int runNumber = entry.openNXData("instrument").getInt("run_number");
  if (runNumber == 0) {
    runNumber = GetRunNumber(entry);
  }
  runDetails.addProperty<int>("run_number", runNumber);

  // now common SANS components
  runDetails.addProperty<std::string>("sample_description", entry.openNXData("sample").getString("description"), true);
  runDetails.addProperty<std::string>("instrument_name", "QUOKKA", true);
  runDetails.addProperty<double>("sample.thickness", 0.1 * entry.openNXData("sample").getFloat("SampleThickness"), "cm",
                                 true);
  runDetails.addProperty<int>("monitor1.data", entry.openNXData("monitor").getInt("bm1_counts"), true);

  double wavelength = entry.openNXData("instrument/velocity_selector").getFloat("wavelength");
  double wvSpread = entry.openNXData("instrument/velocity_selector").getFloat("wavelength_spread");
  double wavelengthRes = 100 * wvSpread / wavelength;
  // round also the wavelength res to avoid unnecessary rebinning during
  // merge runs
  wavelengthRes = std::round(wavelengthRes * 100) / 100.;
  runDetails.addProperty<double>("selector.wavelength", wavelength, "A", true);
  runDetails.addProperty<double>("selector.wavelength_res", wavelengthRes, "%", true);

  // collimation data and attenuation

  // sample aperture
  double apx = entry.openNXData("instrument/parameters").getFloat("EApX");
  double apy = entry.openNXData("instrument/parameters").getFloat("EApZ");
  runDetails.addProperty<double>("Beam.sample_ap_x_or_diam", apx, "mm", true);
  runDetails.addProperty<double>("Beam.sample_ap_y", apy, "mm", true);
}

} // namespace DataHandling
} // namespace Mantid
