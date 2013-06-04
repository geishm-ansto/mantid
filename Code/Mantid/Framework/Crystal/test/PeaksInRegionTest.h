#ifndef MANTID_CRYSTAL_PEAKSINREGIONTEST_H_
#define MANTID_CRYSTAL_PEAKSINREGIONTEST_H_

#include <cxxtest/TestSuite.h>
#include "MantidTestHelpers/WorkspaceCreationHelper.h"
#include "MantidTestHelpers/ComponentCreationHelper.h"
#include "MantidDataObjects/PeaksWorkspace.h"
#include "MantidCrystal/PeaksInRegion.h"
#include <boost/tuple/tuple.hpp>

using Mantid::Crystal::PeaksInRegion;
using namespace Mantid::API;
using namespace Mantid::DataObjects;

/*-------------------------------------------------------------------------------------------------------------------------------------------------------
Functional Tests
-------------------------------------------------------------------------------------------------------------------------------------------------------*/
class PeaksInRegionTest : public CxxTest::TestSuite
{

private:

  typedef boost::tuple<PeaksWorkspace_sptr, std::vector<double> > PeakWorkspaceWithExtents;

  /**
  Helper function. Creates a peaksworkspace with a single peak 
  */
  PeakWorkspaceWithExtents createPeaksWorkspace(const std::string coordFrame, double xMinFromPeak, double xMaxFromPeak, double yMinFromPeak, double yMaxFromPeak, double zMinFromPeak, double zMaxFromPeak)
  {
    PeaksWorkspace_sptr ws = WorkspaceCreationHelper::createPeaksWorkspace(1);
    auto detectorIds = ws->getInstrument()->getDetectorIDs();
    Peak& peak = ws->getPeak(0);
    peak.setDetectorID(detectorIds.front());
    Mantid::Kernel::V3D position;
    if(coordFrame == "Detector space")
    {
      position = peak.getDetector()->getPos();
    }
    else if(coordFrame == "Q (lab frame)")
    {
      position = peak.getQLabFrame();
    }
    else if(coordFrame == "Q (sample frame)")
    {
      position = peak.getQSampleFrame();
    }
    else if(coordFrame == "HKL")
    {
      position = peak.getHKL();
    }
    else
    {
      throw std::runtime_error("Unknown coordinate frame");
    }
    std::vector<double> extents(6);
    extents[0] = position.X() - xMinFromPeak;
    extents[1] = position.X() + xMaxFromPeak;
    extents[2] = position.Y() - yMinFromPeak;
    extents[3] = position.Y() + yMaxFromPeak;
    extents[4] = position.Z() - zMinFromPeak;
    extents[5] = position.Z() + zMaxFromPeak;
    return PeakWorkspaceWithExtents(ws, extents);
  }

public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static PeaksInRegionTest *createSuite() { return new PeaksInRegionTest(); }
  static void destroySuite( PeaksInRegionTest *suite ) { delete suite; }

  void test_setProperties()
  {
    PeaksInRegion alg;
    alg.setRethrows(true);
    alg.initialize() ;
    TS_ASSERT( alg.isInitialized() ) ;
    alg.setProperty("InputWorkspace", WorkspaceCreationHelper::createPeaksWorkspace());
    alg.setPropertyValue("CoordinateFrame", "Q (lab frame)");
    alg.setPropertyValue("Extents", "-1,1,-1,1,-1,1");
    alg.setPropertyValue("OutputWorkspace", "OutWS");
  }

  void test_badExtents()
  {
    PeaksInRegion alg;
    alg.setRethrows(true);
    alg.initialize() ;
    TS_ASSERT( alg.isInitialized() ) ;
    alg.setProperty("InputWorkspace", WorkspaceCreationHelper::createPeaksWorkspace());
    alg.setPropertyValue("CoordinateFrame", "Q (lab frame)");
    alg.setPropertyValue("Extents", "-1,1,-1,1,-1,1, -1"); // Too many
    alg.setPropertyValue("OutputWorkspace", "OutWS");

    TS_ASSERT_THROWS( alg.execute(), std::invalid_argument&);

    alg.setPropertyValue("Extents", "-1"); // Too few

    TS_ASSERT_THROWS( alg.execute(), std::invalid_argument&);
  }

  void do_test_within_bounds(const std::string coordFrame)
  {
    const std::string outName = "OutWS";
    const double xMinFromPeak = 1;
    const double xMaxFromPeak = 1;
    const double yMinFromPeak = 1;
    const double yMaxFromPeak = 1;
    const double zMinFromPeak = 1;
    const double zMaxFromPeak = 1;

    PeakWorkspaceWithExtents tuple = createPeaksWorkspace(coordFrame, xMinFromPeak, xMaxFromPeak, yMinFromPeak, yMaxFromPeak, zMinFromPeak, zMaxFromPeak);

    PeaksInRegion alg;
    alg.setRethrows(true);
    alg.initialize();
    alg.setProperty("InputWorkspace", tuple.get<0>());
    alg.setPropertyValue("CoordinateFrame", coordFrame);
    alg.setProperty("Extents", tuple.get<1>()); 
    alg.setPropertyValue("OutputWorkspace", outName);
    alg.execute();
    
    ITableWorkspace_sptr outWS = AnalysisDataService::Instance().retrieveWS<ITableWorkspace>(outName);

    TS_ASSERT_EQUALS(2, outWS->columnCount());
    TS_ASSERT_EQUALS("PeakIndex", outWS->getColumn(0)->name());
    TS_ASSERT_EQUALS("Intersecting", outWS->getColumn(1)->name());

    TS_ASSERT_EQUALS(1, outWS->rowCount());

    TSM_ASSERT_EQUALS("Peak index should be zero", 0, outWS->cell<int>(0,  0)); 
    TSM_ASSERT_EQUALS("Peak intersect should be true", Boolean(true), outWS->cell<Boolean>(0,  1));
  }

  void do_test_out_of_bounds(const std::string coordFrame, double xMinFromPeak, double xMaxFromPeak, double yMinFromPeak, double yMaxFromPeak, double zMinFromPeak, double zMaxFromPeak)
  {
    const std::string outName = "OutWS";

    PeakWorkspaceWithExtents tuple = createPeaksWorkspace(coordFrame, xMinFromPeak, xMaxFromPeak, yMinFromPeak, yMaxFromPeak, zMinFromPeak, zMaxFromPeak);

    PeaksInRegion alg;
    alg.setRethrows(true);
    alg.initialize();
    alg.setProperty("InputWorkspace", tuple.get<0>());
    alg.setPropertyValue("CoordinateFrame", coordFrame);
    alg.setProperty("Extents", tuple.get<1>()); 
    alg.setPropertyValue("OutputWorkspace", outName);
    alg.execute();
    
    ITableWorkspace_sptr outWS = AnalysisDataService::Instance().retrieveWS<ITableWorkspace>(outName);

    TS_ASSERT_EQUALS(2, outWS->columnCount());
    TS_ASSERT_EQUALS("PeakIndex", outWS->getColumn(0)->name());
    TS_ASSERT_EQUALS("Intersecting", outWS->getColumn(1)->name());

    TS_ASSERT_EQUALS(1, outWS->rowCount());

    TSM_ASSERT_EQUALS("Peak index should be zero", 0, outWS->cell<int>(0,  0)); 
    TSM_ASSERT_EQUALS("Peak intersect should be false", Boolean(false), outWS->cell<Boolean>(0,  1));
  }

  void test_detectorSpace_with_peak_in_bounds()
  {
    do_test_within_bounds("Detector space");
  }

  void test_qLab_with_peak_in_bounds()
  {
    do_test_within_bounds("Q (lab frame)");
  }

  void test_qSample_with_peak_in_bounds()
  {
    do_test_within_bounds("Q (sample frame)");
  }

  void test_HKL_with_peak_in_bounds()
  {
    do_test_within_bounds("HKL");
  }

  void test_detectorSpace_with_peaks_out_of_bounds()
  {
    const std::string coordinateFrame = "Detector space";
    do_test_out_of_bounds(coordinateFrame, -1, 1, 1, 1, 1, 1); // outside xmin
    do_test_out_of_bounds(coordinateFrame, 1, -1, 1, 1, 1, 1); // outside xmax
    do_test_out_of_bounds(coordinateFrame, 1, 1, -1, 1, 1, 1); // outside ymin
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, -1, 1, 1); // outside ymax
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, 1, -1, 1); // outside zmin
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, 1, 1, -1); // outside zmax
  }

  void test_qLab_with_peaks_out_of_bounds()
  {
    const std::string coordinateFrame = "Q (lab frame)";
    do_test_out_of_bounds(coordinateFrame, -1, 1, 1, 1, 1, 1); // outside xmin
    do_test_out_of_bounds(coordinateFrame, 1, -1, 1, 1, 1, 1); // outside xmax
    do_test_out_of_bounds(coordinateFrame, 1, 1, -1, 1, 1, 1); // outside ymin
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, -1, 1, 1); // outside ymax
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, 1, -1, 1); // outside zmin
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, 1, 1, -1); // outside zmax
  }

 void test_qSample_with_peaks_out_of_bounds()
  {
    const std::string coordinateFrame = "Q (sample frame)";
    do_test_out_of_bounds(coordinateFrame, -1, 1, 1, 1, 1, 1); // outside xmin
    do_test_out_of_bounds(coordinateFrame, 1, -1, 1, 1, 1, 1); // outside xmax
    do_test_out_of_bounds(coordinateFrame, 1, 1, -1, 1, 1, 1); // outside ymin
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, -1, 1, 1); // outside ymax
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, 1, -1, 1); // outside zmin
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, 1, 1, -1); // outside zmax
  }

  void test_qHKL_with_peaks_out_of_bounds()
  {
    const std::string coordinateFrame = "HKL";
    do_test_out_of_bounds(coordinateFrame, -1, 1, 1, 1, 1, 1); // outside xmin
    do_test_out_of_bounds(coordinateFrame, 1, -1, 1, 1, 1, 1); // outside xmax
    do_test_out_of_bounds(coordinateFrame, 1, 1, -1, 1, 1, 1); // outside ymin
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, -1, 1, 1); // outside ymax
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, 1, -1, 1); // outside zmin
    do_test_out_of_bounds(coordinateFrame, 1, 1, 1, 1, 1, -1); // outside zmax
  }

};


/*-------------------------------------------------------------------------------------------------------------------------------------------------------
Functional Tests
-------------------------------------------------------------------------------------------------------------------------------------------------------*/
class PeaksInRegionTestPerformance : public CxxTest::TestSuite
{

private:

  Mantid::API::IPeaksWorkspace_sptr inputWS;

public:

  static PeaksInRegionTestPerformance *createSuite() { return new PeaksInRegionTestPerformance(); }
  static void destroySuite( PeaksInRegionTestPerformance *suite ) { delete suite; }

  PeaksInRegionTestPerformance()
  {
    int numPeaks = 4000;
    inputWS = boost::make_shared<PeaksWorkspace>();
    Mantid::Geometry::Instrument_sptr inst = ComponentCreationHelper::createTestInstrumentRectangular2(1, 200);
    inputWS->setInstrument(inst);

    for (int i = 0; i < numPeaks; ++i)
    {
      Peak peak(inst, i, i+0.5);
      inputWS->addPeak(peak);
    }
  }

  void testPerformance()
  {
    const std::string outName = "OutPerfWS";

    PeaksInRegion alg;
    alg.setRethrows(true);
    alg.initialize() ;
    TS_ASSERT( alg.isInitialized() ) ;
    alg.setProperty("InputWorkspace", inputWS);
    alg.setPropertyValue("CoordinateFrame", "Detector space");
    alg.setPropertyValue("Extents", "-1,1,-1,1,-1,1");
    alg.setPropertyValue("OutputWorkspace", outName);
    alg.execute();

    Mantid::API::ITableWorkspace_sptr outWS = AnalysisDataService::Instance().retrieveWS<ITableWorkspace>(outName);

    TS_ASSERT_EQUALS(2, outWS->columnCount());
    TS_ASSERT_EQUALS(inputWS->rowCount(), outWS->rowCount());
  }

};


#endif /* MANTID_CRYSTAL_PEAKSINREGIONTEST_H_ */