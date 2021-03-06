// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidKernel/System.h"
#include "MantidKernel/Timer.h"
#include <cxxtest/TestSuite.h>

#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/Sample.h"
#include "MantidCrystal/LoadIsawUB.h"
#include "MantidCrystal/SelectCellOfType.h"
#include "MantidDataHandling/LoadNexusProcessed.h"
#include "MantidDataObjects/LeanElasticPeaksWorkspace.h"
#include "MantidDataObjects/PeaksWorkspace.h"
#include "MantidGeometry/Crystal/IndexingUtils.h"
#include "MantidGeometry/Crystal/OrientedLattice.h"

using namespace Mantid::Crystal;
using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::DataHandling;
using namespace Mantid::Kernel;
using namespace Mantid::Geometry;

class SelectCellOfTypeTest : public CxxTest::TestSuite {
public:
  void test_Init() {
    SelectCellOfType alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize())
    TS_ASSERT(alg.isInitialized())
  }

  void test_exec() {
    // Name of the loader's output workspace.
    std::string WSName("peaks");
    LoadNexusProcessed loader;
    TS_ASSERT_THROWS_NOTHING(loader.initialize());
    TS_ASSERT(loader.isInitialized());
    loader.setPropertyValue("Filename", "TOPAZ_3007.peaks.nxs");
    loader.setPropertyValue("OutputWorkspace", WSName);

    TS_ASSERT(loader.execute());
    TS_ASSERT(loader.isExecuted());
    PeaksWorkspace_sptr ws;
    TS_ASSERT_THROWS_NOTHING(
        ws = std::dynamic_pointer_cast<PeaksWorkspace>(AnalysisDataService::Instance().retrieve(WSName)));
    TS_ASSERT(ws);
    // set a Niggli UB for run 3007
    // (CuTCA) in the oriented lattice
    V3D row_0(0.0122354, 0.00480056, 0.0860404);
    V3D row_1(-0.1165450, 0.00178145, -0.0045884);
    V3D row_2(-0.0273738, -0.08973560, -0.0252595);

    Matrix<double> UB(3, 3, false);
    UB.setRow(0, row_0);
    UB.setRow(1, row_1);
    UB.setRow(2, row_2);

    auto lattice = std::make_unique<OrientedLattice>();
    lattice->setUB(UB);
    ws->mutableSample().setOrientedLattice(std::move(lattice));

    // now get the UB back from the WS
    UB = ws->sample().getOrientedLattice().getUB();

    SelectCellOfType alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize())
    TS_ASSERT(alg.isInitialized())
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("PeaksWorkspace", WSName));
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("CellType", "Monoclinic"));
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("Centering", "P"));
    TS_ASSERT_THROWS_NOTHING(alg.setProperty("Apply", true));
    TS_ASSERT_THROWS_NOTHING(alg.setProperty("Tolerance", 0.12));
    TS_ASSERT_THROWS_NOTHING(alg.execute(););
    TS_ASSERT(alg.isExecuted());

    int num_indexed = alg.getProperty("NumIndexed");
    TS_ASSERT_EQUALS(num_indexed, 43);
    double average_error = alg.getProperty("AverageError");
    TS_ASSERT_DELTA(average_error, 0.00972862, .0001);

    AnalysisDataService::Instance().remove(WSName);
  }

  void test_exec_LeanElasticPeak() {
    // Name of the loader's output workspace.
    std::string WSName("peaks");
    LoadNexusProcessed loader;
    TS_ASSERT_THROWS_NOTHING(loader.initialize());
    TS_ASSERT(loader.isInitialized());
    loader.setPropertyValue("Filename", "TOPAZ_3007.peaks.nxs");
    loader.setPropertyValue("OutputWorkspace", WSName);

    TS_ASSERT(loader.execute());
    TS_ASSERT(loader.isExecuted());
    PeaksWorkspace_sptr ws;
    TS_ASSERT_THROWS_NOTHING(
        ws = std::dynamic_pointer_cast<PeaksWorkspace>(AnalysisDataService::Instance().retrieve(WSName)));
    TS_ASSERT(ws);

    // Convert PeaksWorkspace to LeanElasticPeaksWorkspace
    auto lpw = std::make_shared<LeanElasticPeaksWorkspace>();
    for (auto peak : ws->getPeaks())
      lpw->addPeak(peak);
    AnalysisDataService::Instance().addOrReplace(WSName, lpw);

    // set a Niggli UB for run 3007
    // (CuTCA) in the oriented lattice
    V3D row_0(0.0122354, 0.00480056, 0.0860404);
    V3D row_1(-0.1165450, 0.00178145, -0.0045884);
    V3D row_2(-0.0273738, -0.08973560, -0.0252595);

    Matrix<double> UB(3, 3, false);
    UB.setRow(0, row_0);
    UB.setRow(1, row_1);
    UB.setRow(2, row_2);

    auto lattice = std::make_unique<OrientedLattice>();
    lattice->setUB(UB);
    lpw->mutableSample().setOrientedLattice(std::move(lattice));

    // now get the UB back from the WS
    UB = lpw->sample().getOrientedLattice().getUB();

    SelectCellOfType alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize())
    TS_ASSERT(alg.isInitialized())
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("PeaksWorkspace", WSName));
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("CellType", "Monoclinic"));
    TS_ASSERT_THROWS_NOTHING(alg.setPropertyValue("Centering", "P"));
    TS_ASSERT_THROWS_NOTHING(alg.setProperty("Apply", true));
    TS_ASSERT_THROWS_NOTHING(alg.setProperty("Tolerance", 0.12));
    TS_ASSERT_THROWS_NOTHING(alg.execute(););
    TS_ASSERT(alg.isExecuted());

    int num_indexed = alg.getProperty("NumIndexed");
    TS_ASSERT_EQUALS(num_indexed, 43);
    double average_error = alg.getProperty("AverageError");
    TS_ASSERT_DELTA(average_error, 0.00972862, .0001);

    AnalysisDataService::Instance().remove(WSName);
  }
};
