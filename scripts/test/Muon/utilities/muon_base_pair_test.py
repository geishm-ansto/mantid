# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
import unittest

from Muon.GUI.Common.muon_base_pair import MuonBasePair
from Muon.GUI.Common.ADSHandler.muon_workspace_wrapper import MuonWorkspaceWrapper
from mantid.simpleapi import CreateWorkspace


class MuonBasePairTest(unittest.TestCase):
    """
    The MuonPair object encapsulates the information the describes a pair:

    - Its name
    - The name of the two groups that are combined to form the pair
    - The value of alpha
    - A workspace (Optional)

    It is intended to be lightweight, but with some type checking and basic logic, intended
    to prevent obvious misuse. The class is used by the MuonAnalysis interface and so these
    tests should be amended with care.
    """

    def test_that_cannot_initialize_MuonBasePair_without_name(self):
        with self.assertRaises(TypeError):
            MuonBasePair()

    def test_that_MuonBasePair_name_is_set_correctly(self):
        pair = MuonBasePair(pair_name="pair1")
        self.assertEqual(pair.name, "pair1")

    def test_that_cannot_set_new_name_on_pair(self):
        pair = MuonBasePair(pair_name="pair1")
        with self.assertRaises(AttributeError):
            pair.name = "new_name"
        self.assertEqual(pair.name, "pair1")

    def test_that_can_only_set_workspace_if_MuonWorkspace_object(self):
        pair = MuonBasePair(pair_name="pair1")
        self.assertEqual(pair.workspace, {})
        dataX = [0, 1, 2, 3, 4, 5]
        dataY = [10, 20, 30, 20, 10]
        input_workspace = CreateWorkspace(dataX, dataY)

        workspace_wrapper = MuonWorkspaceWrapper(input_workspace)
        pair.workspace = workspace_wrapper
        self.assertEqual(pair.workspace, workspace_wrapper)

    def test_that_AttributeError_thrown_if_setting_workspace_to_non_MuonWorkspace_object(self):
        pair = MuonBasePair(pair_name="pair1")

        self.assertEqual(pair.workspace, {})
        with self.assertRaises(AttributeError):
            pair.workspace = [1, 2, 3]
        self.assertEqual(pair.workspace, {})



if __name__ == '__main__':
    unittest.main(buffer=False, verbosity=2)
