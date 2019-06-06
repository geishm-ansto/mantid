# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
""" The table  model contains all the reduction information which is provided via the data table

The main information in the table model are the run numbers and the selected periods. However it also contains
information regarding the custom output name and the information in the options tab.
"""

from __future__ import (absolute_import, division, print_function)

import functools
import os
import re

from sans.common.file_information import find_sans_file
from sans.ansto import table_model as core
from sans.ansto.bilby.file_information import BILBYFileInformation

class TableModel(core.TableModel):

    def create_ansto_file_information(self, entry):
        sample_path = find_sans_file(entry.sample)
        trans_path = find_sans_file(entry.t_sample)
        file_info = BILBYFileInformation(sample_path, trans_path)
        return file_info

    def update_row_from_file_information(self, id, file_information):
        row = self.get_row_from_id(id)
        if file_information:
            # no specific values 
            #rounded_file_thickness = round(file_information.get_thickness(), 2)
            #self._table_entries[row].update_attribute('thickness', rounded_file_thickness)
            pass
        super(TableModel, self).update_row_from_file_information(id, file_information)

class RowModel(core.RowModel):

    @staticmethod
    def column_labels():
        return ["Index","Sample","T Empty Beam","Transmission Sample","Thickness", "Blocked Beam", 
                "Start Time", "End Time","Sample Mask", "Transmission Mask",
                "Gravity Fix","Wide Angle Fix", "Radius Cut","Wave Cut","Suffix","Description"]

    @staticmethod
    def column_keys():
        return ["index", "sample", "t_empty_beam", "t_sample", "thickness", "blocked_beam", 
                "start_time", "end_time", "sample_mask", "transmission_mask", 
                "gravity_correction","wide_angle_correction", "radius_cut", "wave_cut",
                "suffix", "description", ]

    @staticmethod
    def column_options():
        return {'index': [0],
                'time_limits': [6,7],
                'corrections': [10,11],
                'radius_wave': [12,13],
                'masks': [8,9]}
    
    @staticmethod
    def create_empty_row():
        row = [''] * 16
        return RowModel(*row)

    def __init__(self, *argv):
        self.sample = None
        self.t_empty_beam = None
        self.t_sample = None
        self.thickness = None
        self.blocked_beam = None
        self.suffix = None
        self.description = None
        self.start_time = None
        self.end_time = None
        self.gravity_correction = None
        self.wide_angle_correction = None
        self.sample_mask = None
        self.transmission_mask = None
        self.radius_cut = None
        self.wave_cut = None
        super(RowModel, self).__init__(*argv)


