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

from sans.ansto import table_model as core

class RowModel(core.RowModel):

    @staticmethod
    def column_labels():
        return ["Sample","T Empty Beam","Transmission Sample","Thickness", "Blocked Beam", 
                "Suffix","Description","Start Time", "End Time","Gravity Fix","Wide Angle Fix",
                "Sample Mask", "Transmission Mask","Radius Cut","Wave Cut","User File","Output"]

    @staticmethod
    def column_keys():
        return ["sample", "t_empty_beam", "t_sample", "thickness", "blocked_beam", 
                "suffix", "description", "start_time", "end_time", "gravity_correction","wide_angle_correction", 
                "sample_mask", "transmission_mask", "radius_cut", "wave_cut", "user_file", "output_name"]

    @staticmethod
    def column_options():
        return {'corrections': [9,10],
                'radius_wave': [13,14],
                'masks': [11,12]}
    
    @staticmethod
    def create_empty_row():
        row = [''] * 17
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


