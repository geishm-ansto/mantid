# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
""" 
Simple file readers for loading Batch and batch file data.
"""

from __future__ import (absolute_import, division, print_function)

import functools
import os
import re
import sys
import csv

from sans.common.file_information import find_full_file_path
from sans.common.enums import FitType
from sans.user_file.settings_tags import simple_range

numeric_pattern = re.compile('[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?', re.VERBOSE)

def parse_simple_range(value):
    # assumes start step stop as string
    ss = numeric_pattern.findall(value)
    try:
        return simple_range(start=ss[0], stop=ss[2], step=ss[1],step_type='linear')
    except (IndexError, ValueError):
        raise RuntimeError('Invalid range format: {}'.format(value))

valid_fit_types = {'log': FitType.Logarithmic,
                   'logarithmic': FitType.Logarithmic,
                   'polynomial': FitType.Polynomial,   
                   'poly': FitType.Polynomial, 
                   'linear': FitType.Linear,
                   'lin': FitType.Linear,
                   }
def parse_fit_type(value):
    try: 
        return valid_fit_types[value.lower().strip()]
    except KeyError:
        return FitType.Polynomial

def parse_integer(value):
    return int(value)

class UserFileReader(object):
    """
    The reader is responsible for populating the BilbyStateGuiModel from the file.

    Assumes the following format for a csv file.
    index, param1, param2, ...
         , comments follow an empty first cell
        0, 12.3, ...

    The loader creates an index ordered list of dictionaries using the header information as 
    the keys. 
    """
    def __init__(self, file_path):
        super(UserFileReader, self).__init__()
        self._file_path = find_full_file_path(file_path)
        self._header = None

    def read_user_file(self):
        
        # Python 2 and 3 take input in different modes for writing lists to csv files
        if sys.version_info[0] == 2:
            open_type = 'rb'
        else:
            open_type = 'r'

        parsed_values = {}
        csv.register_dialect('dialect', delimiter = ',', skipinitialspace=True)

        with open(self._file_path, open_type) as csv_file:
            header = None
            reader = csv.reader(csv_file, dialect='dialect')
            for _row in reader:
                row = list(_row)
                if header is None:
                    if row[0].lower().strip() != 'index':
                        raise RuntimeError("User file did not contain an index column as the first entry in the header.")
                    header = self.map_header(row[1:])
                    continue
                # assumes comment if the first entry is empty 
                if not row[0] or row[0].lower().strip() == 'end':
                    continue
                try:
                    #index = int(row[0])
                    for tag, value in zip(header, row[1:]):
                        parsed_values[tag] = [self.parse_value(tag, value)]

                except ValueError:
                    raise RuntimeError("Invalid value in user file.")
        csv_file.close()

        return parsed_values

    @staticmethod
    def parse_value(tag, value):
        function_maps = {'wavelength_bins': parse_simple_range,
                         'q1d_bins': parse_simple_range,
                         'transmission_bins': parse_simple_range,
                         'transmission_fit': parse_fit_type,
                         'polynomial_order': parse_integer,
                         }
        try:
            return function_maps[tag](value)
        except KeyError:
            return value


    @staticmethod
    def map_header(tags):
        # maps similar labels to the equivalent Bilby row model column keys
        # the headings are case insenstive
        bilby_map = {'binning_wavelength_ini': 'wavelength_bins',
                     'binning_q': 'q1d_bins',
                     'binning_wavelength_transmission': 'transmission_bins',
                     'radiuscut': 'radius_cut',
                     'wavecut': 'wave_cut',
                     'polynomialorder': 'polynomial_order', 
                     '2d_number_data_points': 'qxy_points',
                    }
        ntags = []
        for tag in tags:
            try:
                ntag = tag.lower()
                ntags.append(bilby_map[ntag])
            except KeyError:
                ntags.append(ntag)
        return ntags

class BatchFileReader(object):
    """
    Assumes the following format for a csv file.
    index, param1, param2, ...
         , comments follow an empty first cell
        0, 12.3, ...

    The loader creates an index ordered list of dictionaries using the header information as 
    the keys. 
    """
    def __init__(self, file_path):
        super(BatchFileReader, self).__init__()
        self._file_path = find_full_file_path(file_path)
        self._header = None

    def parse_batch_file(self):
        
        # Python 2 and 3 take input in different modes for writing lists to csv files
        if sys.version_info[0] == 2:
            open_type = 'rb'
        else:
            open_type = 'r'

        parsed_rows = {}
        csv.register_dialect('dialect', delimiter = ',', skipinitialspace=True)

        with open(self._file_path, open_type) as csv_file:
            header = None
            reader = csv.reader(csv_file, dialect='dialect')
            for _row in reader:
                row = list(_row)
                if header is None:
                    if row[0].lower().strip() != 'index':
                        raise RuntimeError("Batch file did not contain an index column as the first entry in the header.")
                    header = self.map_header(row[1:])
                    continue
                # assumes comment if the first entry is empty 
                if not row[0] or row[0].lower().strip() == 'end':
                    continue
                try:
                    index = int(row[0])
                    parsed_rows[index] = dict(zip(header, self.map_to_run(row[1:])))
                except ValueError:
                    raise RuntimeError("Invalid index value in Batch file.")
        csv_file.close()

        # return the index ordered list of batch parameter values as dictionaries
        skeys = sorted(parsed_rows.keys())
        return [parsed_rows[k] for k in skeys]

    @staticmethod
    def map_to_run(values):
        # any value of the form 'BBY00[nnnn]' is mapped to 'nnnn'
        nvalues = []
        for value in values:
            s = re.search('^bby0*(.+)$', value.lower())
            if s:
                nvalues.append(s.group(1))
            else:
                nvalues.append(value)
        return nvalues

    @staticmethod
    def map_header(tags):
        # maps similar labels to the equivalent Bilby row model column keys
        # the headings are case insenstive
        bilby_map = {'t_emptybeam': 't_empty_beam',
                     'thickness [cm]': 'thickness',
                     'mask_transmission': 'transmission_mask',
                     'blockedbeam': 'blocked_beam',
                     'additional_description': 'description',
                     'mask': 'sample_mask',
                     'starttime': 'start_time',
                     'endtime': 'end_time'}
        ntags = []
        for tag in tags:
            try:
                ntag = tag.lower()
                ntags.append(bilby_map[ntag])
            except KeyError:
                ntags.append(ntag)
        return ntags
