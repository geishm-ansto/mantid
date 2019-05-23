# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
""" 
Simple file readers for loading user and batch file data.
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
    return int(value) if value else None

def parse_float(value):
    return float(value) if value else None

def parse_scale_float(scale, value):
    return scale * float(value) if value else None

boolean_true = ('t', 'true', '1', 'y', 'yes')
boolean_false = ('f', 'false', '0', 'n', 'no')
def parse_boolean(value):
    if value:
        svalue = value.lower().strip()
        if svalue in boolean_true or svalue in boolean_false:
            return svalue in boolean_true
        else:
            raise RuntimeError('Invalid boolean argument: {}'.format(value))
    else:
        return False

class UserFileReader(object):
    """
    The reader is responsible for populating the BilbyStateGuiModel from the file.

    Assumes the following format for a csv file.
    index, param1, param2, ...
         , comments follow an empty first cell
        0, 12.3, ...

    The loader creates an index ordered list of dictionaries using the header information as 
    the keys.
    If the loader follows the csv format then it is only necessary to update the function
    and header maps.
    """
    function_maps = {}
    header_maps = {}

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
                        parsed_values[tag] = self.parse_value(tag, value)

                except ValueError:
                    raise RuntimeError("Invalid value in user file.")
        csv_file.close()

        return parsed_values

    def parse_value(self, tag, value):

        try:
            return self.function_maps[tag](value)
        except KeyError:
            return value

    def map_header(self, tags):
        # the headings are case insenstive
        ntags = []
        for tag in tags:
            try:
                ntag = tag.lower()
                ntags.append(self.header_maps[ntag])
            except KeyError:
                ntags.append(ntag)
        return ntags

def parse_run_number(tag, value):
    # any value of the form 'BBY00[nnnn]' is mapped to 'nnnn'
    pattern = '^{}0*(.+)$'.format(tag)
    s = re.search(pattern, value.lower())
    if s:
        return s.group(1)
    else:
        return value

bby_run_number = functools.partial(parse_run_number, 'bby')

class BatchFileReader(object):
    """
    Assumes the following format for a csv file.
    index, param1, param2, ...
         , comments follow an empty first cell
        0, 12.3, ...

    The loader creates an index ordered list of dictionaries using the header information as 
    the keys. 
    If the loader follows the csv format then it is only necessary to update the function
    and header maps.
    """
    function_maps = {}
    header_maps = {}

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
                    header = self.map_header(row)
                    continue
                # assumes comment if the first entry is empty 
                if not row[0] or row[0].lower().strip() == 'end':
                    continue
                try:
                    parsed_values = {}
                    for tag, value in zip(header, row):
                        parsed_values[tag] = self.parse_value(tag, value)
                    index = int(row[0])
                    parsed_rows[index] = parsed_values
                except ValueError:
                    raise RuntimeError("Invalid index value in Batch file.")
        csv_file.close()

        # return the index ordered list of batch parameter values as dictionaries
        skeys = sorted(parsed_rows.keys())
        return [parsed_rows[k] for k in skeys]

    def parse_value(self, tag, value):

        try:
            return self.function_maps[tag](value)
        except KeyError:
            return value

    def map_header(self, tags):
        # the headings are case insenstive
        ntags = []
        for tag in tags:
            try:
                ntag = tag.lower()
                ntags.append(self.header_maps[ntag])
            except KeyError:
                ntags.append(ntag)
        return ntags
