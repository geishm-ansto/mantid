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
                    parsed_rows[index] = dict(zip(header, row[1:]))
                except ValueError:
                    raise RuntimeError("Invalid index value in Batch file.")
        csv_file.close()

        # return the index ordered list of batch parameter values as dictionaries
        skeys = sorted(parsed_rows.keys())
        return [parsed_rows[k] for k in skeys]

    @staticmethod
    def map_header(tags):
        return tags
