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

from sans.ansto.file_readers import (UserFileReader, BatchFileReader, parse_simple_range, parse_fit_type,
                                     parse_integer, parse_float, parse_scale_float,
                                     parse_boolean, parse_run_number)


class BilbyUserFileReader(UserFileReader):

    function_maps = {'wavelength_bins': parse_simple_range,
                    'q1d_bins': parse_simple_range,
                    'transmission_bins': parse_simple_range,
                    'transmission_fit': parse_fit_type,
                    'polynomial_fit_order': parse_integer,
                    'reduce_2d': parse_boolean,
                    'plot_2d': parse_boolean,
                    'gravity_correction': parse_boolean,
                    'wide_angle_correction': parse_boolean,
                    'blocked_beam_correction': parse_boolean,
                    'radius_cut': parse_float,
                    'wave_cut': parse_float,
                    'qxy_points': parse_integer,
                    'wavelength_interval': parse_float,
                    'slice_wavelength': parse_boolean,
                    'qxy_points': parse_integer,
                    'solid_angle_correction': parse_boolean
                    }

    header_maps = {'binning_wavelength_ini': 'wavelength_bins',
                'binning_q': 'q1d_bins',
                'binning_wavelength_transmission': 'transmission_bins',
                'radiuscut': 'radius_cut',
                'wavecut': 'wave_cut',
                'polynomialorder': 'polynomial_fit_order', 
                '2d_number_data_points': 'qxy_points',
                'reduced_files_folder': 'output_folder',
                'reduce_2D': 'reduce_2d',
                'plot_2D': 'plot_2d',
                'wavelength_intervals': 'slice_wavelength',
                'wav_delta':  'wavelength_interval',
                'sample_mask': 'sample_mask_file',
                'transmission_mask': 'transmission_mask_file',
            }

    def __init__(self, file_path):
        super(BilbyUserFileReader, self).__init__(file_path)

bby_run_number = functools.partial(parse_run_number, 'bby')

class BilbyBatchFileReader(BatchFileReader):

    function_maps = {'sample': bby_run_number,
                    't_empty_beam': bby_run_number,
                    't_sample': bby_run_number,
                    'blocked_beam': bby_run_number,
                    'thickness': functools.partial(parse_scale_float, 1.0),
                    'start_time': parse_float,
                    'end_time': parse_float,
                    }

    header_maps = {'t_emptybeam': 't_empty_beam',
                'thickness [cm]': 'thickness',
                'mask_transmission': 'transmission_mask',
                'blockedbeam': 'blocked_beam',
                'additional_description': 'description',
                'mask': 'sample_mask',
                'starttime': 'start_time',
                'endtime': 'end_time'}

    def __init__(self, file_path):
        super(BilbyBatchFileReader, self).__init__(file_path)
