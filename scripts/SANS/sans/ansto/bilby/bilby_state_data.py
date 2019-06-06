# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
""" The state gui model contains all the reduction parameters for Bilby processing.

This is one of the two models which is used for the data reduction. It contains generally all the settings which
are not available in the model associated with the data table.
"""

from __future__ import (absolute_import, division, print_function)

import json

import sans.common.constants
from sans.state.state_base import (StateBase, StringParameter, PositiveIntegerParameter, BoolParameter,
                                   ClassTypeParameter, PositiveFloatParameter, rename_descriptor_names, 
                                   TypedParameter, validator_sub_state)
from sans.state.state_functions import (is_pure_none_or_not_none, one_is_none, 
                                        is_not_none_and_first_larger_than_second, validation_message)
from sans.common.enums import (RangeStepType, SANSInstrument, SANSFacility)
from sans.state.wavelength import StateWavelength

@rename_descriptor_names
class TimeSliceState(StateBase):
    start_time = PositiveFloatParameter()
    end_time = PositiveFloatParameter()

    def __init__(self):
        super(TimeSliceState, self).__init__()

    def validate(self):
        is_invalid = dict()
        attr_list = [self.start_time, self.end_time]
        if not is_pure_none_or_not_none(attr_list):
            entry = validation_message("An entry has not been set.",
                                       "Make sure all entries are set or cleared.",
                                       {"start_time": self.start_time,
                                        "end_time": self.end_time})
            is_invalid.update(entry)
        elif is_not_none_and_first_larger_than_second(attr_list):
            entry = validation_message("Incorrect bounds.",
                                       "Make sure that lower bound is smaller then upper bound.",
                                       {"start_time": self.start_time,
                                        "end_time": self.end_time})
            is_invalid.update(entry)

        if is_invalid:
            raise ValueError("TimeSliceState: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid, indent=4)))

@rename_descriptor_names
class PositiveRangeState(StateBase):
    low = PositiveFloatParameter()
    high = PositiveFloatParameter()
    step = PositiveFloatParameter()
    step_type = ClassTypeParameter(RangeStepType)

    def __init__(self):
        super(PositiveRangeState, self).__init__()
        self.step_type = RangeStepType.Lin

    def validate(self):
        is_invalid = dict()
        if one_is_none([self.low, self.high, self.step]):
            entry = validation_message("An entry has not been set.",
                                       "Make sure that all entries are set.",
                                       {"low": self.low,
                                        "high": self.high,
                                        "step": self.step})
            is_invalid.update(entry)

        if is_not_none_and_first_larger_than_second([self.low, self.high]):
            entry = validation_message("Incorrect bounds.",
                                       "Make sure that lower bound is smaller then upper bound.",
                                       {"low": self.low,
                                        "high": self.high})
            is_invalid.update(entry)

        if is_invalid:
            raise ValueError("PositiveRangeState: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid, indent=4)))

@rename_descriptor_names
class BilbyStateData(StateBase):
    """
    The state data model is a complete description of the processing to be performed by
    the BilbySANSDataProcessor.
    """
    sample = StringParameter()
    sample_mask = StringParameter()
    blocked_beam = StringParameter()

    transmission_sample = StringParameter()
    transmission_mask = StringParameter()
    transmission_empty = StringParameter()

    empty_beam_spectrum = StringParameter()

    att_pos = PositiveIntegerParameter()
    sample_thickness = PositiveFloatParameter()
    fit_method = StringParameter()
    polynomial_order = PositiveIntegerParameter()
    time_slice = TypedParameter(TimeSliceState, validator_sub_state)

    wavelength = TypedParameter(StateWavelength, validator_sub_state)
    transmission_wavelength = TypedParameter(PositiveRangeState, validator_sub_state)
    backward_q_bin = BoolParameter()
    binning_q = TypedParameter(PositiveRangeState, validator_sub_state)

    external_mode = BoolParameter()
    gravity_correction = BoolParameter()
    solid_angle_weighting = BoolParameter()
    wide_angle_correction = BoolParameter()
    reduce_2D = BoolParameter()
    data_points_2D = PositiveIntegerParameter()

    radius_cut = PositiveFloatParameter()
    wave_cut = PositiveFloatParameter()

    pre_2016_data = BoolParameter()
    pre_may_2016_data = BoolParameter()

    output_folder = StringParameter()
    suffix = StringParameter()
    
    instrument = ClassTypeParameter(SANSInstrument)
    facility = ClassTypeParameter(SANSFacility)
    idf_file_path = StringParameter()
    ipf_file_path = StringParameter()

    def __init__(self):
        super(BilbyStateData, self).__init__()

        self.instrument = SANSInstrument.BILBY
        self.facility = SANSFacility.ANSTO

    def _validate_binning(self, errors):

        # transmission range must be equal or longer than the wavelength binning range for data reduction
        if not self.wavelength or not self.transmission_wavelength:
            entry = validation_message("Wavelengths are not specified.",
                                       {"wavelength": self.wavelength,
                                        "transmission_wavelength": self.transmission_wavelength})
            errors.update(entry)
        elif (self.wavelength.wavelength_low[0] < self.transmission_wavelength.low or
            self.wavelength.wavelength_high[-1] > self.transmission_wavelength.high):
            entry = validation_message("Transmission binning range is less than the sample binning range.",
                                        {"wavelength": self.wavelength,
                                        "transmission_wavelength": self.transmission_wavelength})
            errors.update(entry)

        if not self.external_mode and len(self.wavelength.wavelength_low) > 1:
            entry = validation_message("Attempting multiple wavelength reduction for monochromatic data ",
                                       {})

    def _validate_2D_settings(self, errors):
        if self.reduce_2D:
            if not self.data_points_2D:
                entry = validation_message("Number of data points need to be set for 2D reduction",
                                           {"data_points": self.data_points_2D})
                errors.update(entry)

    def _validate_attenuation(self, errors):
        if self.att_pos < 1 or self.att_pos > 5:
            entry = validation_message("Attenation position outside expected range of 1 .. 5.",
                                       {"att_pos": self.att_pos})
            errors.update(entry)
        elif self.pre_may_2016_data and self.att_pos in [2,4]:
            entry = validation_message("Attenation position 2 or 4 for pre May 2016 data is unexpected.",
                                       {"att_pos": self.att_pos})
            errors.update(entry)

    def validate(self):
        is_invalid = dict()

        if not self.sample:
            entry = validation_message("Sample file was not specified.",
                                       "Make sure that the sample scatter file is specified.",
                                       {"sample": self.sample})
            is_invalid.update(entry)

        if not self.empty_beam_spectrum:
            entry = validation_message("Empty beam spectrum file was not specified.",
                                       "Make sure that the Empty beam spectrum file is specified.",
                                       {"empty_beam_spectrum": self.empty_beam_spectrum})
            is_invalid.update(entry)
    
        if not self.transmission_sample:
            entry = validation_message("Transmission file was not specified.",
                                       "Make sure that the transmission sample file is specified.",
                                       {"transmission_sample": self.transmission_sample})
            is_invalid.update(entry)

        if not self.transmission_mask:
            entry = validation_message("Transmission mask file was not specified.",
                                       "Make sure that the transmission mask file is specified.",
                                       {"transmission_mask": self.transmission_mask})
            is_invalid.update(entry)

        if not self.transmission_empty:
            entry = validation_message("Transmission empty file was not specified.",
                                       "Make sure that the transmission empty file is specified.",
                                       {"transmission_empty": self.transmission_empty})
            is_invalid.update(entry)

        # now more specific checks
        self._validate_binning(is_invalid)
        self._validate_2D_settings(is_invalid)
        self._validate_attenuation(is_invalid)

        if is_invalid:
            raise ValueError("BilbyStateData: The provided inputs are illegal. "
                             "Please see: {0}".format(json.dumps(is_invalid)))

