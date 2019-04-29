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

from sans.user_file.settings_tags import (OtherId, DetectorId, LimitsId, SetId, SampleId, MonId, TransId, GravityId,
                                          QResolutionId, FitId, MaskId, event_binning_string_values, set_scales_entry,
                                          monitor_spectrum, simple_range, monitor_file, det_fit_range,
                                          q_rebin_values, fit_general, mask_angle_entry, range_entry, position_entry)
from sans.common.enums import (ReductionDimensionality, ISISReductionMode, RangeStepType, SaveType,
                               DetectorType, DataType, FitType, SANSInstrument)

from sans.ansto.state_model import AnstoStateGuiModel

def _fname():
    # gets the name of the function in which called
    return sys._getframe(1).f_code.co_name

class BilbyStateGuiModel(AnstoStateGuiModel):
    def __init__(self, user_file_items):
        super(BilbyStateGuiModel, self).__init__(user_file_items)

    # ==================================================================================================================
    # Settings | Reduction Parameters
    # ==================================================================================================================

    @property
    def minimum_wavelength(self):
        return self.get_simple_element_with_attribute("wavelength_bins", 
                                                      default_value="", 
                                                      attribute='start')

    @minimum_wavelength.setter
    def minimum_wavelength(self, value):
        self.update_simple_range("wavelength_bins", start=value)

    @property
    def maximum_wavelength(self):
        return self.get_simple_element_with_attribute("wavelength_bins", 
                                                      default_value="", 
                                                      attribute='stop')

    @maximum_wavelength.setter
    def maximum_wavelength(self, value):
        self.update_simple_range("wavelength_bins", stop=value)

    @property
    def wavelength_step(self):
        return self.get_simple_element_with_attribute("wavelength_bins", 
                                                      default_value="", 
                                                      attribute='step')

    @wavelength_step.setter
    def wavelength_step(self, value):
        self.update_simple_range("wavelength_bins", step=value)

    #-------------------------------------------------------

    @property
    def minimum_q1d(self):
        return self.get_simple_element_with_attribute("q1d_bins", 
                                                      default_value="", 
                                                      attribute='start')

    @minimum_q1d.setter
    def minimum_q1d(self, value):
        self.update_simple_range("q1d_bins", step=value)

    @property
    def maximum_q1d(self):
        return self.get_simple_element_with_attribute("q1d_bins", 
                                                      default_value="", 
                                                      attribute='stop')

    @maximum_q1d.setter
    def maximum_q1d(self, value):
        self.update_simple_range("q1d_bins", stop=value)

    @property
    def q1d_step(self):
        return self.get_simple_element_with_attribute("q1d_bins", 
                                                      default_value="", 
                                                      attribute='step')

    @q1d_step.setter
    def q1d_step(self, value):
        self.update_simple_range("q1d_bins", step=value)
    
    #-------------------------------------------------------

    @property
    def qxy_interval(self):
        return self.get_simple_element("qxy_interval", default_value="")

    @qxy_interval.setter
    def qxy_interval(self, value):
        if not value:
            return
        self._user_file_items["qxy_interval"] = [value]

    @property
    def qxy_points(self):
        return self.get_simple_element("qxy_points", default_value="")

    @qxy_points.setter
    def qxy_points(self, value):
        if not value:
            return
        self._user_file_items["qxy_points"] = [value]

    # ==================================================================================================================
    # Settings | Transmission Parameters
    # ==================================================================================================================

    @property
    def minimum_transmission_wavelength(self):
        return self.get_simple_element_with_attribute("transmission_bins", 
                                                      default_value="", 
                                                      attribute='start')

    @minimum_transmission_wavelength.setter
    def minimum_transmission_wavelength(self, value):
        self.update_simple_range("transmission_bins", start=value)

    @property
    def maximum_transmission_wavelength(self):
        return self.get_simple_element_with_attribute("transmission_bins", 
                                                      default_value="", 
                                                      attribute='stop')

    @maximum_transmission_wavelength.setter
    def maximum_transmission_wavelength(self, value):
        self.update_simple_range("transmission_bins", stop=value)

    @property
    def transmission_wavelength_step(self):
        return self.get_simple_element_with_attribute("transmission_bins", 
                                                      default_value="", 
                                                      attribute='step')

    @transmission_wavelength_step.setter
    def transmission_wavelength_step(self, value):
        self.update_simple_range("transmission_bins", step=value)

    #-------------------------------------------------------

    @property
    def transmission_fit(self):
        return self.get_simple_element("transmission_fit", default_value=FitType.Polynomial)

    @transmission_fit.setter
    def transmission_fit(self, value):
        if not value:
            return
        self._user_file_items["transmission_fit"] = [value]

    @property
    def polynomial_fit_order(self):
        return self.get_simple_element("polynomial_order", default_value=2)

    @polynomial_fit_order.setter
    def polynomial_fit_order(self, value):
        if not value:
            return
        self._user_file_items["polynomial_order"] = [value]

    #-------------------------------------------------------

    @property
    def plot_transmission(self):
        return self.get_bool_element("plot_transmission", False)

    @plot_transmission.setter
    def plot_transmission(self, value):
        self.set_bool_element("plot_transmission", value)

    @property
    def save_transmission(self):
        return self.get_bool_element("save_transmission", False)

    @save_transmission.setter
    def save_transmission(self, value):
        self.set_bool_element("save_transmission", value)

    # ==================================================================================================================
    # Settings | Advanced Parameters
    # ==================================================================================================================


    @property
    def wide_angle_correction(self):
        return self.get_bool_element("wide_angle_correction", False)

    @wide_angle_correction.setter
    def wide_angle_correction(self, value):
        self.set_bool_element("wide_angle_correction", value)

    @property
    def gravity_correction(self):
        return self.get_bool_element("gravity_correction", default_value=False)
        #return self.get_simple_element("gravity_correction", default_value=False)

    @gravity_correction.setter
    def gravity_correction(self, value):
        self.set_bool_element("gravity_correction", value=value)
        #self.set_simple_element("gravity_correction", value=value)

    @property
    def blocked_beam_correction(self):
        return self.get_bool_element("blocked_beam_correction", False)

    @blocked_beam_correction.setter
    def blocked_beam_correction(self, value):
        self.set_bool_element("blocked_beam_correction", value)

    #------------------------------------------------------------------

    @property
    def sample_mask_file(self):
        return self.get_simple_element("sample_mask_file", default_value="")

    @sample_mask_file.setter
    def sample_mask_file(self, value):
        if not value:
            return
        self._user_file_items["sample_mask_file"] = [value]

    @property
    def transmission_mask_file(self):
        return self.get_simple_element("transmission_mask_file", default_value="")

    @transmission_mask_file.setter
    def transmission_mask_file(self, value):
        if not value:
            return
        self._user_file_items["transmission_mask_file"] = [value]

    #------------------------------------------------------------------

    @property
    def radius_cut(self):
        return self.get_simple_element("radius_cut", default_value="")

    @radius_cut.setter
    def radius_cut(self, value):
        if not value:
            return
        self._user_file_items["radius_cut"] = [value]

    @property
    def wave_cut(self):
        return self.get_simple_element("wave_cut", default_value="")

    @wave_cut.setter
    def wave_cut(self, value):
        if not value:
            return
        self._user_file_items["wave_cut"] = [value]

    # ------------------------------------------------------------------------------------------------------------------
    # Reduction dimensionality
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def reduction_dimensionality(self):
        return self.get_simple_element_with_attribute(element_id=OtherId.reduction_dimensionality,
                                                      default_value=ReductionDimensionality.OneDim)

    @reduction_dimensionality.setter
    def reduction_dimensionality(self, value):
        if value is ReductionDimensionality.OneDim or value is ReductionDimensionality.TwoDim:
            if OtherId.reduction_dimensionality in self._user_file_items:
                del self._user_file_items[OtherId.reduction_dimensionality]
            new_state_entries = {OtherId.reduction_dimensionality: [value]}
            self._user_file_items.update(new_state_entries)
        else:
            raise ValueError("A reduction dimensionality was expected, got instead {}".format(value))


    # ------------------------------------------------------------------------------------------------------------------
    # Fit
    # ------------------------------------------------------------------------------------------------------------------
    def _get_transmission_fit(self, data_type, attribute, default_value):
        if FitId.general in self._user_file_items:
            settings = self._user_file_items[FitId.general]
            # Check first if there are data type specific settings, else check if there are general settings
            extracted_settings = [setting for setting in settings if setting.data_type is data_type]
            if not extracted_settings:
                extracted_settings = [setting for setting in settings if setting.data_type is None]
            if extracted_settings:
                setting = extracted_settings[-1]
                return getattr(setting, attribute)
        return default_value

    def _set_transmission_fit(self, data_type, start=None, stop=None, fit_type=None, polynomial_order=None):
        if FitId.general in self._user_file_items:
            # Gather all settings which correspond to the data type and where the data type is none
            settings = self._user_file_items[FitId.general]
            settings_general = [setting for setting in settings if setting.data_type is None]
            settings_for_data_type = [setting for setting in settings if setting.data_type is data_type]
            # We check if there are data-type specific settings.
            # 1. There are data type specific settings. Then we are good.
            # 2. There are no data type specific settings. We create one data type specific setting and populate it
            #    with a general setting if it exists else we create a new entry
            if not settings_for_data_type:
                if settings_general:
                    setting_general = settings_general[-1]
                    settings.append(fit_general(start=setting_general.start, stop=setting_general.stop,
                                                data_type=data_type, fit_type=setting_general.fit_type,
                                                polynomial_order=setting_general.polynomial_order))
                else:
                    settings.append(fit_general(start=None, stop=None, data_type=data_type,
                                                fit_type=FitType.NoFit, polynomial_order=2))
        else:
            settings = [fit_general(start=None, stop=None, data_type=data_type,
                                    fit_type=FitType.NoFit, polynomial_order=2)]

        new_settings = []
        for setting in settings:
            # We only want to modify the settings which are either the data type specific ones or the ones which
            # don't have a specific data type
            if setting.data_type is data_type and setting.data_type is not None:
                new_start = start if start is not None else setting.start
                new_stop = stop if stop is not None else setting.stop
                new_fit_type = fit_type if fit_type is not None else setting.fit_type
                new_polynomial_order = polynomial_order if polynomial_order is not None else setting.polynomial_order
                new_settings.append(fit_general(start=new_start, stop=new_stop, fit_type=new_fit_type,
                                                data_type=setting.data_type, polynomial_order=new_polynomial_order))
            else:
                new_settings.append(setting)
        self._user_file_items.update({FitId.general: new_settings})

    def has_transmission_fit_got_separate_settings_for_sample_and_can(self):
        if FitId.general in self._user_file_items:
            settings = self._user_file_items[FitId.general]
            if settings:
                settings_sample = [setting for setting in settings if setting.data_type is DataType.Sample]
                settings_can = [setting for setting in settings if setting.data_type is DataType.Can]
                # If we have either one or the other
                if settings_sample or settings_can:
                    return True
        return False


    @property
    def show_transmission(self):
        return self.get_simple_element(element_id=OtherId.show_transmission, default_value=True)

    @show_transmission.setter
    def show_transmission(self, value):
        self.set_simple_element(element_id=OtherId.show_transmission, value=value)

    # ------------------------------------------------------------------------------------------------------------------
    # Wavelength- and pixel-adjustment files
    # ------------------------------------------------------------------------------------------------------------------
    def _get_adjustment_file_setting(self, element_id, detector_type, default_value):
        if element_id in self._user_file_items:
            settings = self._user_file_items[element_id]

            # Separate out the correct detector type
            settings = [setting for setting in settings if setting.detector_type is detector_type]
            if settings:
                setting = settings[-1]
                return setting.file_path
        return default_value

    def _set_adjustment_file_setting(self, element_id, detector_type, file_path):
        # Check if we already have items for the particular detector type
        settings_with_correct_detector = []
        settings = []
        if element_id in self._user_file_items:
            settings = self._user_file_items[element_id]
            settings_with_correct_detector = [setting for setting in settings if setting.detector_type is detector_type]

        if not (settings and settings_with_correct_detector):
            settings.append(monitor_file(file_path="", detector_type=detector_type))

        # At this point we have settings with the desired detector type
        new_settings = []
        for setting in settings:
            if setting.detector_type is detector_type:
                new_settings.append(monitor_file(file_path=file_path, detector_type=setting.detector_type))
            else:
                new_settings.append(setting)
        self._user_file_items.update({element_id: new_settings})


    # ------------------------------------------------------------------------------------------------------------------
    # Q Limits
    # ------------------------------------------------------------------------------------------------------------------
    def _set_q_1d_limits(self, min_value=None, max_value=None, rebin_string=None):
        element_id = LimitsId.q
        if element_id in self._user_file_items:
            settings = self._user_file_items[element_id]
        else:
            settings = [q_rebin_values(min=None, max=None, rebin_string=None)]

        # At this point we have settings with the desired detector type
        new_settings = []
        for setting in settings:
            new_min = min_value if min_value is not None else setting.min
            new_max = max_value if max_value is not None else setting.max
            new_rebin_string = rebin_string if rebin_string is not None else setting.rebin_string
            new_settings.append(q_rebin_values(min=new_min, max=new_max, rebin_string=new_rebin_string))
        self._user_file_items.update({element_id: new_settings})

    def _set_q_xy_limits(self, stop_value=None, step_value=None, step_type_value=None):
        element_id = LimitsId.qxy
        if element_id in self._user_file_items:
            settings = self._user_file_items[element_id]
        else:
            settings = [simple_range(start=None, stop=None, step=None, step_type=None)]

        # At this point we have settings with the desired detector type
        new_settings = []
        for setting in settings:
            new_stop = stop_value if stop_value is not None else setting.stop
            new_step = step_value if step_value is not None else setting.step
            new_step_type_value = step_type_value if step_type_value is not None else setting.step_type
            new_settings.append(simple_range(start=None, stop=new_stop, step=new_step, step_type=new_step_type_value))
        self._user_file_items.update({element_id: new_settings})

    @property
    def q_1d_rebin_string(self):
        return self.get_simple_element_with_attribute(element_id=LimitsId.q, default_value="",
                                                      attribute="rebin_string")

    @q_1d_rebin_string.setter
    def q_1d_rebin_string(self, value):
        self._set_q_1d_limits(rebin_string=value)

    @property
    def q_xy_max(self):
        return self.get_simple_element_with_attribute(element_id=LimitsId.qxy, default_value="",
                                                      attribute="stop")

    @q_xy_max.setter
    def q_xy_max(self, value):
        self._set_q_xy_limits(stop_value=value)

    @property
    def q_xy_step(self):
        return self.get_simple_element_with_attribute(element_id=LimitsId.qxy, default_value="",
                                                      attribute="step")

    @q_xy_step.setter
    def q_xy_step(self, value):
        self._set_q_xy_limits(step_value=value)

    @property
    def q_xy_step_type(self):
        return self.get_simple_element_with_attribute(element_id=LimitsId.qxy, default_value=None,
                                                      attribute="step_type")

    @q_xy_step_type.setter
    def q_xy_step_type(self, value):
        self._set_q_xy_limits(step_type_value=value)

    @property
    def r_cut(self):
        return self.get_simple_element(element_id=LimitsId.radius_cut, default_value="")

    @r_cut.setter
    def r_cut(self, value):
        self.set_simple_element(element_id=LimitsId.radius_cut, value=value)

    @property
    def w_cut(self):
        return self.get_simple_element(element_id=LimitsId.wavelength_cut, default_value="")

    @w_cut.setter
    def w_cut(self, value):
        self.set_simple_element(element_id=LimitsId.wavelength_cut, value=value)

    # ------------------------------------------------------------------------------------------------------------------
    # Gravity
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def gravity_on_off(self):
        return self.get_simple_element(element_id=GravityId.on_off, default_value=True)

    @gravity_on_off.setter
    def gravity_on_off(self, value):
        self.set_simple_element(element_id=GravityId.on_off, value=value)

    @property
    def gravity_extra_length(self):
        return self.get_simple_element(element_id=GravityId.extra_length, default_value="")

    @gravity_extra_length.setter
    def gravity_extra_length(self, value):
        self.set_simple_element(element_id=GravityId.extra_length, value=value)

    # ------------------------------------------------------------------------------------------------------------------
    # QResolution
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def use_q_resolution(self):
        return self.get_simple_element(element_id=QResolutionId.on, default_value=False)

    @use_q_resolution.setter
    def use_q_resolution(self, value):
        self.set_simple_element(element_id=QResolutionId.on, value=value)

    @property
    def q_resolution_source_a(self):
        return self.get_simple_element(element_id=QResolutionId.a1, default_value="")

    @q_resolution_source_a.setter
    def q_resolution_source_a(self, value):
        self.set_simple_element(element_id=QResolutionId.a1, value=value)

    @property
    def q_resolution_sample_a(self):
        return self.get_simple_element(element_id=QResolutionId.a2, default_value="")

    @q_resolution_sample_a.setter
    def q_resolution_sample_a(self, value):
        self.set_simple_element(element_id=QResolutionId.a2, value=value)

    @property
    def q_resolution_source_h(self):
        return self.get_simple_element(element_id=QResolutionId.h1, default_value="")

    @q_resolution_source_h.setter
    def q_resolution_source_h(self, value):
        self.set_simple_element(element_id=QResolutionId.h1, value=value)

    @property
    def q_resolution_sample_h(self):
        return self.get_simple_element(element_id=QResolutionId.h2, default_value="")

    @q_resolution_sample_h.setter
    def q_resolution_sample_h(self, value):
        self.set_simple_element(element_id=QResolutionId.h2, value=value)

    @property
    def q_resolution_source_w(self):
        return self.get_simple_element(element_id=QResolutionId.w1, default_value="")

    @q_resolution_source_w.setter
    def q_resolution_source_w(self, value):
        self.set_simple_element(element_id=QResolutionId.w1, value=value)

    @property
    def q_resolution_sample_w(self):
        return self.get_simple_element(element_id=QResolutionId.w2, default_value="")

    @q_resolution_sample_w.setter
    def q_resolution_sample_w(self, value):
        self.set_simple_element(element_id=QResolutionId.w2, value=value)

    @property
    def q_resolution_delta_r(self):
        return self.get_simple_element(element_id=QResolutionId.delta_r, default_value="")

    @q_resolution_delta_r.setter
    def q_resolution_delta_r(self, value):
        self.set_simple_element(element_id=QResolutionId.delta_r, value=value)

    @property
    def q_resolution_moderator_file(self):
        return self.get_simple_element(element_id=QResolutionId.moderator, default_value="")

    @q_resolution_moderator_file.setter
    def q_resolution_moderator_file(self, value):
        self.set_simple_element(element_id=QResolutionId.moderator, value=value)

    @property
    def q_resolution_collimation_length(self):
        return self.get_simple_element(element_id=QResolutionId.collimation_length, default_value="")

    @q_resolution_collimation_length.setter
    def q_resolution_collimation_length(self, value):
        self.set_simple_element(element_id=QResolutionId.collimation_length, value=value)

    # ==================================================================================================================
    # ==================================================================================================================
    # MASK TAB
    # ==================================================================================================================
    # ==================================================================================================================

    # ------------------------------------------------------------------------------------------------------------------
    # Phi limit
    # ------------------------------------------------------------------------------------------------------------------
    def _set_phi_limit(self, min_value=None, max_value=None, use_mirror=None):
        if LimitsId.angle in self._user_file_items:
            settings = self._user_file_items[LimitsId.angle]
        else:
            settings = [mask_angle_entry(min=None, max=None, use_mirror=False)]

        new_settings = []
        for setting in settings:
            new_min = min_value if min_value is not None else setting.min
            new_max = max_value if max_value is not None else setting.max
            new_use_mirror = use_mirror if use_mirror is not None else setting.use_mirror
            new_settings.append(mask_angle_entry(min=new_min, max=new_max, use_mirror=new_use_mirror))
        self._user_file_items.update({LimitsId.angle: new_settings})

    @property
    def phi_limit_min(self):
        return self.get_simple_element_with_attribute(element_id=LimitsId.angle, attribute="min", default_value="-90")

    @phi_limit_min.setter
    def phi_limit_min(self, value):
        self._set_phi_limit(min_value=value)

    @property
    def phi_limit_max(self):
        return self.get_simple_element_with_attribute(element_id=LimitsId.angle, attribute="max", default_value="90")

    @phi_limit_max.setter
    def phi_limit_max(self, value):
        self._set_phi_limit(max_value=value)

    @property
    def phi_limit_use_mirror(self):
        return self.get_simple_element_with_attribute(element_id=LimitsId.angle, attribute="use_mirror", default_value=True)  # noqa

    @phi_limit_use_mirror.setter
    def phi_limit_use_mirror(self, value):
        self._set_phi_limit(use_mirror=value)

    # ------------------------------------------------------------------------------------------------------------------
    # Radius limit
    # ------------------------------------------------------------------------------------------------------------------
    def _set_radius_limit(self, min_value=None, max_value=None):
        if LimitsId.radius in self._user_file_items:
            settings = self._user_file_items[LimitsId.radius]
        else:
            settings = [range_entry(start=None, stop=None)]

        new_settings = []
        for setting in settings:
            new_min = min_value if min_value is not None else setting.start
            new_max = max_value if max_value is not None else setting.stop
            new_settings.append(range_entry(start=new_min, stop=new_max))
        self._user_file_items.update({LimitsId.radius: new_settings})

    @property
    def radius_limit_min(self):
        return self.get_simple_element_with_attribute(element_id=LimitsId.radius, attribute="start", default_value="")

    @radius_limit_min.setter
    def radius_limit_min(self, value):
        self._set_radius_limit(min_value=value)

    @property
    def radius_limit_max(self):
        return self.get_simple_element_with_attribute(element_id=LimitsId.radius, attribute="stop", default_value="")

    @radius_limit_max.setter
    def radius_limit_max(self, value):
        self._set_radius_limit(max_value=value)

    # ------------------------------------------------------------------------------------------------------------------
    # Mask files
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def mask_files(self):
        if MaskId.file in self._user_file_items:
            return self._user_file_items[MaskId.file]
        return []

    @mask_files.setter
    def mask_files(self, value):
        if value is None:
            return
        if MaskId.file in self._user_file_items:
            del self._user_file_items[MaskId.file]
        new_state_entries = {MaskId.file: value}
        self._user_file_items.update(new_state_entries)

    # ------------------------------------------------------------------------------------------------------------------
    # Output name
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def output_name(self):
        return self.get_simple_element(element_id=OtherId.user_specified_output_name, default_value="")

    @output_name.setter
    def output_name(self, value):
        self.set_simple_element(element_id=OtherId.user_specified_output_name, value=value)
