# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
""" The ANSTO run tab presenter.

The core presenter is only responsible for managing the run table for batch processing.
Itis stripped down version of the ISIS run tab presenter that handles all of the
ISIS SANS instruments.
"""

from __future__ import (absolute_import, division, print_function)

import os
import math

from mantid.kernel import DateAndTime
from sans.common.enums import (FitType, ReductionDimensionality)
from sans.state.wavelength import StateWavelength
from sans.ansto.run_tab_presenter import RunTabPresenter
from sans.ansto.bilby.bilby_state_data import (BilbyStateData, TimeSliceState, 
                                               PositiveRangeState)
from ui.ansto.ansto_bilby_gui import BilbyBatchReductionGui

# ----------------------------------------------------------------
# Convert to ascii
#-----------------------------------------------------------------
FitType2Bilby = {FitType.Linear: 'Linear',
                 FitType.Logarithmic: 'Log',
                 FitType.Polynomial: 'Polynomial',
                 FitType.NoFit: None }


def _convert_fit_method(method):
    try:
        return FitType2Bilby[method]
    except KeyError:
        return None

def _create_range_state(min_w, max_w, step_w):
    state = PositiveRangeState()
    state.low = float(min_w)
    state.high = float(max_w)
    state.step = float(step_w)
    #state.validate()
    return state


def _create_wavelength_slices(min_w, max_w, step_w, wave_delta, wavelength_intervals):
    state = StateWavelength()

    if wavelength_intervals and wave_delta > 0:
        n = int(math.floor((max_w - min_w)/wave_delta))
        state.wavelength_low = [min_w + i * wave_delta for i in range(n)]
        state.wavelength_high = [min_w + (i + 1) * wave_delta for i in range(n)]
        state.wavelength_step = step_w
    else:
        state.wavelength_low = [min_w]
        state.wavelength_high = [max_w]
        state.wavelength_step = step_w

    #state.validate()
    return state

def _create_time_slice_state(start_time, end_time):
    state = TimeSliceState()
    fill = (start_time or end_time)
    if fill:
        state.start_time = float(start_time) if start_time else 0.
        state.end_time = float(end_time) if end_time else 0.
    #state.validate()
    return state

class BilbyPresenter(RunTabPresenter):
    class ConcreteBilbyListener(RunTabPresenter.ConcreteRunTabListener, BilbyBatchReductionGui.BilbyListener):
        def __init__(self, presenter):
            super(BilbyPresenter.ConcreteBilbyListener, self).__init__(presenter)

        def on_show_index_selection(self, show):
            self._presenter.on_show_index_selection(show)

        def on_show_times_selection(self, show):
            self._presenter.on_show_times_selection(show)

        def on_show_corrections_selection(self, show):
            self._presenter.on_show_corrections_selection(show)

        def on_show_masks_selection(self, show):
            self._presenter.on_show_masks_selection(show)

        def on_show_radius_wave_selection(self, show):
            self._presenter.on_show_radius_wave_selection(show)

    def __init__(self, facility, view, models):
        super(BilbyPresenter, self).__init__(facility, view, models)

    def _default_gui_setup(self):
        """
        Provides a default setup of the GUI. This is important for the initial start up, when the view is being set.
        """
        super(BilbyPresenter, self)._default_gui_setup()

        fit_types = [FitType.to_string(FitType.Linear),
                     FitType.to_string(FitType.Logarithmic),
                     FitType.to_string(FitType.Polynomial)]
        self._view.transmission_fit = fit_types

    def on_show_index_selection(self, show):
        self._view.show_column_options(['index'], show)

    def on_show_times_selection(self, show):
        self._view.show_column_options(['time_limits'], show)

    def on_show_corrections_selection(self, show):
        self._view.show_column_options(['corrections'], show)

    def on_show_masks_selection(self, show):
        self._view.show_column_options(['masks'], show)

    def on_show_radius_wave_selection(self, show):
        self._view.show_column_options(['radius_wave'], show)

    # ------------------------------------------------------------------------------------------------------------------
    # Table + Actions
    # ------------------------------------------------------------------------------------------------------------------
    def set_view(self, view):
        super(BilbyPresenter, self).set_view(view)

    def _add_listener(self):
        listener = BilbyPresenter.ConcreteBilbyListener(self)
        self._view.add_listener(listener)

    def get_hidden_column_groups(self):
        hidden = []
        if not self._view.show_index:
            hidden.append('index')
        if not self._view.show_times:
            hidden.append('time_limits')
        if not self._view.show_masks:
            hidden.append('masks')
        if not self._view.show_corrections:
            hidden.append('corrections')
        if not self._view.show_radius_wave:
            hidden.append('radius_wave')    
        return hidden

    # ----------------------------------------------------------------------------------------------
    # Processing
    # ----------------------------------------------------------------------------------------------

    def update_view_from_user_setting(self):

        super(BilbyPresenter, self).update_view_from_user_setting()

        self._set_on_view_range("wavelength_bins", 
                                start_attr="minimum_wavelength",
                                step_attr="wavelength_step",
                                stop_attr="maximum_wavelength")

        self._set_on_view_range("q1d_bins", 
                                start_attr="minimum_q1d",
                                step_attr="q1d_step",
                                stop_attr="maximum_q1d")

        self._set_on_view_range("transmission_bins", 
                                start_attr="minimum_transmission_wavelength",
                                step_attr="transmission_wavelength_step",
                                stop_attr="maximum_transmission_wavelength")

        self._set_on_view("wavelength_interval")
        self._set_on_view("slice_wavelength")  
        self._set_on_view("qxy_points")
        self._set_on_view("plot_transmission")
        self._set_on_view("save_transmission")
        self._set_on_view("transmission_fit")
        self._set_on_view("polynomial_fit_order")
        self._set_on_view("solid_angle_correction")
        self._set_on_view("gravity_correction")
        self._set_on_view("wide_angle_correction")
        self._set_on_view("blocked_beam_correction")
        self._set_on_view("radius_cut")
        self._set_on_view("wave_cut")    
        self._set_on_view("sample_mask_file")
        self._set_on_view("transmission_mask_file")

    def _create_row_state(self, row, table_model, facility, instrument, file_lookup):

        # creates and populates the state data model
        state = BilbyStateData()

        view = self._view
        row_model = self._table_model.get_table_entry(row)
        if file_lookup:
            file_information = row_model.file_information
            # att pos, pre 2016 data, ..
            state.att_pos = file_information.get_att_pos()
            sample_date = file_information.get_date()
            state.pre_2016_data = (sample_date < DateAndTime("2016-01-01T00:00:00"))
            state.pre_may_2016_data = (sample_date < DateAndTime("2016-05-01T00:00:00"))
            pre_feb_2017 = (sample_date < DateAndTime("2017-02-01T00:00:00"))
        else:
            state.att_pos = 1
            state.pre_2016_data = False
            state.pre_may_2016_data = False
            pre_feb_2017 = False

        # populate the run numbers for sample, transmission, etc
        state.sample = row_model.sample
        state.sample_mask = row_model.sample_mask
        state.blocked_beam = row_model.blocked_beam
        state.transmission_sample = row_model.t_sample
        state.transmission_mask = row_model.transmission_mask
        state.transmission_empty = row_model.t_empty_beam
        state.empty_beam_spectrum = row_model.t_empty_beam

        state.sample_thickness = row_model.thickness
        state.fit_method = _convert_fit_method(view.transmission_fit)
        state.polynomial_order = view.polynomial_fit_order 
        state.time_slice = _create_time_slice_state(row_model.start_time, row_model.end_time)

        state.external_mode = file_information.is_event_mode() if file_lookup else True
        state.gravity_correction = view.gravity_correction
        state.wide_angle_correction = view.wide_angle_correction
        state.solid_angle_weighting = True  # hard coded in the bilby code reduction
        state.reduce_2D = (view.reduction_dimensionality == ReductionDimensionality.TwoDim)
        state.data_points_2D = view.qxy_points

        # wavelength
        state.wavelength = _create_wavelength_slices(view.minimum_wavelength,
                                                            view.maximum_wavelength,
                                                            view.wavelength_step,
                                                            view.wavelength_interval,
                                                            view.slice_wavelength)

        state.transmission_wavelength = _create_range_state(view.minimum_transmission_wavelength,
                                                                      view.maximum_transmission_wavelength,
                                                                      view.transmission_wavelength_step)

        if state.reduce_2D:
            q_step = (view.minimum_q1d + view.maximum_q1d) / view.qxy_points
        else:
            q_step = float(view.q1d_step)
        state.backward_q_bin = q_step < 0.0
        step_sign = -1.0 if state.backward_q_bin else 1.0   
        state.binning_q = _create_range_state(view.minimum_q1d,
                                                   view.maximum_q1d,
                                                   step_sign * q_step)

        state.radius_cut =float(view.radius_cut) if view.radius_cut else 0.0
        state.wave_cut = float(view.wave_cut) if view.wave_cut else 0.0

        state.output_folder = self._set_output_folder( row_model)
        state.suffix = '_'.join([x for x in [row_model.suffix, row_model.description] if x]).replace(' ','_')

        # display warnings
        if state.wide_angle_correction:
            self.sans_logger.warning("WARNING: Enabling wide_angle_correction will lead to the wrong error bar calculations.")
        if pre_feb_2017 and state.radius_cut:
            self.sans_logger.warning("WARNING: Radius cut is not available for pre Feb 2017 data.")
        if pre_feb_2017 and state.wave_cut:
            self.sans_logger.warning("WARNING: Wave cut is not available for pre Feb 2017 data.")
        if file_lookup and state.time_slice.end_time > 1.1 * file_information.get_duration():
            self.sans_logger.warning(
                "WARNING: End time {:.1f}s exceeds the data collection time {:.1f}s.".format(state.time_slice.end_time,
                                                                                             file_information.get_duration()))
        # check multiple reduction covers the wavelength extent (to 2 decimal places) 
        if (round(view.maximum_wavelength - state.wavelength.wavelength_high[-1], 2)):
            self.sans_logger.warning(
                "WARNING: Maximum for multiple reduction does not cover the maximum wavelength: {:.2f}, ({:.2f})".format(
                    state.wavelength.wavelength_high[-1], view.maximum_wavelength))

        if view.slice_wavelength and not view.wavelength_interval:
            self.sans_logger.warning(
                "WARNING: Use wavelength intervals set but wavelength delta is empty")
                                                            
        return state

    def _set_output_folder(self, row_model):
        # merged the save directory with the output folder in the row
        savedir = self._view.save_directory
        subdir = self._view.output_folder
        return os.path.join(savedir, subdir)

