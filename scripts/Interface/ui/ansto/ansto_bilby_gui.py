# -*- coding: utf-8 -*-
# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
""" Core view for the run tab view to support batch processing.
"""

from __future__ import (absolute_import, division, print_function)

from abc import ABCMeta, abstractmethod
from qtpy.QtWidgets import QListWidgetItem  # noqa

from mantidqt import icons
from mantidqt.utils.qt import load_ui

from sans.common.enums import (FitType, ReductionDimensionality)
from ui.ansto.run_tab_gui import RunTabGui

from qtpy import PYQT4
if PYQT4:
    IN_MANTIDPLOT = False
    try:
        from pymantidplot import proxies
        IN_MANTIDPLOT = True
    except ImportError:
        # We are not in MantidPlot e.g. testing
        pass

Ui_AnstoBilbyWindow, _ = load_ui(__file__, "ansto_bilby_window.ui")

# ----------------------------------------------------------------------------------------------------------------------
# Gui Classes
# ----------------------------------------------------------------------------------------------------------------------
class BilbyBatchReductionGui(RunTabGui):

    data_processor_table = None

    class BilbyListener(RunTabGui.RunTabListener):
        """
        Defines the additional elements which a presenter can listen to in this View
        """

        @abstractmethod
        def on_show_corrections_selection(self, show):
            pass

        @abstractmethod
        def on_show_masks_selection(self, show):
            pass

        @abstractmethod
        def on_show_radius_wave_selection(self, show):
            pass

    def __init__(self, ui_form, parent=None):
        super(BilbyBatchReductionGui, self).__init__(ui_form, parent=parent)

    def _setup_connections(self):
        super(BilbyBatchReductionGui, self)._setup_connections()

    def _setup_page_tabs(self, list_widget):

        super(BilbyBatchReductionGui, self)._setup_page_tabs(list_widget)
        settings_icon = icons.get_icon("fa.cog")
        _ = QListWidgetItem(settings_icon, "Settings", list_widget)  # noqa

    def setup_layout(self, all_columns, column_groups):
        super(BilbyBatchReductionGui, self).setup_layout(all_columns, column_groups)

        # --------------------------------------------------------------------------------------------------------------
        # Runs Page
        # --------------------------------------------------------------------------------------------------------------
        self._ui.show_corrections_checkbox.stateChanged.connect(self._on_show_corrections_selection)
        self._ui.show_masks_checkbox.stateChanged.connect(self._on_show_masks_selection)
        self._ui.show_radius_wave_checkbox.stateChanged.connect(self._on_show_radius_wave_selection)

        # --------------------------------------------------------------------------------------------------------------
        # Settings Page : General Tab
        # --------------------------------------------------------------------------------------------------------------        
        # Set the transmission fit selection
        self._ui.transmission_fit_combo.currentIndexChanged.connect(self._on_transmission_fit_selection_changed)
        self._on_transmission_fit_selection_changed()


    # Bilby specific listener handlers

    def _on_show_corrections_selection(self):
        self._call_settings_listeners(
            lambda listener: listener.on_show_corrections_selection(self.show_corrections))

    def _on_show_masks_selection(self):
        self._call_settings_listeners(
            lambda listener: listener.on_show_masks_selection(self.show_masks))

    def _on_show_radius_wave_selection(self):
        self._call_settings_listeners(
            lambda listener: listener.on_show_radius_wave_selection(self.show_radius_wave))

    def _on_transmission_fit_selection_changed(self):
        pass

    # ==================================================================================================================
    # Runs Page properties
    # ==================================================================================================================
    @property
    def show_corrections(self):
        return self._ui.show_corrections_checkbox.isChecked()

    @show_corrections.setter
    def show_corrections(self, value):
        self._ui.show_corrections_checkbox.setChecked(value)

    @property
    def show_masks(self):
        return self._ui.show_masks_checkbox.isChecked()

    @show_masks.setter
    def show_masks(self, value):
        self._ui.show_masks_checkbox.setChecked(value)

    @property
    def show_radius_wave(self):
        return self._ui.show_radius_wave_checkbox.isChecked()

    @show_radius_wave.setter
    def show_radius_wave(self, value):
        self._ui.show_radius_wave_checkbox.setChecked(value)

    @property
    def reduction_dimensionality(self):
        return ReductionDimensionality.OneDim if self._ui.reduction_dimensionality_1D.isChecked() \
            else ReductionDimensionality.TwoDim

    @reduction_dimensionality.setter
    def reduction_dimensionality(self, value):
        is_1d = value is ReductionDimensionality.OneDim
        self._ui.reduction_dimensionality_1D.setChecked(is_1d)
        self._ui.reduction_dimensionality_2D.setChecked(not is_1d)

    @property
    def plot_results(self):
        return self._ui.plot_results_checkbox.isChecked()

    @plot_results.setter
    def plot_results(self, value):
        self._ui.plot_results_checkbox.setChecked(value)

    # ==================================================================================================================
    # Settings Page properties : Reduction Tab
    # ==================================================================================================================

    @property
    def minimum_wavelength(self):
        return self.get_simple_line_edit_field(line_edit="wavelength_min_edit",
                                               expected_type=float)

    @minimum_wavelength.setter
    def minimum_wavelength(self, value):
        self.update_simple_line_edit_field(line_edit="wavelength_min_edit",
                                           value=value)

    @property
    def maximum_wavelength(self):
        return self.get_simple_line_edit_field(line_edit="wavelength_max_edit",
                                               expected_type=float)

    @maximum_wavelength.setter
    def maximum_wavelength(self, value):
        self.update_simple_line_edit_field(line_edit="wavelength_max_edit",
                                           value=value)

    @property
    def wavelength_step(self):
        return self.get_simple_line_edit_field(line_edit="wavelength_step_edit",
                                               expected_type=float)

    @wavelength_step.setter
    def wavelength_step(self, value):
        self.update_simple_line_edit_field(line_edit="wavelength_step_edit",
                                           value=value)

    #---------------------------------------------------------------------------

    @property
    def minimum_q1d(self):
        return self.get_simple_line_edit_field(line_edit="q1d_min_edit",
                                               expected_type=float)

    @minimum_q1d.setter
    def minimum_q1d(self, value):
        self.update_simple_line_edit_field(line_edit="q1d_min_edit",
                                           value=value)

    @property
    def maximum_q1d(self):
        return self.get_simple_line_edit_field(line_edit="q1d_max_edit",
                                               expected_type=float)

    @maximum_q1d.setter
    def maximum_q1d(self, value):
        self.update_simple_line_edit_field(line_edit="q1d_max_edit",
                                           value=value)

    @property
    def q1d_step(self):
        return self.get_simple_line_edit_field(line_edit="q1d_step_edit",
                                               expected_type=float)

    @q1d_step.setter
    def q1d_step(self, value):
        self.update_simple_line_edit_field(line_edit="q1d_step_edit",
                                           value=value)

    #---------------------------------------------------------------------------

    @property
    def qxy_interval(self):
        return self.get_simple_line_edit_field(line_edit="qxy_interval_edit",
                                               expected_type=float)

    @qxy_interval.setter
    def qxy_interval(self, value):
        self.update_simple_line_edit_field(line_edit="qxy_interval_edit",
                                           value=value)

    @property
    def qxy_points(self):
        return self.get_simple_line_edit_field(line_edit="qxy_data_points_edit",
                                               expected_type=float)

    @qxy_points.setter
    def qxy_points(self, value):
        self.update_simple_line_edit_field(line_edit="qxy_data_points_edit",
                                           value=value)
    
    # ==================================================================================================================
    # Settings Page properties : Transmission Tab
    # ==================================================================================================================    

    @property
    def transmission_fit(self):
        return self._ui.transmission_fit_combo.currentIndex()

    @transmission_fit.setter
    def transmission_fit(self, value):
        self._ui.transmission_fit_combo.setCurrentIndex(value)

    @property
    def transmission_fit(self):
        fit_type_as_string = self._ui.transmission_fit_combo.currentText().encode('utf-8')
        return FitType.from_string(fit_type_as_string)

    @transmission_fit.setter
    def transmission_fit(self, value):
        if value is None:
            self._ui.transmission_fit_combo.setCurrentIndex(0)
        else:
            self.update_gui_combo_box(value=value, expected_type=FitType,
                                      combo_box="transmission_fit_combo")

    @property
    def polynomial_fit_order(self):
        return self._ui.poylnominal_fit_order_spinner.value()

    @polynomial_fit_order.setter
    def polynomial_fit_order(self, value):
        self._set_polynomial_order(self._ui.polynomial_fit_order_spinner, value)

    @staticmethod
    def _set_polynomial_order(spin_box, value):
        minimum = spin_box.minimum()
        maximum = spin_box.maximum()
        if value < minimum or value > maximum:
            raise ValueError("The value for the polynomial order {} has "
                             "to be in the range of {} and {}".format(value, minimum, maximum))
        spin_box.setValue(value)

    #---------------------------------------------------------------------------
  
    @property
    def minimum_transmission_wavelength(self):
        return self.get_simple_line_edit_field(line_edit="t_wavelength_min_edit",
                                               expected_type=float)

    @minimum_transmission_wavelength.setter
    def minimum_transmission_wavelength(self, value):
        self.update_simple_line_edit_field(line_edit="t_wavelength_min_edit",
                                           value=value)

    @property
    def maximum_transmission_wavelength(self):
        return self.get_simple_line_edit_field(line_edit="t_wavelength_max_edit",
                                               expected_type=float)

    @maximum_transmission_wavelength.setter
    def maximum_transmission_wavelength(self, value):
        self.update_simple_line_edit_field(line_edit="t_wavelength_max_edit",
                                           value=value)

    @property
    def transmission_wavelength_step(self):
        return self.get_simple_line_edit_field(line_edit="t_wavelength_step_edit",
                                               expected_type=float)

    @transmission_wavelength_step.setter
    def transmission_wavelength_step(self, value):
        self.update_simple_line_edit_field(line_edit="t_wavelength_step_edit",
                                           value=value)

    #---------------------------------------------------------------------------

    @property
    def plot_transmission(self):
        return self._ui.plot_transmission_checkbox.isChecked()

    @plot_transmission.setter
    def plot_transmission(self, value):
        self._ui.plot_transmission_checkbox.setChecked(value)

    @property
    def save_transmission(self):
        return self._ui.save_transmission_checkbox.isChecked()

    @save_transmission.setter
    def save_transmission(self, value):
        self._ui.save_transmission_checkbox.setChecked(value)

    # ==================================================================================================================
    # Settings Page properties : Advanced Tab
    # ==================================================================================================================

    @property
    def wide_angle_correction(self):
        return self._ui.wide_angle_correction_checkbox.isChecked()

    @wide_angle_correction.setter
    def wide_angle_correction(self, value):
        self._ui.wide_angle_correction_checkbox.setChecked(value)

    @property
    def gravity_correction(self):
        return self._ui.gravity_correction_checkbox.isChecked()

    @gravity_correction.setter
    def gravity_correction(self, value):
        self._ui.gravity_correction_checkbox.setChecked(value)

    @property
    def blocked_beam_correction(self):
        return self._ui.blocked_beam_correction_checkbox.isChecked()

    @blocked_beam_correction.setter
    def blocked_beam_correction(self, value):
        self._ui.blocked_beam_correction_checkbox.setChecked(value)

    @property
    def transmission_mask_file(self):
        return self.get_simple_line_edit_field(line_edit="transmission_mask_edit",
                                               expected_type=str)

    @transmission_mask_file.setter
    def transmission_mask_files(self, value):
        self.update_simple_line_edit_field(line_edit="transmission_mask_edit",
                                           value=value)

    @property
    def sample_mask_file(self):
        return self.get_simple_line_edit_field(line_edit="sample_mask_edit",
                                               expected_type=str)

    @transmission_mask_file.setter
    def sample_mask_files(self, value):
        self.update_simple_line_edit_field(line_edit="sample_mask_edit",
                                           value=value)

    # ==================================================================================================================
    # Settings Page properties : Resolution Tab
    # ==================================================================================================================

    @property
    def radius_cut(self):
        return self.get_simple_line_edit_field(line_edit="radius_cut_edit",
                                               expected_type=float)

    @radius_cut.setter
    def radius_cut(self, value):
        self.update_simple_line_edit_field(line_edit="radius_cut_edit",
                                           value=value)

    #---------------------------------------------------------------------------

    @property
    def wave_cut(self):
        return self.get_simple_line_edit_field(line_edit="wave_cut_edit",
                                               expected_type=float)

    @wave_cut.setter
    def wave_cut(self, value):
        self.update_simple_line_edit_field(line_edit="wave_cut_edit",
                                           value=value)

