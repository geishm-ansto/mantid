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

import copy
import csv
import os
import sys
import time
import traceback

#from mantid.api import (FileFinder)
#from mantid.kernel import Logger, ConfigService

#from sans.common.constants import ALL_PERIODS
#from sans.common.enums import (BatchReductionEntry, RangeStepType, SampleShape, FitType, RowState, SANSInstrument)
#from sans.gui_logic.gui_common import (get_reduction_mode_strings_for_gui, get_string_for_gui_from_instrument,
#                                       add_dir_to_datasearch, remove_dir_from_datasearch)
#from sans.gui_logic.models.batch_process_runner import BatchProcessRunner
#from sans.gui_logic.models.create_state import create_states
#from sans.gui_logic.models.diagnostics_page_model import create_state
#from sans.gui_logic.models.state_gui_model import StateGuiModel
#from sans.gui_logic.presenter.add_runs_presenter import OutputDirectoryObserver as SaveDirectoryObserver
#from sans.user_file.user_file_reader import UserFileReader

#from ui.sans_isis import SANSSaveOtherWindow
#from ui.sans_isis.work_handler import WorkHandler

from sans.ansto.run_tab_presenter import RunTabPresenter
from ui.ansto.ansto_bilby_gui import BilbyBatchReductionGui

class BilbyPresenter(RunTabPresenter):
    class ConcreteBilbyListener(RunTabPresenter.ConcreteRunTabListener, BilbyBatchReductionGui.BilbyListener):
        def __init__(self, presenter):
            super(BilbyPresenter.ConcreteBilbyListener, self).__init__(presenter)

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
        pass

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

    # ----------------------------------------------------------------------------------------------
    # Processing
    # ----------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------
    # Controls
    # ------------------------------------------------------------------------------------------------------------------
    def disable_controls(self):
        """
        Disable all input fields and buttons during the execution of the reduction.
        """
        # TODO: think about enabling and disable some controls during reduction
        pass

    def enable_controls(self):
        """
        Enable all input fields and buttons after the execution has completed.
        """
        # TODO: think about enabling and disable some controls during reduction
        pass


    def _update_view_from_state_model(self):

        self._set_on_view("gravity_correction")
        self._set_on_view("maximum_wavelength")
        self._set_on_view("minimum_wavelength")
        self._set_on_view("wavelength_step")

        self._set_on_view("maximum_q1d")
        self._set_on_view("minimum_q1d")
        self._set_on_view("q1d_step")

        self._set_on_view("qxy_interval")
        self._set_on_view("qxy_points")        

        self._set_on_view("radius_cut")
        self._set_on_view("wave_cut")            
        # TODO Add remaining terms

    def _set_on_view_transmission_fit_sample_settings(self):
        # Set transmission_sample_use_fit
        fit_type = self._state_model.transmission_sample_fit_type
        use_fit = fit_type is not FitType.NoFit
        self._view.transmission_sample_use_fit = use_fit

        # Set the polynomial order for sample
        polynomial_order = self._state_model.transmission_sample_polynomial_order if fit_type is FitType.Polynomial else 2  # noqa
        self._view.transmission_sample_polynomial_order = polynomial_order

        # Set the fit type for the sample
        fit_type = fit_type if fit_type is not FitType.NoFit else FitType.Linear
        self._view.transmission_sample_fit_type = fit_type

        # Set the wavelength
        wavelength_min = self._state_model.transmission_sample_wavelength_min
        wavelength_max = self._state_model.transmission_sample_wavelength_max
        if wavelength_min and wavelength_max:
            self._view.transmission_sample_use_wavelength = True
            self._view.transmission_sample_wavelength_min = wavelength_min
            self._view.transmission_sample_wavelength_max = wavelength_max

    def _set_on_view_transmission_fit(self):
        # Steps for adding the transmission fit to the view
        # 1. Check if individual settings exist. If so then set the view to separate, else set them to both
        # 2. Apply the settings
        separate_settings = self._state_model.has_transmission_fit_got_separate_settings_for_sample_and_can()
        self._view.set_fit_selection(use_separate=separate_settings)

        if separate_settings:
            self._set_on_view_transmission_fit_sample_settings()

            # Set transmission_sample_can_fit
            fit_type_can = self._state_model.transmission_can_fit_type()
            use_can_fit = fit_type_can is FitType.NoFit
            self._view.transmission_can_use_fit = use_can_fit

            # Set the polynomial order for can
            polynomial_order_can = self._state_model.transmission_can_polynomial_order if fit_type_can is FitType.Polynomial else 2  # noqa
            self._view.transmission_can_polynomial_order = polynomial_order_can

            # Set the fit type for the can
            fit_type_can = fit_type_can if fit_type_can is not FitType.NoFit else FitType.Linear
            self.transmission_can_fit_type = fit_type_can

            # Set the wavelength
            wavelength_min = self._state_model.transmission_can_wavelength_min
            wavelength_max = self._state_model.transmission_can_wavelength_max
            if wavelength_min and wavelength_max:
                self._view.transmission_can_use_wavelength = True
                self._view.transmission_can_wavelength_min = wavelength_min
                self._view.transmission_can_wavelength_max = wavelength_max
        else:
            self._set_on_view_transmission_fit_sample_settings()

    def _set_on_view(self, attribute_name):
        attribute = getattr(self._state_model, attribute_name)
        if attribute or isinstance(attribute,
                                   bool):  # We need to be careful here. We don't want to set empty strings, or None, but we want to set boolean values. # noqa
            setattr(self._view, attribute_name, attribute)

    def _set_on_view_with_view(self, attribute_name, view):
        attribute = getattr(self._state_model, attribute_name)
        if attribute or isinstance(attribute,
                                   bool):  # We need to be careful here. We don't want to set empty strings, or None, but we want to set boolean values. # noqa
            setattr(view, attribute_name, attribute)

    def _get_state_model_with_view_update(self):
        """
        Goes through all sub presenters and update the state model based on the views.

        Note that at the moment we have set up the view and the model such that the name of a property must be the same
        in the view and the model. This can be easily changed, but it also provides a good cohesion.
        """
        state_model = copy.deepcopy(self._state_model)

        # If we don't have a state model then return None
        if state_model is None:
            return state_model
        # Run tab view
        self._set_on_state_model("gravity_correction", state_model)
        self._set_on_state_model("maximum_wavelength", state_model)
        self._set_on_state_model("minimum_wavelength", state_model)
        self._set_on_state_model("wavelength_step", state_model)
        self._set_on_state_model("maximum_q1d", state_model)
        self._set_on_state_model("minimum_q1d", state_model)
        self._set_on_state_model("q1d_step", state_model)
        self._set_on_state_model("qxy_interval", state_model)
        self._set_on_state_model("qxy_points", state_model)
        self._set_on_state_model("radius_cut", state_model)
        self._set_on_state_model("wave_cut", state_model)

        return state_model

