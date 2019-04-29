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

from sans.common.enums import FitType
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
        super(BilbyPresenter, self)._default_gui_setup()

        fit_types = [FitType.to_string(FitType.Linear),
                     FitType.to_string(FitType.Logarithmic),
                     FitType.to_string(FitType.Polynomial)]
        self._view.transmission_fit = fit_types

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

        self._set_on_view("maximum_wavelength")
        self._set_on_view("minimum_wavelength")
        self._set_on_view("wavelength_step")

        self._set_on_view("maximum_q1d")
        self._set_on_view("minimum_q1d")
        self._set_on_view("q1d_step")

        self._set_on_view("qxy_interval")
        self._set_on_view("qxy_points")        

        self._set_on_view("maximum_transmission_wavelength")
        self._set_on_view("minimum_transmission_wavelength")
        self._set_on_view("transmission_wavelength_step")
        self._set_on_view("plot_transmission")
        self._set_on_view("save_transmission")
        self._set_on_view("transmission_fit")
        self._set_on_view("polynomial_fit_order")

        self._set_on_view("gravity_correction")
        self._set_on_view("wide_angle_correction")
        self._set_on_view("blocked_beam_correction")
        self._set_on_view("radius_cut")
        self._set_on_view("wave_cut")            

    def _set_on_view_transmission_fit(self):

        # Set transmission_sample_can_fit
        fit_type = self._state_model.transmission_fit
        self._view.transmission_fit = fit_type

        # Set the polynomial order for can
        polynomial_order = self._state_model.polynomial_order if fit_type is FitType.Polynomial else 2  # noqa
        self._view.polynomial_order = polynomial_order

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

