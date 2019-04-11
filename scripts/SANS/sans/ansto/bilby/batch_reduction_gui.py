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
from sans.ansto.run_tab_gui import RunTabGui

from qtpy import PYQT4
if PYQT4:
    IN_MANTIDPLOT = False
    try:
        from pymantidplot import proxies
        IN_MANTIDPLOT = True
    except ImportError:
        # We are not in MantidPlot e.g. testing
        pass

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
        def on_corrections_selection(self, show):
            pass

        @abstractmethod
        def on_masks_selection(self, show):
            pass

        @abstractmethod
        def on_radius_wave_selection(self, show):
            pass

    def __init__(self, ui_form, parent=None):
        super(BilbyBatchReductionGui, self).__init__(ui_form, parent=parent)

    def _setup_connections():
        super(BilbyBatchReductionGui, self)._setup_connections()

    def _setup_page_tabs():

        super(BilbyBatchReductionGui, self)._setup_page_tabs()
        settings_icon = icons.get_icon("fa.cog")
        _ = QListWidgetItem(settings_icon, "Settings", self._ui.tab_choice_list)  # noqa

    def _complete_layout():

        # add the connections
        pass