# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantid workbench
#
#
from __future__ import absolute_import, unicode_literals

from workbench.widgets.settings.newtab.view import NewTabSettingsView


class NewTabSettings(object):

    def __init__(self, parent, view=None):
        self.view = view if view else NewTabSettingsView(parent, self)
        self.parent = parent