# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
""" The state gui model contains all the reduction information which is not explicitly available in the data table.

This is one of the two models which is used for the data reduction. It contains generally all the settings which
are not available in the model associated with the data table.
"""

from __future__ import (absolute_import, division, print_function)

import sys

from sans.user_file.settings_tags import simple_range

class AnstoStateGuiModel(object):
    def __init__(self, user_file_items):
        super(AnstoStateGuiModel, self).__init__()
        self._user_file_items = user_file_items

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    @property
    def settings(self):
        return self._user_file_items

    def get_simple_element(self, element_id, default_value):
        return self.get_simple_element_with_attribute(element_id, default_value)

    def set_simple_element(self, element_id, value):
        if element_id in self._user_file_items:
            del self._user_file_items[element_id]
        new_state_entries = {element_id: [value]}
        self._user_file_items.update(new_state_entries)

    def get_simple_element_with_attribute(self, element_id, default_value, attribute=None):
        if element_id in self._user_file_items:
            element = self._user_file_items[element_id][-1]
            return getattr(element, attribute) if attribute else element
        else:
            return default_value

    def get_bool_element(self, element_id, default_value):
        true_values = ['t', 'true', 'y', 'yes']
        if element_id in self._user_file_items:
            value = self._user_file_items[element_id][-1]
            return value and value.lower() in true_values
        else:
            return default_value

    def set_bool_element(self, element_id, value):
        self._user_file_items[element_id] = ['T'] if value else ['F']

    def update_simple_range(element_id, start=None, stop=None, step=None, step_type=None):
        if not any([start, stop, step, step_type]):
            return
        if element_id in self._user_file_items:
            settings = self._user_file_items[element_id]
        else:
            # If the entry does not already exist, then add it.
            settings = [simple_range(start='', stop='', step='', step_type='')]

        new_settings = []
        for setting in settings:
            new_start = start if start else setting.start
            new_stop = stop if stop else setting.stop
            new_step = step if step else setting.step
            new_step_type = step_type if step_type else setting.step_type

            new_setting = simple_range(start=new_start, stop=new_stop, step=new_step, step_type=new_step_type)
            new_settings.append(new_setting)
        self._user_file_items[element_id] = new_settings



