# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
""" The table  model contains all the reduction information which is provided via the data table

The main information in the table model are the run numbers and the selected periods. However it also contains
information regarding the custom output name and the information in the options tab.
"""

from __future__ import (absolute_import, division, print_function)

import functools
import os
import re

from sans.common.enums import RowState
from sans.common.file_information import SANSFileInformationFactory
from sans.gui_logic.presenter.create_file_information import create_file_information
from ui.sans_isis.work_handler import WorkHandler


class RowModel(object):

    @staticmethod
    def column_labels():
        return ["Index","User File"]

    @staticmethod
    def column_keys():
        return ["index", "user_file"]

    @staticmethod
    def column_options():
        return {}

    @staticmethod
    def create_empty_row():
        row = [''] * 2
        return RowModel(*row)
    
    def __init__(self, *argv):
        super(RowModel, self).__init__()
        self.id = None
        self.index = None
        self.user_file = None
        self.output_name = None
        for tag, value in zip(self.column_keys(), argv):
            setattr(self, tag, value)

        self.row_state = RowState.Unprocessed

        self.tool_tip = ''
        self.file_information = None
        self.file_finding = False

    def update_attribute(self, attribute_name, value):
        setattr(self, attribute_name, value)

    def update_attribute_by_column(self, column, value):
        setattr(self, self.column_keys()[column], value)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return self.__dict__ != other.__dict__

    def to_list(self):
        return [getattr(self, tag) for tag in self.column_keys()]

    def to_batch_list(self):
        """
        :return: a list of data in the order as would typically appear
        in a batch file
        """
        return_list = [self.user_file]
        return_list = list(map(str, return_list))
        return_list = list(map(str.strip, return_list))
        return return_list

    def includes_options(self):
        groups = []
        attrs = self.column_keys()
        for key, cols in self.column_options().items():
            if any(map(functools.partial(getattr, self), map(attrs.__getitem__,cols))):
                groups.append(key)
        return groups

    def is_empty(self):
        return not any(map(functools.partial(getattr, self), self.column_keys()))


class TableModel(object):

    THICKNESS_ROW = 14

    def __init__(self, rowModelClass):
        super(TableModel, self).__init__()
        self._RowModel = rowModelClass
        self._user_file = ""
        self._batch_file = ""
        self._table_entries = []
        self.work_handler = WorkHandler()
        self._subscriber_list = []
        self._id_count = 0

    @staticmethod
    def _validate_file_name(file_name):
        if not file_name:
            return
        if not os.path.exists(file_name):
            raise ValueError("The file {} does not seem to exist.".format(file_name))

    @property
    def user_file(self):
        return self._user_file

    @user_file.setter
    def user_file(self, value):
        self._user_file = value

    def get_row_user_file(self, row_index):
        if row_index < len(self._table_entries):
            return self._table_entries[row_index].user_file
        else:
            raise IndexError("The row {} does not exist.".format(row_index))

    @property
    def batch_file(self):
        return self._batch_file

    @batch_file.setter
    def batch_file(self, value):
        self._batch_file = value

    def column_labels(self):
        return self._RowModel.column_labels()

    def column_keys(self):
        return self._RowModel.column_keys()

    def column_options(self):
        return self._RowModel.column_options()

    def create_empty_row(self):
        return self._RowModel.create_empty_row()

    def create_row(self, *argv):
        return self._RowModel(*argv)

    def get_table_entry(self, index):
        return self._table_entries[index]

    def add_table_entry(self, row, table_index_model):
        table_index_model.id = self._id_count
        self._id_count += 1
        self._table_entries.insert(row, table_index_model)
        if row >= self.get_number_of_rows():
            row = self.get_number_of_rows() - 1
        self.get_thickness_for_rows([row])
        self.notify_subscribers()

    def append_table_entry(self, table_index_model):
        table_index_model.id = self._id_count
        self._id_count += 1
        self._table_entries.append(table_index_model)
        self.get_thickness_for_rows([self.get_number_of_rows() - 1])
        self.notify_subscribers()

    def remove_table_entries(self, rows):
        # For speed rows should be a Set here but don't think it matters for the list sizes involved.
        self._table_entries[:] = [item for i,item in enumerate(self._table_entries) if i not in rows]
        if not self._table_entries:
            row_index_model = self.create_empty_row()
            self.append_table_entry(row_index_model)
        else:
            self.notify_subscribers()

    def replace_table_entries(self, row_to_replace_index, rows_to_insert):
        self.remove_table_entries(row_to_replace_index)
        for row_entry in reversed(rows_to_insert):
            self.add_table_entry(row_to_replace_index[0], row_entry)

    def clear_table_entries(self):
        self._table_entries = []
        row_index_model = self.create_empty_row()
        self.append_table_entry(row_index_model)

    def get_number_of_rows(self):
        return len(self._table_entries)

    def update_table_entry(self, row, column, value):
        self._table_entries[row].update_attribute_by_column(column, value)
        self._table_entries[row].update_attribute('row_state', RowState.Unprocessed)
        self._table_entries[row].update_attribute('tool_tip', '')
        if column == 0:
            self.get_thickness_for_rows([row])
        self.notify_subscribers()

    def is_empty_row(self, row):
        return self._table_entries[row].is_empty()

    def get_non_empty_rows(self, rows):
        return list(filter(lambda x: not self.get_table_entry(x).is_empty(), rows))

    def set_row_to_processed(self, row, tool_tip):
        self._table_entries[row].update_attribute('row_state', RowState.Processed)
        self._table_entries[row].update_attribute('tool_tip', tool_tip)
        self.notify_subscribers()

    def reset_row_state(self, row):
        self._table_entries[row].update_attribute('row_state', RowState.Unprocessed)
        self._table_entries[row].update_attribute('tool_tip', '')
        self.notify_subscribers()

    def set_row_to_error(self, row, tool_tip):
        self._table_entries[row].update_attribute('row_state', RowState.Error)
        self._table_entries[row].update_attribute('tool_tip', tool_tip)
        self.notify_subscribers()

    def get_thickness_for_rows(self, rows=None):
        """
        Read in the sample thickness for the given rows from the file and set it in the table.
        :param rows: list of table rows
        """
        if not rows:
            rows = range(len(self._table_entries))
        for row in rows:
            entry = self._table_entries[row]
            if entry.is_empty():
                continue
            entry.file_finding = True
            success_callback = functools.partial(self.update_thickness_from_file_information, entry.id)

            error_callback = functools.partial(self.failure_handler, entry.id)
            create_file_information(entry.sample_scatter, error_callback, success_callback,
                                    self.work_handler, entry.id)

    def failure_handler(self, id, error):
        row = self.get_row_from_id(id)
        self._table_entries[row].update_attribute('file_information', '')
        self._table_entries[row].file_finding = False
        self.set_row_to_error(row, str(error[1]))

    def update_thickness_from_file_information(self, id, file_information):
        row = self.get_row_from_id(id)
        if file_information:
            rounded_file_thickness = round(file_information.get_thickness(), 2)
            self._table_entries[row].update_attribute('file_information', file_information)
            self._table_entries[row].file_finding = False
            self.reset_row_state(row)

    def subscribe_to_model_changes(self, subscriber):
        self._subscriber_list.append(subscriber)

    def notify_subscribers(self):
        for subscriber in self._subscriber_list:
            subscriber.on_update_rows()

    def get_file_information_for_row(self, row):
        return self._table_entries[row].file_information

    def get_row_from_id(self, id):
        for row, entry in enumerate(self._table_entries):
            if entry.id == id:
                return row
        return None

    def wait_for_file_finding_done(self):
        self.work_handler.wait_for_done()

    def wait_for_file_information(self, row):
        if self._table_entries[row].file_finding:
            self.wait_for_file_finding_done()

    def add_table_entry_no_thread_or_signal(self, row, table_index_model):
        table_index_model.id = self._id_count
        self._id_count += 1
        self._table_entries.insert(row, table_index_model)
        if row >= self.get_number_of_rows():
            row = self.get_number_of_rows() - 1

        entry = self._table_entries[row]
        file_information_factory = SANSFileInformationFactory()
        file_information = file_information_factory.create_sans_file_information(entry.sample_scatter)
        self.update_thickness_from_file_information(entry.id, file_information)

    def set_option(self, row, key, value):
        self._table_entries[row].options_column_model.set_option(key, value)

    def __eq__(self, other):
        return self.equal_dicts(self.__dict__, other.__dict__, ['work_handler'])

    def __ne__(self, other):
        return not self.equal_dicts(self.__dict__, other.__dict__, ['work_handler'])

    @staticmethod
    def equal_dicts(d1, d2, ignore_keys):
        d1_filtered = dict((k, v) for k, v in d1.items() if k not in ignore_keys)
        d2_filtered = dict((k, v) for k, v in d2.items() if k not in ignore_keys)
        return d1_filtered == d2_filtered

def options_column_bool(string):
    """
    Evaluate input string as a bool. Used for UseMirror in Options column,
    as evaluating bool("False") returns True (any string is Truthy).
    :param string: User input string to be evaluated as False or True
    :return: True or False
    """
    truthy_strings = ("true", "1", "yes", "t", "y")  # t short for true, y short for yes
    falsy_strings = ("false", "0", "no", "f", "n")
    if string.lower() in truthy_strings:
        return True
    elif string.lower() in falsy_strings:
        return False
    else:
        raise ValueError("Could not evaluate {} as a boolean value. It should be True or False.".format(string))