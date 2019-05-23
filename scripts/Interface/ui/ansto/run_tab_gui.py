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
from inspect import isclass
from qtpy.QtWidgets import (QListWidgetItem, QMessageBox, QFileDialog, QMainWindow)  # noqa
from qtpy.QtCore import (QRegExp, QSettings)  # noqa
from qtpy.QtGui import (QDoubleValidator, QIcon, QIntValidator, QRegExpValidator)  # noqa
from six import with_metaclass

from reduction_gui.reduction.scripter import execute_script
from mantid.kernel import (Logger)
from mantidqt import icons
from mantidqt.utils.qt import load_ui
from mantidqt.widgets import jobtreeview, manageuserdirectories
from sans.common.enums import (ReductionDimensionality, OutputMode, SaveType, SANSInstrument,
                               RangeStepType, ReductionMode, FitType)
from sans.common.file_information import SANSFileInformationFactory
from sans.gui_logic.gui_common import (GENERIC_SETTINGS,
                                       load_file, load_default_file, set_setting)

from ui.sans_isis.work_handler import WorkHandler
from ui.sans_isis.SANSSaveOtherWindow import SANSSaveOtherDialog

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
class RunTabGui(QMainWindow):

    data_processor_table = None

    class RunTabListener(with_metaclass(ABCMeta, object)):
        """
        Defines the elements which a presenter can listen to in this View
        """

        @abstractmethod
        def on_user_file_load(self):
            pass

        @abstractmethod
        def on_batch_file_load(self):
            pass

        @abstractmethod
        def on_process_selected_clicked(self):
            pass

        @abstractmethod
        def on_process_all_clicked(self):
            pass

        @abstractmethod
        def on_cancel_processing_clicked(self):
            pass

        @abstractmethod
        def on_save_directory_clicked(self):
            pass

        @abstractmethod
        def on_default_directory_clicked(self, show):
            pass

        @abstractmethod
        def on_data_changed(self, row, column, new_value, old_value):
            pass

        @abstractmethod
        def on_row_inserted(self):
            pass

        @abstractmethod
        def on_rows_removed(self):
            pass

        @abstractmethod
        def on_copy_rows_requested(self):
            pass

        @abstractmethod
        def on_paste_rows_requested(self):
            pass

        @abstractmethod
        def on_insert_row(self):
            pass

        @abstractmethod
        def on_erase_rows(self):
            pass

        @abstractmethod
        def on_cut_rows(self):
            pass

    def __init__(self, ui_form, parent=None):
        """
        Initialise the interface
        """
        super(QMainWindow, self).__init__(parent)
        self._ui = ui_form()
        self._ui.setupUi(self)

        # Listeners allow us to to notify all presenters
        self._settings_listeners = []

        # Q Settings
        self.__generic_settings = GENERIC_SETTINGS
        self.__path_key = "sans_path"
        self.__user_file_key = "user_file"
        self.__mask_file_input_path_key = "mask_files"

        # Logger
        self.gui_logger = Logger("ANSTO GUI LOGGER")

        self.instrument = SANSInstrument.NoInstrument

        self._setup_connections()

        # Attach validators
        self._attach_validators()

        self._setup_progress_bar()

    def _setup_connections(self):
        self._ui.paste_button.setIcon(icons.get_icon("fa.paste"))
        self._ui.copy_button.setIcon(icons.get_icon("fa.copy"))
        self._ui.cut_button.setIcon(icons.get_icon("fa.cut"))
        self._ui.erase_button.setIcon(icons.get_icon("fa.eraser"))
        self._ui.delete_row_button.setIcon(icons.get_icon("fa.trash"))
        self._ui.insert_row_button.setIcon(icons.get_icon("fa.table"))

        self._ui.paste_button.clicked.connect(self._paste_rows_requested)
        self._ui.copy_button.clicked.connect(self._copy_rows_requested)
        self._ui.erase_button.clicked.connect(self._erase_rows)
        self._ui.cut_button.clicked.connect(self._cut_rows)

        self._ui.delete_row_button.clicked.connect(self._remove_rows_requested_from_button)
        self._ui.insert_row_button.clicked.connect(self._on_insert_button_pressed)

    def _setup_progress_bar(self):
        self._ui.batch_progress_bar.setMinimum(0)
        self._ui.batch_progress_bar.setMaximum(1)
        self._ui.batch_progress_bar.setValue(0)

    def add_listener(self, listener):
        if not isinstance(listener, RunTabGui.RunTabListener):
            raise ValueError(
                "The listener is not of type RunTabListener but rather {}".format(type(listener)))
        self._settings_listeners.append(listener)

    def clear_listeners(self):
        self._settings_listeners = []

    def _call_settings_listeners(self, target):
        for listener in self._settings_listeners:
            target(listener)

    def set_current_page(self, index):
        self._ui.main_stacked_widget.setCurrentIndex(index)

    def setup_layout(self, all_columns, column_options, hidden_groups):
        """
        Do further setup that could not be done in the designer.
        So far only two menus have been added, we need to add the processing table manually.
        """
        # --------------------------------------------------------------------------------------------------------------
        # Tab selection
        # --------------------------------------------------------------------------------------------------------------
        self._setup_page_tabs(self._ui.tab_choice_list)
        # Set the 0th row enabled
        self._ui.tab_choice_list.setCurrentRow(0)

        # --------------------------------------------------------------------------------------------------------------
        # Main Tab
        # --------------------------------------------------------------------------------------------------------------
		# Initial setup hides all the columns listed in the groups
        self._column_options = column_options
        self._create_data_table(all_columns, hidden_groups)

        self._setup_main_tab()
        self.reset_all_fields_to_default()

        # Add the ui connections
        self._ui.process_selected_button.clicked.connect(self._process_selected_clicked)
        self._ui.process_all_button.clicked.connect(self._process_all_clicked)
        self._ui.cancel_processing_button.clicked.connect(self._cancel_processing_clicked)
        self._ui.save_directory_button.clicked.connect(self._save_directory_clicked)
        self._ui.export_table_button.clicked.connect(self._export_table_clicked)
        self._ui.help_button.clicked.connect(self._on_help_button_clicked)
        self._ui.default_directory_checkbox.stateChanged.connect(self._default_directory_clicked)

        return True

    def _setup_page_tabs(self, list_widget):

        list_widget.setAlternatingRowColors(True)
        list_widget.setSpacing(10)
        list_widget.currentRowChanged.connect(self.set_current_page)
        self.set_current_page(0)

        runs_icon = icons.get_icon("fa.play-circle-o")
        _ = QListWidgetItem(runs_icon, "Runs", list_widget)  # noqa


    def _create_data_table(self, all_columns, hide_groups):

        # Delete an already existing table
        if self.data_processor_table:
            self.data_processor_table.setParent(None)

        self.data_processor_table = jobtreeview.JobTreeView(all_columns, 
															self.cell(""), self)

        self.data_processor_table.setRootIsDecorated(False)

        row_entry = [''] * len(all_columns)
        self.add_row(row_entry)
        self._call_settings_listeners(lambda listener: listener.on_row_inserted(0, row_entry))

        self.table_signals = \
            jobtreeview.JobTreeViewSignalAdapter(self.data_processor_table, self)
        # The signal adapter subscribes to events from the table
        # and emits signals whenever it is notified.

        self.show_column_options(hide_groups, show=False)

        self._ui.data_processor_widget_layout.addWidget(self.data_processor_table)
        self.table_signals.cellTextChanged.connect(self._data_changed)
        self.table_signals.rowInserted.connect(self._row_inserted)
        self.table_signals.removeRowsRequested.connect(self._remove_rows_requested)
        self.table_signals.copyRowsRequested.connect(self._copy_rows_requested)
        self.table_signals.pasteRowsRequested.connect(self._paste_rows_requested)

    def show_column_options(self, key_list, show=True):
        for key in key_list:
            try:
                for col in self._column_options[key]:
                    if show:
                        self.data_processor_table.showColumn(col)
                    else:
                        self.data_processor_table.hideColumn(col)
            except KeyError:
                pass

    def hide_column_options(self):
        hide_list = list(self._column_options.keys())
        self.show_column_options(hide_list, show=False)

    def cell(self, text):
        background_color = 'white'
        border_thickness = 1
        border_color = "black"
        border_opacity = 255
        is_editable = True
        return jobtreeview.Cell(text, background_color, border_thickness,
                                border_color, border_opacity, is_editable)

    def row(self, path):
        return jobtreeview.RowLocation(path)

    def _setup_main_tab(self):
        self._ui.user_file_button.clicked.connect(self._on_user_file_load)
        self._ui.batch_file_button.clicked.connect(self._on_batch_file_load)

        # Disable the line edit fields. The user should not edit the paths manually.
        # They have to use the button.
        self._ui.user_file_line_edit.setDisabled(True)
        self._ui.batch_line_edit.setDisabled(True)

    def _process_selected_clicked(self):
        """
        Process runs
        """
        self._call_settings_listeners(lambda listener: listener.on_process_selected_clicked())

    def _process_all_clicked(self):
        """
        Process All button clicked
        """
        self._call_settings_listeners(lambda listener: listener.on_process_all_clicked())

    def _cancel_processing_clicked(self):
        """
        Cancel processing of remaining states clicked
        """
        self._call_settings_listeners(lambda listener: listener.on_cancel_processing_clicked())

    def _save_directory_clicked(self):
        self._call_settings_listeners(lambda listener: listener.on_save_directory_clicked())

    def _export_table_clicked(self):
        self._call_settings_listeners(lambda listener: listener.on_export_table_clicked())

    def _processing_finished(self):
        """
        Clean up
        """
        self._call_settings_listeners(lambda listener: listener.on_processing_finished())

    def _data_changed(self, row_location, column, old_value, new_value):
        row = row_location.rowRelativeToParent()
        self._call_settings_listeners(
            lambda listener: listener.on_data_changed(row, column, str(new_value), (old_value)))

    def _row_inserted(self, row_location):
        if row_location.depth() > 1:
            self.data_processor_table.removeRowAt(row_location)
        else:
            index = row_location.rowRelativeToParent()
            row = self.get_row(row_location)
            self._call_settings_listeners(lambda listener: listener.on_row_inserted(index, row))

    def _remove_rows_requested(self, rows):
        rows = [item.rowRelativeToParent() for item in rows]
        self._call_settings_listeners(lambda listener: listener.on_rows_removed(rows))

    def _remove_rows_requested_from_button(self):
        rows = self.get_selected_rows()
        self._call_settings_listeners(lambda listener: listener.on_rows_removed(rows))

    def _copy_rows_requested(self):
        self._call_settings_listeners(lambda listener: listener.on_copy_rows_requested())

    def _erase_rows(self):
        self._call_settings_listeners(lambda listener: listener.on_erase_rows())

    def _cut_rows(self):
        self._call_settings_listeners(lambda listener: listener.on_cut_rows())

    def _paste_rows_requested(self):
        self._call_settings_listeners(lambda listener: listener.on_paste_rows_requested())

    def _on_insert_button_pressed(self):
        self._call_settings_listeners(lambda listener: listener.on_insert_row())

    def _on_help_button_clicked(self):
        if PYQT4:
            proxies.showCustomInterfaceHelp('ANSTO Bilby')

    def _on_user_file_load(self):
        """
        Load the user file
        """
        # Load the user file
        load_file(self._ui.user_file_line_edit, "*.*", self.__generic_settings, self.__path_key,
                  self.get_user_file_path)

        # Set full user file path for default loading
        set_setting(self.__generic_settings, self.__user_file_key, self.get_user_file_path())

        # Notify presenters
        self._call_settings_listeners(lambda listener: listener.on_user_file_load())

    def on_user_file_load_failure(self):
        set_setting(self.__generic_settings, self.__user_file_key, "")
        self._ui.user_file_line_edit.setText("")

    def set_out_default_user_file(self):
        """
        Load a default user file, called on view set-up
        """
        load_default_file(self._ui.user_file_line_edit, self.__generic_settings, self.__user_file_key)

        if self.get_user_file_path() != "":
            self._call_settings_listeners(lambda listener: listener.on_user_file_load())

    def _on_batch_file_load(self):
        """
        Load the batch file
        """
        load_file(self._ui.batch_line_edit, "*.*", self.__generic_settings, self.__path_key,
                  self.get_batch_file_path)
        self._call_settings_listeners(lambda listener: listener.on_batch_file_load())

    def enable_buttons(self, enabled):
        self._ui.process_selected_button.setEnabled(enabled)
        self._ui.process_all_button.setEnabled(enabled)
        self._ui.cancel_processing_button.setEnabled(not enabled)
        self._ui.user_file_button.setEnabled(enabled)
        self._ui.batch_file_button.setEnabled(enabled)
        self._ui.save_directory_button.setEnabled(enabled)
        self._ui.export_table_button.setEnabled(enabled)
        self._ui.default_directory_checkbox.setEnabled(enabled)

    def enable_process_buttons(self, enabled):
        self._ui.process_selected_button.setEnabled(enabled)
        self._ui.process_all_button.setEnabled(enabled)
        self._ui.cancel_processing_button.setEnabled(not enabled)

    def enable_cancel_button(self, enabled):
        self._ui.cancel_processing_button.setEnabled(enabled)

    def enable_save_directory_button(self, enabled):
        self._ui.save_directory_button.setEnabled(enabled)

    def display_message_box(self, title, message, details):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)

        message_length = len(message)

        # This is to ensure that the QMessage box if wide enough to display nicely.
        msg.setText(10 * ' ' + message + ' ' * (30 - message_length))
        msg.setWindowTitle(title)
        msg.setDetailedText(details)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.setDefaultButton(QMessageBox.Ok)
        msg.setEscapeButton(QMessageBox.Ok)
        msg.exec_()

    def display_save_file_box(self, title, default_path, file_filter):
        filename = QFileDialog.getSaveFileName(self, title, default_path, filter=file_filter)
        return filename

    def get_user_file_path(self):
        return str(self._ui.user_file_line_edit.text())

    def get_batch_file_path(self):
        return str(self._ui.batch_line_edit.text())

    def set_out_file_directory(self, out_file_directory):
        self._ui.output_directory_location.setText("{}".format(out_file_directory))

    def _default_directory_clicked(self):
        self._call_settings_listeners(
            lambda listener: listener.on_default_directory_clicked(self.use_default_directory))

    def display_save_directory_box(self, title, default_path):
        directory = QFileDialog.getExistingDirectory(self, title, default_path, 
                                                     QFileDialog.ShowDirsOnly)
        return directory

    # ------------------------------------------------------------------------------------------------------------------
    # Elements which can be set and read by the model
    # ------------------------------------------------------------------------------------------------------------------
    def update_gui_combo_box(self, value, expected_type, combo_box):
        # There are three types of values that can be passed:
        # Lists: we set the combo box to the values in the list
        # expected_type: we set the expected type
        # str (in the case of "Variable" Q rebin): We set the combo box to the text if it is an option
        gui_element = getattr(self._ui, combo_box)
        if isinstance(value, list):
            gui_element.clear()
            for element in value:
                self._add_list_element_to_combo_box(gui_element=gui_element, element=element,
                                                    expected_type=expected_type)
        elif expected_type.has_member(value):
            self._set_enum_as_element_in_combo_box(gui_element=gui_element, element=value,
                                                   expected_type=expected_type)
        elif isinstance(value, str):
            index = gui_element.findText(value)
            if index != -1:
                gui_element.setCurrentIndex(index)
        else:
            raise RuntimeError("Expected an input of type {}, but got {}".format(expected_type, type(value)))

    def _add_list_element_to_combo_box(self, gui_element, element, expected_type=None):
        if expected_type is not None and isclass(element) and issubclass(element, expected_type):
            self._add_enum_as_element_in_combo_box(gui_element=gui_element, element=element,
                                                   expected_type=expected_type)
        else:
            gui_element.addItem(element)

    @staticmethod
    def _set_enum_as_element_in_combo_box(gui_element, element, expected_type):
        value_as_string = expected_type.to_string(element)
        index = gui_element.findText(value_as_string)
        if index != -1:
            gui_element.setCurrentIndex(index)

    def _add_enum_as_element_in_combo_box(self, gui_element, element, expected_type):
        value_as_string = expected_type.to_string(element)
        gui_element.addItem(value_as_string)

    def get_simple_line_edit_field(self, expected_type, line_edit):
        gui_element = getattr(self._ui, line_edit)
        value_as_string = gui_element.text()
        return expected_type(value_as_string) if value_as_string else None

    def update_simple_line_edit_field(self, line_edit, value):
        if value:
            gui_element = getattr(self._ui, line_edit)
            gui_element.setText(str(value))

    @staticmethod
    def _set_limited_spinner_value(spin_box, value):
        minimum = spin_box.minimum()
        maximum = spin_box.maximum()
        if value < minimum or value > maximum:
            raise ValueError("The value for the spinner box {} has "
                             "to be in the range of {} and {}".format(value, minimum, maximum))
        spin_box.setValue(value)

    # ==================================================================================================================
    # START PROPERTIES
    # ==================================================================================================================

    # -------------------------------------------------------
    # FRONT TAB
    # -------------------------------------------------------
    @property
    def progress_bar_minimum(self):
        return self._ui.batch_progress_bar.minimum()

    @progress_bar_minimum.setter
    def progress_bar_minimum(self, value):
        self._ui.batch_progress_bar.setMinimum(value)

    @property
    def progress_bar_maximum(self):
        return self._ui.batch_progress_bar.maximum()

    @progress_bar_maximum.setter
    def progress_bar_maximum(self, value):
        self._ui.batch_progress_bar.setMaximum(value)

    @property
    def progress_bar_value(self):
        return self._ui.batch_progress_bar.value()

    @progress_bar_value.setter
    def progress_bar_value(self, progress):
        self._ui.batch_progress_bar.setValue(progress)

    @property
    def output_folder(self):
        return self.get_simple_line_edit_field(line_edit="output_folder_edit",
                                               expected_type=str)

    @output_folder.setter
    def output_folder(self, value):
        self.update_simple_line_edit_field(line_edit="output_folder_edit",
                                           value=value)

    @property
    def save_directory(self):
        return str(self._ui.output_directory_location.text())

    # -----------------------------------------------------------------
    # Global options
    # -----------------------------------------------------------------
    
    @property
    def use_default_directory(self):
        return self._ui.default_directory_checkbox.isChecked()

    @use_default_directory.setter
    def use_default_directory(self, value):
        self._ui.default_directory_checkbox.setChecked(value)


    # ==================================================================================================================
    # END PROPERTIES
    # ==================================================================================================================

    def _attach_validators(self):
        # Setup the list of validators
        double_validator = QDoubleValidator()
        positive_double_validator = QDoubleValidator()
        positive_double_validator.setBottom(0.0)
        positive_integer_validator = QIntValidator()
        positive_integer_validator.setBottom(1)

        # -------------------------------
        # General tab
        # -------------------------------

    def reset_all_fields_to_default(self):
        # ------------------------------
        # General tab
        # ------------------------------
        self.use_default_directory = True

    # ----------------------------------------------------------------------------------------------
    # Table interaction
    # ----------------------------------------------------------------------------------------------

    def get_cell(self, row, column, convert_to=None):
        row_location = self.row([row])
        value = self.data_processor_table.cellAt(row_location, column).contentText()
        return value if convert_to is None else convert_to(value)

    def set_cell(self, value, row, column):
        value_as_str = str(value)
        cell = self.data_processor_table.cellAt(row, column)
        cell.setContentText(value_as_str)
        self.data_processor_table.setCellAt(row, column, cell)

    def change_row_color(self, color, row):
        row_location = self.row([row])
        cell_data = self.data_processor_table.cellsAt(row_location)
        for index, cell in enumerate(cell_data):
            cell.setBackgroundColor(color)
            self.data_processor_table.setCellAt(row_location, index, cell)

    def set_row_tooltip(self, tool_tip, row):
        row_location = self.row([row])
        cell_data = self.data_processor_table.cellsAt(row_location)
        for index, cell in enumerate(cell_data):
            cell.setToolTip(tool_tip)
            self.data_processor_table.setCellAt(row_location, index, cell)

    def get_selected_rows(self):
        row_locations = self.data_processor_table.selectedRowLocations()
        rows = [x.rowRelativeToParent() for x in row_locations]
        return rows

    def get_row(self, row_location):
        cell_data = self.data_processor_table.cellsAt(row_location)
        return [str(x.contentText()) for x in cell_data]

    def clear_table(self):
        self.data_processor_table.removeAllRows()

    def clear_selection(self):
        self.data_processor_table.clearSelection()

    def update_table_selection(self, row_locations):
        row_locations = [self.row([x]) for x in row_locations]
        self.data_processor_table.setSelectedRowLocations(row_locations)

    def add_row(self, value):
        value = [self.cell(x) for x in value]
        self.data_processor_table.appendChildRowOf(self.row([]), value)

    def remove_rows(self, rows):
        rows = [self.row([item]) for item in rows]
        self.data_processor_table.removeRows(rows)

    def insert_empty_row(self, row_index):
        self.data_processor_table.insertChildRowOf(self.row([]), row_index)

    def set_hinting_line_edit_for_column(self, column, hint_strategy):
        self.data_processor_table.setHintsForColumn(column, hint_strategy)

    def _run_python_code(self, text):
        """
        Re-emits 'runPytonScript' signal
        """
        execute_script(text)

    def closeEvent(self, event):
        for child in self.children():
            if isinstance(child, SANSSaveOtherDialog):
                child.done(0)
        super(QMainWindow, self).closeEvent(event)
