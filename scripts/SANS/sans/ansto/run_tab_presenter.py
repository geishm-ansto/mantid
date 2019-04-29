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

from mantid.api import (FileFinder)
from mantid.kernel import Logger, ConfigService

from sans.common.constants import ALL_PERIODS
from sans.common.enums import (BatchReductionEntry, RangeStepType, SampleShape, FitType, RowState, SANSInstrument)
from sans.gui_logic.gui_common import (get_reduction_mode_strings_for_gui, get_string_for_gui_from_instrument,
                                       add_dir_to_datasearch, remove_dir_from_datasearch)
from sans.gui_logic.models.batch_process_runner import BatchProcessRunner
from sans.gui_logic.models.create_state import create_states
from sans.gui_logic.models.diagnostics_page_model import create_state
from sans.gui_logic.presenter.add_runs_presenter import OutputDirectoryObserver as SaveDirectoryObserver
from sans.user_file.user_file_reader import UserFileReader

from ui.sans_isis import SANSSaveOtherWindow
from ui.sans_isis.work_handler import WorkHandler

from ui.ansto.run_tab_gui import RunTabGui

from qtpy import PYQT4
IN_MANTIDPLOT = False
if PYQT4:
    try:
        from mantidplot import graph, newGraph
        IN_MANTIDPLOT = True
    except ImportError:
        pass
else:
    from mantidqt.plotting.functions import get_plot_fig

row_state_to_colour_mapping = {RowState.Unprocessed: '#FFFFFF', RowState.Processed: '#d0f4d0',
                               RowState.Error: '#accbff'}


def log_times(func):
    """
    Generic decorator to time the execution of the function and
    print it to the logger.
    """

    def run(*args, **kwargs):
        t0 = time.time()
        result = func(*args, **kwargs)
        t1 = time.time()
        time_taken = t1 - t0
        # args[0] is the self parameter
        args[0].sans_logger.information("The generation of all states took {}s".format(time_taken))
        return result

    return run


class RunTabPresenter(object):
    class ConcreteRunTabListener(RunTabGui.RunTabListener):
        def __init__(self, presenter):
            super(RunTabPresenter.ConcreteRunTabListener, self).__init__()
            self._presenter = presenter

        def on_user_file_load(self):
            self._presenter.on_user_file_load()

        def on_batch_file_load(self):
            self._presenter.on_batch_file_load()

        def on_process_selected_clicked(self):
            self._presenter.on_process_selected_clicked()

        def on_process_all_clicked(self):
            self._presenter.on_process_all_clicked()

        def on_load_clicked(self):
            self._presenter.on_load_clicked()

        def on_export_table_clicked(self):
            self._presenter.on_export_table_clicked()

        def on_data_changed(self, row, column, new_value, old_value):
            self._presenter.on_data_changed(row, column, new_value, old_value)

        def on_row_inserted(self, index, row):
            self._presenter.on_row_inserted(index, row)

        def on_rows_removed(self, rows):
            self._presenter.on_rows_removed(rows)

        def on_copy_rows_requested(self):
            self._presenter.on_copy_rows_requested()

        def on_paste_rows_requested(self):
            self._presenter.on_paste_rows_requested()

        def on_insert_row(self):
            self._presenter.on_insert_row()

        def on_erase_rows(self):
            self._presenter.on_erase_rows()

        def on_cut_rows(self):
            self._presenter.on_cut_rows_requested()

    class ProcessListener(WorkHandler.WorkListener):
        def __init__(self, presenter):
            super(RunTabPresenter.ProcessListener, self).__init__()
            self._presenter = presenter

        def on_processing_finished(self, result):
            self._presenter.on_processing_finished(result)

        def on_processing_error(self, error):
            self._presenter.on_processing_error(error)

    def __init__(self, facility, view, models):
        super(RunTabPresenter, self).__init__()
        self._facility = facility
        self._models = models
        # Logger
        self.sans_logger = Logger("SANS")
        # Name of graph to output to
        self.output_graph = 'SANS-Latest'
        # For matplotlib continuous plotting
        self.output_fig = None
        self.progress = 0

        # Models that are being used by the presenter
        self._state_model = self._models.StateModel({})
        self._table_model = self._models.TableModel(self._models.TableRowModel)
        self._table_model.subscribe_to_model_changes(self)

        # Presenter needs to have a handle on the view since it delegates it
        self._view = None
        self.set_view(view)
        self._processing = False
        self.work_handler = WorkHandler()
        self.batch_process_runner = BatchProcessRunner(self._models.BatchReduction(), 
													   self.notify_progress,
                                                       self.on_processing_finished,
                                                       self.on_processing_error)

        # File information for the first input
        self._file_information = None
        self._clipboard = []

        # Check save dir for display
        self._save_directory_observer = \
            SaveDirectoryObserver(self._handle_output_directory_changed)

    def _default_gui_setup(self):
        """
        Provides a default setup of the GUI. This is important for the initial start up, when the view is being set.
        """
        pass

    def _handle_output_directory_changed(self, new_directory):
        """
        Update the gui to display the new save location for workspaces
        :param new_directory: string. Current save directory for files
        :return:
        """
        self._view.set_out_file_directory(new_directory)

    # ------------------------------------------------------------------------------------------------------------------
    # Table + Actions
    # ------------------------------------------------------------------------------------------------------------------
    def set_view(self, view):
        """
        Sets the view
        :param view: the view is a derivative of the RunTabGui. The presenter needs to access some of the API
        """
        if view is not None:
            self._view = view

            # Add a listener to the view
            self._add_listener()

            # Default gui setup
            self._default_gui_setup()
            self._view.disable_process_buttons()

            all_columns = self._table_model.column_labels()
            col_groups = self._table_model.column_options()
            self._view.setup_layout(all_columns, col_groups)

            self._view.set_out_file_directory(ConfigService.Instance().getString("defaultsave.directory"))

            self._view.set_out_default_user_file()

    def _add_listener(self):

        listener = RunTabPresenter.ConcreteRunTabListener(self)
        self._view.add_listener(listener)

    def on_user_file_load(self):
        """
        Loads the user file. Populates the models and the view.
        """
        error_msg = "Loading of the user file failed"
        try:
            # 1. Get the user file path from the view
            user_file_path = self._view.get_user_file_path()

            if not user_file_path:
                return
            # 2. Get the full file path
            user_file_path = FileFinder.getFullPath(user_file_path)
            if not os.path.exists(user_file_path):
                raise RuntimeError(
                    "The user path {} does not exist. Make sure a valid user file path"
                    " has been specified.".format(user_file_path))
        except RuntimeError as path_error:
            # This exception block runs if user file does not exist
            self._view.on_user_file_load_failure()
            self.display_errors(path_error, error_msg + " when finding file.")
        else:
            try:
                self._table_model.user_file = user_file_path
                # Clear out the current view
                self._view.reset_all_fields_to_default()

                # 3. Read and parse the user file
                user_file_reader = self._models.UserFileReader(user_file_path)
                user_file_items = user_file_reader.read_user_file()
            except (RuntimeError, ValueError) as e:
                # It is in this exception block that loading fails if the file is invalid (e.g. a csv)
                self._view.on_user_file_load_failure()
                self.display_errors(e, error_msg + " when reading file.", use_error_name=True)
            else:
                try:
                    # 4. Populate the model
                    self._state_model = self._models.StateModel(user_file_items)
                    # 5. Update the views.
                    self._update_view_from_state_model()

                    # 6. Warning if user file did not contain a recognised instrument
                    #if self._view.instrument == SANSInstrument.NoInstrument:
                    #    raise RuntimeError("User file did not contain a SANS Instrument.")

                except RuntimeError as instrument_e:
                    # This exception block runs if the user file does not contain an parsable instrument
                    self._view.on_user_file_load_failure()
                    self.display_errors(instrument_e, error_msg + " when reading instrument.")
                except Exception as other_error:
                    # If we don't catch all exceptions, SANS can fail to open if last loaded
                    # user file contains an error that would not otherwise be caught
                    traceback.print_exc()
                    self._view.on_user_file_load_failure()
                    self.display_errors(other_error, "Unknown error in loading user file.", use_error_name=True)

    def on_batch_file_load(self):
        """
        Loads a batch file and populates the batch table based on that.
        """
        try:
            # 1. Get the batch file from the view
            batch_file_path = self._view.get_batch_file_path()

            if not batch_file_path:
                return

            datasearch_dirs = ConfigService["datasearch.directories"]
            batch_file_directory, datasearch_dirs = add_dir_to_datasearch(batch_file_path, datasearch_dirs)
            ConfigService["datasearch.directories"] = datasearch_dirs

            if not os.path.exists(batch_file_path):
                raise RuntimeError(
                    "The batch file path {} does not exist. Make sure a valid batch file path"
                    " has been specified.".format(batch_file_path))

            self._table_model.batch_file = batch_file_path

            # 2. Read the batch file
            batch_file_parser = self._models.BatchFileReader(batch_file_path)
            parsed_rows = batch_file_parser.parse_batch_file()

            # 3. Populate the table
            self._table_model.clear_table_entries()
            for index, row in enumerate(parsed_rows):
                self._add_row_to_table_model(row, index)
            self._table_model.remove_table_entries([len(parsed_rows)])
        except RuntimeError as e:
            if batch_file_directory:
                # Remove added directory from datasearch.directories
                ConfigService["datasearch.directories"] = remove_dir_from_datasearch(batch_file_directory, datasearch_dirs)

            self.sans_logger.error("Loading of the batch file failed. {}".format(str(e)))
            self.display_warning_box('Warning', 'Loading of the batch file failed', str(e))

    def _add_row_to_table_model(self, row, index):
        """
        Adds a row to the table
        """

        def get_string_entry(_tag, _row):
            _element = ""
            if _tag in _row:
                _element = _row[_tag]
            return _element

        # 1. Pull out the entries
        row_entry = [get_string_entry(tag, row) for tag in self._table_model.column_keys()]
        table_index_model = self._table_model.create_row(*row_entry)

        self._table_model.add_table_entry(index, table_index_model)

    def on_update_rows(self):
        self.update_view_from_table_model()

    def update_view_from_table_model(self):
        self._view.clear_table()
        self._view.hide_column_options()
        for row_index, row in enumerate(self._table_model._table_entries):
            row_entry = [str(x) for x in row.to_list()]
            self._view.add_row(row_entry)
            self._view.change_row_color(row_state_to_colour_mapping[row.row_state], row_index + 1)
            self._view.set_row_tooltip(row.tool_tip, row_index + 1)
            show_groups = row.includes_options()
            if show_groups:
                self._view.show_column_options(show_groups)
        self._view.remove_rows([0])
        self._view.clear_selection()

    def on_data_changed(self, row, column, new_value, old_value):
        self._table_model.update_table_entry(row, column, new_value)
        self._view.change_row_color(row_state_to_colour_mapping[RowState.Unprocessed], row)
        self._view.set_row_tooltip('', row)

    # ----------------------------------------------------------------------------------------------
    # Processing
    # ----------------------------------------------------------------------------------------------

    def _handle_get_states(self, rows):
        """
        Return the states for the supplied rows, calling on_processing_error for any errors
        which occur.
        """
        states, errors = self.get_states(row_index=rows)
        for row, error in errors.items():
            self.on_processing_error(row, error)
        return states

    def _plot_graph(self):
        """
        Plot a graph if continuous output specified.
        """
        if self._view.plot_results:
            if IN_MANTIDPLOT:
                if not graph(self.output_graph):
                    newGraph(self.output_graph)
            elif not PYQT4:
                ax_properties = {'yscale': 'log',
                                 'xscale': 'log'}
                fig, _ = get_plot_fig(ax_properties=ax_properties, window_title=self.output_graph)
                fig.show()
                self.output_fig = fig

    def _set_progress_bar_min_max(self, min, max):
        """
        The progress of the progress bar is given by min / max
        :param min: Current value of the progress bar.
        :param max: The value at which the progress bar is full
        """
        setattr(self._view, 'progress_bar_value', min)
        setattr(self._view, 'progress_bar_maximum', max)

    def _process_rows(self, rows):
        """
        Processes a list of rows. Any errors cause the row to be coloured red.
        """
        try:
            for row in rows:
                self._table_model.reset_row_state(row)
            self.update_view_from_table_model()

            self._view.disable_buttons()
            self._processing = True
            self.sans_logger.information("Starting processing of batch table.")

            states = self._handle_get_states(rows)
            if not states:
                raise Exception("No states found")

            self._plot_graph()
            self.progress = 0
            self._set_progress_bar_min_max(self.progress, len(states))
            save_can = self._view.save_can

            # MantidPlot and Workbench have different approaches to plotting
            output_graph = self.output_graph if PYQT4 else self.output_fig
            self.batch_process_runner.process_states(states,
                                                     self._view.use_optimizations,
                                                     self._view.output_mode,
                                                     self._view.plot_results,
                                                     output_graph,
                                                     save_can)

        except Exception as e:
            self.on_processing_finished(None)
            self.sans_logger.error("Process halted due to: {}".format(str(e)))
            self.display_warning_box('Warning', 'Process halted', str(e))

    def on_process_all_clicked(self):
        """
        Process all entries in the table, regardless of selection.
        """
        all_rows = range(self._table_model.get_number_of_rows())
        all_rows = self._table_model.get_non_empty_rows(all_rows)
        if all_rows:
            self._process_rows(all_rows)

    def on_process_selected_clicked(self):
        """
        Process selected table entries.
        """
        selected_rows = self._view.get_selected_rows()
        selected_rows = self._table_model.get_non_empty_rows(selected_rows)
        if selected_rows:
            self._process_rows(selected_rows)

    def on_processing_error(self, row, error_msg):
        """
        An error occurs while processing the row with index row, error_msg is displayed as a
        tooltip on the row.
        """
        self.increment_progress()
        self._table_model.set_row_to_error(row, error_msg)
        self.update_view_from_table_model()

    def on_processing_finished(self, result):
        self._view.enable_buttons()
        self._processing = False

    def on_load_clicked(self):
        try:
            self._view.disable_buttons()
            self._processing = True
            self.sans_logger.information("Starting load of batch table.")

            selected_rows = self._get_selected_rows()
            selected_rows = self._table_model.get_non_empty_rows(selected_rows)
            states, errors = self.get_states(row_index=selected_rows)

            for row, error in errors.items():
                self.on_processing_error(row, error)

            if not states:
                self.on_processing_finished(None)
                return

            self.progress = 0
            setattr(self._view, 'progress_bar_value', self.progress)
            setattr(self._view, 'progress_bar_maximum', len(states))
            self.batch_process_runner.load_workspaces(states)
        except Exception as e:
            self._view.enable_buttons()
            self.sans_logger.error("Process halted due to: {}".format(str(e)))
            self.display_warning_box("Warning", "Process halted", str(e))

    def on_export_table_clicked(self):
        non_empty_rows = self.get_row_indices()
        if len(non_empty_rows) == 0:
            self.sans_logger.notice("Cannot export table as it is empty.")
            return

        # Python 2 and 3 take input in different modes for writing lists to csv files
        if sys.version_info[0] == 2:
            open_type = 'wb'
        else:
            open_type = 'w'

        try:
            self._view.disable_buttons()

            default_filename = self._table_model.batch_file
            filename = self.display_save_file_box("Save table as", default_filename, "*.csv")

            if filename:
                self.sans_logger.notice("Starting export of table.")
                if filename[-4:] != '.csv':
                    filename += '.csv'

                with open(filename, open_type) as outfile:
                    # Pass filewriting object rather than filename to make testing easier
                    writer = csv.writer(outfile)
                    self._export_table(writer, non_empty_rows)
                    self.sans_logger.notice("Table exporting finished.")

            self._view.enable_buttons()
        except Exception as e:
            self._view.enable_buttons()
            self.sans_logger.error("Export halted due to : {}".format(str(e)))
            self.display_warning_box("Warning", "Export halted", str(e))

    def display_errors(self, error, context_msg, use_error_name=False):
        """
        Code for alerting the user to a caught error
        :param error: a caught exception
        :param context_msg: string. Text to explain what SANS was trying to do
                            when the error occurred. e.g. 'Loading of the user file failed'.
        :param use_error_name: bool. If True, append type of error (e.g. RuntimeError) to context_msg
        :return:
        """
        logger_msg = context_msg
        if use_error_name:
            logger_msg += " {}:".format(type(error).__name__)
        logger_msg += " {}"
        self.sans_logger.error(logger_msg.format(str(error)))
        self.display_warning_box('Warning', context_msg, str(error))

    def display_warning_box(self, title, text, detailed_text):
        self._view.display_message_box(title, text, detailed_text)

    def display_save_file_box(self, title, default_path, file_filter):
        filename = self._view.display_save_file_box(title, default_path, file_filter)
        return filename

    def notify_progress(self, row, out_shift_factors, out_scale_factors):
        self.increment_progress()
        if out_scale_factors and out_shift_factors:
            self._table_model.set_option(row, 'MergeScale', round(out_scale_factors[0], 3))
            self._table_model.set_option(row, 'MergeShift', round(out_shift_factors[0], 3))

        self._table_model.set_row_to_processed(row, '')

    def increment_progress(self):
        self.progress = self.progress + 1
        setattr(self._view, 'progress_bar_value', self.progress)

    # ----------------------------------------------------------------------------------------------
    # Row manipulation
    # ----------------------------------------------------------------------------------------------

    def num_rows(self):
        return self._table_model.get_number_of_rows()

    def on_row_inserted(self, index, row):
        """
        Insert a row at a selected point
        """
        row_table_index = self._table_model.create_row(*row)
        self._table_model.add_table_entry(index, row_table_index)

    def on_insert_row(self):
        """
        Add an empty row to the table after the first selected row (or at the end of the table
        if nothing is selected).
        """
        selected_rows = self._view.get_selected_rows()

        selected_row = selected_rows[0] + 1 if selected_rows else self.num_rows()
        empty_row = self._table_model.create_empty_row()
        self._table_model.add_table_entry(selected_row, empty_row)

    def on_erase_rows(self):
        """
        Make all selected rows empty.
        """
        selected_rows = self._view.get_selected_rows()
        empty_row = self._table_model.create_empty_row()
        for row in selected_rows:
            empty_row = self._table_model.create_empty_row()
            self._table_model.replace_table_entries([row], [empty_row])

    def on_rows_removed(self, rows):
        """
        Remove rows from the table
        """
        self._table_model.remove_table_entries(rows)

    def on_copy_rows_requested(self):
        selected_rows = self._view.get_selected_rows()
        self._clipboard = []
        for row in selected_rows:
            data_from_table_model = self._table_model.get_table_entry(row).to_list()
            self._clipboard.append(data_from_table_model)

    def on_cut_rows_requested(self):
        self.on_copy_rows_requested()
        rows = self._view.get_selected_rows()
        self.on_rows_removed(rows)

    def on_paste_rows_requested(self):
        if self._clipboard:
            selected_rows = self._view.get_selected_rows()
            selected_rows = selected_rows if selected_rows else [self.num_rows()]
            replacement_table_index_models = [self._table_model.create_row(*x) for x in self._clipboard]
            self._table_model.replace_table_entries(selected_rows, replacement_table_index_models)

    def on_manage_directories(self):
        self._view.show_directory_manager()

    def on_sample_geometry_view_changed(self, show_geometry):
        if show_geometry:
            self._view.show_geometry()
        else:
            self._view.hide_geometry()

    def on_compatibility_unchecked(self):
        self.display_warning_box('Warning', 'Are you sure you want to uncheck compatibility mode?',
                                 'Non-compatibility mode has known issues. DO NOT USE if applying bin masking'
                                 ' to event workspaces.')

    def get_row_indices(self):
        """
        Gets the indices of row which are not empty.
        :return: a list of row indices.
        """
        row_indices_which_are_not_empty = []
        number_of_rows = self._table_model.get_number_of_rows()
        for row in range(number_of_rows):
            if not self.is_empty_row(row):
                row_indices_which_are_not_empty.append(row)
        return row_indices_which_are_not_empty

    def on_mask_file_add(self):
        """
        We get the added mask file name and add it to the list of masks
        """
        new_mask_file = self._view.get_mask_file()
        if not new_mask_file:
            return
        new_mask_file_full_path = FileFinder.getFullPath(new_mask_file)
        if not new_mask_file_full_path:
            return

        # Add the new mask file to state model
        mask_files = self._state_model.mask_files

        mask_files.append(new_mask_file)
        self._state_model.mask_files = mask_files

        # Make sure that the sub-presenters are up to date with this change
        self._masking_table_presenter.on_update_rows()
        self._settings_diagnostic_tab_presenter.on_update_rows()
        self._beam_centre_presenter.on_update_rows()

    def is_empty_row(self, row):
        """
        Checks if a row has no entries. These rows will be ignored.
        :param row: the row index
        :return: True if the row is empty.
        """
        return self._table_model.is_empty_row(row)

    # def _validate_rows(self):
    #     """
    #     Validation of the rows. A minimal setup requires that ScatterSample is set.
    #     """
    #     # If SampleScatter is empty, then don't run the reduction.
    #     # We allow empty rows for now, since we cannot remove them from Python.
    #     number_of_rows = self._table_model.get_number_of_rows()
    #     for row in range(number_of_rows):
    #         if not self.is_empty_row(row):
    #             sample_scatter = self._view.get_cell(row, 0)
    #             if not sample_scatter:
    #                 raise RuntimeError("Row {} has not SampleScatter specified. Please correct this.".format(row))

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

    # ----------------------------------------------------------------------------------------------
    # Table Model and state population
    # ------------------------------------------------------------------------------------------------------------------
    def _get_selected_rows(self):
        selected_rows = self._view.get_selected_rows()
        selected_rows = selected_rows if selected_rows else range(self._table_model.get_number_of_rows())
        for row in selected_rows:
            self._table_model.reset_row_state(row)
        self.update_view_from_table_model()

        return selected_rows

    @log_times
    def get_states(self, row_index=None, file_lookup=True):
        """
        Gathers the state information for all rows.
        :param row_index: if a single row is selected, then only this row is returned,
                          else all the state for all rows is returned.
        :return: a list of states.
        """
        # 1. Update the state model
        state_model_with_view_update = self._get_state_model_with_view_update()
        # 2. Update the table model
        table_model = self._table_model
        # 3. Go through each row and construct a state object
        states, errors = None, None
        if table_model and state_model_with_view_update:
            states, errors = create_states(state_model_with_view_update, table_model,
                                           self._view.instrument,
                                           self._facility,
                                           row_index=row_index,
                                           file_lookup=file_lookup)

        if errors:
            self.sans_logger.warning("Errors in getting states...")
            for _, v in errors.items():
                self.sans_logger.warning("{}".format(v))

        return states, errors

    def get_state_for_row(self, row_index, file_lookup=True):
        """
        Creates the state for a particular row.
        :param row_index: the row index
        :return: a state if the index is valid and there is a state else None
        """
        states, errors = self.get_states(row_index=[row_index], file_lookup=file_lookup)
        if states is None:
            self.sans_logger.warning(
                "There does not seem to be data for a row {}.".format(row_index))
            return None

        if row_index in list(states.keys()):
            if states:
                return states[row_index]
        return None

    def _update_view_from_state_model(self):
        pass

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
        pass

    def _set_on_state_model(self, attribute_name, state_model):
        attribute = getattr(self._view, attribute_name)
        if attribute is not None and attribute != '':
            setattr(state_model, attribute_name, attribute)

    def get_cell_value(self, row, column):
        return self._view.get_cell(row=row, column=self.table_index[column], convert_to=str)

    def _export_table(self, filewriter, rows):
        """
        Take the current table model, and create a comma delimited csv file
        :param filewriter: File object to be written to
        :param rows: list of indices for non-empty rows
        :return: Nothing
        """
        filewriter.writerow(self._table_model.column_keys())
        for row in rows:
                table_row = self._table_model.get_table_entry(row).to_list()
                #batch_file_row = self._create_batch_entry_from_row(table_row)
                filewriter.writerow(table_row)

    @staticmethod
    def _create_batch_entry_from_row(row):
        batch_file_keywords = ["sample_sans",
                               "output_as",
                               "sample_trans",
                               "sample_direct_beam",
                               "can_sans",
                               "can_trans",
                               "can_direct_beam",
                               "user_file"]

        loop_range = min(len(row), len(batch_file_keywords))
        new_row = [''] * (2 * loop_range)

        for i in range(loop_range):
            key = batch_file_keywords[i]
            value = row[i]
            new_row[2*i] = key
            new_row[2*i + 1] = value

        return new_row

