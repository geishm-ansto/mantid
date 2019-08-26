# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#pylint: disable=invalid-name
"""
    Script used to start the Test Interface from MantidPlot
"""
from sans.common.enums import SANSFacility
from sans.ansto import run_tab_models, run_tab_presenter, table_model
from sans.ansto.bilby import (bilby_presenter, bilby_state_data,
                              bilby_batch_reduction, file_information)
from sans.ansto.bilby import table_model as bilby_table_model
from sans.ansto.bilby import file_readers as bilby_file_readers
from ui.ansto import run_tab_gui, ansto_bilby_gui
from sans.ansto.table_model import TableModel
from sans.ansto import file_readers

import imp

# reload modules during development to avoid relaunching the append
# TODO remove before mantid integration

imp.reload(file_information)
imp.reload(bilby_state_data)
imp.reload(file_readers)
imp.reload(bilby_file_readers)
imp.reload(bilby_batch_reduction)
imp.reload(table_model)
imp.reload(bilby_table_model)
imp.reload(run_tab_gui)
imp.reload(ansto_bilby_gui)
imp.reload(run_tab_presenter)
imp.reload(bilby_presenter)
imp.reload(run_tab_models)

#--------------------------------------------------
# Create the models needed
#--------------------------------------------------
models = run_tab_models.RunTabModels(tableModel=bilby_table_model.TableModel,
                                     tableRowModel=bilby_table_model.RowModel,
                                     batchReduction=bilby_batch_reduction.BilbyBatchReduction,
                                     userFileReader=bilby_file_readers.BilbyUserFileReader,
                                     batchFileReader=bilby_file_readers.BilbyBatchFileReader)

# -------------------------------------------------
# Create view
# ------------------------------------------------
ui = ansto_bilby_gui.BilbyBatchReductionGui(ansto_bilby_gui.Ui_AnstoBilbyWindow)

# -----------------------------------------------
# Create presenter
# -----------------------------------------------
presenter = bilby_presenter.BilbyPresenter(SANSFacility.ANSTO, ui, models)

# Show
ui.show()
