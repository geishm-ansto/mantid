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
from sans.ansto import run_tab_models, run_tab_presenter, table_model, state_model
from sans.ansto.bilby import bilby_presenter, bilby_state_model
from ui.ansto import run_tab_gui, ansto_bilby_gui

import imp

imp.reload(run_tab_gui)
imp.reload(table_model)
imp.reload(ansto_bilby_gui)
imp.reload(run_tab_presenter)
imp.reload(bilby_presenter)
imp.reload(state_model)
imp.reload(bilby_state_model)
imp.reload(run_tab_models)


#--------------------------------------------------
# Create the models needed
#--------------------------------------------------
models = run_tab_models.RunTabModels()

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
