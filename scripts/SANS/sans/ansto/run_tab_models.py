# -*- coding: utf-8 -*-
# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
Primary factory for creating the table models, reduction code, etc.
"""
from __future__ import (absolute_import, division, print_function)

from sans.sans_batch import SANSBatchReduction
from sans.gui_logic.models.state_gui_model import StateGuiModel
from sans.ansto.table_model import TableModel
from sans.ansto.bilby import bilby_state_model
from sans.ansto.bilby import table_model as bilby_models
from sans.user_file import user_file_reader
from sans.ansto import file_readers

import imp

imp.reload(bilby_models)
imp.reload(file_readers)

class RunTabModels(object):
	
	def __init__(self, tableModel=TableModel, 
			  tableRowModel=bilby_models.RowModel, 
			  stateModel=bilby_state_model.BilbyStateGuiModel, batchReduction=SANSBatchReduction,
			  userFileReader=file_readers.UserFileReader, 
			  batchFileReader=file_readers.BatchFileReader):
		self._table_model = tableModel
		self._table_row_model = tableRowModel
		self._state_model = stateModel
		self._batch_reduction = batchReduction
		self._user_file_reader = userFileReader
		self._batch_file_reader = batchFileReader

	@property
	def TableModel(self):
		return self._table_model

	@property
	def TableRowModel(self):
		return self._table_row_model

	@property
	def StateModel(self):
		return self._state_model

	@property
	def BatchReduction(self):
		return self._batch_reduction

	@property
	def UserFileReader(self):
		return self._user_file_reader

	@property
	def BatchFileReader(self):
		return self._batch_file_reader