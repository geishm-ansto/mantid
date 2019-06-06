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

class RunTabModels(object):
	
	def __init__(self, tableModel, tableRowModel, batchReduction,
              userFileReader, batchFileReader):
		self._table_model = tableModel
		self._table_row_model = tableRowModel
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
	def BatchReduction(self):
		return self._batch_reduction

	@property
	def UserFileReader(self):
		return self._user_file_reader

	@property
	def BatchFileReader(self):
		return self._batch_file_reader