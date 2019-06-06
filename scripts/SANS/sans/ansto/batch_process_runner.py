# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
from qtpy.QtCore import Slot, QThreadPool, Signal, QObject
from ui.sans_isis.worker import Worker

class BatchProcessRunner(QObject):
    row_processed_signal = Signal(int, list)
    row_failed_signal = Signal(int, str)

    def __init__(self, batch_reduction, notify_progress, notify_done, notify_error):
        super(BatchProcessRunner, self).__init__()
        self.row_processed_signal.connect(notify_progress)
        self.row_failed_signal.connect(notify_error)
        self.notify_done = notify_done
        self.batch_processor = batch_reduction
        self._worker = None
        self._finish_up = False

    def finish_up_processing(self):
        self._finish_up = True

    @Slot()
    def on_finished(self):
        result = self._worker.result if self._worker else None
        self._worker = None
        self.notify_done(result)

    @Slot()
    def on_error(self, error):
        self._worker = None

    def process_states(self, states, use_optimizations, plot_results, save_results):
        self._finish_up = False
        self._worker = Worker(self._process_states_on_thread, states, use_optimizations, 
                              plot_results, save_results)
        self._worker.signals.finished.connect(self.on_finished)
        self._worker.signals.error.connect(self.on_error)

        QThreadPool.globalInstance().start(self._worker)

    def _process_states_on_thread(self, states, use_optimizations, plot_results, save_results):
        for key, state in states.items():
            try:
                results = []
                if not self._finish_up:
                    results = \
                        self.batch_processor([state], use_optimizations, plot_results, save_results)
                self.row_processed_signal.emit(key, results)

            except Exception as e:
                self.row_failed_signal.emit(key, str(e))

