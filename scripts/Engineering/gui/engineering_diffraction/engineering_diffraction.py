# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
# pylint: disable=invalid-name
from qtpy import QtCore, QtWidgets

from mantidqt.icons import get_icon
from mantidqt.utils.observer_pattern import GenericObserverWithArgPassing, GenericObservable
from mantidqt.utils.qt import load_ui
from Engineering.gui.engineering_diffraction.presenter import EngineeringDiffractionPresenter
from .tabs.common import SavedirObserver
from .settings.settings_model import SettingsModel
from .settings.settings_view import SettingsView
from .settings.settings_presenter import SettingsPresenter

Ui_main_window, _ = load_ui(__file__, "main_window.ui")


class EngineeringDiffractionGui(QtWidgets.QMainWindow, Ui_main_window):
    """
    The engineering diffraction interface
    """

    # TODO 2021
    # Currently project save testing cannot work with the state of play in engdiffui
    # Need to seperate the ui into presenter and view, however the view will probably need to own the presenter


    status_savdirMaxwidth = 300

    def __init__(self, parent=None, window_flags=None):
        if window_flags is not None:
            super(EngineeringDiffractionGui, self).__init__(parent, window_flags)
        else:
            super(EngineeringDiffractionGui, self).__init__(parent)

        # Main Window
        self.setupUi(self)
        self.tabs = self.tab_main
        self.settings_presenter = None
        self.setFocusPolicy(QtCore.Qt.StrongFocus)

        self.btn_settings.setIcon(get_icon("mdi.settings", "black", 1.2))

        # Create status bar widgets
        self.status_label = QtWidgets.QLabel()
        self.savedir_label = QtWidgets.QLabel()
        self.savedir_label.setMaximumWidth(self.status_savdirMaxwidth)

        # observables
        self.close_event_observable = GenericObservable()

        # observers
        self.update_statusbar_text_observable = GenericObserverWithArgPassing(self.set_statusbar_text)
        self.update_savedir_observable = GenericObserverWithArgPassing(self.update_savedir)
        self.savedir_observer = SavedirObserver(self)

        # this presenter needs to be accessible to this view so that it can be accessed by project save
        self.presenter = self.setup_presenter()

        # setup that can only happen with presenter created
        self.setup_settings()
        self.setup_statusbar()

        # Usage Reporting
        try:
            import mantid

            # register startup
            mantid.UsageService.registerFeatureUsage(mantid.kernel.FeatureType.Interface,
                                                     "Engineering Diffraction", False)
        except ImportError:
            pass

    def setup_settings(self):
        model = SettingsModel()
        view = SettingsView(self)
        self.settings_presenter = SettingsPresenter(model, view)
        self.settings_presenter.load_settings_from_file_or_default()
        self.setup_savedir_notifier()

    def setup_presenter(self):
        presenter = EngineeringDiffractionPresenter(self)
        self.close_event_observable.add_subscriber(presenter.close_event_observer)
        presenter.statusbar_observable.add_subscriber(self.update_statusbar_text_observable)
        presenter.savedir_observable.add_subscriber(self.update_savedir_observable)
        self.set_on_settings_clicked(presenter.open_settings)
        self.set_on_help_clicked(presenter.open_help_window)
        self.set_on_instrument_changed(presenter.calibration_presenter.set_instrument_override)
        self.set_on_rb_num_changed(presenter.calibration_presenter.set_rb_num)
        self.set_on_instrument_changed(presenter.focus_presenter.set_instrument_override)
        self.set_on_rb_num_changed(presenter.focus_presenter.set_rb_num)
        return presenter

    def setup_savedir_notifier(self):
        self.settings_presenter.savedir_notifier.add_subscriber(self.savedir_observer)

    def closeEvent(self, _):
        self.close_event_observable.notify_subscribers()

    def setup_statusbar(self):
        self.statusbar.addWidget(self.status_label)
        self.set_statusbar_text("No Calibration Loaded.")
        self.statusbar.addWidget(self.savedir_label)
        self.update_savedir(self.settings_presenter.settings["save_location"])

    def set_statusbar_text(self, text):
        self.status_label.setText(text)

    def update_savedir(self, savedir):
        savedir_text = "SaveDir: " + savedir
        self.savedir_label.setToolTip(savedir_text)
        self.savedir_label.setText(savedir_text)

    def get_rb_no(self):
        return self.lineEdit_RBNumber.text()

    def set_on_help_clicked(self, slot):
        self.pushButton_help.clicked.connect(slot)

    def set_on_settings_clicked(self, slot):
        self.btn_settings.clicked.connect(slot)

    def set_on_rb_num_changed(self, slot):
        self.lineEdit_RBNumber.textChanged.connect(slot)

    def set_on_instrument_changed(self, slot):
        self.comboBox_instrument.currentIndexChanged.connect(slot)
