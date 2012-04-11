from PyQt4 import QtGui, uic, QtCore
import reduction_gui.widgets.util as util
import math
import os
import time
import sys
from functools import partial
from reduction_gui.reduction.reflectometer.refl_sf_calculator_data_script import DataSets as REFLDataSets
from reduction_gui.reduction.reflectometer.refl_sf_calculator_data_series import DataSeries
from reduction_gui.settings.application_settings import GeneralSettings
#from base_ref_reduction import BaseRefWidget
from reduction_gui.widgets.base_widget import BaseWidget as BaseRefWidget
import ui.reflectometer.ui_refl_sf_calculator

IS_IN_MANTIDPLOT = False
try:
    import mantidplot
    from MantidFramework import *
    mtd.initialise(False)
    from mantidsimple import *
    import _qti
    from reduction.instruments.reflectometer import data_manipulation

    IS_IN_MANTIDPLOT = True
except:
    pass

class DataReflSFCalculatorWidget(BaseRefWidget):    
    """
        Widget that present instrument details to the user
    """
    ## Widget name
    name = "Data"      
    instrument_name = 'REF_L'
    short_name = 'REFL'
    peak_pixel_range = []
    background_pixel_range = []

    def __init__(self, parent=None, state=None, settings=None, name="REFL", data_proxy=None):      
        super(DataReflSFCalculatorWidget, self).__init__(parent, state, settings, data_proxy=data_proxy) 

        class SummaryFrame(QtGui.QFrame, ui.reflectometer.ui_refl_sf_calculator.Ui_Frame):
            def __init__(self, parent=None):
                QtGui.QFrame.__init__(self, parent)
                self.setupUi(self)
                
        self.short_name = name
        self._settings.instrument_name = name
            
        self._summary = SummaryFrame(self)
        self.initialize_content()
        self._layout.addWidget(self._summary)

        if state is not None:
            self.set_state(state)
        else:
            self.set_state(DataSeries(data_class=REFLDataSets))

    def initialize_content(self):
        """
        Make sure the text fields accept only the right format of data
        """

        #hide labels
        self._summary.waiting_label.hide()
        self._summary.data_run_number_processing.hide()
                
        self._summary.data_run_number_edit.setValidator(QtGui.QIntValidator(self._summary.data_run_number_edit))
        self._summary.number_of_attenuator.setValidator(QtGui.QIntValidator(self._summary.number_of_attenuator))
        self._summary.data_peak_from_pixel.setValidator(QtGui.QIntValidator(self._summary.data_peak_from_pixel))
        self._summary.data_peak_to_pixel.setValidator(QtGui.QIntValidator(self._summary.data_peak_to_pixel))
        self._summary.data_background_from_pixel.setValidator(QtGui.QIntValidator(self._summary.data_background_from_pixel))
        self._summary.data_background_to_pixel.setValidator(QtGui.QIntValidator(self._summary.data_background_to_pixel))
        
        #Event connections
        self.connect(self._summary.data_run_number_edit, QtCore.SIGNAL("returnPressed()"), self.data_run_number_validated)
        self.connect(self._summary.add_dataset_btn, QtCore.SIGNAL("clicked()"), self._add_data)
        self.connect(self._summary.data_background_switch, QtCore.SIGNAL("clicked(bool)"), self._data_background_clicked)
        self.connect(self._summary.remove_btn, QtCore.SIGNAL("clicked()"), self._remove_item)
        self.connect(self._summary.plot_count_vs_y_btn, QtCore.SIGNAL("clicked()"), self._plot_count_vs_y)
        self.connect(self._summary.plot_count_vs_y_bck_btn, QtCore.SIGNAL("clicked()"), self._plot_count_vs_y_bck)
        self.connect(self._summary.angle_list, QtCore.SIGNAL("itemSelectionChanged()"), self._angle_changed)
        
        #Catch edited controls        
        #Incident medium (selection or text changed)
        call_back = partial(self._edit_event, ctrl=self._summary.incident_medium_combobox)
        self.connect(self._summary.incident_medium_combobox, QtCore.SIGNAL("editTextChanged(QString)"), call_back)
        #Number of attenuator value changed
        call_back = partial(self._edit_event, ctrl=self._summary.number_of_attenuator)
        self.connect(self._summary.number_of_attenuator, QtCore.SIGNAL("textChanged(QString)"), call_back)
        #peak selection (from and to) changed
        call_back = partial(self._edit_event, ctrl=self._summary.data_peak_from_pixel)
        self.connect(self._summary.data_peak_from_pixel, QtCore.SIGNAL("textChanged(QString)"), call_back)
        call_back = partial(self._edit_event, ctrl=self._summary.data_peak_to_pixel)
        self.connect(self._summary.data_peak_to_pixel, QtCore.SIGNAL("textChanged(QString)"), call_back)
        #background flag and from/to textEdit changed
        call_back = partial(self._edit_event, ctrl=self._summary.data_background_switch)
        self.connect(self._summary.data_background_switch, QtCore.SIGNAL("clicked()"), call_back)
        call_back = partial(self._edit_event, ctrl=self._summary.data_background_from_pixel)
        self.connect(self._summary.data_background_from_pixel, QtCore.SIGNAL("textChanged(QString)"), call_back)
        call_back = partial(self._edit_event, ctrl=self._summary.data_background_to_pixel)
        self.connect(self._summary.data_background_to_pixel, QtCore.SIGNAL("textChanged(QString)"), call_back)
                
    def _ref_instrument_selected(self):
        self.instrument_name = "REF_L"
        self._summary.center_pix_radio.hide()
        self._summary.center_pix_edit.hide()
        self._summary.angle_radio.hide()
        self._summary.angle_edit.hide()
        self._summary.angle_unit_label.hide()
        self._summary.angle_offset_label.show()
        self._summary.angle_offset_edit.show()
        self._summary.angle_offset_pm_label.show()
        self._summary.angle_offset_error_edit.show()
        self._summary.angle_offset_unit_label.show()
        self._summary.det_angle_offset_check.hide()
        self._summary.det_angle_offset_edit.hide()
        self._summary.det_angle_offset_unit_label.hide()
        self._summary.det_angle_check.hide()
        self._summary.det_angle_edit.hide()
        self._summary.det_angle_unit_label.hide()
        self._summary.direct_pixel_check.hide()
        self._summary.direct_pixel_edit.hide()
        self._summary.q_bins_label.hide()
        self._summary.ref_pix_estimate.hide()
        
        # Output directory
        self._summary.outdir_label.hide()
        self._summary.outdir_edit.hide()
        self._summary.outdir_browse_button.hide()

        #TODO: allow log binning
        self._summary.log_scale_chk.hide()
                 
    def _plot_count_vs_y(self, is_peak=True):
        """
            Plot counts as a function of high-resolution pixels
            and select peak range
            For REFM, this is X
            For REFL, this is Y
        """
        min, max = self._integrated_plot(True,
                                         self._summary.data_run_number_edit,
                                         self._summary.data_peak_from_pixel,
                                         self._summary.data_peak_to_pixel)
#        self._summary.data_peak_from_pixel_estimate.setText(str(int(math.ceil(min))))
#        self._summary.data_peak_to_pixel_estimate.setText(str(int(math.ceil(max))))
#        util.set_tiny(self._summary.data_peak_from_pixel_estimate)
#        util.set_tiny(self._summary.data_peak_to_pixel_estimate)

    def _plot_count_vs_y_bck(self):
        """
            Plot counts as a function of high-resolution pixels
            and select background range
            For REFM, this is X
            For REFL, this is Y
        """
        self._integrated_plot(True,
                              self._summary.data_run_number_edit,
                              self._summary.data_background_from_pixel,
                              self._summary.data_background_to_pixel)

    def _integrated_plot(self, is_high_res, file_ctrl, min_ctrl, max_ctrl):
        """
            Plot counts as a function of:
            
            Low-resolution pixels
                For REFM, this is Y
                For REFL, this is X
                
            High-resolution pixels
                For REFM, this is X
                For REFL, this is Y  
                
            @param is_high_res: True if we are plotting the high-res pixel distribution
            @param file_ctrl: control widget containing the data file name
            @param min_ctrl: control widget containing the range minimum
            @param max_ctrl: control widget containing the range maximum
        """
        if not IS_IN_MANTIDPLOT:
            return
        
        f = FileFinder.findRuns("%s%s" % (self.instrument_name, str(file_ctrl.text())))

        range_min = int(min_ctrl.text())
        range_max = int(max_ctrl.text())

        if len(f)>0 and os.path.isfile(f[0]):
            def call_back(xmin, xmax):
                min_ctrl.setText("%-d" % int(xmin))
                max_ctrl.setText("%-d" % int(xmax))
            
            # For REFL, Y is high-res
            is_pixel_y = is_high_res
            # For REFM it's the other way around
            if self.short_name == "REFM":
                is_pixel_y = not is_pixel_y
                
            min, max = data_manipulation.counts_vs_pixel_distribution(f[0], is_pixel_y=is_pixel_y,
                                                                      callback=call_back,
                                                                      range_min=range_min,
                                                                      range_max=range_max,
                                                                      high_res=is_high_res,
                                                                      instrument=self.short_name)
            return min, max

    def _remove_item(self):
        if self._summary.angle_list.count()==0:
            return
        self._summary.angle_list.setEnabled(False)        
        self._summary.remove_btn.setEnabled(False)  
        row = self._summary.angle_list.currentRow()
        if row>=0:
            self._summary.angle_list.takeItem(row)
        self._summary.angle_list.setEnabled(True)        
        self._summary.remove_btn.setEnabled(True)  

    def _edit_event(self, text=None, ctrl=None):
        self._summary.edited_warning_label.show()
        util.set_edited(ctrl,True)

    def _run_number_changed(self):
        self._edit_event(ctrl=self._summary.data_run_number_edit)

    def data_run_number_validated(self):
        self._summary.data_run_number_processing.show()
        run_number = self._summary.data_run_number_edit.text()
        _file = FileFinder.findRuns("REF_L%d"%int(run_number))
        lambdaRequest = ''
        S1H, S2H, S1W, S2W, lambdaRequest = self.getSlitsValueAndLambda(_file[0])
        self._summary.s1h.setText(str(S1H))
        self._summary.s2h.setText(str(S2H))
        self._summary.s1w.setText(str(S1W))
        self._summary.s2w.setText(str(S2W))
        self._summary.lambda_request.setText(str(lambdaRequest))
        self._summary.data_run_number_processing.hide()
        
    def _add_data(self):
        state = self.get_editing_state()
        in_list = False
        # Check whether it's already in the list
        run_numbers = self._summary.data_run_number_edit.text()
        if (run_numbers == ''):
            return
        list_items = self._summary.angle_list.findItems(run_numbers, QtCore.Qt.MatchFixedString)
        if len(list_items)>0:
            list_items[0].setData(QtCore.Qt.UserRole, state)
            in_list = True
        else:
            item_widget = QtGui.QListWidgetItem(run_numbers, self._summary.angle_list)
            item_widget.setData(QtCore.Qt.UserRole, state)
        
#        # Read logs
#        if not in_list and self.short_name == "REFM":
#            self._read_logs()
        
        self._reset_warnings()
    
    def _data_background_clicked(self, is_checked):
        """
            This is reached when the user clicks the Background switch and will enabled or not
            the widgets that follow that button
        """
        self._summary.data_background_from_pixel.setEnabled(is_checked)
        self._summary.data_background_from_pixel_label.setEnabled(is_checked)
        self._summary.data_background_to_pixel.setEnabled(is_checked)
        self._summary.data_background_to_pixel_label.setEnabled(is_checked)
        self._summary.plot_count_vs_y_bck_btn.setEnabled(is_checked)
        self._edit_event(None, self._summary.data_background_switch)
    
    def _reset_warnings(self):
        self._summary.edited_warning_label.hide()
        util.set_edited(self._summary.data_run_number_edit, False)
        util.set_edited(self._summary.incident_medium_combobox, False)
        util.set_edited(self._summary.number_of_attenuator, False)
        util.set_edited(self._summary.data_peak_from_pixel, False)
        util.set_edited(self._summary.data_peak_to_pixel, False)
        util.set_edited(self._summary.data_background_switch, False)
        util.set_edited(self._summary.data_background_from_pixel, False)
        util.set_edited(self._summary.data_background_to_pixel, False)
    
    def _angle_changed(self):
        if self._summary.angle_list.count()==0:
            return
        self._summary.angle_list.setEnabled(False)  
        self._summary.remove_btn.setEnabled(False)  
        current_item =  self._summary.angle_list.currentItem()
        if current_item is not None:
            state = current_item.data(QtCore.Qt.UserRole).toPyObject()
            self.set_editing_state(state)
            self._reset_warnings()
        self._summary.angle_list.setEnabled(True)
        self._summary.remove_btn.setEnabled(True)  

    def set_state(self, state):
        """
            Populate the UI elements with the data from the given state. 
            @param state: data object    
        """
        super(DataReflSFCalculatorWidget, self).set_state(state)
        
#        if len(state.data_sets)>0:
#            self._summary.q_step_edit.setText(str(math.fabs(state.data_sets[0].q_step)))
            
#        self._reset_warnings()

    def set_editing_state(self, state):
    #    super(DataReflSFCalculatorWidget, self).set_editing_state(state)
                
        self._summary.data_run_number_edit.setText(state.data_file)
        self._summary.incident_medium_combobox.clear()
        self._summary.incident_medium_combobox.addItems(state.incident_medium_list)
        self._summary.incident_medium_combobox.setCurrentIndex(state.incident_medium_index_selected)
        self._summary.number_of_attenuator.setText(str(state.number_attenuator))
        self._summary.data_peak_from_pixel.setText(str(state.peak_selection[0]))
        self._summary.data_peak_to_pixel.setText(str(state.peak_selection[1]))
        self._summary.data_background_from_pixel.setText(str(state.back_selection[0]))
        self._summary.data_background_to_pixel.setText(str(state.back_selection[1]))
        self._summary.data_background_switch.setChecked(state.back_flag)
        self._summary.lambda_request.setText(state.lambda_requested)
        self._summary.s1h.setText(state.s1h)
        self._summary.s2h.setText(state.s2h)
        self._summary.s1w.setText(state.s1w)
        self._summary.s2w.setText(state.s2w)
    
    def get_state(self):
        """
            Returns an object with the state of the interface
        """
        m = self.get_editing_state()
        state = DataSeries(data_class=REFLDataSets)
        state_list = []
        
#        # Common Q binning
#        q_min = float(self._summary.q_min_edit.text())
#        q_step = float(self._summary.q_step_edit.text())
#        if self._summary.log_scale_chk.isChecked():
#            q_step = -q_step
#            
#        # Angle offset
#        if hasattr(m, "angle_offset"):
#            angle_offset = float(self._summary.angle_offset_edit.text())
#            angle_offset_error = float(self._summary.angle_offset_error_edit.text())
#                
#        for i in range(self._summary.angle_list.count()):
#            data = self._summary.angle_list.item(i).data(QtCore.Qt.UserRole).toPyObject()
#            # Over-write Q binning with common binning
#            data.q_min = q_min
#            data.q_step = q_step
#        
#            # Over-write angle offset
#            if hasattr(data, "angle_offset"):
#                data.angle_offset = angle_offset
#                data.angle_offset_error = angle_offset_error
#
#            state_list.append(data)
        state.data_sets = state_list
        
        return state
    
    def get_editing_state(self):
        
        m = REFLDataSets()

        #run number
        m.data_file = str(self._summary.data_run_number_edit.text())
        
        #incident medium
        m.incident_medium_list = [self._summary.incident_medium_combobox.itemText(i) 
                                  for i in range(self._summary.incident_medium_combobox.count())]
        m.incident_medium_index_selected = self._summary.incident_medium_combobox.currentIndex()
        
        #number of attenuator
        m.number_attenuator = int(self._summary.number_of_attenuator.text())

        #peak selection
        m.peak_selection = [int(self._summary.data_peak_from_pixel.text()),
                            int(self._summary.data_peak_to_pixel.text())]
        
        #background flag
        m.back_flag = self._summary.data_background_switch.isChecked()

        
        #background
        m.back_selection = [int(self._summary.data_background_from_pixel.text()),
                            int(self._summary.data_background_to_pixel.text())]
        
        #lambda request
        m.lambda_requested = self._summary.lambda_request.text()
        
        #s1h, s2h, s1w and s2w
        m.s1h = self._summary.s1h.text()
        m.s2h = self._summary.s2h.text()
        m.s1w = self._summary.s1w.text()
        m.s2w = self._summary.s2w.text()
        
        return m
    
    def getLambdaValue(self,mt):
        """
        return the lambdaRequest value
        """
        mt_run = mt.getRun()
        _lambda = mt_run.getProperty('LambdaRequest').value
        return _lambda
    
    def getSh(self,mt, top_tag, bottom_tag):
        """
            returns the height and units of the given slits
        """
        mt_run = mt.getRun()
        st = mt_run.getProperty(top_tag).value
        sb = mt_run.getProperty(bottom_tag).value
        sh = float(sb[0]) - float(st[0])
        units = mt_run.getProperty(top_tag).units
        return sh, units
    
    def getS1h(self,mt=None):
        """    
            returns the height and units of the slit #1 
        """
        if mt != None:
            _h, units = self.getSh(mt, 's1t', 's1b')
            return _h, units
        return None, ''
    
    def getS2h(self,mt=None):
        """    
            returns the height and units of the slit #2 
        """
        if mt != None:
            _h, units = self.getSh(mt, 's2t', 's2b') 
            return _h, units
        return None, None

    def getSw(self,mt, left_tag, right_tag):
        """
            returns the width and units of the given slits
        """
        mt_run = mt.getRun()
        sl = mt_run.getProperty(left_tag).value
        sr = mt_run.getProperty(right_tag).value
        sw = float(sl[0]) - float(sr[0])
        units = mt_run.getProperty(left_tag).units
        return sw, units

    def getS1w(self,mt=None):
        """    
            returns the width and units of the slit #1 
        """
        if mt != None:
            _w, units = self.getSw(mt, 's1l', 's1r') 
            return _w, units
        return None, ''
    
    def getS2w(self,mt=None):
        """    
            returns the width and units of the slit #2 
        """
        if mt != None:
            _w, units = self.getSh(mt, 's2l', 's2r') 
            return _w, units
        return None, None

    def getSlitsValueAndLambda(self,file):
        """
        Retrieve the S1H (slit 1 height), 
                     S2H (slit 2 height), 
                     S1W (slit 1 width), 
                     S2W (slit 2 width) and 
                     lambda requested values
        """
        _full_file_name = file
        LoadEventNexus(Filename=_full_file_name,
                       OutputWorkspace='tmpWks',
                       MetaDataOnly='1')
        mt1 = mtd['tmpWks']
        _s1h_value, _s1h_units = self.getS1h(mt=mt1)
        _s2h_value, _s2h_units = self.getS2h(mt=mt1)
        S1H = "%2.4f" %(_s1h_value)
        S2H = "%2.4f" %(_s2h_value)
        
        _s1w_value, _s1w_units = self.getS1w(mt=mt1)
        _s2w_value, _s2w_units = self.getS2w(mt=mt1)
        S1W = "%2.4f" %(_s1w_value)
        S2W = "%2.4f" %(_s2w_value)
        
        _lambda_value = self.getLambdaValue(mt=mt1)
        lambdaRequest = "%2.2f" %(_lambda_value[0])
        
        return S1H, S2H, S1W, S2W, lambdaRequest
        
        
    