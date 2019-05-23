# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
# pylint: disable=invalid-name
""" SANBatchReduction algorithm is the starting point for any new type reduction, event single reduction"""
from __future__ import (absolute_import, division, print_function)

import os

from mantid.api import (FileFinder)
from mantid.simpleapi import (LoadBBY, LoadMask, MaskDetectors, ConvertUnits, Rebin, 
                              SumSpectra, BilbySANSDataProcessor, SaveNISTDAT, SaveAscii)
from sans.state.state_base import StateBase
from mantidplot import plotSpectrum, DistrFlag, mergePlots, plot2D
import BilbyCustomFunctions_Reduction as bby

class BilbyBatchReduction(object):
    def __init__(self):
        super(BilbyBatchReduction, self).__init__()

    def __call__(self, states, use_optimizations=True, plot_results = False, save_results=False):
        """
        This is the start of any reduction.

        :param states: This is a list of sans states. Each state in the list corresponds to a single reduction.
        :param use_optimizations: if True then the optimizations for file reloading are used.
        :param output_mode: The output mode defines how the reduced data should be published. This can be
                            1. PublishToADS
                            2. SaveToFile
                            3. Both
        """
        self.validate_inputs(states, use_optimizations, plot_results)

        return self._execute(states, use_optimizations, plot_results, save_results)

    def _execute(self, states, use_optimizations, plot_results, save_results):

        # Iterate over each state, load the data and perform the reduction
        results = []
        for state in states:
            result = single_reduction_for_batch(state, use_optimizations, plot_results, save_results)
            if result:
                results.append(result)
        return results

    def validate_inputs(self, states, use_optimizations, plot_results):
        # We are strict about the types here.
        # 1. states has to be a list of sans state objects
        # 2. use_optimizations has to be bool
        if not isinstance(states, list):
            raise RuntimeError("The provided states are not in a list. They have to be in a list.")

        for state in states:
            if not isinstance(state, StateBase):
                raise RuntimeError("The entries have to be sans state objects. "
                                   "The provided type is {0}".format(type(state)))

        if not isinstance(use_optimizations, bool):
            raise RuntimeError("The optimization has to be a boolean. The provided type is"
                               " {0}".format(type(use_optimizations)))

        if not isinstance(plot_results, bool):
            raise RuntimeError("The plot_result has to be a boolean. The provided type is"
                               " {0}".format(type(plot_results)))

        errors = self._validate_inputs(states)
        if errors:
            raise RuntimeError("The provided states are not valid: {}".format(errors))

    def _validate_inputs(self, states):
        errors = dict()
        # Check that the input can be converted into the right state object
        try:
            for state in states:
                state.validate()
        except ValueError as err:
            errors.update({"BilbyBatchReduction": str(err)})
        return errors

# -------------------------------------------------------------------------------------
# Bilby specific - stuff above to be generalized
# -------------------------------------------------------------------------------------

# values for att_pos 2 and 4 shall not make sense; those attenuators have not been in use that time
attenuation_pre_may_2016 = {1: 0.007655, 2: -1.0, 3: 1.0, 4: -1.0, 5: 0.005886}
attenuation_post_may_2016 = {1: 1.0, 2: 0.00955, 3: 0.005886, 4: 0.00290, 5: 0.00062}

def load_workspaces(state):
    # returns a dictionary map to the workspaces
    workspaces = {'blocked': None}

    ws_sample = LoadBBY(Filename=state.sample,
                        FilterByTimeStart=state.time_slice.start_time,
                        FilterByTimeStop=state.time_slice.end_time)
    workspaces['sample'] = ws_sample

    # apply tube shift correction to sample workspace
    # ** this should move to the Bilby loader 
    tube_shift_fpath = FileFinder.getFullPath('shift_assembled.csv')
    bby.correction_tubes_shift(ws_sample, tube_shift_fpath)

    ws_empty = LoadBBY(Filename=state.transmission_empty)
    workspaces['t_empty'] = ws_empty

    ws_trans = LoadBBY(Filename=state.transmission_sample)
    workspaces['t_sample'] = ws_trans

    if state.blocked_beam:
        ws_blocked = LoadBBY(Filename=state.blocked_beam)
        workspaces['blocked'] = ws_blocked

        # TODO tube correction should be applied to blocked beam
        bby.correction_tubes_shift(ws_blocked, tube_shift_fpath)

    masks = {'sample': None}
    mk_trans = LoadMask('BILBY', state.transmission_mask + '.xml')
    masks['trans'] = mk_trans
    if  state.sample_mask:
        mk_sample = LoadMask('BILBY', state.sample_mask + '.xml')
        masks['sample'] = mk_sample

    return workspaces, masks

def _build_output_name(state, low_wave, high_wave):

    tag2D = '2D' if state.reduce_2D else ''
    ts = state.time_slice
    start = '{:.1f}'.format(ts.start_time) if ts.start_time else ''
    end = '{:.1f}'.format(ts.end_time) if ts.end_time else ''
    wtag = '{:.2f}_{:.2f}'.format(low_wave, high_wave)
    tags = [state.sample, tag2D, wtag, start, end, state.suffix]
    return '_'.join([x for x in tags if x])

def _build_output_header(external_mode, used_wl_range, ws_sample, sample_thickness,
                         sample_transmission, empty_beam_transmission, blocked_beam, 
                         sample_mask, transmission_mask):
    # log number value as a string 
    def log_value(run, tag, ndigits, scale = 1.0):
        value = scale * run.getProperty(tag).value
        return str(round(value, ndigits))

    srun = ws_sample.run()
    header = []
    header.append('Velocity selector set wavelength: ' + log_value(srun, 'wavelength', 3) + ' Angstrom')

    if external_mode:
        header.append('Double choppers pair: ' + log_value(srun, 'master1_chopper_id', 0) + \
                      ' and ' + log_value(srun, 'master2_chopper_id', 0))
        frequency = 'Data defining pulse frequency (equal or slower than the Double pair frequency): ' \
                    + str(format(1e6/float(ws_sample.run().getProperty("period").value), '.2f')) + ' Hz'
        header.append(frequency)
        wavelength_range = 'Wavelength range used for the data reduction: ' + str(format(float(used_wl_range[0]), '.2f')) + ' to ' \
                                                                            + str(format(float(used_wl_range[2]), '.2f')) + ' Angstrom'
        header.append(wavelength_range)
        resolution_value = float(used_wl_range[1])
        if resolution_value < 0:
            resolution = 'Resolution used for calculation of dQ: ' + str(format((-100 * resolution_value), '.2f')) + '%'
        else:
            resolution = 'Resolution taken as wavelength binning;' + '\n' + 'the value is set to ' + \
                str(format(resolution_value, '.2f')) + '%'  # on linear scale, hence the dQ calculation is meaningless'
        header.append(resolution)
    else:
        resolution = "Nominal resolution: 10%"
        header.append(resolution)

    header.append('L1: ' + log_value(srun, 'L1', 3) + ' m')
    header.append('L2 to rear detector: ' + log_value(srun, 'L2_det_value', 3) + ' m')
    header.append('L2 to vertical curtains: ' + log_value(srun, 'L2_curtainr_value', 3) + ' m')
    header.append('Left curtain separation: ' + log_value(srun, 'D_curtainl_value', 3) + ' m')
    header.append('Right curtain separation: ' + log_value(srun, 'D_curtainr_value', 3) + ' m')
    header.append('Top curtain separation: ' + log_value(srun, 'D_curtainu_value', 3) + ' m')
    header.append('Bottom curtain separation: ' + log_value(srun, 'D_curtaind_value', 3) + ' m')

    header.append('Source and sample apertures diameters: ' + log_value(srun, 'source_aperture', 1) + \
                  ' mm and ' + log_value(srun, 'sample_aperture', 1) + ' mm')

    header.append('Sample thickness and transmission: ' + \
                  format(float(sample_thickness), '.2f') + ' cm and ' + sample_transmission)
    header.append('Empty beam transmission and blocked beam scattering: ' + \
                  empty_beam_transmission + ' and ' + blocked_beam)
    header.append('Sample and trasmission masks: ' + sample_mask + ' and ' + transmission_mask + '\n')

    return header
    

def single_reduction_for_batch(state, use_optimizations, plot_results, save_results):
    # this implements the equivalent bilby reduction code
    print ('single_reduction_for_batch')

    # load the workspaces and masks  
    wks, mks = load_workspaces(state)  
    
    # create the output directory if it does not exist
    if not os.path.exists(state.output_folder):
        os.makedirs(state.output_folder)

    # set attenuation factor
    if state.pre_may_2016_data:
        attenuation = attenuation_pre_may_2016[state.att_pos]
    else:
        attenuation = attenuation_post_may_2016[state.att_pos]

    # empty beam normalisation
    # does not have to be ws_tranMsk, can be a specific mask
    ws_emp = wks['t_empty']
    t_mask = mks['trans'] 
    MaskDetectors(ws_emp, MaskedWorkspace=t_mask)
    ws_emp = ConvertUnits(ws_emp, Target="Wavelength")

    # run through the slices in the wavelength
    plot1Dgraph = None
    for (lo_wave, hi_wave) in zip(state.wavelength.wavelength_low,
                                  state.wavelength.wavelength_high):

        wave_step = state.wavelength.wavelength_step
        ws_emp_partial = Rebin(ws_emp, Params=[lo_wave, wave_step, hi_wave])
        ws_emp_partial = SumSpectra(ws_emp_partial, IncludeMonitors=False)  

        base_output_name = _build_output_name(state, lo_wave, hi_wave)        
        
        def format_binning(low, step, high, ndigits): 
            s_low = str(round(low, ndigits))
            s_high = str(round(high, ndigits))
            s_step = str(round(step, ndigits+1))
            return ','.join([s_low, s_step, s_high])

        wave_bins = format_binning(lo_wave, wave_step, hi_wave, 3)
        ts = state.transmission_wavelength
        trans_bins = format_binning(ts.low, ts.step, ts.high, 3)
        qs = state.binning_q
        factor = -1. if state.backward_q_bin else 1.
        q_bins = format_binning(qs.low, factor * qs.step, qs.high, 5)

        output_workspace, transmission_fit = \
            BilbySANSDataProcessor(
                InputWorkspace=wks['sample'], InputMaskingWorkspace=mks['sample'],
                BlockedBeamWorkspace=wks['blocked'], EmptyBeamSpectrumShapeWorkspace=ws_emp_partial, SensitivityCorrectionMatrix=None,
                TransmissionWorkspace=wks['t_sample'], TransmissionEmptyBeamWorkspace=wks['t_empty'], TransmissionMaskingWorkspace=mks['trans'],
                ScalingFactor=attenuation, SampleThickness=state.sample_thickness, 
                FitMethod=state.fit_method, PolynomialOrder=str(state.polynomial_order),
                BinningWavelength=wave_bins, BinningWavelengthTransm=trans_bins, BinningQ=q_bins,
                TimeMode=state.external_mode, AccountForGravity=state.gravity_correction, SolidAngleWeighting=state.solid_angle_weighting,
                RadiusCut = state.radius_cut, WaveCut = state.wave_cut, 
                WideAngleCorrection=state.wide_angle_correction,
                Reduce2D=state.reduce_2D,
                OutputWorkspace=base_output_name)

        # if plot results
        if state.reduce_2D:
            plot2Dgraph = plot2D(base_output_name)
            n_2D = output_workspace.name() + '.png'
            save_plot = os.path.join(os.path.expanduser(state.output_folder), n_2D)
            plot2Dgraph.export(save_plot)
            print("2D File Exists: {}".format(os.path.exists(save_plot)))
            save_nxs = os.path.join(os.path.expanduser(state.output_folder), output_workspace.name() + '.nxs')
            SaveNISTDAT(output_workspace.name(), save_nxs)
            if not plot_results: 
                plot2Dgraph.close() 
                # is there more elegant way to do it? Problem is that plot2Dgraph creates and plot the graph file at the same time...
        else:
            bby.strip_NaNs(output_workspace, base_output_name)              
            if plot_results:
                if plot1Dgraph is None:
                    plot1Dgraph = plotSpectrum(base_output_name, 0, distribution=DistrFlag.DistrFalse,  clearWindow=False) # to create first graph to stick all the rest to it; perhaps there is more elegant way of creating initial empty handler, but I am not aware ... yet
                    plot1Dgraph.activeLayer().logLogAxes()
                else:             
                    mergePlots(plot1Dgraph, plotSpectrum(base_output_name, 0, distribution=DistrFlag.DistrFalse, clearWindow=False))
        
        #Section for file saving
        if save_results:
            n_1D = base_output_name +".dat"      
            savefile = os.path.join(os.path.expanduser(state.output_folder), n_1D) 

            # not clear what the objective is in AS code - 'something new 26 March 2019'
            header = _build_output_header(state.external_mode, [lo_wave, wave_step, hi_wave], 
                                          wks['sample'], state.sample_thickness, 
                                          state.transmission_sample, state.transmission_empty, 
                                          state.blocked_beam, state.sample_mask, 
                                          state.transmission_mask)
            with open(savefile, 'w') as f_out:
                for line in header:
                    f_out.write(line + '\n')

            # now save the file
            SaveAscii(InputWorkspace = base_output_name, Filename = savefile, WriteXError = True, 
                      WriteSpectrumID = False, Separator = "CSV", AppendToFile = True) 

    # finally just return the state processed 
    return state
 