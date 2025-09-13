# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 12:02:58 2020

@author: Kollarlab
"""
import copy 
import time
import os
from utility.measurement_helpers import estimate_time
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
import utility.plotting_tools as plots

def get_default_settings():
    '''
    get_default_settings _summary_
    '''    
    fullsettings = {}
    settings = {}
    autoscan_settings = {}
    #Save location
    settings['scanname']    = 'beta_calibration_Scan'
    settings['meas_type']   = 'spec_beta_calibration_scan'
    #settings['project_dir'] = r'Z:\Data'
    
    #Sweep parameter
#    settings['CAV_Attenuation'] = 30
    settings['amplitude_start'] = 0.001
    settings['amplitude_stop']  = 0.001
    settings['amplitude_pts']   = 11
    settings['qubit_num'] = 1
    
    settings['RFport'] = 3
    settings['Mport'] = 2
    settings['CAVport'] = 1
    
    settings['channel'] = 1
    settings['avg_time'] = 10
    settings['measurement'] = 'S21'
    settings['start_freq'] = 6.5e9
    settings['stop_freq'] = 9.5e9
    settings['freq_points'] = 3001
    settings['CAVpower'] = -55
    settings['RFpower'] = -15
    settings['ifBW'] = .2e3
    
    settings['high_power_spec'] = False
    settings['unwrap_phase'] = True
    
#    settings['CAV_Attenuation'] = 30
#    settings['Qbit_Attenuation'] = 10
    
    autoscan_settings['channel'] = 1
    autoscan_settings['measurement'] = 'S21'
    autoscan_settings['freq_points'] = 501
    autoscan_settings['ifBW'] = settings['ifBW']
    autoscan_settings['avg_time'] = 15
    autoscan_settings['start_freq'] = 7.6e9
    autoscan_settings['stop_freq'] = 7.7e9
    autoscan_settings['RFpower'] = settings['CAVpower']
    autoscan_settings['background_subtract'] = False
    autoscan_settings['unwrap_phase'] = False
    
    fullsettings['spec'] = settings
    fullsettings['autoscan'] = autoscan_settings
    
    return fullsettings

def vna_spec_beta_calibration_scan(instruments, settings):
    '''
    vna_spec_beta_calibration_scan _summary_

    :param instruments: _description_
    :type instruments: _type_
    :param settings: _description_
    :type settings: _type_
    :raises ValueError: _description_
    :return: _description_
    :rtype: _type_
    '''    
    ##Instruments used
    vna = instruments['VNA']
    Dual_gen = instruments['DCsupply']

    vna.reset()

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    spec_set = exp_settings['spec']
    autoscan_set = exp_settings['autoscan']
    
    background_subtract = autoscan_set['background_subtract']

    ##Data saving and naming
    #saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = spec_set['scanname'] + '_' + stamp

    CAV_Attenuation  = exp_globals['CAV_Attenuation']
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']

    #Input a flux pair
    flux_pair = exp_settings['spec']['flux_pair']

    #Calculate corresponding voltage pair
    v2f = exp_globals['v2f']
    f2v = np.linalg.inv(v2f)
    flux_offsets = exp_globals['flux_offsets']

    intermediate_flux_pair = flux_pair - flux_offsets
    voltage_pair = f2v@intermediate_flux_pair

    qubit_num = exp_settings['spec']['qubit_num']

    #Create a modulation amplitude array
    start_amplitude  = exp_settings['spec']['amplitude_start']
    stop_amplitude   = exp_settings['spec']['amplitude_stop']
    amplitude_points = exp_settings['spec']['amplitude_pts']
    amplitudes = np.linspace(start_amplitude, stop_amplitude, amplitude_points)
    
    spec_set['CAVpower'] = spec_set['CAVpower'] + CAV_Attenuation
    spec_set['RFpower']  = spec_set['RFpower']  + Qbit_Attenuation
    if not spec_set['high_power_spec']:
        print('Using CAVpower from general settings, ignoring autoscan power')
        autoscan_set['RFpower'] = spec_set['CAVpower']
    else:
        autoscan_set['RFpower'] = autoscan_set['RFpower'] + CAV_Attenuation
    
    #set amplitude sweep
    max_voltage = 4
    if np.amax(np.abs(voltage_pair)) > max_voltage:
        raise ValueError('max voltage too large!')
    else:
        exp_settings['voltages'] = voltage_pair
    
    trans_mags   = np.zeros((len(amplitudes), autoscan_set['freq_points']))
    trans_phases = np.zeros((len(amplitudes), autoscan_set['freq_points']))
    
    mags   = np.zeros((len(amplitudes), spec_set['freq_points']))
    phases = np.zeros((len(amplitudes), spec_set['freq_points']))
    
    if Dual_gen.Ch1_output == 0:
        Dual_gen.Ch1_dc_voltage_ramp(0)
        Dual_gen.Ch1_output = 1

    if Dual_gen.Ch2_output == 0:
        Dual_gen.Ch2_dc_voltage_ramp(0)
        Dual_gen.Ch2_output = 1

    Dual_gen.Ch1_dc_voltage_ramp(-np.round(voltage_pair[0], 4)) #minus sign here is because FBL_1 has the opposite polarity with the magnet and FBL_2
    Dual_gen.Ch2_dc_voltage_ramp(np.round(voltage_pair[1], 4))
    
    tstart = time.time()
    
    identifier = 'Cav Power : ' + str(spec_set['CAVpower'] - CAV_Attenuation) + ' dB'
    
    if background_subtract:
        vna.reset()  
        print('Collecting background ripple, turning cavity power to 0 dBm (on vna)')
        back_settings = copy.deepcopy(autoscan_set)
        back_settings['RFpower'] = 0
        back_data = vna.trans_meas(back_settings)
        
    for aind in range(len(amplitudes)):

        if aind == 0:
            print('dc_offset_1: {} V'.format(-np.round(voltage_pair[0], 4)))
            print('dc_offset_2: {} V'.format(np.round(voltage_pair[1], 4)))
            print('you are modulating qubit {}'.format(qubit_num))

        print('')
        print('Progress: ' + str(aind+1) + '/' + str(len(amplitudes)))                            
        print('Current amplitude {}, Ending amplitude {}'.format(str(np.round(amplitudes[aind],3)),str(amplitudes[-1])))
        print('')

        time.sleep(0.1)
            
        
        vna.reset()  
        vna.output = 'on'
        
        print('trans')
        trans_data  = vna.trans_meas(autoscan_set)
        trans_freqs = trans_data['xaxis']
        trans_mags[aind]   = trans_data['mag']
        trans_phases[aind] = trans_data['phase']
        
        if background_subtract:
            trans_mags[aind] = trans_mags[aind] - back_data['mag']
        else:
            trans_mags[aind] = trans_mags[aind]

        hanger = exp_globals['hanger']
        if hanger:
            spec_set['CAVfreq'] = trans_freqs[np.argmin(trans_mags[aind])] 
        else:
            spec_set['CAVfreq'] = trans_freqs[np.argmax(trans_mags[aind])]

        print('spec, CAV power: {}, cav freq: {}'.format(spec_set['CAVpower'], spec_set['CAVfreq']))

        #Turn on modulation after finding the cavity (only for the first time, afterwards it leaves the modulation on)

        if qubit_num == 1:
            Dual_gen.Ch1_sin_gen(amplitudes[aind], exp_settings['spec']['modulation_frequency'], phase=0, offset=-voltage_pair[0]) #notice the minus sign here
        else:
            Dual_gen.Ch2_sin_gen(amplitudes[aind], exp_settings['spec']['modulation_frequency'], phase=0, offset=voltage_pair[1])
        

        data = vna.spec_meas(spec_set)
        
        vna.autoscale()
        
        mags[aind]   = data['mag']
        phases[aind] = data['phase']

        freqs = data['xaxis']  
        
        if aind==0:
            tstop=time.time()
            estimate_time(tstart, tstop, len(amplitudes))
            
        transdata = {}
        transdata['xaxis'] = trans_freqs/1e9
        transdata['mags'] = trans_mags[0:aind+1,:]
        transdata['phases'] = trans_phases[0:aind+1,:]
        
        specdata = {}
        specdata['xaxis'] = freqs/1e9
        specdata['mags'] = mags[0:aind+1,:]
        specdata['phases'] = phases[0:aind+1,:]
        
        singledata = {}
        singledata['xaxis'] = freqs/1e9
        singledata['mag'] = data['mag'] 
        singledata['phase'] = data['phase']
        
        trans_labels = ['Freq (GHz)','Amplitude (Vpp)']
        spec_labels  = ['Freq (GHz)','Amplitude (Vpp)']
        
        #modify the spec data to subtract the offset in amp and phase
        #and then plot the modified version
        specplotdata = {}
        specplotdata['xaxis']  = specdata['xaxis']
        specplotdata['mags']   = specdata['mags']
        specplotdata['phases'] = specdata['phases']
        
        mat = np.copy(specplotdata['mags'])
        for ind in range(0, mat.shape[0]):
            mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
        specplotdata['mags'] = mat
        
        mat = np.copy(specplotdata['phases'])
        for ind in range(0, mat.shape[0]):
            mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
        specplotdata['phases'] = mat
        
        plots.autoscan_plot(transdata, specplotdata, singledata, amplitudes[0:aind+1], filename, trans_labels, spec_labels, identifier, fig_num = 1)
            
        userfuncs.SaveFull(saveDir, filename, ['transdata', 'specdata', 'singledata', 'flux_pair', 'voltage_pair', 'amplitudes',
                                       'filename', 'trans_labels', 'spec_labels'], 
                                       locals(), expsettings=settings, instruments=instruments)
        plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
    
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))

    #return to zero voltage
    Dual_gen.Ch1_dc_voltage_ramp(0)
    Dual_gen.Ch2_dc_voltage_ramp(0)
    Dual_gen.Ch1_output = 0
    Dual_gen.Ch2_output = 0
    
    data = {'saveDir': saveDir, 'filename': filename, 'transdata':transdata, 'specdata':specdata, 'amplitudes':amplitudes}

    return data

    