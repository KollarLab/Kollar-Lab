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
    settings['scanname']    = 'flux_Scan'
    settings['meas_type']   = 'spec_flux_scan'
    #settings['project_dir'] = r'Z:\Data'
    
    #Sweep parameter
#    settings['CAV_Attenuation'] = 30
    settings['flux_start'] = np.array([0,0,0])
    settings['flux_stop']  = np.array([0,0,0])
    settings['flux_pts']   = 11
    settings['qubit_num'] = False
    
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

def vna_spec_flux_scan(instruments, settings):
    '''
    vna_spec_flux_scan _summary_

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

    #Make array of fluxes
    flux_start = exp_settings['flux_start']
    flux_stop  = exp_settings['flux_stop']
    flux_pts   = exp_settings['flux_pts']
    flux_holder = []

    for iind in range(2):
        flux_range = np.linspace(flux_start[iind],flux_stop[iind],flux_pts)
        flux_holder.append(flux_range)
    
    flux_holder = tuple(flux_holder)

    full_fluxes = np.stack(flux_holder,axis=1)

    #Calculate corresponding voltage array
    v2f = exp_globals['v2f']
    f2v = np.linalg.inv(v2f)
    flux_offsets = exp_globals['flux_offsets']

    full_voltages = np.zeros(full_fluxes.shape)

    for i in range(flux_pts):
        intermediate_fluxes = full_fluxes[i] - flux_offsets
        full_voltages[i] = f2v@intermediate_fluxes

    qubit_num = exp_settings['qubit_num']

    if qubit_num:
        fluxes = np.transpose(full_fluxes)[qubit_num-1]
    else:
        fluxes = np.linspace(1,len(full_voltages),len(full_voltages))

    if fluxes[0] == fluxes[-1]:
        print('')
        print('Warning! Selected qubit is not being swept, reverting to general numbering system')
        print('')
        exp_settings['qubit_num'] = False
        fluxes = np.linspace(1,len(full_voltages),len(full_voltages))
    
    
    spec_set['CAVpower'] = spec_set['CAVpower'] + CAV_Attenuation
    spec_set['RFpower']  = spec_set['RFpower']  + Qbit_Attenuation
    if not spec_set['high_power_spec']:
        print('Using CAVpower from general settings, ignoring autoscan power')
        autoscan_set['RFpower'] = spec_set['CAVpower']
    else:
        autoscan_set['RFpower'] = autoscan_set['RFpower'] + CAV_Attenuation
    
    #set voltage sweep
    max_voltage = 4
    if np.amax(np.abs(full_voltages)) > max_voltage:
        raise ValueError('max voltage too large!')
    else:
        exp_settings['voltages'] = full_voltages
    
    trans_mags   = np.zeros((len(full_voltages), autoscan_set['freq_points']))
    trans_phases = np.zeros((len(full_voltages), autoscan_set['freq_points']))
    
    mags   = np.zeros((len(full_voltages), spec_set['freq_points']))
    phases = np.zeros((len(full_voltages), spec_set['freq_points']))
    
    if Dual_gen.Ch1_output == 0:
        Dual_gen.Ch1_dc_voltage_ramp(0)
        Dual_gen.Ch1_output = 1

    if Dual_gen.Ch2_output == 0:
        Dual_gen.Ch2_dc_voltage_ramp(0)
        Dual_gen.Ch2_output = 1
    
    tstart = time.time()
    
    identifier = 'Cav Power : ' + str(spec_set['CAVpower'] - CAV_Attenuation) + ' dB'
    
    if background_subtract:
        vna.reset()  
        print('Collecting background ripple, turning cavity power to 0 dBm (on vna)')
        back_settings = copy.deepcopy(autoscan_set)
        back_settings['RFpower'] = 0
        back_data = vna.trans_meas(back_settings)
        
    for vind in range(len(full_voltages)):
        Dual_gen.Ch1_dc_voltage_ramp(-np.round(full_voltages[vind][0], 4))
        Dual_gen.Ch2_dc_voltage_ramp(np.round(full_voltages[vind][1], 4))
        time.sleep(0.1)

        if vind == 0:
            print('Voltage {}: initial {}, final {}'.format(str(1),-np.round(full_voltages[vind][0], 4),-np.round(full_voltages[-1][0], 4)))
            print('Voltage {}: initial {}, final {}'.format(str(2),np.round(full_voltages[vind][1], 4),np.round(full_voltages[-1][1], 4)))

        print('')
        print('Progress: ' + str(vind+1) + '/' + str(len(full_voltages)))                            
        print('Current Flux {}, Ending Flux {}'.format(str(np.round(fluxes[vind],4)),str(fluxes[-1])))
        print('')

        time.sleep(0.1)
            
        
        vna.reset()  
        vna.output = 'on'
        
        print('trans')
        trans_data  = vna.trans_meas(autoscan_set)
        trans_freqs = trans_data['xaxis']
        trans_mags[vind]   = trans_data['mag']
        trans_phases[vind] = trans_data['phase']
        
        if background_subtract:
            trans_mags[vind] = trans_mags[vind] - back_data['mag']
        else:
            trans_mags[vind] = trans_mags[vind]

        hanger = exp_globals['hanger']
        if hanger:
            spec_set['CAVfreq'] = trans_freqs[np.argmin(trans_mags[vind])] 
        else:
            spec_set['CAVfreq'] = trans_freqs[np.argmax(trans_mags[vind])]

        print('spec, CAV power: {}, cav freq: {}'.format(spec_set['CAVpower'], spec_set['CAVfreq']))

        data = vna.spec_meas(spec_set)
        
        vna.autoscale()
        
        mags[vind]   = data['mag']
        phases[vind] = data['phase']

        freqs = data['xaxis']  
        
        if vind==0:
            tstop=time.time()
            estimate_time(tstart, tstop, len(full_voltages))
            
        transdata = {}
        transdata['xaxis'] = trans_freqs/1e9
        transdata['mags'] = trans_mags[0:vind+1,:]
        transdata['phases'] = trans_phases[0:vind+1,:]
        
        specdata = {}
        specdata['xaxis'] = freqs/1e9
        specdata['mags'] = mags[0:vind+1,:]
        specdata['phases'] = phases[0:vind+1,:]
        
        singledata = {}
        singledata['xaxis'] = freqs/1e9
        singledata['mag'] = data['mag'] 
        singledata['phase'] = data['phase']
        
        trans_labels = ['Freq (GHz)','Flux']
        spec_labels  = ['Freq (GHz)','Flux']
        
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
        
        plots.autoscan_plot(transdata, specplotdata, singledata, fluxes[0:vind+1], filename, trans_labels, spec_labels, identifier, fig_num = 1)
            
        userfuncs.SaveFull(saveDir, filename, ['transdata', 'specdata', 'singledata', 'full_fluxes', 'fluxes' 
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
    
    data = {'saveDir': saveDir, 'filename': filename, 'transdata':transdata, 'specdata':specdata, 'voltages':voltages}

    return data

    