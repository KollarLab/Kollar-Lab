# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 12:02:58 2020

@author: Kollarlab
"""
import copy 
import time
import os
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
import utility.plotting_tools as plots

def get_default_settings():
    fullsettings = {}
    settings = {}
    autoscan_settings = {}
    #Save location
    settings['scanname']    = 'flux_Scan'
    settings['meas_type']   = 'spec_flux_scan'
    settings['project_dir'] = r'Z:\Data'
    
    settings['start_voltage'] = -2
    settings['stop_voltage'] =  2
    settings['voltage_points'] = 251
    
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
    
    settings['CAV_Attenuation'] = 30
    settings['Qbit_Attenuation'] = 10
    
    autoscan_settings['channel'] = 1
    autoscan_settings['measurement'] = 'S21'
    autoscan_settings['freq_points'] = 501
    autoscan_settings['ifBW'] = settings['ifBW']
    autoscan_settings['avg_time'] = 15
    autoscan_settings['start_freq'] = 7.6e9
    autoscan_settings['stop_freq'] = 7.7e9
    autoscan_settings['RFpower'] = settings['CAVpower']
    autoscan_settings['background_subtract'] = False
    
    fullsettings['spec'] = settings
    fullsettings['autoscan'] = autoscan_settings
    
    return fullsettings

def vna_spec_flux_scan(instruments, fullsettings):
    settings = fullsettings['spec']
    autoscan_set = fullsettings['autoscan']
    
    background_subtract = autoscan_set['background_subtract']
    ##Instruments used
    vna = instruments['VNA']
    SRS = instruments['DCsupply']

    ##Data saving and naming
    saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp = userfuncs.timestamp()
    filename = settings['scanname'] + '_' + stamp

    CAV_Attenuation = settings['CAV_Attenuation']
    Qbit_Attenuation = settings['Qbit_Attenuation']
    
    
    settings['CAVpower'] = settings['CAVpower'] + CAV_Attenuation
    settings['RFpower'] = settings['RFpower'] + Qbit_Attenuation
    autoscan_set['RFpower'] = settings['CAVpower']
    
    #set voltage sweep
    start_voltage = settings['start_voltage']
    stop_voltage  = settings['stop_voltage']
    voltage_points = settings['voltage_points']
    voltages = np.round(np.linspace(settings['start_voltage'], settings['stop_voltage'], settings['voltage_points']),6)
    max_voltage = 3.5
    SRS.range = '10 V'
    if np.max(voltages) > max_voltage:
        raise ValueError('max voltage too large!')
    else:
        settings['voltages'] = voltages
    
    trans_mags = np.zeros((settings['voltage_points'], autoscan_set['freq_points']))
    trans_phases = np.zeros((settings['voltage_points'], autoscan_set['freq_points']))
    
    mags = np.zeros((settings['voltage_points'], settings['freq_points']))
    phases = np.zeros((settings['voltage_points'], settings['freq_points']))
    
    SRS.Volt = 0
    SRS.Output = 'On'
    
    t0 = time.time()
    
    identifier = 'Cav Power : ' + str(settings['CAVpower'] - CAV_Attenuation) + ' dB'
    
    if background_subtract:
        vna.reset()
        print('Collecting background ripple, turning cavity power to 0 dBm (on vna)')
        back_settings = copy.deepcopy(autoscan_set)
        back_settings['RFpower'] = 0
        back_data = vna.trans_meas(back_settings)
        
    for vind in range(len(voltages)):
        voltage = voltages[vind]
        print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
        
        SRS.voltage_ramp(voltage)
        time.sleep(0.1)
        
        vna.reset()
        vna.output = 'on'
        
        print('trans')
        trans_data = vna.trans_meas(autoscan_set)
        trans_freqs = trans_data['xaxis']
        trans_mags[vind] = trans_data['mag']
        trans_phases[vind] = trans_data['phase']
        
        if background_subtract:
            trans_mags[vind] = trans_mags[vind] - back_data['mag']
        else:
            trans_mags[vind] = trans_mags[vind]
        
        settings['CAVfreq'] = trans_freqs[np.argmax(trans_mags[vind])]
        print('spec, CAV power: {}, cav freq: {}'.format(settings['CAVpower'], settings['CAVfreq']))
        data = vna.spec_meas(settings)
        
        vna.autoscale()
        
        ##took out this normalization
        #it was here from transmission power scans (and used to be devision.)
        #I've switched it to subtraction, but basically I do the real offset subtraction that I want lower down.
        
    #    normalization = np.mean(data['mag'][-10:]) #I think that this is only the -10th data point. And for log data we might want to add or subtract it, not devide by
    #    mags[vind] = data['mag'] - normalization
    #    
        
        mags[vind] = data['mag']
        phases[vind] = data['phase']
        freqs = data['xaxis']  
        
        if vind==0:
            t1=time.time()
            tdiff = t1-t0
            ttotal = tdiff*len(voltages)
            
            estimatedTime = tdiff*len(voltages)
            
            print('    ')
            print('Single run time: {}, estimated total time: {}'.format(tdiff, ttotal))
            print('    ')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60, 1)) + ' minutes')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60/60, 2)) + ' hours')
            print('    ')
            
        transdata = {}
        transdata['xaxis'] = trans_freqs
        transdata['mags'] = trans_mags[0:vind+1,:]
        transdata['phases'] = trans_phases[0:vind+1,:]
        
        specdata = {}
        specdata['xaxis'] = freqs
        specdata['mags'] = mags[0:vind+1,:]
        specdata['phases'] = phases[0:vind+1,:]
        
        singledata = {}
        singledata['xaxis'] = freqs
        singledata['mag'] = data['mag'] 
        singledata['phase'] = data['phase']
        
        trans_labels = ['Freq (GHz)','Voltage (V)']
        spec_labels = ['Freq (GHz)','Voltage (V)']
        
        #modify the spec data to subtract the offset in amp and phase
        #and then plot the modified version
        specplotdata = {}
        specplotdata['xaxis'] = specdata['xaxis']
        specplotdata['mags'] = specdata['mags']
        specplotdata['phases'] = specdata['phases']
        
        mat = np.copy(specplotdata['mags'])
        for ind in range(0, mat.shape[0]):
            mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
        specplotdata['mags'] = mat
        
        mat = np.copy(specplotdata['phases'])
        for ind in range(0, mat.shape[0]):
            mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
        specplotdata['phases'] = mat
        
        plots.autoscan_plot(transdata, specplotdata, singledata, voltages[0:vind+1], filename, trans_labels, spec_labels, identifier, fig_num = 1)
            
        if np.mod(vind,5) == 1:
            userfuncs.SaveFull(saveDir, filename, ['transdata', 'specdata', 'singledata', 'voltages', 
                                           'filename', 'trans_labels', 'spec_labels'], 
                                           locals(), expsettings=settings)
            plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
    
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-t0))
    
    #return to zero voltage
    SRS.voltage_ramp(0)
    SRS.Output = 'Off'
    vna.Output = 'Off'
    
    transdata = {}
    transdata['xaxis'] = trans_freqs
    transdata['mags'] = trans_mags
    transdata['phases'] = trans_phases
    
    specdata = {}
    specdata['xaxis'] = freqs
    specdata['mags'] = mags
    specdata['phases'] = phases
    
    singledata = {}
    singledata['xaxis'] = freqs
    singledata['mag'] = data['mag']
    singledata['phase'] = data['phase']
    
    trans_labels = ['Freq (GHz)','Voltage (V)']
    spec_labels = ['Freq (GHz)','Voltage (V)']
    
    
    plots.autoscan_plot(transdata, specplotdata, singledata, voltages, filename, trans_labels, spec_labels, identifier, fig_num = 1)
    
    userfuncs.SaveFull(saveDir, filename, ['transdata', 'specdata', 'singledata', 'voltages', 
                                           'filename', 'trans_labels', 'spec_labels'], 
                                           locals(), expsettings=settings)
    
    plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)