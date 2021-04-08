# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 18:31:45 2020

@author: Kollarlab
"""
import time
import os
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
import utility.plotting_tools as plots

def get_default_settings():
    settings = {}
    
    #Save location
    settings['scanname']    = 'scanname'
    settings['meas_type']   = 'Autler_Townes'
    settings['project_dir'] = r'Z:\Data\defaultdir'

    #Sweep parameters
    settings['CAV_Attenuation'] = 30
    settings['Qbit_Attenuation'] = 10
    settings['Autler_Attenuation'] = 10

    settings['ext_flux'] = 0
    settings['autler_power'] = -20
    settings['start_autler_freq']  = 3.5e9
    settings['stop_autler_freq']   = 4.5e9
    settings['autler_points'] = 31

    #VNA settings
    settings['channel'] = 1
    settings['avg_time'] = 30
    settings['measurement'] = 'S21'
    settings['start_freq'] = 3.5e9
    settings['stop_freq'] = 4.5e9
    settings['freq_points'] = 501
    settings['RFpower'] = -25
    settings['RFport'] = 3
    settings['Mport'] = 2
    settings['CAVport'] = 1
    settings['CAVpower'] = -55
    settings['CAVfreq'] = 8.12555e9
    settings['ifBW'] = 1e3

    return settings

def vna_autler_townes(instruments, settings):
    #Instruments used
    vna = instruments['VNA']
    autlergen = instruments['RFsource']
    SRS = instruments['SRS']

    vna.reset()
    
    #Data saving and naming
    saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp = userfuncs.timestamp()
    filename = settings['scanname'] + '_' + stamp
    scanname = settings['scanname']

    CAV_Attenuation = settings['CAV_Attenuation']
    Qbit_Attenuation = settings['Qbit_Attenuation']
    Autler_Attenuation = settings['Autler_Attenuation']
    settings['CAVpower'] = settings['CAVpower'] + CAV_Attenuation
    settings['RFpower'] = settings['RFpower'] + Qbit_Attenuation
    settings['autler_power'] = settings['autler_power'] + Autler_Attenuation
    
    autlergen.power = settings['autler_power']
    autlergen.output = 'On'
    SRS.output = 'On'
    SRS.voltage_ramp(settings['ext_flux'])

    start_autler_freq = settings['start_autler_freq']
    stop_autler_freq = settings['stop_autler_freq']
    autler_points = settings['autler_points']
    autler_freqs = np.round(np.linspace(start_autler_freq, stop_autler_freq, autler_points),-3)

    findices = np.array(list(range(len(autler_freqs))))
    if settings['reverse']:
        findices = np.flipud(findices)
    
    if settings['random']:
        np.random.shuffle(findices)
        
    mags = np.zeros((len(autler_freqs), settings['freq_points']))
    phases = np.zeros((len(autler_freqs), settings['freq_points']))

    tstart = time.time()
    for freqind in findices:
        autler_freq = autler_freqs[freqind]
        print('Freq: {}, final freq: {}'.format(autler_freq, autler_freqs[-1]))
        
        autlergen.freq = autler_freq

        data = vna.spec_meas(settings)

        vna.autoscale()

        mags[freqind] = data['mag']
        phases[freqind] = data['phase']

        if freqind==0:
            tstop = time.time()
            singlePointTime = tstop-tstart
                
            estimatedTime = singlePointTime*len(autler_freqs)
            
            print('    ')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60, 1)) + ' minutes')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60/60, 2)) + ' hours')
            print('    ')

        freqs = data['xaxis']
        labels = ['Freq (GHz)', 'Autler freq (GHz)']
        full_data = {}
        single_data = {}
        if not settings['random']:
            if settings['reverse']:
                full_data = {}
                full_data['xaxis'] = freqs
                full_data['mags'] = mags[freqind:]
                full_data['phases'] = phases[freqind:]
            
                single_data = data
                yaxis = autler_freqs[freqind:]
            else:
                full_data = {}
                full_data['xaxis'] = freqs
                full_data['mags'] = mags[0:freqind+1]
                full_data['phases'] = phases[0:freqind+1]
            
                single_data = data
                yaxis = autler_freqs[0:freqind+1]
            
        
            plots.simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1)
            
        userfuncs.SaveFull(saveDir, filename, ['full_data', 'single_data', 'autler_freqs', 'labels', 'filename'], locals(), expsettings=settings)
        plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))
    if settings['random']:
        full_data = {}
        full_data['xaxis'] = freqs
        full_data['mags'] = mags
        full_data['phases'] = phases
    
        single_data = data
        yaxis = autler_freqs
        plots.simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1)
                
#    SRS.voltage_ramp(0.)
#    SRS.output = 'Off'
    autlergen.output = 'Off'
    
    userfuncs.SaveFull(saveDir, filename, ['full_data', 'single_data', 'autler_freqs', 'labels', 'filename'], locals(), expsettings=settings)
    plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)