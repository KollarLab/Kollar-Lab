# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 14:07:57 2020

@author: Kollarlab
"""
import time
import os
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
import plotting_tools as plots

def get_default_settings():

    settings = {}
    
    #Save location
    settings['scanname']    = 'initial_power_scan_q4'
    settings['meas_type']   = 'Trans'
    settings['project_dir'] = r'Z:\Data\HouckQuadTransmon'
    
    #Sweep parameter
    settings['CAV_Attenuation'] = 30

    settings['start_power']  = -20
    settings['stop_power']   = 10
    settings['power_points'] = 31

    #VNA settings
    settings['channel']  = 1
    settings['avg_time'] = 1
    settings['measurement'] = 'S21'
    settings['start_freq']  = 8.7*1e9  
    settings['stop_freq']   = 8.8*1e9 
    settings['freq_points'] = 1001
    settings['ifBW'] = 1e3

    return settings

def vna_trans(instruments, settings):
    ##Instruments used
    vna = instruments['VNA']

    ##Data saving and naming
    saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp = userfuncs.timestamp()
    filename = settings['scanname'] + '_' + stamp

    CAV_Attenuation = settings['CAV_attenuation']
    scanname = settings['scanname']

    start_power = settings['start_power'] + CAV_Attenuation
    stop_power = settings['stop_power'] + CAV_Attenuation
    power_points = settings['power_points']
    powers = np.linspace(start_power, stop_power, power_points)

    mags = np.zeros((len(powers), settings['freq_points']))
    phases = np.zeros((len(powers), settings['freq_points']))

    t0 = time.time()
    for powerind in range(len(powers)):
        power = powers[powerind]
        print('Power: {}, final power: {}'.format(power-CAV_Attenuation, powers[-1]-CAV_Attenuation))
        settings['RFpower'] = power

        data = vna.trans_meas(settings)

        vna.autoscale()

        mags[powerind] = data['mag']
        phases[powerind] = data['phase']

        if powerind==0:
            tstop = time.time()
            singlePointTime = tstop-tstart
                
            estimatedTime = singlePointTime*len(powers)
            
            print('    ')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60, 1)) + ' minutes')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60/60, 2)) + ' hours')
            print('    ')
        
        freqs = data['xaxis']   

        full_data = {}
        full_data['xaxis'] = freqs
        full_data['mags'] = mags[0:powerind+1]
        full_data['phases'] = phases[0:powerind+1]
    
        single_data = data
    
        labels = ['Freq (GHz)', 'Power (dBm)']
        yaxis = powers[0:powerind+1]-CAV_Attenuation
        plots.simplescan_plot(full_data, single_data, yaxis, scanname, labels, identifier='', fig_num=2)
    
        userfuncs.SaveFull(saveDir, filename, ['mags', 'phases', 'freqs', 'powers'], locals(), expsettings=settings)
        
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-t0))

    plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)