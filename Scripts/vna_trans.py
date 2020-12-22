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

    start_power = settings['start_power']
    stop_power = settings['stop_power']
    power_points = settings['power_points']
    powers = np.linspace(start_power, stop_power, power_points)

    mags = np.zeros((len(powers), settings['freq_points']))
    phases = np.zeros((len(powers), settings['freq_points']))

    t0 = time.time()
    for powerind in range(len(powers)):
        power = powers[powerind]
        print('Power: {}, final power: {}'.format(power, powers[-1]))
        settings['RFpower'] = power

        data = vna.trans_meas(settings)

        vna.autoscale()

        mags[powerind] = data['mag']
        phases[powerind] = data['phase']

        if powerind==0:
            t1=time.time()
            tdiff = t1-t0
            ttotal = tdiff*len(powers)
            print('Single run time: {}, estimated total time: {}'.format(tdiff, ttotal))

    t2 = time.time()
    print('Elapsed time: {}'.format(t2-t0))

    freqs = data['xaxis']   

    full_data = {}
    full_data['xaxis'] = freqs
    full_data['mags'] = mags
    full_data['phases'] = phases

    single_data = {}
    single_data['xaxis'] = freqs
    single_data['mag'] = mags[-1]

    labels = ['Freq (GHz)', 'Power (dBm)']
    yaxis = powers
    plots.simplescan_plot(full_data, single_data, yaxis, scanname, labels, identifier='', fig_num=2)
    ###Martin's power plot version
    #VNAplots.power_plot(freqs, mags, phases, powers, scanname, -CAV_Attenuation)

    ####Alicia's general plot version
    #VNAplots.general_VNAplot(freqs, mags, phases, powers, scanname, 
    #                                 xlabel = 'Frequency (GHz)', ylabel = 'Power (dB)', identifier = '',
    #                                 HWattenuation = CAV_Attenuation,
    #                                 fig_num = 2)

    ####manual plotting example
    #fig = plt.figure(1,figsize=(13,8))
    #fig.clf()

    #ax = plt.subplot(1,2,1)
    #VNAplots.general_colormap_subplot(ax,freqs, powers-CAV_Attenuation, mags)
    #plt.xlabel('Frequency (GHz)')
    #plt.ylabel('Power (dB)')
    #plt.title('S21 mag')

    #ax = plt.subplot(1,2,2)
    #VNAplots.general_colormap_subplot(ax,freqs, powers-CAV_Attenuation, phases)
    #plt.xlabel('Frequency (GHz)')
    #plt.ylabel('Power (dB)')
    #plt.title('S21 phase')

    #plt.suptitle('Filename: {}'.format(scanname))
    #plt.show()

    userfuncs.SaveFull(saveDir, filename, ['mags', 'phases', 'freqs', 'powers'], locals(), expsettings=settings)
    plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)