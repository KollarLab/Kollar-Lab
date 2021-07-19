# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 17:37:06 2021

@author: Kollarlab
"""

import time
import os
from utility.measurement_helpers import estimate_time
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
import utility.plotting_tools as plots
from utility.userfits import fit_lorentzian

def get_default_settings():
    settings = {}
    
    #Save location
    settings['scanname']    = 'initial_power_scan_q4'
    settings['meas_type']   = 'stability'
    #settings['project_dir'] = r'Z:\Data'

    #Sweep parameters
    #settings['CAV_Attenuation'] = 30
    #settings['Qbit_Attenuation'] = 10

    settings['num_points'] = 10

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
    settings['ifBW'] = 2e2
    settings['mode'] = 'MOV'
    
    return settings

def vna_spec_stability(instruments, settings):
    #Instruments used
    vna = instruments['VNA']
    
    vna.reset()

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    
    #Data saving and naming
    #saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = settings['scanname'] + '_' + stamp

    CAV_Attenuation  = exp_globals['CAV_Attenuation']
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']

    exp_settings['CAVpower'] = exp_settings['CAVpower'] + CAV_Attenuation
    exp_settings['RFpower']  = exp_settings['RFpower'] + Qbit_Attenuation
    
    timing = np.linspace(0, exp_settings['num_points']*exp_settings['avg_time'], exp_settings['num_points'])

    mags   = np.zeros((exp_settings['num_points'], exp_settings['freq_points']))
    phases = np.zeros((exp_settings['num_points'], exp_settings['freq_points']))
    
    centers = np.zeros(exp_settings['num_points'])
    widths  = np.zeros(exp_settings['num_points'])
    amps    = np.zeros(exp_settings['num_points'])
    
    tstart = time.time()
    for tind in range(len(timing)):
        print('Measurement {} out of {}'.format(tind, exp_settings['num_points']))

        data = vna.spec_meas(exp_settings)

        vna.autoscale()

        mags[tind]   = data['mag']
        phases[tind] = data['phase']

        if tind==0:
            tstop = time.time()
            estimate_time(tstart, tstop, len(timing))

        freqs = data['xaxis']
    
        full_data = {}
        full_data['xaxis'] = freqs
        full_data['mags'] = mags[0:tind+1]
        full_data['phases'] = phases[0:tind+1]
    
        single_data = data
        yaxis = timing[0:tind+1]
        labels = ['Freq (GHz)', 'Time elapsed (s)']
    
        plots.simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1)
        plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
        center, width, amp = fit_lorentzian(freqs, data['mag'], plot=True, fig_num=4)
        
        centers[tind] = center
        widths[tind]  = width
        amps[tind]    = amp
        
        fig = plt.figure(2, figsize=(13,8))
        plt.clf()
        plt.subplot(311)
        plt.title('Center vs time (GHz)')
        plt.plot(timing[0:tind+1], centers[0:tind+1]/1e9)
        plt.subplot(312)
        plt.title('FWHM vs time (MHz)')
        plt.plot(timing[0:tind+1], widths[0:tind+1]/1e6)
        plt.subplot(313)
        plt.title('Contrast (dB)')
        plt.plot(timing[0:tind+1], amps[0:tind+1])
        
        plt.suptitle('Filename:{}'.format(filename))
        
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'fits.png'), dpi = 150)
        
        userfuncs.SaveFull(saveDir, filename, 
                           ['full_data', 'single_data', 'timing', 'labels', 'filename', 
                            'amps', 'centers', 'widths'], 
                           locals(), expsettings=settings, instruments=instruments)
        
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))

    