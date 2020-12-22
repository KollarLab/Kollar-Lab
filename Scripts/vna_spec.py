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
import VNAplottingTools as VNAplots

def get_default_settings():
    settings = {}
    
    #Save location
    settings['scanname']    = 'initial_power_scan_q4'
    settings['meas_type']   = 'Spec'
    settings['project_dir'] = r'Z:\Data\HouckQuadTransmon'

    #Sweep parameters
    settings['CAV_Attenuation'] = 30
    settings['Qbit_Attenuation'] = 10

    settings['start_power']  = -20
    settings['stop_power']   = 10
    settings['power_points'] = 31

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

    return settings

def vna_spec(instruments, settings):
    #Instruments used
    vna = instruments['VNA']

    #Data saving and naming
    saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp = userfuncs.timestamp()
    filename = settings['scanname'] + '_' + stamp

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

        data = vna.spec_meas(settings)

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

    VNAplots.power_plot(freqs, mags, phases, powers, filename, settings['CAV_Attenuation'])

    userfuncs.SaveFull(saveDir, filename, ['mags', 'phases', 'freqs', 'powers'], locals(), expsettings=settings)
    plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)