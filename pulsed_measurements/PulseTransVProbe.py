# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 19:12:43 2020

@author: Kollarlab
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, estimate_time, read_and_process

def get_default_settings():
    settings = {}
    
    settings['scanname']  = 'scanname'
    settings['meas_type'] = 'PulsedTransVProbe'
    
    #Sweep parameters
    settings['CAV_Power']   = -60
    settings['Q_Freq']      = 3.2*1e9

    settings['start_freq']  = 4.15*1e9  
    settings['stop_freq']   = 4.25*1e9 
    settings['freq_points'] = 50

    settings['start_power']  = -20
    settings['stop_power']   = 10
    settings['power_points'] = 31

    #Card settings
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 5e3
    
    return settings

def pulsed_trans_v_probe(instruments, settings):
    ##Instruments used
    cavitygen = instruments['cavitygen']
    qubitgen  = instruments['qubitgen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 

    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    CAV_Attenuation  = exp_globals['CAV_Attenuation']

    ##Data saving and naming
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = settings['scanname'] + '_' + stamp
    
    ##Sweep settings
    start_freq  = exp_settings['start_freq'] 
    stop_freq   = exp_settings['stop_freq'] 
    freq_points = exp_settings['freq_points']
    freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
 
    start_power  = exp_settings['start_power'] + Qbit_Attenuation
    stop_power   = exp_settings['stop_power']  + Qbit_Attenuation
    power_points = exp_settings['power_points']
    powers = np.round(np.linspace(start_power,stop_power,power_points),2)
    
    ## Generator settings
    qubitgen.freq   = settings['Q_Freq']
    qubitgen.power  = powers[0]
    qubitgen.IQ.Mod = 'Off'
    
    
    cavitygen.Freq   = freqs[0]
    cavitygen.Power  = exp_settings['cavity_power'] + CAV_Attenuation
    cavitygen.IQ.Mod = 'On'

    cavitygen.Output = 'On'
    qubitgen.Output  = 'On'
    
    LO.Power  = 12
    LO.Freq   = freqs[0] - exp_globals['IF']
    LO.Output = 'On'
    
    ##Card settings
    configure_card(card, settings)

    ##HDAWG settings
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='True')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\Control\HDAWG_sequencer_codes\pulsedtrans.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse  = exp_globals['measurement_pulse'] 
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))
    
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)

    powerdat = np.zeros((len(powers), len(freqs)))
    phasedat = np.zeros((len(powers), len(freqs)))
    
    tstart   = time.time()
    first_it = True
    
    for powerind in range(len(powers)):
        qubitgen.Power = powers[powerind]
        time.sleep(0.2)
        
        amps   = np.zeros((len(freqs), card.samples))
        phases = np.zeros((len(freqs), card.samples))
        
        print('Current power:{}, max:{}'.format(powers[powerind] - Qbit_Attenuation, powers[-1]))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]
            
            cavitygen.freq = freq

            LO.freq   = freq - exp_globals['IF']
            LO.output = 'On'
            time.sleep(0.2)
            
            amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot=first_it)

            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, len(powers)*len(freqs))
            
            amps[find,:]   = amp_full
            phases[find,:] = phase_full
            powerdat[powerind, find] = np.mean(amp)
            phasedat[powerind, find] = np.mean(phase)
            
        full_data = {}
        full_data['xaxis']  = freqs/1e9
        full_data['mags']   = powerdat[0:powerind+1]
        full_data['phases'] = phasedat[0:powerind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mag']   = powerdat[powerind]
        single_data['phase'] = phasedat[powerind]

        yaxis  = powers[0:powerind+1] - Qbit_Attenuation
        labels = ['Freq (GHz)', 'Power (dBm)']
        simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) 
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

        full_time = {}
        full_time['xaxis']  = xaxis*1e6
        full_time['mags']   = amps
        full_time['phases'] = phases

        single_time = {}
        single_time['xaxis'] = xaxis*1e6
        single_time['mag']   = amp
        single_time['phase'] = phase

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier  = 'Power: {}dBm'.format(powers[powerind]-Qbit_Attenuation)
        
        simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces\n'+filename, time_labels, identifier, fig_num=2)
        plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)

        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis', 'full_data', 'single_data', 'full_time', 'single_time'], 
                            locals(), expsettings=settings, instruments=instruments)
        
    t2 = time.time()
    print('elapsed time = ' + str(t2-tstart))