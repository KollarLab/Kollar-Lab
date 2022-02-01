# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 18:36:21 2021

@author: Kollarlab
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
from utility.userfits import fit_T2
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process

import scipy.signal as signal

def GetDefaultSettings():
    settings = {}
    
    settings['scanname'] = 'pulse_calibration'
    settings['meas_type'] = 'Tmeas'

    settings['Q_Freq']  = 4.20431e9
    settings['Q_Power'] = -11

    settings['CAV_Freq']  = 8.126e9
    settings['CAV_Power'] = -18
    
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 25e3
    
    settings['pulse_count'] = 10
    
    return settings

def pulse_calibration(instruments, settings):
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 

    Q_Freq    = exp_settings['Q_Freq']
    Q_Power   = exp_settings['Q_Power']
    CAV_Freq  = exp_settings['CAV_Freq']
    CAV_Power = exp_settings['CAV_Power']
    
    CAV_Attenuation  = exp_globals['CAV_Attenuation']
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp
    
    ## Configure generators
    cavitygen.Freq   = CAV_Freq
    cavitygen.Power  = CAV_Power + CAV_Attenuation
    cavitygen.IQ.Mod = 'On'
    
    qubitgen.Freq   = Q_Freq
    qubitgen.Power  = Q_Power + Qbit_Attenuation
    qubitgen.IQ.Mod = 'On'
    
    LO.power = 12
    LO.freq = '{} GHz'.format((cavitygen.Freq-exp_globals['IF'])/1e9)
    LO.output = 'On'
    
    cavitygen.Output = 'On'
    qubitgen.Output  = 'On'
    
    ## Configure card
    configure_card(card, settings)
    
    ## Configure HDAWG
    configure_hdawg(hdawg, settings)
    #hard coded the delay between the two AWG cores
    hdawg.Channels[0].configureChannel(amp=0.5,marker_out='Marker', hold='False', delay=30e-9)
    hdawg.Channels[1].configureChannel(amp=0.5,marker_out='Marker', hold='False', delay=30e-9)
    
    HDAWG_dir = r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes"
    progFile = open(os.path.join(HDAWG_dir,'pulse_calibration.cpp'),'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()

    q_pulse = exp_globals['qubit_pulse']
    m_pulse = exp_globals['measurement_pulse']
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']+2e-6))
    loadprog = loadprog.replace('_wait_time_',str(q_pulse['delay']))
    loadprog = loadprog.replace('_qsigma_',str(q_pulse['sigma']))
    loadprog = loadprog.replace('_num_sigma_',str(q_pulse['num_sigma']))
    loadprog = loadprog.replace('_piAmp_',str(q_pulse['piAmp']))
    
    ## Start main measurement loop
    pi_count = exp_settings['pulse_count']
    amp_int = np.zeros(pi_count)
    ang_int = np.zeros(pi_count)
    amps    = np.zeros((pi_count,card.samples))
    
    tstart = time.time()
    first_it = True
    
    if exp_globals['IF'] != 0:
        #create Chebychev type II digital filter
        filter_N = exp_globals['ddc_config']['order']
        filter_rs = exp_globals['ddc_config']['stop_atten']
        filter_cutoff = np.abs(exp_globals['ddc_config']['cutoff'])
        LPF = signal.cheby2(filter_N, filter_rs, filter_cutoff, btype='low', analog=False, output='sos', fs=card.sampleRate)
        
        xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
        digLO_sin = np.sin(2*np.pi*exp_globals['IF']*xaxis)
        digLO_cos = np.cos(2*np.pi*exp_globals['IF']*xaxis)
        
        #store in settings so that the processing functions can get to them
        settings['digLO_sin'] = digLO_sin 
        settings['digLO_cos'] = digLO_cos
        settings['LPF'] = LPF
        
    for num_pulses in range(pi_count):
        
        print('{} pulses'.format(2*num_pulses+1))
        finalprog = loadprog
        finalprog = finalprog.replace('_pi_count_',str(2*num_pulses+1))
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)

        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, plot=first_it, IQstorage=True)    
  
        amps[num_pulses] = np.sqrt((I_full**2+Q_full**2))
        amp_int[num_pulses] = np.sqrt(np.mean(I_window)**2+np.mean(Q_window)**2)
        ang_int[num_pulses] = np.arctan2(np.mean(Q_window),np.mean(I_window))

        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, pi_count)
            first_it = False
        
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(amp_int)
        plt.suptitle('Live data')
        plt.xlabel('Pulse count (odd repeats)')
        plt.ylabel('Amplitude')
        plt.subplot(122)
        plt.plot(ang_int)
        plt.ylabel('Angle')
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)

        fig2 = plt.figure(2,figsize=(13,8))
        plt.clf()

        ax = plt.subplot(1,1,1)
        general_colormap_subplot(ax, xaxis*1e6, list(range(pi_count)), amps, ['Time (us)', 'Pulse count (odd repeats)'], 'Raw data\n'+filename)

        plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
        userfuncs.SaveFull(saveDir, filename, ['xaxis', 'amps', 'amp_int'], locals(), expsettings=settings, instruments=instruments)
    
    t2 = time.time()
    print('Elapsed time:{}'.format(t2-tstart))
    return amp_int, ang_int