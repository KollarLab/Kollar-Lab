# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:52:18 2020

@author: Kollarlab
"""
import os
import time
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
from utility.userfits import fit_T1
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import configure_card, generate_filter, estimate_time, read_and_process
from utility.scheduler import scheduler
import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'T1_meas'
    settings['meas_type'] = 'Tmeas'
    
    settings['Q_Freq']  = 4.21109e9
    settings['Q_Power'] = -18

    settings['CAV_Freq']  = 8.12357e9
    settings['CAV_Power'] = -42.9
    
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 25e3
    
    settings['Tau_min']    = 200e-9
    settings['Tau_max']    = 30e-6
    settings['Tau_points'] = 5
    settings['spacing']    = 'Linear'
    
    settings['T1_guess'] = 10e-6
    
    return settings
    
def hold_time_sweep(instruments, settings):
    ## Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    ## Bookkeeping and setting up the save directory
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
    LO.power  = 12
    LO.freq = CAV_Freq-exp_globals['IF']
    LO.output = 'On'
    
    cavitygen.freq   = CAV_Freq
    cavitygen.power  = CAV_Power + CAV_Attenuation
    cavitygen.enable_pulse()
    
    qubitgen.Freq   = Q_Freq
    qubitgen.Power  = Q_Power + Qbit_Attenuation
    qubitgen.enable_IQ()
    qubitgen.enable_pulse()
    
    cavitygen.Output = 'On'
    qubitgen.Output  = 'On'
    
    cavitygen.output = 'On'
    qubitgen.output  = 'On'
    
    ## Configure card
    configure_card(card, settings)
    generate_filter(card, settings)

    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()

    m_pulse = exp_globals['measurement_pulse']
    q_pulse = exp_globals['qubit_pulse']
    
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    
    awg_sched = scheduler(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='Qubit_I')
    awg_sched.add_analog_channel(2, name='Qubit_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    
    qubit_I       = awg_sched.analog_channels['Qubit_I']
    qubit_Q       = awg_sched.analog_channels['Qubit_Q']
    qubit_marker  = awg_sched.digital_channels['Qubit_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']

    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    
    cavity_marker.add_window(start_time, start_time+window_time)

    ## Set up array of taus and randomize it
    if exp_settings['spacing']=='Log':
        tau_list = np.logspace(np.log10(exp_settings['Tau_min']), np.log10(exp_settings['Tau_max']), exp_settings['Tau_points'])
    else:
         tau_list = np.linspace(exp_settings['Tau_min'], exp_settings['Tau_max'], exp_settings['Tau_points'])
    taus = np.round(tau_list, 9)
    
    indices = list(range(len(taus)))

    ## Start the main measurement loop 
    total_samples = card.samples
    
    amp_int = np.zeros(len(taus))
    ang_int = np.zeros(len(taus))
    amps    = np.zeros((len(taus),total_samples))
    angles  = np.zeros((len(taus),total_samples))
        
    tstart = time.time()
    first_it = True

    for tind in indices:
        
        tau = taus[tind]
        print('Hold: {}'.format(tau))
        
        hdawg.AWGs[0].stop()
        qubit_I.reset()
        qubit_marker.reset()
        
        position = start_time-delay-num_sigma*sigma
        
        qubit_I.add_pulse('gaussian_square', position=position-tau, 
                          amplitude=q_pulse['piAmp'], length = tau, 
                          ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
        qubit_marker.add_window(position-tau-150e-9, position+num_sigma*sigma+150e-9)
        if tind==0 or tind==len(indices)-1:
            awg_sched.plot_waveforms()
        
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        
        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[0].run_loop()
        qubitgen.output='On'
        
        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
            I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
            qubitgen.freq=Q_Freq
        else:
            I_window_b, Q_window_b, I_full_b, Q_full_b = 0,0,0,0
        
        ##Useful handles for variables
        I_sig, Q_sig   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
        I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)] #<I>, <Q> for background trace
        theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
        theta_back = np.arctan2(Q_back, I_back)*180/np.pi #angle relative to x axis in IQ plane 
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back #compute <Q_net> in the data window
        
        amps[tind] = np.sqrt((I_full-I_full_b)**2+(Q_full-Q_full_b)**2)
        angles[tind] = np.arctan2((Q_full-Q_full_b), (I_full-I_full_b))*180/np.pi
        
        amp_int[tind] = np.sqrt(I_final**2+Q_final**2)
        ang_int[tind] = np.arctan2(Q_final, I_final)*180/np.pi
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(taus))
            first_it = False      

        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(taus*1e6, amp_int, 'x')
        plt.xlabel('Hold time (us)')
        plt.ylabel('Amplitude')  
        plt.subplot(122)
        plt.plot(taus*1e6, ang_int, 'x')
        plt.xlabel('Hold time (us)')
        plt.ylabel('Phase')  
        plt.title('Live data (no fit)\n'+filename)
        fig.canvas.draw()
        fig.canvas.flush_events()
        

    plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)
    fig2 = plt.figure(2,figsize=(13,8))
    plt.clf()

    ax = plt.subplot(1,1,1)
    general_colormap_subplot(ax, xaxis*1e6, taus*1e6, amps, ['Time (us)', 'Hold time (us)'], 'Raw data\n'+filename)

    plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
    userfuncs.SaveFull(saveDir, filename, ['taus','xaxis', 'amps', 'amp_int'],
                         locals(), expsettings=settings, instruments=instruments)
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))

    return taus, amp_int