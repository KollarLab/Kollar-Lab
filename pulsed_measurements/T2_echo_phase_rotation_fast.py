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
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import configure_card, generate_filter, estimate_time, read_and_process
from utility.scheduler import scheduler

import scipy.signal as signal

def GetDefaultSettings():
    settings = {}
    
    settings['scanname'] = 'T2_meas'
    settings['meas_type'] = 'Tmeas'

    settings['Q_Freq']  = 4.20431e9
    settings['Q_Power'] = -11

    settings['CAV_Freq']  = 8.126e9
    settings['CAV_Power'] = -18
    
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 25e3
    
    settings['Tau_min'] = 200e-9
    settings['Tau_max'] = 30e-6
    settings['Tau_points'] = 5
    settings['pulse_count'] = 1
    settings['phase_rotation_f'] = 1e6
    settings['detuning'] = 1e6
    settings['T2_mode'] = 'phase_rotation'
    
    settings['verbose'] = True
    
    return settings

def meas_T2_phase_rotation(instruments, settings):

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
    cavitygen.freq   = CAV_Freq
    cavitygen.power  = CAV_Power + CAV_Attenuation
    cavitygen.enable_pulse()
    
    if exp_settings['T2_mode'] == 'detuning':
        qubitgen.freq   = Q_Freq + exp_settings['detuning']
    elif exp_settings['T2_mode'] == 'phase_rotation':
        qubitgen.freq   = Q_Freq
    else:
        raise ValueError('Invalid T2_mode')
    qubitgen.power  = Q_Power + Qbit_Attenuation
    qubitgen.enable_IQ()
    qubitgen.enable_pulse()
    
    LO.power = 12
    LO.freq = cavitygen.freq-exp_globals['IF']
    LO.output = 'On'
    
    cavitygen.Output = 'On'
    qubitgen.Output  = 'On'
    
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
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    
    qubit_I       = awg_sched.analog_channels['Qubit_I']
    qubit_Q       = awg_sched.analog_channels['Qubit_Q']
    qubit_marker  = awg_sched.digital_channels['Qubit_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    hold_time = q_pulse['hold_time']

    cavity_marker.add_window(start_time, start_time+window_time+1e-6)

    taus = np.linspace(exp_settings['Tau_min'],exp_settings['Tau_max'],exp_settings['Tau_points'])
    taus = np.round(taus, 9)
    
    ## Start main measurement loop
    amp_int = np.zeros(len(taus))
    ang_int = np.zeros(len(taus))
    amps    = np.zeros((len(taus),card.samples))
    angles  = np.zeros(amps.shape)

    fig1 = plt.figure(1,figsize=(13,8))
    fig1.clf()
    fig1.suptitle('Live T2 data (no fit), {} pi pulses'.format(exp_settings['pulse_count']))
    ax11 = fig1.add_subplot(121)
    ax12 = fig1.add_subplot(122)

    fig2 = plt.figure(2,figsize=(13,8))
    fig2.clf()
    ax21 = fig2.add_subplot(111)

    t_start = time.time() 
    first_it = True
    for tind in range(len(taus)):
        #Generating pulse sequence and upload
        
        tau = taus[tind]
        if exp_settings['verbose']:
            print('Tau: {}'.format(tau))
        hdawg.AWGs[0].stop()
        qubit_I.reset()
        qubit_Q.reset()
        qubit_marker.reset()
        
        position = start_time-delay-num_sigma*sigma
        qubit_time = num_sigma*sigma
        
        #the main pulses    
        qubit_I.add_pulse('gaussian_square', position=position-tau-hold_time, 
                          amplitude=q_pulse['piAmp']/2,
                          length = hold_time,
                          ramp_sigma=q_pulse['sigma'], 
                          num_sigma=q_pulse['num_sigma'])

        if exp_settings['T2_mode'] == 'detuning':
            if exp_settings['basis'] =='X':
                qubit_I.add_pulse('gaussian_square', position=position-hold_time, 
                                  amplitude=q_pulse['piAmp']/2,
                                  length = hold_time,
                                  ramp_sigma=q_pulse['sigma'], 
                                  num_sigma=q_pulse['num_sigma'])
            elif exp_settings['basis'] =='Y':
                qubit_Q.add_pulse('gaussian_square', position=position-hold_time, 
                                  amplitude=q_pulse['piAmp']/2,
                                  length = hold_time,
                                  ramp_sigma=q_pulse['sigma'], 
                                  num_sigma=q_pulse['num_sigma'])
            else:
                print(exp_settings['basis'])
                print('Invalid basis selected for second pulse ("X" or "Y")')
        elif exp_settings['T2_mode'] == 'phase_rotation':
            qubit_I.add_pulse('gaussian_square', position=position-hold_time,
                              length = hold_time,
                              amplitude= np.cos(2*np.pi * tau* exp_settings['phase_rotation_f'])*q_pulse['piAmp']/2, 
                              ramp_sigma=q_pulse['sigma'], 
                              num_sigma=q_pulse['num_sigma'])
            qubit_Q.add_pulse('gaussian_square', position=position-hold_time, 
                              length = hold_time, 
                              amplitude= np.sin(2*np.pi * tau* exp_settings['phase_rotation_f'])*q_pulse['piAmp']/2, 
                              ramp_sigma=q_pulse['sigma'], 
                              num_sigma=q_pulse['num_sigma'])
        else:
            raise ValueError('Invalid T2_mode')
        
        qubit_marker.add_window(position-tau-100e-9, position-tau+100e-9)
        qubit_marker.add_window(position-100e-9, position+100e-9)
    
        
        if exp_settings['pulse_count'] > 0:
            numPulses = exp_settings['pulse_count']
            CPMG_times = tau*(np.linspace(1, numPulses, numPulses)-0.5)/numPulses
            temp = np.linspace(0,tau, numPulses + 2)
            pulseTimes = CPMG_times#temp[1:-1]
            
            for tp in pulseTimes:
                qubit_I.add_pulse('gaussian_square', position=position-tp-hold_time, 
                                  length = hold_time,
                                  amplitude=q_pulse['piAmp'], 
                                  ramp_sigma=q_pulse['sigma'], 
                                  num_sigma=q_pulse['num_sigma'])
                qubit_marker.add_window(position-tp-100e-9, position-tp+100e-9)
      
        awg_sched.plot_waveforms()
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[0].run_loop()
        qubitgen.output='On'

        #Acquiring data
        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
        
        #Background collection
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
            I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
            qubitgen.freq=Q_Freq + exp_settings['detuning']
        else:
            I_window_b, Q_window_b, I_full_b, Q_full_b = 0,0,0,0
        
        #Plotting and saving 
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
        
        ax11.cla()
        ax12.cla()
        
        ax11.plot(taus*1e6, amp_int)
        ax11.set_xlabel('Tau (us)')
        ax11.set_ylabel('Amplitude')

        ax12.plot(taus*1e6, ang_int)
        ax12.set_ylabel('Angle')
        ax12.set_xlabel('Tau (us)')
             
        fig1.canvas.draw()
        fig1.canvas.flush_events()            

        first_it = False
    

    


    # fig2 = plt.figure(48,figsize=(13,8))
    # plt.clf()
    # ax21 = fig2.add_subplot(111)
    # # # plt.figure(fig2.number)

    # plt.show()  

    # #hail mary plot of everything at the end
    # fig47 = plt.figure(47)
    # ax11 = fig47.add_subplot(121)
    # ax12 = fig47.add_subplot(122)
    # ax11.plot(taus*1e6, amp_int)
    # ax11.set_xlabel('Tau (us)')
    # ax11.set_ylabel('Amplitude')

    # ax12.plot(taus*1e6, ang_int)
    # ax12.set_ylabel('Angle')
    # ax12.set_xlabel('Tau (us)')

    # plt.suptitle('Hail MAry Figure')

    # fig47.canvas.draw()
    # fig47.canvas.flush_events()   
    # #########




    ax21.cla()

    cbar = general_colormap_subplot(ax21, xaxis*1e6, 
                                    taus*1e6, amps, 
                                    ['Time (us)', 'Tau (us)'], 
                                    'Raw data\n'+filename)
    fig1.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 250)
    fig2.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 250)
    userfuncs.SaveFull(saveDir, filename, ['taus','xaxis', 'amps', 'amp_int'],
                         locals(), expsettings=settings, instruments=instruments)
    
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-t_start))
    return taus, amp_int