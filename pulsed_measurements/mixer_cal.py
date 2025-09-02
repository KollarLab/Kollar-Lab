'''
8-25-21 AK modifying to normalize the amplitudes to the drive power.

8-25-21 AK modifying to return the raw data.

'''


import os
import time
import numpy as np 
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, generate_filter, remove_IQ_ellipse
from utility.ellipse_fitting import fitEllipse, make_elipse
from utility.scheduler import scheduler

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'mixer_cal'
    settings['meas_type'] = 'mixer_cal'
    
    #
    settings['CAV_Freq'] = 5e9
    settings['CAV_Power'] = -60
    
    settings['plot'] = False
    #Card settings
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 5e3
    
    return settings

def configure_sequence(hdawg, exp_globals):
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    
    awg_sched = scheduler(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='Qubit_I')
    awg_sched.add_analog_channel(2, name='Qubit_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=55e-9, HW_offset_off=0e-9)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    cavity_marker.add_window(start_time, start_time+window_time)
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)

def extract_data_local(raw, card, xaxis, measurement_pulse):
    init_buffer = measurement_pulse['init_buffer']
    emp_delay   = measurement_pulse['emp_delay']
    meas_window = measurement_pulse['meas_window']
    
    timestep   = 1/card.sampleRate
    data_start = int((init_buffer + emp_delay+0.1*meas_window)/timestep)
    window_width = int(0.85*meas_window/timestep)
    data    = raw[data_start:data_start+window_width]
    return data, xaxis[data_start:data_start+window_width]

def mixer_cal(instruments, settings, ax=None):
    ##Instruments used
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 

    CAV_Attenuation = exp_globals['CAV_Attenuation']
    
    ## Generator settings
    cavitygen.freq   = exp_settings['CAV_Freq']
    cavitygen.power  = exp_settings['CAV_Power']+CAV_Attenuation
    
    cavitygen.enable_pulse()
    cavitygen.output = 'On'

    LO.power  = 12
    LO.freq   = exp_settings['CAV_Freq'] - exp_globals['IF']
    LO.output = 'On'
    
    ##Card settings
    configure_card(card, settings)
    generate_filter(card, settings)

    configure_sequence(hdawg, exp_globals)
    
    card.ArmAndWait()
    xaxis = np.linspace(0, card.samples, card.samples)/card.sampleRate
    
    I, Q = card.ReadAllData()
    I_net = np.mean(I,0)*1e3
    Q_net = np.mean(Q,0)*1e3
    
    I_net = I_net-np.mean(I_net)
    Q_net = Q_net-np.mean(Q_net)
    m_pulse = exp_globals['measurement_pulse']
    I_win, xI = extract_data_local(I_net, card, xaxis, m_pulse)
    Q_win, xQ = extract_data_local(Q_net, card, xaxis, m_pulse)
    
    try:
        axes, center, phi = fitEllipse(I_win, Q_win)
    except:
        print('Fit failed, returning neutral')
        axes = [1,1]
        center = [0,0]
        phi = 0
        
    mixer_config = {}
    mixer_config['axes'] = axes
    mixer_config['center'] = center
    mixer_config['phi'] = phi
    I_corr, Q_corr = remove_IQ_ellipse(I_win, Q_win, mixer_config)
    
    if exp_settings['plot']:
        plt.figure(211, figsize=(15,5))
        plt.clf()
        plt.subplot(131)
        xx,yy = make_elipse(axes, center, phi, 101)
        plt.plot(I_win, Q_win, 'b',label='Data')
        plt.plot(xx,yy, 'orange',label='Fit')
        plt.gca().set_aspect('equal')
        plt.title('Ellipse fit')
        plt.legend(loc='upper right')
        plt.subplot(132)
        plt.plot(xaxis*1e6, I_net, label='I')
        plt.plot(xaxis*1e6, Q_net, label='Q')
        plt.plot(xI*1e6, I_win, label='window')
        plt.title('Raw data')
        plt.xlabel('Time (us)')
        plt.ylabel('Amp (mV)')
        plt.legend()
        plt.subplot(133)
        plt.plot(xaxis, np.sqrt(I_net**2+Q_net**2))
        plt.plot(xI, np.sqrt(I_corr**2+Q_corr**2))
    if ax:
        xx,yy = make_elipse(axes, center, phi, 101)
        ax.plot(I_win, Q_win, 'b',label='Data')
        ax.plot(xx,yy, 'orange',label='Fit')
        ax.set_aspect('equal')
        ax.set_title('Ellipse fit')
        ax.legend(loc='upper right')
    return mixer_config