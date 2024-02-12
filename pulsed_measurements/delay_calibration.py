# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 18:02:57 2024

@author: kollarlab
"""
import os
import userfuncs
import numpy as np
import matplotlib.pyplot as plt
import time

from utility.measurement_helpers import read_and_process, configure_card, generate_filter
from utility.scheduler import scheduler

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'delay_calibration'
    settings['meas_type'] = 'topo_tuneup'
    
    #Cavity parameters
    settings['CAV_Power']    = -45
    settings['CAV_Freq']    = 8e9
    settings['Q_Power']     = -20
    settings['Q_Freq']      = 5e9

    #Hold time 
    settings['start_time']  = 10e-9
    settings['stop_time']   = 500e-9
    settings['time_points'] = 40

    #Card settings
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 5e3
    
    #Pulse settings
    settings['Z_pulse_buffer'] = 0
    settings['ramp_sigma'] = 20e-9
    settings['flux_length'] = 200e-9
    settings['flux_offset'] = 0
    settings['flux_amp'] = 0
    
    settings['X_pulse_buffer'] = 0
    
    return settings
    
def X_Z_delay_cal(instruments, settings):
    
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp
    
    ##Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = exp_settings['CAV_Power'] + CAV_Attenuation
    CAV_freq  = exp_settings['CAV_Freq']
    
    ##Qubit settings
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power   = exp_settings['Q_Power'] + Qbit_Attenuation
    Qbit_freq    = exp_settings['Q_Freq']
    
    start_time   = exp_settings['start_time']
    stop_time    = exp_settings['stop_time']
    time_points  = exp_settings['time_points']
    times = np.round(np.linspace(start_time,stop_time,time_points), 9)
    
    ## Generator settings
    cavitygen.freq   = CAV_freq
    cavitygen.power  = CAV_power
    cavitygen.enable_pulse()

    qubitgen.freq   = Qbit_freq
    qubitgen.power  = Qbit_power
    
    qubitgen.enable_pulse()
    qubitgen.enable_IQ()

    cavitygen.output = 'On'
    qubitgen.output  = 'On'
    
    LO.power  = 12
    LO.freq   = CAV_freq - exp_globals['IF']
    LO.output = 'On'
    
    ##Card settings
    configure_card(card, settings)
    
    ##Schedule building
    progFile = open(r"C:\Users\kollarlab\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
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
    awg_sched.add_analog_channel(3, name='Flux_pulse')
    awg_sched.add_analog_channel(4, name='blank')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=50e-9, HW_offset_off=50e-9)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(3, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    
    qubit_I       = awg_sched.analog_channels['Qubit_I']
    qubit_Q       = awg_sched.analog_channels['Qubit_Q']
    Flux_pulse    = awg_sched.analog_channels['Flux_pulse']
    qubit_marker  = awg_sched.digital_channels['Qubit_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']

    cavity_marker.add_window(start_time, start_time+window_time+1e-6)
    
    #Generate flux pulse that we're going to try to find with the x pulse
    buffer = exp_settings['Z_pulse_buffer']
    ramp_len = 2*exp_settings['ramp_sigma']
    Bz_pos = buffer+2*ramp_len+exp_settings['flux_length']
    Flux_pulse.add_pulse('gaussian_square', 
                         position=start_time-Bz_pos-exp_settings['flux_offset'], 
                         amplitude=exp_settings['flux_amp'],
                         length = exp_settings['flux_length'], 
                         ramp_sigma=exp_settings['ramp_sigma'], 
                         num_sigma=4) 

    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
        
    ## Starting main measurement loop 

    generate_filter(card, settings)

    first_it = True
    timedat = np.zeros(len(times))
    phasedat = np.zeros(len(times))
    
    t1 = time.time()
    for timeind in range(len(times)):
        hold_time = times[timeind]
        
        #Update schedule
        hdawg.AWGs[0].stop()
        qubit_I.reset()
        qubit_marker.reset()
        
        qubit_time = num_sigma*sigma+hold_time
        position = start_time-qubit_time-q_pulse['delay']-exp_settings['X_pulse_buffer']
        
        qubit_I.add_pulse('gaussian_square', 
                          position=position, 
                          amplitude=q_pulse['piAmp'], 
                          length = hold_time, 
                          ramp_sigma=q_pulse['sigma'], 
                          num_sigma=q_pulse['num_sigma'])
        
        qubit_marker.add_window(position, position+qubit_time)
        awg_sched.plot_waveforms()
        
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Flux_pulse', 'blank'], ['blank1', 'blank2'])

        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
        hdawg.AWGs[0].run_loop()
        
            
        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                     plot=first_it, 
                                                                     IQstorage = True)
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
            I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
            qubitgen.output='On'
        else:
            I_window_b, Q_window_b, I_full_b, Q_full_b = 0,0,0,0
        
        if first_it:
            first_it=False
        ##Useful handles for variables
        I_sig, Q_sig   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
        I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)] #<I>, <Q> for background trace
        theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
        theta_back = np.arctan2(Q_back, I_back)*180/np.pi #angle relative to x axis in IQ plane 
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back #compute <Q_net> in the data window
        I_full_net = I_full-I_full_b #full I data with background subtracted
        Q_full_net = Q_full-Q_full_b #full Q data with background subtracted
        
        amp = np.sqrt(I_final**2 + Q_final**2)
        phase = np.arctan2(Q_final, I_final)*180/np.pi
        amp_full = np.sqrt(I_full_net**2+Q_full_net**2)  
        phase_full = np.arctan2(Q_full_net, I_full_net)*180/np.pi             

        timedat[timeind]  = np.mean(amp)
        phasedat[timeind] = np.mean(phase)
        
        fig = plt.figure(62)
        plt.clf()
        plt.subplot(121)
        plt.plot(times[:timeind], timedat[:timeind])
        plt.xlabel('Hold time (ns)')
        plt.ylabel('Mag (arb.)')
        plt.title('Mag')
        plt.subplot(122)
        plt.plot(times[:timeind], phasedat[:timeind])
        plt.xlabel('Hold time (ns)')
        plt.ylabel('Ang (deg.)')
        plt.title('Ang')
        fig.canvas.draw()
        fig.canvas.flush_events()
        userfuncs.SaveFull(saveDir, filename, ['times','timedat', 'phasedat','xaxis'], 
                        locals(), expsettings=settings, instruments=instruments, saveHWsettings=first_it)
    t2 = time.time()
    plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
    print('elapsed time = ' + str(t2-t1))

    cavitygen.Output = 'Off'
    qubitgen.Output  = 'Off'
    LO.Output        = 'Off'
    full_data = {}
    full_data['amps'] = timedat
    full_data['phases'] = phasedat
    full_data['time_ax'] = xaxis
    
    return full_data