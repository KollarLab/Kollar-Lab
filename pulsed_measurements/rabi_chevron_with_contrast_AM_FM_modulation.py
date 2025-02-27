# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:28:57 2021

@author: Kollarlab
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process
from utility.scheduler import scheduler

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'rabi_chevron'
    settings['meas_type'] = 'rabi_chevron'
    
    #Measurement params (from qubit calibration)
    settings['CAV_Power']  = -45
    settings['CAV_Freq']   = 5e9
    settings['Q_Power'] = 0
    settings['Q_Freq']  = 4e9
    
    settings['subtract_background'] = True
    settings['back_rate']  = 1 #retake background every n measurements
    #Sweep params
    settings['Drive_Power'] = -20
    
    settings['start_freq']  = 4*1e9  
    settings['stop_freq']   = 5*1e9 
    settings['freq_points'] = 50

    settings['start_time']  = 10e-9
    settings['stop_time']   = 500e-9
    settings['time_points'] = 40

    #Card settings
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 5e3
    
    settings['FM_range'] = 0
    settings['FM_freq'] = 0
    return settings

def rabi_chevron_AM_FM(instruments, settings):
    
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
    start_freq  = exp_settings['start_freq']
    stop_freq   = exp_settings['stop_freq']
    freq_points = exp_settings['freq_points']
    freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
    
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    drive_power   = exp_settings['Drive_Power'] + Qbit_Attenuation
    qubit_power   = exp_settings['Q_Power'] + Qbit_Attenuation
    start_time   = exp_settings['start_time']
    stop_time    = exp_settings['stop_time']
    time_points  = exp_settings['time_points']
    times = np.round(np.linspace(start_time,stop_time,time_points), 9)
    
    ## Generator settings
    cavitygen.freq   = CAV_freq
    cavitygen.power  = CAV_power
    cavitygen.enable_pulse()

    qubitgen.freq   = 4e9
    qubitgen.power  = drive_power
    
    qubitgen.enable_pulse()
    qubitgen.enable_IQ()

    cavitygen.output = 'On'
    qubitgen.output  = 'On'
    
    LO.power  = 12
    LO.freq   = CAV_freq - exp_globals['IF']
    LO.output = 'On'
    
    ##Card settings
    configure_card(card, settings)

    ##HDAWG settings
    configure_hdawg(hdawg, settings)
    
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

    cavity_marker.add_window(start_time, start_time+window_time+1e-6)

    ## Starting main measurement loop 
    timedat  = np.zeros((len(times), len(freqs)))
    phasedat = np.zeros((len(times), len(freqs)))
    
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
        
    tstart = time.time()
    first_it = True
    

    
    for timeind in range(len(times)):
        hold_time = times[timeind]
        
        hdawg.AWGs[0].stop()
        qubit_I.reset()
        qubit_Q.reset()
        qubit_marker.reset()
        
        position = start_time-delay-num_sigma*sigma
        qubit_time = num_sigma*sigma
        
        AWG_x = np.linspace(0, qubit_I.samples, qubit_I.samples)/2.4e9-position
        f_mod = exp_settings['FM_freq']
        f_det = exp_settings['FM_range']
        phi = f_det*np.cos(2*np.pi*f_mod*AWG_x)/f_mod
        Imod = np.cos(phi)
        Qmod = np.sin(phi)
        mod_depth = exp_settings['AM_depth']
        mod_sig = 1+mod_depth*(np.cos(2*np.pi*f_mod*AWG_x)-1)
        
        qubit_I.add_pulse('gaussian_square', position=position-hold_time, amplitude=q_pulse['piAmp'], length = hold_time, ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
        qubit_I.wave_array = qubit_I.wave_array*Imod*mod_sig
        qubit_Q.add_pulse('gaussian_square', position=position-hold_time, amplitude=q_pulse['piAmp'], length = hold_time, ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
        qubit_Q.wave_array = qubit_Q.wave_array*Qmod*mod_sig
        
        qubit_marker.add_window(position-qubit_time-hold_time, position+2*qubit_time+hold_time)
        awg_sched.plot_waveforms()
        
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        
        loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
        hdawg.AWGs[0].load_program(loadprog)
        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[0].run_loop()
#        time.sleep(0.1)
        
        total_samples = card.samples
        amps   = np.zeros((len(freqs), card.samples))
        phases = np.zeros((len(freqs), card.samples))
#        Is_full  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
#        Qs_full  = np.zeros((len(freqs), total_samples))
        
        print('Current hold time:{}, max:{}'.format(hold_time, times[-1]))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]
            
            qubitgen.freq  = freq
            qubitgen.power = drive_power
#            time.sleep(0.1)
            
            I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            # Acquire g and e cavity traces to convert voltage to population 
            if exp_settings['subtract_background']:
                if find%exp_settings['back_rate']==0:
                    #Acquire g trace
                    qubitgen.output='Off'
#                    time.sleep(0.1)
                    I_window_g, Q_window_g, I_full_g, Q_full_g, xaxis_g = read_and_process(card, settings, 
                                                                     plot=first_it, 
                                                                     IQstorage = True)
                    #Acquire the e trace
                    qubitgen.power = qubit_power
                    qubitgen.freq  = exp_settings['Q_Freq']
                    qubitgen.output='On'
                    hdawg.AWGs[0].stop()
                    qubit_I.reset()
                    qubit_marker.reset()
                    
                    position = start_time-delay-num_sigma*sigma
                    qubit_time = num_sigma*sigma
                    
                    qubit_I.add_pulse('gaussian', position=position, amplitude=q_pulse['piAmp'], sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
                    
                    qubit_marker.add_window(position-qubit_time, position+2*qubit_time+hold_time)
                    awg_sched.plot_waveforms()
                    
                    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
                    
                    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
                    hdawg.AWGs[0].load_program(loadprog)
                    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
                    hdawg.AWGs[0].run_loop()
#                    time.sleep(0.1)
                    I_window_e, Q_window_e, I_full_e, Q_full_e, xaxis_e = read_and_process(card, settings, 
                                                                     plot=first_it, 
                                                                     IQstorage = True)
                    
                    #Reset the qubit pulse to the right length
                    hdawg.AWGs[0].stop()
                    qubit_I.reset()
                    qubit_marker.reset()
                    
                    position = start_time-delay-num_sigma*sigma
                    qubit_time = num_sigma*sigma
                    
                    qubit_I.add_pulse('gaussian_square', position=position-hold_time, amplitude=q_pulse['piAmp'], length = hold_time, ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
                    
                    qubit_marker.add_window(position-qubit_time-hold_time, position+2*qubit_time+hold_time)
                    awg_sched.plot_waveforms()
                    
                    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
                    
                    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
                    hdawg.AWGs[0].load_program(loadprog)
                    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
                    hdawg.AWGs[0].run_loop()
#                    time.sleep(0.1)
                    
                    I_ground, Q_ground = [np.mean(I_window_g), np.mean(Q_window_g)] #<I>, <Q> for ground trace
                    I_excited, Q_excited = [np.mean(I_window_e), np.mean(Q_window_e)] #<I>, <Q> for ground trace
            
                    contrast = np.sqrt((I_excited-I_ground)**2+(Q_excited-Q_ground)**2)
                
            else:
                I_window_g, Q_window_g, I_full_g, Q_full_g = 0,0,0,0
                contrast = 1
            
            if first_it:
                first_it=False
            ##Useful handles for variables
            I_sig, Q_sig   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
            I_ground, Q_ground = [np.mean(I_window_g), np.mean(Q_window_g)] #<I>, <Q> for ground trace
            theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
            theta_ground = np.arctan2(Q_ground, I_ground)*180/np.pi #angle relative to x axis in IQ plane 
            
            I_final = (I_sig-I_ground)/contrast #compute <I_net> in the data window
            Q_final = (Q_sig-Q_ground)/contrast #compute <Q_net> in the data window
            I_full_net = (I_full-I_full_g)/contrast #full I data with background subtracted
            Q_full_net = (Q_full-Q_full_g)/contrast #full Q data with background subtracted
            
            amp = np.sqrt(I_final**2 + Q_final**2)
            phase = np.arctan2(Q_final, I_final)*180/np.pi
            amp_full = np.sqrt(I_full_net**2+Q_full_net**2)  
            phase_full = np.arctan2(Q_full_net, I_full_net)*180/np.pi             
            amps[find,:]   = amp_full
            phases[find,:] = phase_full

            timedat[timeind, find]  = np.mean(amp)
            phasedat[timeind, find] = np.mean(phase)
            

        full_data = {}
        full_data['xaxis']  = freqs/1e9
        full_data['mags']   = timedat[0:timeind+1]
        full_data['phases'] = phasedat[0:timeind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mag']   = timedat[timeind]
        single_data['phase'] = phasedat[timeind]

        yaxis = times[0:timeind+1]
        labels = ['Freq (GHz)', 'Hold time (us)']
        simplescan_plot(full_data, single_data, yaxis*1e6, filename, labels, identifier='', fig_num=1) 
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

        full_time = {}
        full_time['xaxis']  = xaxis*1e6
        full_time['mags']   = amps
        full_time['phases'] = phases

        single_time = {}
        single_time['xaxis'] = xaxis*1e6
        single_time['mag']   = amp_full
        single_time['phase'] = phase_full
        
#        full_time = {}
#        full_time['xaxis']  = xaxis*1e6
#        full_time['Is']   = Is_full
#        full_time['Qs'] = Qs_full
#
#        single_time = {}
#        single_time['xaxis'] = xaxis*1e6
#        single_time['I']   = I_full
#        single_time['Q'] = Q_full

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Hold time: {}us'.format(hold_time*1e6)
        simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces\n'+filename, time_labels, identifier, fig_num=2, IQdata=False)
        plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)
        
        userfuncs.SaveFull(saveDir, filename, ['times','freqs', 'timedat', 'phasedat','xaxis', 'full_data', 'single_data', 'full_time', 'single_time'], 
                            locals(), expsettings=settings, instruments=instruments)
        

    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))

    cavitygen.Output = 'Off'
    qubitgen.Output  = 'Off'
    LO.Output        = 'Off'
    return t2-tstart