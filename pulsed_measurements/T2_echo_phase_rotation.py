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

    settings['T2_guess'] = 10e-6
    
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
    
    qubitgen.freq   = Q_Freq + exp_settings['detuning']
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
    
    ## Configure HDAWG
    configure_hdawg(hdawg, settings)
    #hard coded the delay between the two AWG cores
#    hdawg.Channels[0].configureChannel(amp=0.5,marker_out='Marker', hold='False', delay=30e-9)
#    hdawg.Channels[1].configureChannel(amp=0.5,marker_out='Marker', hold='False', delay=30e-9)
#    
#    HDAWG_dir = r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes"
#    progFile = open(os.path.join(HDAWG_dir,'T2_echo_phase_rotation.cpp'),'r')
#    rawprog  = progFile.read()
#    loadprog = rawprog
#    progFile.close()
#
#    q_pulse = exp_globals['qubit_pulse']
#    m_pulse = exp_globals['measurement_pulse']
#    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
#    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))
#    loadprog = loadprog.replace('_wait_time_',str(q_pulse['delay']))
#    loadprog = loadprog.replace('_qsigma_',str(q_pulse['sigma']))
#    loadprog = loadprog.replace('_num_sigma_',str(q_pulse['num_sigma']))
#    loadprog = loadprog.replace('_piAmp_',str(q_pulse['piAmp']))
#    
#    loadprog = loadprog.replace('_pi_count_', str(exp_settings['pulse_count']))
#    if exp_settings['T2_mode'] == 'phase_rotation':
#        loadprog = loadprog.replace('_IF_', str(exp_settings['phase_rotation_f']))
#    elif exp_settings['T2_mode'] == 'detuning':
#        loadprog = loadprog.replace('_IF_', str(0))
#        qubitgen.Freq   = Q_Freq + exp_settings['detuning']
#    else:
#        print('Invalid T2 mode set, valid options are: "phase_rotation" or "detuning"')
#        return
    
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

    taus = np.linspace(exp_settings['Tau_min'],exp_settings['Tau_max'],exp_settings['Tau_points'])
    taus = np.round(taus, 9)
    
    ## Start main measurement loop
    amp_int = np.zeros(len(taus))
    ang_int = np.zeros(len(taus))
    amps    = np.zeros((len(taus),card.samples))
    
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
        
    for tind in range(len(taus)):
            
        tau = taus[tind]
        print('Tau: {}'.format(tau))
#        finalprog = loadprog
#        finalprog = finalprog.replace('_tau_',str(tau))
#        hdawg.AWGs[0].load_program(finalprog)
#        hdawg.AWGs[0].run_loop()
#        time.sleep(0.1)
        hdawg.AWGs[0].stop()
        qubit_I.reset()
        qubit_marker.reset()
        
        position = start_time-delay-num_sigma*sigma
        qubit_time = num_sigma*sigma
        
        qubit_I.add_pulse('gaussian', position=position-tau, amplitude=q_pulse['piAmp']/2, sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
        qubit_I.add_pulse('gaussian', position=position, amplitude=q_pulse['piAmp']/2, sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
        
        qubit_marker.add_window(position-qubit_time-tau, position+2*qubit_time-tau)
        qubit_marker.add_window(position-qubit_time, position+2*qubit_time)
        awg_sched.plot_waveforms()
        
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        
        loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
        hdawg.AWGs[0].load_program(loadprog)
        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)

#        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, plot=first_it, IQstorage=True)    
#  
#        amps[tind] = np.sqrt((I_full**2+Q_full**2))
#        amp_int[tind] = np.sqrt(np.mean(I_window)**2+np.mean(Q_window)**2)
#        ang_int[tind] = np.arctan2(np.mean(Q_window),np.mean(I_window))
        amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot=first_it, IQstorage=False)
        amps[tind] = amp_full
        amp_int[tind] = np.mean(amp)
        ang_int[tind] = np.mean(phase)
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(taus))
            first_it = False
        
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(taus*1e6, amp_int)
        plt.suptitle('Live T2 data (no fit), {} pi pulses'.format(exp_settings['pulse_count']))
        plt.xlabel('Tau (us)')
        plt.ylabel('Amplitude')
        plt.subplot(122)
        plt.plot(taus*1e6, ang_int)
        plt.ylabel('Angle')
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)

        fig2 = plt.figure(2,figsize=(13,8))
        plt.clf()

        ax = plt.subplot(1,1,1)
        general_colormap_subplot(ax, xaxis*1e6, taus*1e6, amps, ['Time (us)', 'Tau (us)'], 'Raw data\n'+filename)

        plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
        userfuncs.SaveFull(saveDir, filename, ['taus','xaxis', 'amps', 'amp_int'], locals(), expsettings=settings, instruments=instruments)
    
    t2 = time.time()
    
    print('Elapsed time: {}'.format(t2-tstart))
    
    T2_guess     = exp_settings['T2_guess']
    amp_guess    = max(amp_int)-min(amp_int)
    offset_guess = np.mean(amp_int[-10:])
    if exp_settings['T2_mode']=='detuning':
        freq_guess = exp_settings['detuning']
    else:
        freq_guess   = exp_settings['phase_rotation_f']
    phi_guess    = 0

    fit_guess = [T2_guess, amp_guess, offset_guess, freq_guess, phi_guess]
    T2, amp, offset, freq, phi, fit_xvals, fit_yvals = fit_T2(taus, amp_int, fit_guess)
    fig3 = plt.figure(3)
    plt.clf()
    plt.plot(taus*1e6, amp_int)
    plt.plot(fit_xvals*1e6, fit_yvals)
    plt.title('T2:{}us freq:{}MHz. {} pi pulses \n {}'.format(np.round(T2*1e6,3), np.round(freq/1e6, 3), exp_settings['pulse_count'], filename))
    plt.xlabel('Time (us)')
    plt.ylabel('Amplitude')
    fig3.canvas.draw()
    fig3.canvas.flush_events()
    plt.savefig(os.path.join(saveDir, filename+'_fit.png'), dpi=150)

    userfuncs.SaveFull(saveDir, filename, ['taus','xaxis', 'amps', 'amp_int', 'tau', 'amp', 'offset', 'freq', 'phi', 'fit_guess'],
                         locals(), expsettings=settings, instruments=instruments)

    return T2, freq, taus, amp_int