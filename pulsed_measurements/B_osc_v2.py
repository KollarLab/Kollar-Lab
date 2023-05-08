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
from utility.scheduler import scheduler, gaussian

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'rabi_chevron'
    settings['meas_type'] = 'rabi_chevron'
    
    #Cavity parameters
    settings['CAVpower']    = -45
    settings['CAV_freq']    = 5e9
    settings['Q_power']     = -20
    
    #Qubit parameters
    settings['start_freq']  = 4*1e9  
    settings['stop_freq']   = 5*1e9 
    settings['freq_points'] = 50

    #Rabi chevron parameters
    settings['start_time']  = 10e-9
    settings['stop_time']   = 500e-9
    settings['time_points'] = 40

    #Card settings
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 5e3
    
    return settings

def Bx_WF(settings):
    osc_freq    = settings['osc_freq']
    osc_amp     = settings['Bx_osc_amp']
    osc_time    = settings['osc_time']
    osc_points = int(2.4e9*osc_time)
    t_ax = np.linspace(0, osc_points, osc_points)/2.4e9
    osc = np.sin(2*np.pi*osc_freq*t_ax)*osc_amp 
    
    final = np.concatenate([osc])
    return final

settings = {}
settings['osc_freq'] = 5e6
settings['Bx_osc_amp'] = 0.2
settings['osc_time'] = 1e-6
settings['ramp_sigma'] = 100e-9
settings['ramp_hold'] = 200e-9
settings['mod_depth'] = 0.5
settings['Bx_hold'] = 100e-9

def Bx_WF_v2(settings):
    osc_freq    = settings['osc_freq']
    osc_amp     = settings['Bx_osc_amp']
    osc_time    = settings['osc_time']
    ramp_sigma  = settings['ramp_sigma']
    ramp_hold   = settings['Bx_hold']
    mod_depth   = settings['mod_depth']
    
    t = np.linspace(0, 4*ramp_sigma, int(4*ramp_sigma*2.4e9))
    gauss = np.exp(-(t-2*ramp_sigma)**2/(2*ramp_sigma**2))
    gauss_ramp = gauss-gauss[0]
    gauss_ramp = osc_amp*gauss_ramp/max(gauss_ramp)
    hold = np.ones(int(ramp_hold*2.4e9))*osc_amp
    osc_points = int(2.4e9*osc_time)
    t_ax = np.linspace(0, osc_points, osc_points)/2.4e9
#    osc = np.cos(2*np.pi*osc_freq*t_ax)*osc_amp 
    osc = osc_amp*(1+mod_depth*(np.cos(2*np.pi*osc_freq*t_ax)-1))
    
    final = np.concatenate([gauss_ramp[:len(gauss_ramp)//2], hold, osc])
    return final

def Bz_WF(settings):
    #Creating the ramp
    ramp_length = settings['ramp_length']
    ramp_amp    = settings['ramp_amp']
    ramp_hold   = settings['ramp_hold']
    
    osc_freq    = settings['osc_freq']
    osc_amp     = settings['Bz_osc_amp']
    osc_time    = settings['osc_time']
    
    ramp_points = int(2.4e9*ramp_length)
    hold_points = int(2.4e9*ramp_hold)
    ramp = np.linspace(0,1,ramp_points)*ramp_amp
    ramp_down = np.linspace(ramp_amp,-osc_amp,ramp_points)
    hold = np.ones(hold_points)*ramp_amp
    
    osc_points = int(2.4e9*osc_time)
    t_ax = np.linspace(0, osc_points, osc_points)/2.4e9
    osc = -np.cos(2*np.pi*osc_freq*t_ax)*osc_amp
    
    final = np.concatenate([ramp, hold, ramp_down, osc])
    return final

def Bz_WF_v2(settings):
    #Creating the ramp
    ramp_length = settings['ramp_length']
    ramp_amp    = settings['ramp_amp']
    ramp_hold   = settings['ramp_hold']
    
    osc_freq    = settings['osc_freq']
    osc_amp     = settings['Bz_osc_amp']
    osc_time    = settings['osc_time']
    
    ramp_points = int(2.4e9*ramp_length)
    hold_points = int(2.4e9*ramp_hold)
    ramp = np.linspace(0,1,ramp_points)*ramp_amp
    ramp_down = np.flip(ramp)
    hold = np.ones(hold_points)*ramp_amp
    
    osc_points = int(2.4e9*osc_time)
    t_ax = np.linspace(0, osc_points, osc_points)/2.4e9
    osc = np.sin(2*np.pi*osc_freq*t_ax)*osc_amp
    
    final = np.concatenate([ramp, hold, ramp_down, osc])
    return final
    
def Bz_osc(instruments, settings):
    
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
    Q_freq = exp_settings['Q_Freq']
    Q_power = exp_settings['Q_Power']+exp_globals['Qbit_Attenuation']
    
    start_time   = exp_settings['start_time']
    stop_time    = exp_settings['stop_time']
    time_points  = exp_settings['time_points']
    times = np.round(np.linspace(start_time,stop_time,time_points), 9)
    
    ## Generator settings
    cavitygen.freq   = CAV_freq
    cavitygen.power  = CAV_power
    cavitygen.enable_pulse()

    qubitgen.freq   = Q_freq
    qubitgen.power  = Q_power
    
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
#    configure_hdawg(hdawg, settings)
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
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
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(3, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    
    qubit_I       = awg_sched.analog_channels['Qubit_I']
    qubit_Q       = awg_sched.analog_channels['Qubit_Q']
    Flux_pulse    = awg_sched.analog_channels['Flux_pulse']
    qubit_marker  = awg_sched.digital_channels['Qubit_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)

    cavity_marker.add_window(start_time, start_time+window_time+1e-6)

    ## Starting main measurement loop 
    timedat  = np.zeros(len(times))
    phasedat = np.zeros(len(times))
    
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
    
    #Generate the XZ field. this is just a long array which we'll pick stuff from in the loop
    #The "hold time" really becomes an evolution time here
    
    Bx = Bx_WF_v2(exp_settings)
    Bz = Bz_WF_v2(exp_settings)
    ramp_length = exp_settings['ramp_length']*2+exp_settings['ramp_hold']
    Bx_offset = exp_settings['ramp_sigma']*2+exp_settings['Bx_hold']
    #0 should correspond to the end of the ramping stuff
    delay = q_pulse['delay']
    for timeind in range(len(times)):
        hold_time = times[timeind]
        
        hdawg.AWGs[0].stop()
        qubit_I.reset()
        qubit_Q.reset()
        qubit_marker.reset()
        Flux_pulse.reset()
        q_pulse_len = q_pulse['sigma']*q_pulse['num_sigma']+q_pulse['hold_time']
        Bx_ind = int((hold_time+Bx_offset)*2.4e9)
        ##adding a 50 ns buffer on Bz so that it definitely turns off after omega turns off
        Bz_ind = int((hold_time+ramp_length)*2.4e9) 
        Bx_evolv = Bx[0:Bx_ind]
        Bz_evolv = Bz[0:Bz_ind]
        
        if exp_settings['sigma']=='X':
#            print('X!')
            qubit_I.add_pulse('gaussian_square', position=start_time-delay-q_pulse_len, 
                                  amplitude=q_pulse['piAmp']/2,
                                  length = q_pulse['hold_time'],
                                  ramp_sigma=q_pulse['sigma'], 
                                  num_sigma=q_pulse['num_sigma'])
            qubit_marker.add_window(start_time-delay-100e-9, start_time+100e-9)
            Bx_hold  = Bx_evolv[-1]*np.ones(int(exp_settings['buffer_z']*2.4e9))
            Bx_evolv = np.concatenate([Bx_evolv, Bx_hold])
            position = start_time-delay-hold_time-2*q_pulse_len-exp_settings['buffer_z']
        elif exp_settings['sigma']=='Y':
#            print('X!')
            qubit_Q.add_pulse('gaussian_square', position=start_time-delay-q_pulse_len, 
                                  amplitude=q_pulse['piAmp']/2,
                                  length = q_pulse['hold_time'],
                                  ramp_sigma=q_pulse['sigma'], 
                                  num_sigma=q_pulse['num_sigma'])
            qubit_marker.add_window(start_time-delay-100e-9, start_time+100e-9)
            Bx_hold  = Bx_evolv[-1]*np.ones(int(exp_settings['buffer_z']*2.4e9))
            Bx_evolv = np.concatenate([Bx_evolv, Bx_hold])
            position = start_time-delay-hold_time-2*q_pulse_len-exp_settings['buffer_z']
        else:
            position = start_time-delay-hold_time-2*q_pulse_len
            Bz_hold  = Bz_evolv[-1]*np.ones(int(exp_settings['buffer_z']*2.4e9))
            Bz_evolv = np.concatenate([Bz_evolv, Bz_hold])
            
        qubit_I.add_pulse('custom', position-Bx_offset, wave=Bx_evolv)
        qubit_marker.add_window(position-Bx_offset-100e-9, position+hold_time+150e-9)
        Flux_pulse.add_pulse('custom', position-ramp_length-exp_settings['FBL_offset'], wave=Bz_evolv)
        
        awg_sched.plot_waveforms()
        
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Flux_pulse', 'blank'], ['blank1', 'blank2'])


        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
        hdawg.AWGs[0].run_loop()
#        time.sleep(0.1)
        
        total_samples = card.samples
        amps   = np.zeros((len(times), card.samples))
        phases = np.zeros((len(times), card.samples))
        
        print('Current hold time:{}, max:{}'.format(hold_time, times[-1]))
            
        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                     plot=first_it, 
                                                                     IQstorage = True)
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
            time.sleep(0.1)
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
        
        amps[timeind,:]   = amp_full
        phases[timeind,:] = phase_full

        timedat[timeind]  = np.mean(amp)
        phasedat[timeind] = np.mean(phase)
        
        fig = plt.figure(1)
        plt.clf()
        plt.subplot(121)
        plt.plot(times[:timeind+1], timedat[:timeind+1])
        plt.xlabel('Hold time (us)')
        plt.ylabel('Mag')
        plt.subplot(122)
        plt.plot(times[:timeind+1], phasedat[:timeind+1])
        plt.xlabel('Hold time (us)')
        plt.ylabel('Phase')
        plt.suptitle('Bz oscillation\n'+filename)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
        
        userfuncs.SaveFull(saveDir, filename, ['times', 'amps', 'phases',
                                               'timedat', 'phasedat'], 
                            locals(), expsettings=settings, instruments=instruments, 
                            saveHWsettings=first_it)
        
    userfuncs.SaveFull(saveDir, filename, ['times', 'amps', 'phases',
                                           'timedat', 'phasedat'], 
                        locals(), expsettings=settings, instruments=instruments)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))

    cavitygen.Output = 'Off'
    qubitgen.Output  = 'Off'
    LO.Output        = 'Off'
    return timedat, phasedat, times