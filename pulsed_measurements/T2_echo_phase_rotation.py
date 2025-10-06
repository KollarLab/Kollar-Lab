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
#from utility.userfits import fit_T2
from utility.FitT import fit_T2, fit_T1
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process
from utility.scheduler import scheduler

import scipy.signal as signal
from matplotlib.cm import ScalarMappable

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
    settings['fit_data'] = True
    
    settings['verbose'] = True
    
    return settings

def meas_T2_phase_rotation(instruments, settings):
    tstart = time.time()
    t_init_start = time.time()
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
#    qubitgen.freq   = Q_Freq + exp_settings['detuning']
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
#    configure_hdawg(hdawg, settings)
    
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

    
    time_array = np.zeros((len(taus),5))
    t_pulse_array = np.zeros((len(taus),5))
    save_plot = np.zeros((len(taus), 6))
    
    fig1 = plt.figure(1,figsize=(13,8))
    fig1.clf()
    fig1.suptitle('Live T2 data (no fit), {} pi pulses'.format(exp_settings['pulse_count']))
    ax11 = fig1.add_subplot(121)
    ax12 = fig1.add_subplot(122)
    fig2 = plt.figure(2,figsize=(13,8))
    fig2.clf()
    ax21 = fig2.add_subplot(111)
#    cax = fig2.add_subplot(1,2,2)
    t_init_stop = time.time() 
    for tind in range(len(taus)):
        t_loop_start  = time.time() 
        #Generating pulse sequence and upload
        t_pulse_start = time.time()
        
        tau = taus[tind]
        if exp_settings['verbose']:
            print('Tau: {}'.format(tau))
        hdawg.AWGs[0].stop()
        qubit_I.reset()
        qubit_Q.reset()
        qubit_marker.reset()
        t_reset_t = time.time()
        t_reset = t_reset_t-t_loop_start
        position = start_time-delay
        qubit_time = num_sigma*sigma + q_pulse['hold_time']
        
# =============================================================================
#         if qubit_time >= np.min(taus):
#             print('qubit_time is longer than tau_min')
#             print('Two qubit pulses will overlap!')
#             print('qubit_time = ' + str(qubit_time))
#             
#             break
# =============================================================================
        
        #the main pulses    
        qubit_I.add_pulse('gaussian_square', position=position-tau-2*qubit_time, #this factor of 2 is to make sure the separation between the end of this pulse
                          amplitude=q_pulse['piAmp']/2,                          #and the beginning of 2nd pulse is tau.
                          length = hold_time,
                          ramp_sigma=q_pulse['sigma'], 
                          num_sigma=q_pulse['num_sigma'])

        if exp_settings['T2_mode'] == 'detuning':
            if exp_settings['basis'] =='X':
                qubit_I.add_pulse('gaussian_square', position=position-qubit_time, 
                                  amplitude=q_pulse['piAmp']/2,
                                  length = hold_time,
                                  ramp_sigma=q_pulse['sigma'], 
                                  num_sigma=q_pulse['num_sigma'])
            elif exp_settings['basis'] =='Y':
                qubit_Q.add_pulse('gaussian_square', position=position-qubit_time, 
                                  amplitude=q_pulse['piAmp']/2,
                                  length = hold_time,
                                  ramp_sigma=q_pulse['sigma'], 
                                  num_sigma=q_pulse['num_sigma'])
            else:
                print(exp_settings['basis'])
                print('Invalid basis selected for second pulse ("X" or "Y")')
                
        elif exp_settings['T2_mode'] == 'phase_rotation':
            qubit_I.add_pulse('gaussian_square', position=position-qubit_time,
                              length = hold_time,
                              amplitude= np.cos(2*np.pi * tau* exp_settings['phase_rotation_f'])*q_pulse['piAmp']/2, 
                              ramp_sigma=q_pulse['sigma'], 
                              num_sigma=q_pulse['num_sigma'])
            qubit_Q.add_pulse('gaussian_square', position=position-qubit_time, 
                              length = hold_time, 
                              amplitude= np.sin(2*np.pi * tau* exp_settings['phase_rotation_f'])*q_pulse['piAmp']/2, 
                              ramp_sigma=q_pulse['sigma'], 
                              num_sigma=q_pulse['num_sigma'])
        else:
            raise ValueError('Invalid T2_mode')
            
        
        qubit_marker.add_window(position-tau-2*qubit_time-250e-9, position-tau-qubit_time+250e-9) # the extra 250e-9 is to include the qubit pulse fully
        qubit_marker.add_window(position-qubit_time-250e-9, position+250e-9) # same here
        #qubit_marker.add_window(position-tau-2*qubit_time-250e-9, position+250e-9) # this is a temporary solution just to check if the marker
        
        if exp_settings['pulse_count'] > 0:
            numPulses = exp_settings['pulse_count']
            CPMG_times = tau*(np.linspace(1, numPulses, numPulses)-0.5)/numPulses
            temp = np.linspace(0,tau, numPulses + 2)
            pulseTimes = CPMG_times#temp[1:-1]
            
            if tau <= qubit_time:
                print('qubit_time is longer than pulse separation tau')
                print('Qubit pulses will overlap!')
                print('qubit_time = ' + str(qubit_time))
                print('Separation = ' + str(tau))
            
                break
            
            for tp in pulseTimes:
                qubit_I.add_pulse('gaussian_square', position=position-tp-qubit_time, 
                                  length = hold_time,
                                  amplitude=q_pulse['piAmp'], 
                                  ramp_sigma=q_pulse['sigma'], 
                                  num_sigma=q_pulse['num_sigma'])
                qubit_marker.add_window(position-tp-qubit_time-10e-9, position-tp+10e-9) # same here for including the qubit pulses
                print('Separation = ' + str(tau-qubit_time))
      
        t_pulse_gen_t = time.time()
        t_pulse_gen = t_pulse_gen_t-t_reset_t
        awg_sched.plot_waveforms()
        t_plot_wf_t = time.time()
        t_plot_wf = t_plot_wf_t-t_pulse_gen_t
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        t_compile_t = time.time()
        t_compile = t_compile_t-t_plot_wf_t
        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[0].run_loop()
        qubitgen.output='On'
#        time.sleep(0.1)
        t_upload_t = time.time()
        t_upload = t_upload_t-t_compile_t
        t_pulse_stop = time.time()
        t_pulse_array[tind] = [t_reset, t_pulse_gen, t_plot_wf, t_compile, t_upload]
        #Acquiring data
        t_data_start = time.time()
        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
        t_data_stop = time.time()
        
        #Background collection
        t_back_start = time.time()
        if exp_settings['subtract_background']:
            #Acquire background trace
#            qubitgen.freq=3.8e9
            qubitgen.output='Off'
#            time.sleep(0.1)
            I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
            qubitgen.freq=Q_Freq + exp_settings['detuning']
        else:
            I_window_b, Q_window_b, I_full_b, Q_full_b = 0,0,0,0
        t_back_stop = time.time()
        
        #Plotting and saving 
        t_plot_start = time.time()
        t1 = time.time()
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
        t2 = time.time()
        
        ax11.cla()
        ax12.cla()
        
        ax11.plot(taus*1e6, amp_int)
        ax11.set_xlabel('Tau (us)')
        ax11.set_ylabel('Amplitude')

        ax12.plot(taus*1e6, ang_int)
        ax12.set_ylabel('Angle')
        ax12.set_xlabel('Tau (us)')
             
        ax21.cla()
        if first_it:
            cbar = general_colormap_subplot(ax21, xaxis*1e6, 
                                            taus*1e6, amps, 
                                            ['Time (us)', 'Tau (us)'], 
                                            'Raw data\n'+filename)
        else:
            general_colormap_subplot(ax21, xaxis*1e6, 
                                            taus*1e6, amps, 
                                            ['Time (us)', 'Tau (us)'], 
                                            'Raw data\n'+filename, cbar=cbar)
        fig1.canvas.draw()
        fig1.canvas.flush_events()            
        if tind%exp_settings['num_save']==0:    
            t3 = time.time()
            fig1.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 100)
            ta = time.time()
            fig2.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 100)
            tb = time.time()
            userfuncs.SaveFull(saveDir, filename, 
                               ['taus','xaxis', 'amps', 'amp_int'], 
                               locals(), 
                               expsettings=settings, 
                               instruments=instruments, 
                               saveHWsettings=first_it)
            t4 = time.time()
            
            data_manip = t2-t1
            data_plot = t3-t2
            data_save = t4-t3
            save_plot[tind] = [data_manip, data_plot, data_save, ta-t3, tb-ta, t4-tb]
        
        t_plot_stop = time.time()
        
        t_pulse = t_pulse_stop-t_pulse_start
        t_data = t_data_stop-t_data_start
        t_back = t_back_stop-t_back_start
        t_plot = t_plot_stop-t_plot_start
        
        t_loop_stop = time.time()
        t_loop = t_loop_stop-t_loop_start
        time_array[tind] = np.array([t_pulse, t_data, t_back,t_plot, t_loop])
        first_it = False

    t_finalize_start = time.time()
    
    if exp_settings['fit_data']:
        T2_guess     = exp_settings['T2_guess']
#        amp_guess    = max(amp_int)-min(amp_int)
#        offset_guess = np.mean(amp_int[-10:])
#        if exp_settings['T2_mode']=='detuning':
#            freq_guess = exp_settings['detuning']
#        else:
#            freq_guess   = exp_settings['phase_rotation_f']
#        phi_guess    = 0
#    
#        fit_guess = [T2_guess, amp_guess, offset_guess, freq_guess, phi_guess]
#        T2, amp, offset, freq, phi, fit_xvals, fit_yvals = fit_T2(taus, amp_int, fit_guess)
        if exp_settings['pulse_count'] == 0:
            datafit = fit_T2(taus, amp_int, T2_guess)
            
        else:
            datafit = fit_T1(taus, amp_int)
            datafit['freq'] = 0
            datafit['phi'] = 0
        T2, amp, offset, freq, phi = datafit['tau'], datafit['amp'], datafit['offset'], datafit['freq'], datafit['phi']
        fig3 = plt.figure(3)
        plt.clf()
        plt.plot(taus*1e6, amp_int)
#        plt.plot(fit_xvals*1e6, fit_yvals)
        plt.plot(datafit['ts']*1e6, datafit['fit_curve'])
        plt.title('T2:{}us freq:{}MHz. {} pi pulses \n {}'.format(np.round(T2*1e6,3), 
                  np.round(freq/1e6, 3), exp_settings['pulse_count'], filename))
        plt.xlabel('Time (us)')
        plt.ylabel('Amplitude')
        fig3.canvas.draw()
        fig3.canvas.flush_events()
        fig3.savefig(os.path.join(saveDir, filename+'_fit.png'), dpi=250)
    else:
        T2, freq, amp, offset, phi = 0,0,0,0,0
#    plt.savefig(os.path.join(saveDir, filename+'_fit.png'), dpi=250)
    fig1.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 250)
    fig2.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 250)
    userfuncs.SaveFull(saveDir, filename, ['taus','xaxis', 'amps', 'amp_int', 
                                           'tau', 'amp', 'offset', 'freq', 'phi'],
                         locals(), expsettings=settings, instruments=instruments)
    t_finalize_stop = time.time()
    
    t_init = t_init_stop-t_init_start
    t_loop_avg = np.mean(time_array[:,-1])
    t_finalize = t_finalize_stop-t_finalize_start
    
    if exp_settings['verbose']:
        print('Timing info: \n initialization time:{}, \n avg loop time:{}, \n wrap up time:{}'.format(t_init, t_loop_avg, t_finalize))
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))
    return T2, freq, taus, amp_int, time_array, [t_pulse_array, save_plot]




