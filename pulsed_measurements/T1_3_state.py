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
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process
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
    
def meas_T1_3_state(instruments, settings):
    ## Instruments used
    qubitgen  = instruments['qubitgen']
    extragen  = instruments['extragen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    ## Bookkeeping and setting up the save directory
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    Q_Freq    = exp_settings['Q_Freq']
    Q_Power   = exp_settings['Q_Power']
    ef_Freq    = exp_settings['ef_Freq']
    ef_Power   = exp_settings['ef_Power']
    CAV_Freq  = exp_settings['CAV_Freq']
    CAV_Power = exp_settings['CAV_Power']
    
    CAV_Attenuation  = exp_globals['CAV_Attenuation']
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Extra_Attenuation = exp_globals['Extragen_Attenuation']
    
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
    extragen.Freq   = ef_Freq
    extragen.Power  = ef_Power + Extra_Attenuation
    qubitgen.enable_IQ()
    qubitgen.enable_pulse()
    extragen.enable_IQ()
    extragen.enable_pulse()
    
    cavitygen.Output = 'On'
    qubitgen.Output  = 'On'
    
    cavitygen.output = 'On'
    qubitgen.output  = 'On'
    
    ## Configure card
    configure_card(card, settings)
   
    ## Configure HDAWG
    configure_hdawg(hdawg, settings)
    
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
    awg_sched.add_analog_channel(3, name='Extra_I')
    awg_sched.add_analog_channel(4, name='Extra_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(3, name='Extra_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(4, name='blank', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    qubit_I       = awg_sched.analog_channels['Qubit_I']
    extra_I       = awg_sched.analog_channels['Extra_I']
    qubit_marker  = awg_sched.digital_channels['Qubit_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    extra_marker  = awg_sched.digital_channels['Extra_enable']
    
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']

    cavity_marker.add_window(start_time, start_time+window_time)

    ## Set up array of taus and randomize it
    if exp_settings['spacing']=='Log':
        tau_list = np.logspace(np.log10(exp_settings['Tau_min']), np.log10(exp_settings['Tau_max']), exp_settings['Tau_points'])
    else:
         tau_list = np.linspace(exp_settings['Tau_min'], exp_settings['Tau_max'], exp_settings['Tau_points'])
    taus = np.round(tau_list, 9)
    
    indices = list(range(len(taus)))
#    np.random.shuffle(indices)

    ## Start the main measurement loop 
    total_samples = card.samples
    
    amp_int = np.zeros(len(taus))
    ang_int = np.zeros(len(taus))
    amps    = np.zeros((len(taus),total_samples))
    angles  = np.zeros((len(taus),total_samples))
    #amps    = np.zeros((len(taus),card.samples))
    #angles  = np.zeros((len(taus),card.samples))
    
        ##create the digital down conversion filter if needed.
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
    
    f_pop = np.zeros(len(taus))
    e_pop = np.zeros(len(taus))
    g_pop = np.zeros(len(taus))
    alpha = np.zeros(len(taus))
    beta  = np.zeros(len(taus))
    IQ_data = np.zeros((len(taus),2))
    IQ_f = np.zeros((len(taus),2))
    IQ_e = np.zeros((len(taus),2))
    IQ_g = np.zeros((len(taus),2))
    for tind in indices:
        
        tau = taus[tind]
        print('Tau: {}'.format(tau))
        
        hdawg.AWGs[0].stop()
        qubit_I.reset()
        extra_I.reset()
        qubit_marker.reset()
        extra_marker.reset()
        
        position = start_time-delay-num_sigma*sigma
        qubit_time = num_sigma*sigma
        
        #First pulse
        qubit_I.add_pulse('gaussian', position=position-tau-1*qubit_time, amplitude=q_pulse['piAmp'], sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
        qubit_marker.add_window(position-2*qubit_time-tau, position+qubit_time-tau)
        extra_I.add_pulse('gaussian', position=position-tau-0*qubit_time, amplitude=q_pulse['piAmp'], sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
        extra_marker.add_window(position-1*qubit_time-tau, position-tau+2*qubit_time)
#        #Second pi series (to bring the population back down again)
#        qubit_I.add_pulse('gaussian', position=position+qubit_time, amplitude=q_pulse['piAmp'], sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
#        qubit_marker.add_window(position-qubit_time, position+2*qubit_time)
#        extra_I.add_pulse('gaussian', position=position-qubit_time, amplitude=q_pulse['piAmp'], sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
#        extra_marker.add_window(position-qubit_time, position+2*qubit_time)
#        awg_sched.plot_waveforms()
        awg_sched.plot_waveforms()
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        
        [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Extra_I', 'Extra_Q'], ['Extra_enable', 'blank'])
        
        loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
        hdawg.AWGs[0].load_program(loadprog)
        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)
        
        qubitgen.output = 'On'
        extragen.output = 'On'
        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
        if exp_settings['subtract_background']:
            ##Configure schedule to acquire F, E, G traces
            tau = 0
            qubit_I.reset()
            extra_I.reset()
            qubit_marker.reset()
            extra_marker.reset()
            qubit_I.add_pulse('gaussian', position=position-tau-1*qubit_time, amplitude=q_pulse['piAmp'], sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
            qubit_marker.add_window(position-2*qubit_time-tau, position+qubit_time-tau)
            extra_I.add_pulse('gaussian', position=position-tau-0*qubit_time, amplitude=q_pulse['piAmp'], sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
            extra_marker.add_window(position-1*qubit_time-tau, position-tau+2*qubit_time)
            awg_sched.plot_waveforms()
            [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
            
            [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Extra_I', 'Extra_Q'], ['Extra_enable', 'blank'])
            
            loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
            hdawg.AWGs[0].load_program(loadprog)
            hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
            hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
            hdawg.AWGs[0].run_loop()
            time.sleep(0.1)

            ## F trace
            time.sleep(0.1)
            I_window_f, Q_window_f, I_full_f, Q_full_f, xaxis_b = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
            
            ## E trace
            extragen.output='Off'
            time.sleep(0.1)
            I_window_e, Q_window_e, I_full_e, Q_full_e, xaxis_b = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
            ## G trace
            qubitgen.output='Off'
            time.sleep(0.1)
            I_window_g, Q_window_g, I_full_g, Q_full_g, xaxis_b = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
        else:
            I_window_b, Q_window_b, I_full_b, Q_full_b = 0,0,0,0
        
        [I_sig, Q_sig] = [np.mean(I_window), np.mean(Q_window)]
        [I_f, Q_f] = [np.mean(I_window_f), np.mean(Q_window_f)]
        [I_e, Q_e] = [np.mean(I_window_e), np.mean(Q_window_e)]
        [I_g, Q_g] = [np.mean(I_window_g), np.mean(Q_window_g)]
        
        ##Store this good stuff
        IQ_data[tind] = [I_sig, Q_sig]
        IQ_f[tind] = [I_f, Q_f]
        IQ_e[tind] = [I_e, Q_e]
        IQ_g[tind] = [I_g, Q_g]
        ## Since everything is centered around the |e> state (that's where we parked the cavity), 
        ## we should compute vectors relative to the e state, right? 
        ge_vec = np.array([I_g-I_e, Q_g-Q_e])
        fe_vec = np.array([I_f-I_e, Q_f-Q_e])
        se_vec = np.array([I_sig-I_e, I_sig-I_e])
        sg_vec = np.array([I_sig-I_g, I_sig-I_g])
        
        ge_contrast = np.dot(ge_vec, ge_vec)
        fe_contrast = np.dot(fe_vec, fe_vec)
        alpha[tind] = np.dot(sg_vec, -ge_vec)/ge_contrast
        beta[tind]  = np.dot(se_vec, fe_vec)/fe_contrast
        
        ##Nominally we should be able to plot the g,e,f populations now
        total = 1+alpha[tind]+alpha[tind]*beta[tind]
        f_pop[tind] = alpha[tind]*beta[tind]/(total)
        e_pop[tind] = alpha[tind]/total
        g_pop[tind] = 1/total
        
        fig = plt.figure(1)
        plt.clf()
        plt.plot(IQ_data[:tind,0], IQ_data[:tind,1], label='Signal')
        plt.plot(IQ_g[:tind,0], IQ_g[:tind,1], label='G')
        plt.plot(IQ_e[:tind,0], IQ_e[:tind,1], label='E')
        plt.plot(IQ_f[:tind,0], IQ_f[:tind,1], label='F')
        plt.legend()
        fig.canvas.draw()
        fig.canvas.flush_events()
        
        fig2 = plt.figure(2)
        plt.clf()
        plt.plot(taus[:tind]*1e6, g_pop[:tind], label = 'G')
        plt.plot(taus[:tind]*1e6, e_pop[:tind], label = 'E')
        plt.plot(taus[:tind]*1e6, f_pop[:tind], label = 'F')
        plt.legend()
        fig2.canvas.draw()
        fig2.canvas.flush_events()
        
#        ##Useful handles for variables
#        I_sig, Q_sig   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
#        I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)] #<I>, <Q> for background trace
#        theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
#        theta_back = np.arctan2(Q_back, I_back)*180/np.pi #angle relative to x axis in IQ plane 
#        
#        I_final = I_sig-I_back #compute <I_net> in the data window
#        Q_final = Q_sig-Q_back #compute <Q_net> in the data window
#        
#        amps[tind] = np.sqrt((I_full-I_full_b)**2+(Q_full-Q_full_b)**2)
#        angles[tind] = np.arctan2((Q_full-Q_full_b), (I_full-I_full_b))*180/np.pi
#        
#        amp_int[tind] = np.sqrt(I_final**2+Q_final**2)
#        ang_int[tind] = np.arctan2(Q_final, I_final)*180/np.pi
#        if first_it:
#            tstop = time.time()
#            estimate_time(tstart, tstop, len(taus))
#            first_it = False      
#
#        fig = plt.figure(1, figsize=(13,8))
#        plt.clf()
#        plt.subplot(121)
#        plt.plot(taus*1e6, amp_int, 'x')
#        plt.xlabel('Tau (us)')
#        plt.ylabel('Amplitude')  
#        plt.subplot(122)
#        plt.plot(taus*1e6, ang_int, 'x')
#        plt.xlabel('Tau (us)')
#        plt.ylabel('Phase')  
#        plt.title('Live T1 data (no fit)\n'+filename)
#        fig.canvas.draw()
#        fig.canvas.flush_events()
#        plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)
#
#        fig2 = plt.figure(2,figsize=(13,8))
#        plt.clf()
#
#        ax = plt.subplot(1,1,1)
#        general_colormap_subplot(ax, xaxis*1e6, taus*1e6, amps, ['Time (us)', 'Tau (us)'], 'Raw data\n'+filename)
#
#        plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
        userfuncs.SaveFull(saveDir, filename, ['taus','xaxis', 
                                               'IQ_data', 'IQ_g', 'IQ_e', 'IQ_f', 
                                               'alpha', 'beta', 'g_pop','e_pop', 'f_pop'], 
                           locals(), expsettings=settings, instruments=instruments)

    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))

#    T1_guess = exp_settings['T1_guess']
#    amp_guess = max(amp_int)-min(amp_int)
#    offset_guess = np.mean(amp_int[-10:])
#
#    fit_guess = [T1_guess, amp_guess, offset_guess]
#    T1, amp, offset, fit_xvals, fit_yvals = fit_T1(taus, amp_int, fit_guess)
#    fig3 = plt.figure(3)
#    plt.clf()
#    plt.plot(taus*1e6, amp_int)
#    plt.plot(fit_xvals*1e6, fit_yvals)
#    plt.title('T1:{}us \n {}'.format(np.round(T1*1e6,3), filename))
#    plt.xlabel('Time (us)')
#    plt.ylabel('Amplitude')
#    fig3.canvas.draw()
#    fig3.canvas.flush_events()
#    plt.savefig(os.path.join(saveDir, filename+'_fit.png'), dpi=150)
#
#    userfuncs.SaveFull(saveDir, filename, ['taus','xaxis', 'amps', 'amp_int', 'tau', 'amp', 'offset', 'fit_guess'],
#                         locals(), expsettings=settings, instruments=instruments)
    
    IQ_full = {}
    IQ_full['s'] = IQ_data
    IQ_full['g'] = IQ_g
    IQ_full['e'] = IQ_e
    IQ_full['f'] = IQ_f
    pops = {}
    pops['g'] = g_pop
    pops['e'] = e_pop
    pops['f'] = f_pop
    greeks = {}
    greeks['alpha'] = alpha
    greeks['beta'] = beta
    return taus, IQ_full, pops, greeks