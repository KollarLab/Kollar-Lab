# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:52:18 2020

@author: Kollarlab
"""
import time
import os
from utility.measurement_helpers import configure_card, estimate_time, read_and_process
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import base_power_plot_imshow, general_colormap_subplot
from utility.userfits import fit_T2

def GetDefaultSettings():
    settings = {}
    
    settings['scanname']  = 'T2_meas'
    settings['meas_type'] = 'Tmeas'

    
    settings['Q_Freq']    = 4.20431e9
    settings['Q_Power']   = -11
    settings['CAV_Freq']  = 8.126e9
    settings['CAV_Power'] = -18
    
    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 25e3
    
    settings['Tau_min'] = 200e-9
    settings['Tau_max'] = 30e-6
    settings['Tau_points'] = 5
    
    return settings

def meas_T2(instruments, settings):
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
    
    ## Generator settings
    LO.Power = 12
    LO.Freq = CAV_Freq - exp_globals['IF']
    LO.Output = 'On'

    cavitygen.Freq   = CAV_Freq
    cavitygen.Power  = CAV_Power + CAV_Attenuation
    cavitygen.IQ.Mod = 'On'
    
    qubitgen.Freq   = Q_Freq
    qubitgen.Power  = Q_Power + Qbit_Attenuation
    qubitgen.IQ.Mod = 'On'
    
    cavitygen.Output = 'On'
    qubitgen.Output = 'On'
    
    ## Card Settings
    configure_card(card, settings)

    ## HDAWG
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\T2.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
        
    q_pulse = exp_globals['qubit_pulse']
    m_pulse = exp_globals['measurement_pulse']
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))
    loadprog = loadprog.replace('_wait_time_',str(q_pulse['delay']))
    loadprog = loadprog.replace('_qsigma_',str(q_pulse['sigma']))
    loadprog = loadprog.replace('_num_sigma_',str(q_pulse['num_sigma']))
        
    taus = np.linspace(settings['Tau_min'],settings['Tau_max'] , settings['Tau_points'] )
    
    amp_int = np.zeros(len(taus))
    amps    = np.zeros((len(taus),card.samples))
    
    tstart = time.time()
    first_it = True

    for tind in range(len(taus)):
            
        tau = taus[tind]
        print('Tau: {}'.format(tau))
        finalprog = loadprog
        finalprog = finalprog.replace('_tau_',str(tau))
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)
        
        amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings)

        amps[tind]    = amp_full
        amp_int[tind] = np.mean(amp)
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop)
    
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.plot(taus*1e6, amps)
        plt.title('Live T2 data (no fit)')
        plt.xlabel('Tau (us)')
        plt.ylabel('Amplitude')  
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
    amp_guess    = max(amps)-min(amps)
    offset_guess = np.mean(amps[-10:])
    freq_guess = exp_settings['detuning']
    phi_guess    = 0

    fit_guess = [T2_guess, amp_guess, offset_guess, freq_guess, phi_guess]
    T2, amp, offset, freq, phi, fit_xvals, fit_yvals = fit_T2(taus, amp_int, fit_guess)
    fig3 = plt.figure(3)
    plt.plot(fit_xvals*1e6, amps)
    plt.plot(fit_xvals*1e6, fit_yvals)
    plt.title('T2:{}us freq:{}MHz \n {}'.format(np.round(T2*1e6,3), np.round(freq/1e6, 3), exp_settings['pulse_count'], filename))
    plt.xlabel('Time (us)')
    plt.ylabel('Amplitude')
    fig3.canvas.draw()
    fig3.canvas.flush_events()
    plt.savefig(os.path.join(saveDir, filename+'_fit.png'), dpi=150)

    userfuncs.SaveFull(saveDir, filename, ['taus','xaxis', 'amps', 'amp_int', 'tau', 'amp', 'offset', 'freq', 'phi', 'fit_guess'],
                         locals(), expsettings=settings, instruments=instruments)

    return T2, freq, taus, amp_int