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
    LO.Power = 12
    LO.Freq = CAV_Freq - 1e6
    LO.Output = 'On'

    cavitygen.Freq   = CAV_Freq
    cavitygen.Power  = CAV_Power + CAV_Attenuation
    cavitygen.IQ.Mod = 'On'
    
    qubitgen.Freq   = Q_Freq
    qubitgen.Power  = Q_Power + Qbit_Attenuation
    qubitgen.IQ.Mod = 'On'
    
    cavitygen.Output = 'On'
    qubitgen.Output  = 'On'
    
    ## Configure card
    configure_card(card, settings)
    
    ## Configure HDAWG
    configure_hdawg(hdawg, settings)
    
    HDAWG_dir = r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes"
    progFile = open(os.path.join(HDAWG_dir,'T2_echo_phase_rotation.cpp'),'r')
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

    loadprog = loadprog.replace('_pi_count_', str(exp_settings['pulse_count']))
    if exp_settings['T2_mode'] == 'phase_rotation':
        loadprog = loadprog.replace('_IF_', str(exp_settings['phase_rotation_f']))
    elif exp_settings['T2_mode'] == 'detuning':
        loadprog = loadprog.replace('_IF_', str(0))
        qubitgen.Freq   = Q_Freq + exp_settings['detuning']
    else:
        print('Invalid T2 mode set, valid options are: "phase_rotation" or "detuning"')
        return

    taus = np.linspace(exp_settings['Tau_min'],exp_settings['Tau_max'],exp_settings['Tau_points'])
    taus = np.round(taus, 9)
    
    ## Start main measurement loop
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

        amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot=first_it)    
  
        amps[tind] = amp_full
        amp_int[tind] = amp
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(taus))
            first_it = False
        
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.plot(taus*1e6, amps)
        plt.title('Live T2 data (no fit), {} pi pulses'.format(exp_settings['pulse_count']))
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
    if exp_settings['T2_mode']=='detuning':
        freq_guess = exp_settings['detuning']
    else:
        freq_guess   = exp_settings['phase_rotation_f']
    phi_guess    = 0

    fit_guess = [T2_guess, amp_guess, offset_guess, freq_guess, phi_guess]
    T2, amp, offset, freq, phi, fit_xvals, fit_yvals = fit_T2(taus, amp_int, fit_guess)
    fig3 = plt.figure(3)
    plt.plot(fit_xvals*1e6, amps)
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