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

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'T1_drainage_meas'
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
    
def meas_T1(instruments, settings):
    ## Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    ## Bookkeeping and setting up the save directory
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
    LO.Power  = 12
    LO.Freq   = CAV_Freq - 1e6
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

    # We'll be using the "hold" utility to keep the qubit drive high by default and pulse it low when we are measuring
    hdawg.Channels[1].configureChannel(amp=exp_globals['hdawg_config']['amplitude'], marker_out='Marker', hold=True) 
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\T1_drainage.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))

    ## Set up array of taus and randomize it
    if settings['spacing']=='Log':
        tau_list = np.logspace(np.log10(settings['Tau_min']), np.log10(settings['Tau_max']), settings['Tau_points'])
    else:
         tau_list = np.linspace(settings['Tau_min'], settings['Tau_max'], settings['Tau_points'])
    taus = np.round(tau_list, 9)
    
    indices = list(range(len(taus)))
    np.random.shuffle(indices)

    ## Start the main measurement loop 
    amp_int = np.zeros(len(taus))
    amps    = np.zeros((len(taus),card.samples))
    
    tstart = time.time()
    first_it = True

    for tind in indices:
        
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
        plt.title('Live T1 data (no fit)\n'+filename)
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

    T1_guess = exp_settings['T1_guess']
    amp_guess = max(amps)-min(amps)
    offset_guess = np.mean(amps[-10:])

    fit_guess = [T1_guess, amp_guess, offset_guess]
    T1, amp, offset, fit_xvals, fit_yvals = fit_T1(taus, amp_int, fit_guess)
    fig3 = plt.figure(3)
    plt.plot(fit_xvals*1e6, amps)
    plt.plot(fit_xvals*1e6, fit_yvals)
    plt.title('T1:{}us \n {}'.format(np.round(T1*1e6,3), filename))
    plt.xlabel('Time (us)')
    plt.ylabel('Amplitude')
    fig3.canvas.draw()
    fig3.canvas.flush_events()
    plt.savefig(os.path.join(saveDir, filename+'_fit.png'), dpi=150)

    userfuncs.SaveFull(saveDir, filename, ['taus','xaxis', 'amps', 'amp_int', 'tau', 'amp', 'offset', 'fit_guess'],
                         locals(), expsettings=settings, instruments=instruments)

    return T1, taus, amp_int