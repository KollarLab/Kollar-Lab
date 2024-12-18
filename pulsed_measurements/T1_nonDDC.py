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
    
def meas_T1_nonDDC(instruments, settings):
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
    LO.power  = 12
    LO.freq = '{} GHz'.format((cavitygen.Freq-exp_globals['IF'])/1e9)
    LO.output = 'On'
    
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
    #hard coded the delay between the two AWG cores
    hdawg.Channels[0].configureChannel(amp=0.5,marker_out='Marker', hold='False', delay=30e-9)
    hdawg.Channels[1].configureChannel(amp=0.5,marker_out='Marker', hold='False', delay=30e-9)
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\T1.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    q_pulse = exp_globals['qubit_pulse']
    m_pulse = exp_globals['measurement_pulse']
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))
    loadprog = loadprog.replace('_qsigma_',str(q_pulse['sigma']))
    loadprog = loadprog.replace('_num_sigma_',str(q_pulse['num_sigma']))
    loadprog = loadprog.replace('_piAmp_',str(q_pulse['piAmp']))
    loadprog = loadprog.replace('_piAngle_',str(q_pulse['piAngle']))

    ## Set up array of taus and randomize it
    if exp_settings['spacing']=='Log':
        tau_list = np.logspace(np.log10(exp_settings['Tau_min']), np.log10(exp_settings['Tau_max']), exp_settings['Tau_points'])
    else:
         tau_list = np.linspace(exp_settings['Tau_min'], exp_settings['Tau_max'], exp_settings['Tau_points'])
    taus = np.round(tau_list, 9)
    
    indices = list(range(len(taus)))
#    np.random.shuffle(indices)

    ## Start the main measurement loop 
    amp_int = np.zeros(len(taus))
    ang_int = np.zeros(len(taus))
    amps    = np.zeros((len(taus),card.samples))
    angles  = np.zeros((len(taus),card.samples))
    
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

    for tind in indices:
        
        tau = taus[tind]
        print('Tau: {}'.format(tau))
        finalprog = loadprog
        finalprog = finalprog.replace('_tau_',str(tau))
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)
        
        
        amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot=first_it,
                                                                   IQstorage = False) 
        
        amps[tind] = amp_full
        angles[tind] = phase_full
        amp_int[tind] = np.mean(amp) 
        ang_int[tind] = np.mean(phase) #I think that this is the right phase to store
        
        
#        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
#                                                                         plot=first_it, 
#                                                                         IQstorage = True)
#            
#        I_final = np.mean(I_window) #compute <I> in the data window
#        Q_final = np.mean(Q_window) #compute <Q> in the data window
#        
#        amps[tind] = np.sqrt(I_full**2+Q_full**2)
#        angles[tind] = np.arctan2(Q_full, I_full)*180/np.pi
#        
#        amp_int[tind] = np.sqrt(I_final**2+Q_final**2)
#        ang_int[tind] = np.arctan2(Q_final, I_final)*180/np.pi
        
        
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(taus))
            first_it = False      

        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(taus*1e6, amp_int, 'x')
        plt.xlabel('Tau (us)')
        plt.ylabel('Amplitude')  
        plt.subplot(122)
        plt.plot(taus*1e6, ang_int, 'x')
        plt.xlabel('Tau (us)')
        plt.ylabel('Phase')  
        plt.title('Live T1 data (no fit)\n'+filename)
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
    amp_guess = max(amp_int)-min(amp_int)
    offset_guess = np.mean(amp_int[-10:])

    fit_guess = [T1_guess, amp_guess, offset_guess]
    T1, amp, offset, fit_xvals, fit_yvals = fit_T1(taus, amp_int, fit_guess)
    fig3 = plt.figure(3)
    plt.clf()
    plt.plot(taus*1e6, amp_int)
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