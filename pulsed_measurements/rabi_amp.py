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
    
    settings['scanname'] = 'rabi'
    settings['meas_type'] = 'rabi'
    
    settings['Q_Freq']  = 4.21109e9
    settings['Q_Power'] = -18

    settings['CAV_Freq']  = 8.12357e9
    settings['CAV_Power'] = -42.9
    
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 25e3
    
    settings['amp_min']    = 200e-9
    settings['amp_max']    = 30e-6
    settings['amp_points'] = 5
    
    return settings
    
def pi_amp_sweep(instruments, settings):
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
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\rabi.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    q_pulse = exp_globals['qubit_pulse']
    m_pulse = exp_globals['measurement_pulse']
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']+2e-6))
    loadprog = loadprog.replace('_qsigma_',str(q_pulse['sigma']))
    loadprog = loadprog.replace('_num_sigma_',str(q_pulse['num_sigma']))
    loadprog = loadprog.replace('_tau_',str(q_pulse['delay']))

    amp_vals = np.linspace(exp_settings['amp_min'], exp_settings['amp_max'], exp_settings['amp_points'])

    ## Start the main measurement loop 
    amp_int = np.zeros(len(amp_vals))
    ang_int = np.zeros(len(amp_vals))
    amps    = np.zeros((len(amp_vals),card.samples))
    angles  = np.zeros((len(amp_vals),card.samples))
    
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

    for ind, amp in enumerate(amp_vals):
        
        print('Amp: {}'.format(amp))
        finalprog = loadprog
        finalprog = finalprog.replace('_piAmp_',str(amp))
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)
        
#        amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot=first_it) 
#        
#        amps[tind] = amp_full
#        amp_int[tind] = amp 
        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            
        I_final = np.mean(I_window) #compute <I> in the data window
        Q_final = np.mean(Q_window) #compute <Q> in the data window
        
        amps[ind] = np.sqrt(I_full**2+Q_full**2)
        angles[ind] = np.arctan2(Q_full, I_full)*180/np.pi
        
        amp_int[ind] = np.sqrt(I_final**2+Q_final**2)
        ang_int[ind] = np.arctan2(Q_final, I_final)*180/np.pi
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(amp_vals))
            first_it = False      

        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(amp_vals[:ind], amp_int[:ind], 'x')
        plt.xlabel('pi pulse amplitude')
        plt.ylabel('Amplitude')  
        plt.subplot(122)
        plt.plot(amp_vals[:ind], ang_int[:ind], 'x')
        plt.xlabel('pi pulse amplitude')
        plt.ylabel('Phase')  
        plt.title('Rabi rate vs pi pulse amp\n'+filename)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)

        fig2 = plt.figure(2,figsize=(13,8))
        plt.clf()

        ax = plt.subplot(1,1,1)
        general_colormap_subplot(ax, xaxis*1e6, amp_vals, amps, ['Time (us)', 'pi pulse amp'], 'Raw data\n'+filename)

        plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
        userfuncs.SaveFull(saveDir, filename, ['amp_vals','xaxis', 'amps', 'amp_int'], locals(), expsettings=settings, instruments=instruments)

    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))

    return amp_vals, amp_int