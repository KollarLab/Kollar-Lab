# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:52:18 2020


08-30-20 AK making a lot of changes. Differentiating settinds and exp_settings
-It also looks like the data extraction is not quite right. The integration of the trace in time 
is missing.
-switched to plotting amp_int instead of amps, because amps is the raw time traces.
-And fixed offset subtraction since we were cutting the 2D array amps in a funny way.
-switched built-in max and min to np.max, np.min
- fixed the live plot so that it only plots the tau points that have been done 
instead of plotting all of them and getting lost of fake zeros. Also switched to 
plotting dots instead of lines.
-updated IF to follow expt globals.
-switched hold in HDAWG config command to 1 instead of True. That still doesn't fix it.
-In the old version, there was no counter to keep track of where in the loop over randomized taus you are.
This makes autoplotting very hard.
-Had to change the randomized tau list to draw from exp_settings instead of settings, because the higher-level
settings dictionary doesn't have the right fields.

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
#    LO.Freq   = CAV_Freq - 1e6 #old hard coded version
    LO.Freq   = CAV_Freq- exp_globals['IF']
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
#    hdawg.Channels[1].configureChannel(amp=exp_globals['hdawg_config']['amplitude'], marker_out='Marker', hold=True)#hold needs to a string 'True'
    hdawg.Channels[1].configureChannel(amp=exp_globals['hdawg_config']['amplitude'], marker_out='Marker', hold='True')
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\T1_drainage_generalized.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))

    ## Set up array of taus and randomize it
    if exp_settings['spacing']=='Log':
        tau_list = np.logspace(np.log10(exp_settings['Tau_min']), np.log10(exp_settings['Tau_max']), exp_settings['Tau_points'])
    else:
         tau_list = np.linspace(exp_settings['Tau_min'], exp_settings['Tau_max'], exp_settings['Tau_points'])
    taus = np.round(tau_list, 9)
    
    indices = list(range(len(taus)))
    np.random.shuffle(indices)

    ## Start the main measurement loop 
    amp_int = np.zeros(len(taus))
    amps    = np.zeros((len(taus),card.samples))
    
    tstart = time.time()
    first_it = True

#    for tind in indices: #old version that has no explicit counter for the loop
    for lind in range(0,len(indices)):
        tind = indices[lind] #tind is the index in the tau array. lind is the index in the loop
        
        tau = taus[tind]
        print('Tau: {}'.format(tau))
        finalprog = loadprog
        finalprog = finalprog.replace('_tau_',str(tau))
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)
        
        amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot=first_it) 
        
#        #old version that errors because amp is a vector
#        amps[tind] = amp_full
#        amp_int[tind] = amp 
        
        #new version
        amps[tind] = amp_full
        amp_int[tind] = np.mean(amp)

        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(taus))
            first_it = False      

        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
#        plt.plot(taus*1e6, amps)#old version. Plotted a ton of raw data. amps is a matrix of the full time traces
#        plt.plot(taus*1e6, amp_int) #this plots all taus, including ones that aren't done yet
#        plt.plot(taus[indices[0:tind]]*1e6, amp_int[indices[0:tind]], 'o') #hopefully this does only the ones that have been done
        #nope. it takes the wrong cut. 
        plt.plot(taus[indices[0:lind+1]]*1e6, amp_int[indices[0:lind+1]], 'o') # this should take the right one
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
#    amp_guess = max(amps)-min(amps) # this crashes because amps is 2D
    amp_guess = np.max(amps)-np.min(amps)
#    offset_guess = np.mean(amps[-10:]) #this is also wierd. amps is 2D
    offset_guess = np.mean(amps[:,-50:]) #now this will take the end of each row of faw data

    fit_guess = [T1_guess, amp_guess, offset_guess]
    T1, amp, offset, fit_xvals, fit_yvals = fit_T1(taus, amp_int, fit_guess)
    fig3 = plt.figure(3)
    plt.clf() #adding a clear figure
#    plt.plot(fit_xvals*1e6, amps) #old version. amps is a matrix of the full time traces
#    plt.plot(fit_xvals*1e6, amp_int) #next problems. fit_xvals is not the same shape as taus
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