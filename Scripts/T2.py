# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:52:18 2020

@author: Kollarlab
"""
import time
import os
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
from plotting_tools import base_power_plot_imshow
from UserFits import fit_T2

def GetDefaultSettings():
    settings = {}
    
    settings['scanname'] = 'T2_meas'
    settings['project_dir'] = r'Z:\Data\HouckQuadTransmon'
    settings['meas_type'] = 'Tmeas'

    
    settings['Q_Freq'] = 4.20431e9
    settings['Q_Power'] = -11
    settings['CAV_Freq'] = 8.126e9
    settings['CAV_Power'] = -18
    
    Meas_pos = 80e-6
    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 25e3
    settings['activeChannels']   = [1,2]
    settings['sampleRate']       = 2e9/8
    settings['trigger_buffer']   = Meas_pos
    settings['meas_window']      = 20e-6
    
    settings['Tau_min'] = 200e-9
    settings['Tau_max'] = 30e-6
    settings['Tau_points'] = 5
    
    #Pulse settings
    settings['Measurement_pos'] = Meas_pos
    settings['wait_time'] = 200e-9
    settings['pulse_width'] = 80e-9
    
    return settings

def meas_T2(instruments, settings):
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    Q_Freq    = settings['Q_Freq'] + settings['detuning']
    Q_Power   = settings['Q_Power']
    CAV_Freq  = settings['CAV_Freq']
    CAV_Power = settings['CAV_Power']
    
    CAV_Attenuation  = settings['CAV_Attenuation']
    Qbit_Attenuation = settings['Qbit_Attenuation']
#    extra_shift = 500e3
    detuning = settings['detuning']
#    Q_Freq = 4.20431e9 + detuning + extra_shift
    
    stamp = userfuncs.timestamp()
    filename = settings['scanname'] + '_' + stamp
    saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    
    ## Generator settings
#    LO.Ref.Source = 'EXT'
#    LO.Power = 12
#    LO.Freq = CAV_Freq
#    LO.Output = 'On'
    LO.inst.write('ROSC EXT')
    LO.inst.write('SENS:SWE:TYPE CW')
    LO.inst.write('SOUR:FREQ:CW {}'.format(CAV_Freq))
    LO.inst.write('SOUR:POW 12')
    LO.output = 'On' 
    
    cavitygen.Freq   = CAV_Freq
    cavitygen.Power  = CAV_Power + CAV_Attenuation
    cavitygen.IQ.Mod = 'On'
    
    qubitgen.Freq   = Q_Freq
    qubitgen.Power  = Q_Power + Qbit_Attenuation
    qubitgen.IQ.Mod = 'On'
    
    cavitygen.Output = 'On'
    qubitgen.Output = 'On'
    
    ## Card Settings
    meas_samples = settings['sampleRate']*settings['meas_window']
    
    card.averages       = settings['averages']
    card.segments       = settings['segments']
    card.sampleRate     = settings['sampleRate']
    card.activeChannels = settings['activeChannels']
    card.triggerDelay   = settings['trigger_buffer']+1e-6
    card.timeout        = 30
    card.samples        = int(meas_samples*2.5)
    card.channelRange   = 2.5
    card.SetParams()
    
    data_window = int(meas_samples) #for 100us measurement pulses 
    
    xaxis = (np.array(range(card.samples))/card.sampleRate)
    xaxis_us = xaxis*1e6 + settings['trigger_buffer']
    
    ## HDAWG
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='True')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\Control\HDAWG_sequencer_codes\T2.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    loadprog = loadprog.replace('_max_time_', str(settings['Measurement_pos']))
    loadprog = loadprog.replace('_meas_window_', str(settings['meas_window']))
    loadprog = loadprog.replace('_wait_time_', str(settings['wait_time']))
    loadprog = loadprog.replace('_qwidth_', str(settings['pulse_width']))
    
    taus = np.linspace(settings['Tau_min'],settings['Tau_max'] , settings['Tau_points'] )
    
    plt.figure(1)
    plt.clf()
    amp_int = np.zeros(len(taus))
    amps    = np.zeros((len(taus),card.samples))
    
    tstart = time.time()
    
    for tind in range(len(taus)):
            
        tau = taus[tind]
        print('Tau: {}'.format(tau))
        finalprog = loadprog
        finalprog = finalprog.replace('_tau_',str(tau))
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)
        
        card.ArmAndWait()
        I,Q = card.ReadAllData()
        
        DC_I = np.mean(I[0][-data_window:])
        DC_Q = np.mean(Q[0][-data_window:])
        
        Idat = I[0]-DC_I
        Qdat = Q[0]-DC_Q
        
        amp = np.sqrt(Idat**2+Qdat**2)
        
        amps[tind] = amp
        amp_int[tind] = np.mean(amp[0:int(data_window)])
        if tind == 0:
            tstop = time.time()
            singlePointTime = tstop-tstart
            
            estimatedTime = singlePointTime*len(taus)
            print('    ')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60, 1)) + ' minutes')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60/60, 2)) + ' hours')
            print('    ')
    
    userfuncs.SaveFull(saveDir, filename, ['taus','xaxis_us', 'amps', 'amp_int'], locals(), expsettings=settings)
    
    t2 = time.time()
    
    print('Elapsed time: {}'.format(t2-tstart))
    
    fig = plt.figure(1,figsize=(13,8))
    plt.clf()
        
    ax = plt.subplot(1,1,1)
    base_power_plot_imshow(fig, ax, xaxis_us, taus*1e6, amps, ['Time (us)', 'Tau (us)', 'Amp'], attenuation=0)
    
    plt.suptitle(filename)
    plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
    
    fig = plt.figure(2,figsize=(13,8))
    plt.clf()
    tau0 = 5e-9
    offset0 = np.mean(amp_int[-10])
    amp0 = np.mean(max(amp_int)-offset0)
    freq0 = detuning
    phi0 = 0
    
    fit_guess = [tau0, amp0, offset0, freq0, phi0]
    try:
        T2, detuning = fit_T2(taus, amp_int, fit_guess)
    except:
        print('T2 fit did not converge')
        T2 = 0
        detuning = 0
    
    plt.suptitle(filename)
    plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
    
    
    return T2, detuning