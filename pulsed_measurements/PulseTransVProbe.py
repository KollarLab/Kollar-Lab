# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 19:12:43 2020

@author: Kollarlab
"""

import time
import numpy 
import os
import matplotlib.pyplot as plt

import userfuncs
from plotting_tools import simplescan_plot


pi = numpy.pi

def remove_IQ_ellipse(Is, Qs, axes, center, phi):
   
    Isprime = numpy.cos(phi)*(Is-center[0]) + numpy.sin(phi)*(Qs-center[1]) + center[0]
    Qsprime = -numpy.sin(phi)*(Is-center[0]) + numpy.cos(phi)*(Qs-center[1]) 
    Qsprime = Qsprime*axes[0]/axes[1] + center[1]
    return Isprime, Qsprime

center = [0.027, -0.034]
phi = 0.044*2*pi/180
axes = [0.024, 0.018]

def GetDefaultSettings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'PulsedTransVProbe'
    settings['project_dir'] = r'Z:\Data\defaultdir'
    
    #Sweep parameters
    settings['start_freq']      = 4.15*1e9  
    settings['stop_freq']       = 4.25*1e9 
    settings['freq_points']     = 50

    settings['CAV_Attenuation'] = 30
    settings['Qbit_Attenuation'] = 10

    settings['cavity_power']    = -60
    
    settings['probe_freq']      = 3.2*1e9
    settings['start_power']     = -20
    settings['stop_power']      = 10
    settings['power_points']    = 31

    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 5e3
    settings['activeChannels']   = [1,2]
    settings['channelRange']     = 0.5
    settings['sampleRate']       = 2e9/8
    settings['timeout']          = 30
    
    ##Pulse settings
    settings['meas_pos'] = 80e-6
    settings['meas_window'] = 10e-6
    return settings

def PulsedTransVProbe(instruments, settings):
    ##Instruments used
    cavitygen = instruments['cavitygen']
    qubitgen = instruments['qubitgen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    ##Data saving and naming
    saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp = userfuncs.timestamp()
    filename = settings['scanname'] + '_' + stamp
    
    ##Sweep settings
    start_freq  = settings['start_freq'] 
    stop_freq   = settings['stop_freq'] 
    freq_points = settings['freq_points']
    freqs  = numpy.round(numpy.linspace(start_freq,stop_freq,freq_points),-3)
    
    start_power  = settings['start_power'] + settings['Qbit_Attenuation']
    stop_power   = settings['stop_power'] + settings['Qbit_Attenuation']
    power_points = settings['power_points']
    powers = numpy.round(numpy.linspace(start_power,stop_power,power_points),2)
    
    ## Generator settings
    qubitgen.Freq   = settings['probe_freq']
    qubitgen.Power  = powers[0]
    qubitgen.IQ.Mod = 'Off'
    
    
    cavitygen.Freq   = freqs[0]
    cavitygen.Power  = settings['cavity_power'] + settings['CAV_Attenuation']
    cavitygen.IQ.Mod = 'On'

    cavitygen.Output = 'On'
    qubitgen.Output = 'On'
    
    LO.Ref.Source = 'EXT'
    LO.Power = 12
    LO.Freq = freqs[0]
    LO.Output = 'On'
    
    ##Card settings
    meas_samples = settings['sampleRate']*settings['meas_window']

    card.averages       = settings['averages']
    card.segments       = settings['segments']
    card.sampleRate     = settings['sampleRate']
    card.activeChannels = settings['activeChannels']
    card.triggerDelay   = settings['meas_pos']
    card.timeout        = settings['timeout']
    card.samples        = int(meas_samples*1.2)
    card.channelRange   = settings['channelRange']
    card.SetParams()
    
    ##HDAWG settings
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='True')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\Control\HDAWG_sequencer_codes\pulsedtrans.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    loadprog = loadprog.replace('_max_time_', str(settings['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(settings['meas_window']))
    
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)

    data_window = int(settings['meas_window']*card.sampleRate)
    buffer = int(1e-6*card.sampleRate)
    
    xaxis = (numpy.array(range(card.samples))/card.sampleRate)
    xaxis_us = xaxis*1e6
    
    powerdat = numpy.zeros((len(powers), len(freqs)))
    phasedat = numpy.zeros((len(powers), len(freqs)))
    
    t1 = time.time()
    
    for powerind in range(len(powers)):
        power = powers[powerind]
        qubitgen.Power = powers[powerind]
        time.sleep(0.2)
        
        amps = numpy.zeros((len(freqs), len(xaxis) ))

        phases = numpy.zeros((len(freqs),len(xaxis) ))
        
        Is = numpy.zeros((len(freqs), len(xaxis) ))
        Qs = numpy.zeros((len(freqs), len(xaxis) ))
        
        print('Current power:{}, max:{}'.format(powers[powerind] - settings['Qbit_Attenuation'], powers[-1]))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]
            
            if powerind == 0 and find == 0:
                tstart = time.time()
            cavitygen.Freq = freq
            LO.freq = freq
            LO.output = 'On'
            time.sleep(0.2)
            
            card.ArmAndWait()
            
            I,Q = card.ReadAllData()
            
            if powerind == 0 and find == 0:
                tstop = time.time()
                singlePointTime = tstop-tstart
                
                estimatedTime = singlePointTime*len(freqs)*len(powers)
                print('    ')
                print('estimated time for this scan : ' + str(numpy.round(estimatedTime/60, 1)) + ' minutes')
                print('estimated time for this scan : ' + str(numpy.round(estimatedTime/60/60, 2)) + ' hours')
                print('    ')
            
            # mixer correction
#            Ip, Qp = remove_IQ_ellipse(I[0], Q[0], axes, center, phi)
            
            # No mixer correction
            Ip, Qp = I[0], Q[0]
            
            DC_I = numpy.mean(Ip[-50:])
            DC_Q = numpy.mean(Qp[-50:])
            Idat = Ip-DC_I
            Qdat = Qp-DC_Q
            
            amp = numpy.sqrt(Idat**2+Qdat**2)
            phase = numpy.arctan2(Qdat, Idat)*180/numpy.pi
            
            amps[find,:] = amp
            phases[find,:] = phase
            Is[find,:] = Idat 
            Qs[find,:] = Qdat
        
#        powerslice = numpy.mean(amps[:,buffer:int(data_window+buffer)], axis=1)/(10**(power/20))
        powerslice = numpy.mean(amps[:,buffer:int(data_window+buffer)], axis=1) #not sweeping cavity power. No need to devide.
        phaseslice = numpy.mean(phases[:,buffer:int(data_window+buffer)], axis=1)
        
        powerdat[powerind,:] = powerslice
        phasedat[powerind,:] = phaseslice
        
        full_data = {}
        full_data['xaxis'] = freqs/1e9
        full_data['mags'] = powerdat[0:powerind+1]
        full_data['phases'] = phasedat[0:powerind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mag'] = powerslice
        single_data['phase'] = phaseslice

        yaxis = powers[0:powerind+1] - settings['Qbit_Attenuation']
        labels = ['Freq (GHz)', 'Power (dBm)']
        simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) 

        full_time = {}
        full_time['xaxis'] = xaxis_us
        full_time['mags'] = amps
        full_time['phases'] = phases

        single_time = {}
        single_time['xaxis'] = xaxis_us
        single_time['mag'] = amp
        single_time['phase'] = phase

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Power: {}dBm'.format(power-settings['Qbit_Attenuation'])
        
        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','xaxis_us', 'full_time'], locals(), expsettings=settings)
        
        simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces', time_labels, identifier, fig_num=2)
        
        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','xaxis_us'], locals(), expsettings=settings)
        
#        ###plot the full power dependent data 
#        if len(powers) > 1:
##    
##            plt.suptitle(filename)
#            fig = plt.figure(1, figsize=(13,8))
#            plt.clf()
#            pulsed_debug(fig, freqs, powers[0:powerind+2], powerdat[0:powerind+1], phasedat[0:powerind+1], powerslice, phaseslice, filename, power)
#            fig.canvas.draw()
#            fig.canvas.flush_events()
#            plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
#        else:
#            if powerind == 1:
#                fig = plt.figure(1)
#                plt.clf()
        
    
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-t1))
    
    simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) 
    plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
    simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces', time_labels, identifier=identifier, fig_num=2)
    plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)
    
#    ###plot the full power dependent data 
#    if len(powers) > 1:
#
#        fig = plt.figure(1, figsize=(13,8))
#        plt.clf()
#        pulsed_debug(fig, freqs, powers[0:powerind+1], powerdat[0:powerind+1], phasedat[0:powerind+1], powerslice, phaseslice, filename, power)
#        fig.canvas.draw()
#        fig.canvas.flush_events()
#            
#        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
#    else:
#        fig = plt.figure(1)
#        plt.clf()
#    
#    
#    #plot the raw data at (the last?) power
#    fig = plt.figure(2, figsize=(13,8))
#    plt.clf()
#    ax = plt.subplot(1,2,1)
#    base_raw_time_plot_spec(fig, ax, times = xaxis_us, ydata = amps, ys = freqs/1e9, ylabel = 'Freq (GHz)', zlabel = 'Volts', scanname = 'Raw trace amps', scanformat = '')
#    
#    
#    
#    ax = plt.subplot(1,2,2)
#    base_raw_time_plot_spec(fig, ax, times = xaxis_us, ydata = phases, ys = freqs/1e9, ylabel = 'Freq (GHz)', zlabel = 'Phase (deg)', scanname = 'Raw trace amps', scanformat = '')
    
    
    
#    #plot the last row as a trace.
#    ax = plt.subplot(2,2,3)
#    plt.plot(freqs/1e9, powerslice)
#    plt.xlabel('freqs (GHz)')
#    plt.ylabel('spec voltage')
#    
#    #plot the last row as a trace.
#    ax = plt.subplot(2,2,4)
#    plt.plot(freqs/1e9, phaseslice)
#    plt.xlabel('freqs (GHz)')
#    plt.ylabel('spec phase')
    
    #datatip()
    
#    plt.suptitle(filename)
#    plt.savefig(os.path.join(saveDir, filename+'_singleRawData.png'), dpi = 150)
