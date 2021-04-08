# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:28:57 2021

@author: Kollarlab
"""

import time
import numpy
import os
import matplotlib.pyplot as plt

import userfuncs
from plotting_tools import simplescan_plot

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'rabi_chevron'
    settings['meas_type'] = 'rabi_chevron'
    settings['project_dir'] = r'Z:\Data\defaultdir'
    
    #Cavity parameters
    settings['CAVpower']        = -45
    settings['CAV_freq']        = 5e9
    settings['Q_power']         = -20
    
    #Qubit parameters
    settings['start_freq']      = 4*1e9  
    settings['stop_freq']       = 5*1e9 
    settings['freq_points']     = 50

    settings['start_time']     = 10e-9
    settings['stop_time']      = 500e-9
    settings['time_points']    = 40


    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 5e3
    settings['activeChannels']   = [1,2]
    settings['channelRange']     = 0.5
    settings['sampleRate']       = 2e9/8
    #settings['trigger_buffer']   = 0e-6
    settings['meas_window']      = 10e-6
    settings['timeout']          = 30
    
    ##Pulse settings
    settings['meas_pos'] = 80e-6
    settings['empirical_delay']  = 1e-6
    settings['pulse_delay'] = 200e-9
    settings['pulse_width'] = 80e-9
    
    return settings

def rabi_chevron(instruments, settings):
    
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    ##Data saving and naming
    saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp = userfuncs.timestamp()
    filename = settings['scanname'] + '_' + stamp
    
    ##Cavity settings
    CAV_Attenuation = settings['CAV_Attenuation']
    CAV_power = settings['CAVpower'] + CAV_Attenuation
    CAV_freq  = settings['CAV_freq']
    
    ##Qubit settings
    start_freq  = settings['start_freq']
    stop_freq   = settings['stop_freq']
    freq_points = settings['freq_points']
    freqs  = numpy.round(numpy.linspace(start_freq,stop_freq,freq_points),-3)
    
    Qbit_Attenuation = settings['Qbit_Attenuation']
    Qbit_power   = settings['Q_power'] + Qbit_Attenuation
    start_time   = settings['start_time']
    stop_time    = settings['stop_time']
    time_points  = settings['time_points']
    times = numpy.linspace(start_time,stop_time,time_points)
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.IQ.Mod = 'On'

    qubitgen.Freq   = 4e9
    qubitgen.Power  = Qbit_power
    
    qubitgen.IQ.Mod = 'On'

    cavitygen.Output = 'On'
    qubitgen.Output = 'On'
    
    LO.Ref.Source = 'EXT'
    LO.Power = 12
    LO.Freq = CAV_freq
    LO.Output = 'On'
    
    ##Card settings
    meas_samples = settings['sampleRate']*settings['meas_window']

    card.averages       = settings['averages']
    card.segments       = settings['segments']
    card.sampleRate     = settings['sampleRate']
    card.activeChannels = settings['activeChannels']
    card.triggerDelay   = settings['meas_pos'] + settings['empirical_delay']
    card.timeout        = settings['timeout']
    card.samples        = int(meas_samples*1.2)
    card.channelRange   = settings['channelRange']
    card.SetParams()

    ##HDAWG settings
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\Control\HDAWG_sequencer_codes\chevron.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    loadprog = loadprog.replace('_max_time_', str(settings['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(settings['meas_window']))
    loadprog = loadprog.replace('_meas_delay_',str(settings['pulse_delay']))
    loadprog = loadprog.replace('_min_q_time_',str(settings['pulse_width']))
    
    ###########################################
    
    data_window = int(meas_samples)
    xaxis = (numpy.array(range(card.samples))/card.sampleRate)
    xaxis_us = xaxis*1e6
    
    timedat  = numpy.zeros((len(times), len(freqs)))
    phasedat = numpy.zeros((len(times), len(freqs)))
    
    t1 = time.time()
    
    for timeind in range(len(times)):
        hold_time = times[timeind]
        finalprog = loadprog.replace('_hold_time_',str(hold_time))
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)
        
        amps = numpy.zeros((len(freqs), len(xaxis) ))
        phases = numpy.zeros((len(freqs),len(xaxis) ))
        
        Is = numpy.zeros((len(freqs), len(xaxis) ))
        Qs = numpy.zeros((len(freqs), len(xaxis) ))
        
        print('Current hold time:{}, max:{}'.format(hold_time, times[-1]))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]
            
            if timeind == 0 and find == 0:
                tstart = time.time()
            qubitgen.Freq = freq
            time.sleep(0.2)
            
            card.ArmAndWait()
            
            I,Q = card.ReadAllData()
            
            if timeind == 0 and find == 0:
                tstop = time.time()
                singlePointTime = tstop-tstart
                
                estimatedTime = singlePointTime*len(freqs)*len(times)
                print('    ')
                print('estimated time for this scan : ' + str(numpy.round(estimatedTime/60, 1)) + ' minutes')
                print('estimated time for this scan : ' + str(numpy.round(estimatedTime/60/60, 2)) + ' hours')
                print('    ')
            
           
            # no mixer correction, but average different segments together.
            Ip = numpy.mean(I, 0)
            Qp = numpy.mean(Q, 0)
            
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
        
        timeslice = numpy.mean(amps[:,:int(data_window)], axis=1)
        phaseslice = numpy.mean(phases[:,:int(data_window)], axis=1)
        
        timedat[timeind,:] = timeslice
        phasedat[timeind,:] = phaseslice

        full_data = {}
        full_data['xaxis'] = freqs/1e9
        full_data['mags'] = timedat[0:timeind+1]
        full_data['phases'] = phasedat[0:timeind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mag'] = timeslice
        single_data['phase'] = phaseslice

        yaxis = times[0:timeind+1]
        labels = ['Freq (GHz)', 'Hold time (us)']
        simplescan_plot(full_data, single_data, yaxis*1e6, filename, labels, identifier='', fig_num=1) 

        full_time = {}
        full_time['xaxis'] = xaxis_us
        full_time['mags'] = amps
        full_time['phases'] = phases

        single_time = {}
        single_time['xaxis'] = xaxis_us
        single_time['mag'] = amp
        single_time['phase'] = phase

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Hold time: {}us'.format(hold_time*1e6)
        
        userfuncs.SaveFull(saveDir, filename, ['times','freqs', 'timedat', 'phasedat','xaxis','xaxis_us', 'full_time'], locals(), expsettings=settings)
        
        simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces', time_labels, identifier, fig_num=2)

    t2 = time.time()
    
    print('elapsed time = ' + str(t2-t1))

    simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) 
    plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
    simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces', time_labels, identifier=identifier, fig_num=2)
    plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)
       
    userfuncs.SaveFull(saveDir, filename, ['times','freqs', 'timedat', 'phasedat','xaxis','xaxis_us', 'full_time'], locals(), expsettings=settings)
    
    cavitygen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.Output = 'Off'