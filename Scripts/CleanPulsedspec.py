import time
from mplcursors import cursor as datatip
import numpy
import os
from VNAplottingTools import base_power_plot, base_raw_time_plot_spec
import matplotlib.pyplot as plt
import userfuncs

def GetDefaultSettings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'PulsedSpec'
    settings['project_dir'] = r'Z:\Data\HouckQuadTransmon'
#    settings['saveDir']  = r'Z:\Data\HouckQuadTransmon\PulsedSpec\20201206'
    
    #Cavity parameters
    settings['CAVpower']        = -18
    settings['CAV_freq']        = 8.126e9
    
    #Qubit parameters
    settings['start_freq']      = 4.15*1e9  
    settings['stop_freq']       = 4.25*1e9 
    settings['freq_points']     = 50

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
    #settings['trigger_buffer']   = 0e-6
    settings['meas_window']      = 10e-6
    settings['timeout']          = 30
    
    ##Pulse settings
    settings['meas_pos'] = 80e-6
    settings['pulse_delay'] = 200e-9
    settings['pulse_width'] = 80e-9
    
    return settings

def PulsedSpec(instruments, settings):
    
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
    CAV_power = settings['CAVpower']
    CAV_freq  = settings['CAV_freq']
    
    ##Qubit settings
    start_freq  = settings['start_freq']
    stop_freq   = settings['stop_freq']
    freq_points = settings['freq_points']
    freqs  = numpy.round(numpy.linspace(start_freq,stop_freq,freq_points),-3)
    
    start_power  = settings['start_power']
    stop_power   = settings['stop_power']
    power_points = settings['power_points']
    powers = numpy.round(numpy.linspace(start_power,stop_power,power_points),2)
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.IQ.Mod = 'On'

    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20
    qubitgen.IQ.Mod = 'On'

    cavitygen.Output = 'On'
    qubitgen.Output = 'On'
    
    LO.Ref.Source = 'EXT'
    LO.Power = 12
    LO.Freq = CAV_freq
    LO.Output = 'On'
    #vna.inst.write('ROSC EXT')
    #vna.inst.write('SENS:SWE:TYPE CW')
    #vna.inst.write('SOUR:POW 12')
    #vna.inst.write('SOUR:FREQ:CW {}'.format(CAV_freq))
    #vna.output = 'On' 
    
    ##Card settings
    meas_samples = settings['sampleRate']*settings['meas_window']

    card.averages       = settings['averages']
    card.segments       = settings['segments']
    card.sampleRate     = settings['sampleRate']
    card.activeChannels = settings['activeChannels']
    card.triggerDelay   = settings['meas_pos']
    card.timeout        = settings['timeout']
    card.samples        = int(meas_samples*1.5)
    card.channelRange   = settings['channelRange']
    card.SetParams()

    ##HDAWG settings
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='True')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\Control\HDAWG_sequencer_codes\T1.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    loadprog = loadprog.replace('_max_time_', str(settings['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(settings['meas_window']))
    loadprog = loadprog.replace('_tau_',str(settings['pulse_delay']))
    loadprog = loadprog.replace('_qwidth_',str(settings['pulse_width']))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)

    ###########################################
    
    data_window = int(meas_samples)
    
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
        
        print('Current power:{}, max:{}'.format(powers[powerind], powers[-1]))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]
            
            if powerind == 0 and find == 0:
                tstart = time.time()
            qubitgen.Freq = freq
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
            
            DC_I = numpy.mean(I[0][-50:])
            DC_Q = numpy.mean(Q[0][-50:])
            
            Idat = I[0]-DC_I
            Qdat = Q[0]-DC_Q
            
            amp = numpy.sqrt(Idat**2+Qdat**2)
            phase = numpy.arctan2(Qdat, Idat)*180/numpy.pi
            
            amps[find,:] = amp
            phases[find,:] = phase
            Is[find,:] = Idat 
            Qs[find,:] = Qdat
        
        powerslice = numpy.mean(amps[:,0:int(data_window)], axis=1)
        phaseslice = numpy.mean(phases[:,0:int(data_window)], axis=1)
        
        powerdat[powerind,:] = powerslice
        phasedat[powerind,:] = phaseslice
        
        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','xaxis_us'], locals(), expsettings=settings)
        
        ###plot the full power dependent data 
        if len(powers) > 1:
            fig = plt.figure(1)
            plt.clf()
            ax = plt.subplot(1,2,1)
            base_power_plot(fig, ax, freqs, powerdat[0:powerind+1,:], powers[0:powerind+1], 'Freq sweep', 'mag', HWattenuation = -10)  
            
            ax = plt.subplot(1,2,2)
            base_power_plot(fig, ax, freqs, phasedat[0:powerind+1,:], powers[0:powerind+1], 'Freq sweep', 'phase', HWattenuation = -10)
            
            plt.suptitle(filename)
            fig.canvas.draw()
            fig.canvas.flush_events()
    #        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
        else:
            if powerind == 1:
                fig = plt.figure(1)
                plt.clf()
        
    
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-t1))
    
    ###plot the full power dependent data 
    if len(powers) > 1:
        fig = plt.figure(1)
        plt.clf()
        ax = plt.subplot(1,2,1)
        base_power_plot(fig, ax, freqs, powerdat[0:powerind+1,:], powers[0:powerind+1], 'Freq sweep', 'mag', HWattenuation = -10)  
        
        ax = plt.subplot(1,2,2)
        base_power_plot(fig, ax, freqs, phasedat[0:powerind+1,:], powers[0:powerind+1], 'Freq sweep', 'phase', HWattenuation = -10)
        
        plt.suptitle(filename)
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
    else:
        fig = plt.figure(1)
        plt.clf()
    
    
    #plot the raw data at (the last?) power
    fig = plt.figure(2)
    plt.clf()
    ax = plt.subplot(2,2,1)
    #base_power_plot_spec(fig, ax, xaxis, amps, freqs, 'Raw trace amps','',0)
    base_raw_time_plot_spec(fig, ax, times = xaxis_us, ydata = amps, ys = freqs/1e9, scanname = 'Raw trace amps', scanformat = '')
    plt.xlabel('time (us)')
    plt.ylabel('frequency (GHz)')
    
    
    ax = plt.subplot(2,2,2)
    #base_power_plot_spec(fig, ax, xaxis, phases, freqs, 'Raw trace phases','',0)
    base_raw_time_plot_spec(fig, ax, times = xaxis_us, ydata = phases, ys = freqs/1e9, scanname = 'Raw trace amps', scanformat = '')
    plt.xlabel('time (us)')
    plt.ylabel('frequency (GHz)')
    
    
    #plot the last row as a trace.
    ax = plt.subplot(2,2,3)
    plt.plot(freqs/1e9, powerslice)
    plt.xlabel('freqs (GHz)')
    plt.ylabel('spec voltage')
    
    #plot the last row as a trace.
    ax = plt.subplot(2,2,4)
    plt.plot(freqs/1e9, phaseslice)
    plt.xlabel('freqs (GHz)')
    plt.ylabel('spec phase')
    
    plt.suptitle(filename)
    plt.savefig(os.path.join(saveDir, filename+'_singleRawData.png'), dpi = 150)