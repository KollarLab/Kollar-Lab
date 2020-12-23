import time
import numpy 
import os
import matplotlib.pyplot as plt

import userfuncs
from ellipse_fitting import remove_IQ_ellipse
from plotting_tools import simplescan_plot


pi = numpy.pi

center = [0.027, -0.034]
phi = 0.044*2*pi/180
axes = [0.024, 0.018]

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'PulsedTrans'
    settings['project_dir'] = r'Z:\Data\HouckQuadTransmon'
    
    #Sweep parameters
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
    settings['timeout']          = 30
    
    ##Pulse settings
    settings['meas_pos'] = 80e-6
    settings['empirical_delay']  = 1e-6
    settings['meas_window'] = 10e-6
    return settings

def pulsed_trans(instruments, settings):
    ##Instruments used
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    ##Data saving and naming
    saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp = userfuncs.timestamp()
    filename = settings['scanname'] + '_' + stamp

    scanname = settings['scanname']

    ##Sweep settings
    start_freq  = settings['start_freq']
    stop_freq   = settings['stop_freq']
    freq_points = settings['freq_points']
    freqs  = numpy.round(numpy.linspace(start_freq,stop_freq,freq_points),-3)
    
    start_power  = settings['start_power']
    stop_power   = settings['stop_power']
    power_points = settings['power_points']
    powers = numpy.round(numpy.linspace(start_power,stop_power,power_points),2)
    
    ## Generator settings
    cavitygen.Freq   = freqs[0]
    cavitygen.Power  = powers[0]
    cavitygen.IQ.Mod = 'On'

    cavitygen.Output = 'On'
    
#    LO.Ref.Source = 'EXT'
#    LO.Power = 12
#    LO.Freq = CAV_freq
#    LO.Output = 'On'
    LO.inst.write('ROSC EXT')
    LO.inst.write('SENS:SWE:TYPE CW')
    LO.inst.write('SOUR:POW 12')
    LO.Freq = freqs[0]
    LO.output = 'On' 
    
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
    
    xaxis = (numpy.array(range(card.samples))/card.sampleRate)
    xaxis_us = xaxis*1e6
    
    LO.inst.write('SENS:SWE:TYPE CW')
    LO.inst.write('SOUR:POW 12') 
    
    powerdat = numpy.zeros((len(powers), len(freqs)))
    phasedat = numpy.zeros((len(powers), len(freqs)))
    
    t1 = time.time()
    
    for powerind in range(len(powers)):
        power = powers[powerind]
        cavitygen.Power = powers[powerind]
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
            Ip, Qp = remove_IQ_ellipse(I[0], Q[0], axes, center, phi)
            
            # No mixer correction
#            Ip, Qp = I[0], Q[0]
            
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
        
        powerslice = numpy.mean(amps[:,:int(data_window)], axis=1)/(10**(power/20))
        phaseslice = numpy.mean(phases[:,:int(data_window)], axis=1)
        
        powerdat[powerind,:] = powerslice
        phasedat[powerind,:] = phaseslice
        
        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','xaxis_us'], locals(), expsettings=settings)

        full_data = {}
        full_data['xaxis'] = freqs/1e9
        full_data['mags'] = powerdat[0:powerind+1]
        full_data['phases'] = phasedat[0:powerind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mag'] = powerslice
        single_data['phase'] = phaseslice

        yaxis = powers[0:powerind+1]
        labels = ['Freq (GHz)', 'Power (dBm)']
        simplescan_plot(full_data, single_data, yaxis, scanname, labels, identifier='', fig_num=1) 

        full_time = {}
        full_time['xaxis'] = xaxis_us
        full_time['mags'] = amps
        full_time['phases'] = phases

        single_time = {}
        single_time['xaxis'] = xaxis_us
        single_time['mag'] = amp
        single_time['phase'] = phase

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Power: {}dBm'.format(power)
        simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces', time_labels, identifier, fig_num=2)
    
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-t1))

    simplescan_plot(full_data, single_data, yaxis, scanname, labels, identifier='', fig_num=1) 
    plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
    simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces', time_labels, identifier='', fig_num=2)
    plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)

    userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','xaxis_us'], locals(), expsettings=settings)