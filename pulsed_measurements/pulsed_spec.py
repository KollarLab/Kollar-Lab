import time
import numpy
import os
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import check_inputs, extract_data, remove_IQ_ellipse, configure_card, estimate_time, read_and_process

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'pulsed_spec'
    
    #Cavity parameters
    settings['CAVpower']        = -60
    settings['CAV_freq']        = 7e9
    
    #Qubit parameters
    settings['start_freq']      = 4.15*1e9  
    settings['stop_freq']       = 4.25*1e9 
    settings['freq_points']     = 50

    settings['start_power']     = -30
    settings['stop_power']      = -20
    settings['power_points']    = 5

    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 5e3
    
    #Measurement settings
    settings['Quasi_CW']    = False
    
    settings['subtract_background'] = True
    
    return settings

def pulsed_spec(instruments, settings):
    
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    exp_globals = settings['exp_globals']
    exp_settings = settings['exp_settings']

    ##Data saving and naming
    root_dir = exp_globals['root_folder']
    project_name = exp_globals['project_name']
    device_name = exp_globals['device_name']
    device_dir = os.path.join(root_dir, project_name, device_name)
    saveDir = userfuncs.saveDir(exp_globals, settings['meas_type'])
    stamp = userfuncs.timestamp()
    filename = exp_settings['scanname'] + '_' + stamp
    
    ##Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = exp_settings['CAVpower'] + CAV_Attenuation
    CAV_freq  = exp_settings['CAV_freq']
    
    ##Qubit settings
    start_freq  = exp_settings['start_freq']
    stop_freq   = exp_settings['stop_freq']
    freq_points = exp_settings['freq_points']
    freqs  = numpy.round(numpy.linspace(start_freq,stop_freq,freq_points),-3)
    
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    start_power  = exp_settings['start_power'] + Qbit_Attenuation
    stop_power   = exp_settings['stop_power'] + Qbit_Attenuation
    power_points = exp_settings['power_points']
    powers = numpy.round(numpy.linspace(start_power,stop_power,power_points),2)
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.IQ.Mod = 'On'

    #setting qubit generator to some safe starting point before we turn it on
    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20
    
    if settings['Quasi_CW']:
        qubitgen.IQ.Mod = 'Off'
    else:
        qubitgen.IQ.Mod = 'On'

    cavitygen.Output = 'On'
    qubitgen.Output = 'On'
    
    LO.Ref.Source = 'EXT'
    LO.Power = 12
    LO.Freq = CAV_freq - exp_globals['IF']
    LO.Output = 'On'
    
    configure_card(card, fullsettings)

    ##HDAWG settings
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\T1.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    q_pulse = exp_globals['qubit_pulse']
    m_pulse = exp_globals['measurement_pulse']
    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))
    loadprog = loadprog.replace('_tau_',str(q_pulse['delay']))
    loadprog = loadprog.replace('_qsigma_',str(q_pulse['sigma']))
    loadprog = loadprog.replace('_num_sigma_',str(q_pulse['sigma']))

    hdawg.AWGs[0].load_program(loadprog)

    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)
    
    ###########################################
    
    #Replace with extract data function
    data_window = int(settings['meas_window']*card.sampleRate)
    data_start  = int((settings['empirical_delay']+settings['init_buffer'])*card.sampleRate)
    back_start  = int((settings['empirical_delay']+settings['init_buffer']+settings['meas_window']+settings['pulse_buffer'])*card.sampleRate)
#    start_points = int((settings['meas_pos'] - card.settings['triggerDelay'] + settings['empirical_delay'])*card.sampleRate)
#    print(start_points)
#    start_points = int(1.2e-6*card.sampleRate)
    
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

        amps_full = numpy.zeros((len(freqs), len(xaxis) ))
        phases_full = numpy.zeros((len(freqs),len(xaxis) ))
        
        print('Current power:{}, max:{}'.format(powers[powerind]-Qbit_Attenuation, powers[-1]-Qbit_Attenuation))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]
            
            if powerind == 0 and find == 0:
                tstart = time.time()
                plot_single_IQ = True

            qubitgen.Freq = freq
            time.sleep(0.2)

            amp, phase, amp_full, phase_full = read_and_process(card, settings, plot_single_IQ)

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
            
#            # No mixer correction
#            Ip, Qp = I[0], Q[0]
            
            # no mixer correction, but average different segments together.
            Ip = numpy.mean(I, 0)
            Qp = numpy.mean(Q, 0)
            
            DC_I = numpy.mean(Ip[back_start:])
            DC_Q = numpy.mean(Qp[back_start:])
            Idat = Ip-DC_I
            Qdat = Qp-DC_Q
#            Idat = Ip
#            Qdat = Qp
            
            amps[find,:]   = amp
            phases[find,:] = phase

            amps_full[find,:]   = amp_full
            phases_full[find,:] = phase_full

        powerslice = numpy.mean(amps[:,:], axis=1)
        phaseslice = numpy.mean(phases[:,:], axis=1)
        
        
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

        yaxis = powers[0:powerind+1] - Qbit_Attenuation
        labels = ['Freq (GHz)', 'Power (dBm)']
        simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) 

        full_time = {}
        full_time['xaxis'] = xaxis_us
        full_time['mags'] = amps_full
        full_time['phases'] = phases_full

        single_time = {}
        single_time['xaxis'] = xaxis_us
        single_time['mag'] = amp_full
        single_time['phase'] = phase_full

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Power: {}dBm'.format(power-Qbit_Attenuation)
        
        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','xaxis_us', 'full_time'], locals(), expsettings=settings)
        
        simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces\n'+filename, time_labels, identifier, fig_num=2)

    t2 = time.time()
    
    print('elapsed time = ' + str(t2-t1))

    simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) 
    plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
    simplescan_plot(full_time, single_time, freqs/1e9, 'Raw_time_traces', time_labels, identifier=identifier, fig_num=2)
    plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)
       
    userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','xaxis_us', 'full_time'], locals(), expsettings=settings)
    
    cavitygen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.Output = 'Off'
