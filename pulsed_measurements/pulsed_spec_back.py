'''
8-25-21 AK modifying to normalize the amplitudes to the drive power. Undid that.

9-2-21 AK made it return the data

'''


import os
import time
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process
from utility.scheduler import scheduler

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'pulsed_spec'
    
    #Cavity parameters
    settings['CAV_Power']        = -60
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
    
    return settings

def pulsed_spec(instruments, settings):
    
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    ##Data saving and naming
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp
    
    ##Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = exp_settings['CAV_Power'] + CAV_Attenuation
    CAV_freq  = exp_settings['CAV_freq']
    
    ##Qubit settings
    start_freq  = exp_settings['start_freq']
    stop_freq   = exp_settings['stop_freq']
    freq_points = exp_settings['freq_points']
    if exp_settings['reverse']:
        freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points)[::-1],-3)
    else:
        freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
    
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    start_power  = exp_settings['start_power'] + Qbit_Attenuation
    stop_power   = exp_settings['stop_power']  + Qbit_Attenuation
    power_points = exp_settings['power_points']
    powers = np.round(np.linspace(start_power,stop_power,power_points),2)
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power

    #setting qubit generator to some safe starting point before we turn it on
    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20
    
    if exp_settings['Quasi_CW']:
        qubitgen.disable_pulse()
        qubitgen.disable_IQ()
    else:
        qubitgen.enable_pulse()
        qubitgen.enable_IQ()

    cavitygen.Output = 'On'
    qubitgen.Output  = 'On'
    
    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_pulse()
    cavitygen.output = 'On'
    
    configure_card(card, settings)

    ##HDAWG settings
    configure_hdawg(hdawg, settings)
    
#    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\T1.cpp",'r')
#    rawprog  = progFile.read()
#    loadprog = rawprog
#    progFile.close()
#    
#    q_pulse = exp_globals['qubit_pulse']
#    m_pulse = exp_globals['measurement_pulse']
#    loadprog = loadprog.replace('_max_time_', str(m_pulse['meas_pos']))
#    loadprog = loadprog.replace('_meas_window_', str(m_pulse['meas_window']))
#    loadprog = loadprog.replace('_tau_',str(q_pulse['delay']))
#    loadprog = loadprog.replace('_qsigma_',str(q_pulse['sigma']))
#    loadprog = loadprog.replace('_num_sigma_',str(q_pulse['num_sigma']))
#    loadprog = loadprog.replace('_piAmp_',str(q_pulse['piAmp']))
    
#    hdawg.AWGs[0].load_program(loadprog)
#
#    hdawg.AWGs[0].run_loop()
#    time.sleep(0.1)
    
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    q_pulse = exp_globals['qubit_pulse']
    
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    
    awg_sched = scheduler(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='Qubit_I')
    awg_sched.add_analog_channel(2, name='Qubit_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    
    qubit_I       = awg_sched.analog_channels['Qubit_I']
    qubit_marker  = awg_sched.digital_channels['Qubit_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']

    position = start_time-delay-num_sigma*sigma
    qubit_I.add_pulse('gaussian', position=position, amplitude=q_pulse['piAmp'], sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    
    qubit_marker.add_window(position-num_sigma*sigma, position+2*num_sigma*sigma)
    cavity_marker.add_window(start_time, start_time+window_time)
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)
    
    
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
    
    powerdat = np.zeros((len(powers), len(freqs)))
    phasedat = np.zeros((len(powers), len(freqs)))
    
    tstart = time.time()
    first_it = True
    
#    drive_powers_lin = 10**(powers/10)
#    drive_amps_lin = np.sqrt(drive_powers_lin)
        
    for powerind in range(len(powers)):
        qubitgen.Power = powers[powerind]
        time.sleep(0.2)
        
        total_samples = card.samples

#        amps_full = np.zeros((len(freqs), total_samples))
#        phases_full = np.zeros((len(freqs), total_samples))
        Is_full  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_full  = np.zeros((len(freqs), total_samples))
        
        print('Current power:{}, max:{}'.format(powers[powerind]-Qbit_Attenuation, powers[-1]-Qbit_Attenuation))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]
            
            qubitgen.Freq = freq
            qubitgen.output='On'
            time.sleep(0.2)

#            amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot = first_it, IQstorage=False)
            I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            qubitgen.output='Off'
            time.sleep(0.2)
            I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process(card, settings, 
                                                             plot=first_it, 
                                                             IQstorage = True)
            
            I_final = np.mean(I_window-I_window_b) #compute <I> in the data window
            Q_final = np.mean(Q_window-Q_window_b) #compute <Q> in the data window
#            I_final = np.mean(I_window)
#            Q_final = np.mean(Q_window)

            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, len(powers)*len(freqs))
                first_it = False
            
#            powerdat[powerind, find] = np.mean(amp)
#            phasedat[powerind, find] = np.mean(phase)
#            amps_full[find,:]   = amp_full
#            phases_full[find,:] = phase_full
            Is_full[find,:]   = I_full-I_full_b
            Qs_full[find,:] = Q_full-Q_full_b
            powerdat[powerind, find] = np.sqrt(I_final**2 + Q_final**2)
            #phasedat[powerind, find] = np.arctan2(Q_final, I_final)*180/np.pi
            I_sig, Q_sig = [np.mean(I_window), np.mean(Q_window)]
            I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)]
            theta_sig = np.arctan2(Q_sig,I_sig)*180/np.pi
            theta_back = np.arctan2(Q_back, I_back)*180/np.pi
            phasedat[powerind, find] = theta_sig-theta_back

        full_data = {}
        full_data['xaxis'] = freqs/1e9
        full_data['mags'] = powerdat[0:powerind+1]
        full_data['phases'] = phasedat[0:powerind+1]
        
#        #rescale to fractional amplitude for the plots.
#        plot_data = {}
#        plot_data['xaxis']  = freqs/1e9
#        plot_data['mags']   = np.transpose(   np.transpose(powerdat[0:powerind+1])/drive_amps_lin[0:powerind+1])
#        plot_data['phases'] = phasedat[0:powerind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mag'] = powerdat[powerind, :]
        single_data['phase'] = phasedat[powerind, :]

        yaxis = powers[0:powerind+1] - Qbit_Attenuation
        labels = ['Freq (GHz)', 'Power (dBm)']
        simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) #unscaled amplitudes
#        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=1)#scaled amplitudes
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

#        full_time = {}
#        full_time['xaxis'] = xaxis*1e6
#        full_time['mags'] = amps_full
#        full_time['phases'] = phases_full
        full_time = {}
        full_time['xaxis']  = xaxis*1e6
        full_time['Is']   = Is_full
        full_time['Qs'] = Qs_full

#        single_time = {}
#        single_time['xaxis'] = xaxis*1e6
#        single_time['mag'] = amp_full
#        single_time['phase'] = phase_full
        single_time = {}
        single_time['xaxis'] = xaxis*1e6
        single_time['I']   = I_full
        single_time['Q'] = Q_full

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Power: {}dBm'.format(powers[powerind]-Qbit_Attenuation)
        
#        simplescan_plot(full_time, 
#                        single_time, 
#                        freqs/1e9, 
#                        'Raw_time_traces\n'+filename, 
#                        time_labels, 
#                        identifier, 
#                        fig_num=2)
        simplescan_plot(full_time, single_time, freqs/1e9, 
                        'Raw_time_traces\n'+filename, 
                        time_labels, 
                        identifier, 
                        fig_num=2,
                        IQdata = True)
        plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)

        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','full_data', 'single_data', 'full_time', 'single_time'],
                             locals(), expsettings=settings, instruments=instruments)

    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
       
    cavitygen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.output = 'Off'
    
    return full_data