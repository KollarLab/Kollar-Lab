'''
8-25-21 AK modifying to normalize the amplitudes to the drive power.

8-25-21 AK modifying to return the raw data.

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
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'PulsedTransStability'
    
    #Sweep parameters
    settings['start_freq']   = 4.15*1e9  
    settings['stop_freq']    = 4.25*1e9 
    settings['freq_points']  = 50

    settings['start_power']  = -20
    settings['stop_power']   = 10
    settings['power_points'] = 31

    #Card settings
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 5e3
    
    return settings

def pulsed_trans_stability(instruments, settings):
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    ##Sweep settings
    cav_start_freq  = exp_settings['cav_start_freq']
    cav_stop_freq   = exp_settings['cav_stop_freq']
    cav_freq_points = exp_settings['cav_freq_points']
    freqs  = np.round(np.linspace(cav_start_freq,cav_stop_freq,cav_freq_points),-3)
    
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = exp_settings['CAV_power']  

    # Qubit settings
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power  = exp_settings['Qbit_power'] + Qbit_Attenuation
    Qbit_freq = exp_settings['Qbit_freq']
    avgs = exp_settings['averages']
    qubitgen.freq = Qbit_freq
    qubitgen.power = Qbit_power
    if exp_settings['Quasi_CW']:
        qubitgen.disable_pulse()
        qubitgen.disable_IQ()
    else:
        qubitgen.enable_pulse()
        qubitgen.enable_IQ()
    
    stab_points = exp_settings['stability_points']
    
    ## Generator settings
    cavitygen.freq   = freqs[0]
    cavitygen.power  = CAV_power
    
    cavitygen.enable_pulse()
    cavitygen.disable_IQ()
    cavitygen.output = 'On'

    LO.power  = 12
    LO.freq   = freqs[0] - exp_globals['IF']
    LO.output = 'On'
    
    ##Card settings
    configure_card(card, settings)

    ##HDAWG settings
    ## Sequencer program
    #####################################################################################################################################
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

    position = start_time-delay-num_sigma*sigma-q_pulse['hold_time']
    qubit_I.add_pulse('gaussian_square', position=position, 
                              amplitude=q_pulse['piAmp'], length = q_pulse['hold_time'], 
                              ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    
    qubit_marker.add_window(position-160e-9, position+2*160e-9+q_pulse['hold_time'])
    cavity_marker.add_window(start_time, start_time+window_time)
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)

    #############################################################################################################################
    
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

    ## Start main measurement loop 
    powerdat = np.zeros((len(stab_points), len(freqs)))
    phasedat = np.zeros((len(stab_points), len(freqs)))
    
    tstart = time.time()
    first_it = True

    drive_stab_points_lin = 10**(stab_points/10)
    drive_amps_lin = np.sqrt(drive_stab_points_lin)

    print('asdf')
    for powerind in range(len(stab_points)):
#        cavitygen.Power = stab_points[powerind]
        cavitygen.power = stab_points[powerind]
        
        # time.sleep(0.2)

        total_samples = card.samples 
#        amps   = np.zeros((len(freqs), total_samples))
#        phases = np.zeros((len(freqs), total_samples))
        Is  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs  = np.zeros((len(freqs), total_samples))
        
        print('Current power:{}, max:{}'.format(stab_points[powerind] - CAV_Attenuation, stab_points[-1] - CAV_Attenuation))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]

            cavitygen.freq = freq

            #LO.freq   = freq-exp_globals['IF']
            LO.freq = freq - exp_globals['IF']
            LO.output = 'On'
            
            cavitygen.phase = 0
            LO.phase = 0
            # time.sleep(0.2)

#            amp, phase, amp_full, phase_full, xaxis = read_and_process(card, settings, plot=first_it, IQstorage = False)
            I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            
            I_final = np.mean(I_window) #compute <I> in the data window
            Q_final = np.mean(Q_window) #compute <Q> in the data window
            
            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, len(stab_points)*len(freqs))
                first_it = False
                
            
#            amps[find,:]   = amp_full
#            phases[find,:] = phase_full
#            powerdat[powerind, find] = np.mean(amp)
#            phasedat[powerind, find] = np.mean(phase)   ####!!!! this does not work with heterodyne. Because the phase is wrapping.
            Is[find,:]   = I_full
            Qs[find,:] = Q_full
            powerdat[powerind, find] = np.sqrt(I_final**2 + Q_final**2)
            phasedat[powerind, find] = np.arctan2(Q_final, I_final)*180/np.pi
            
        full_data = {}
        full_data['xaxis']  = freqs/1e9
        full_data['mags']   = powerdat[0:powerind+1]
        full_data['phases'] = phasedat[0:powerind+1]
        
        #rescale to fractional amplitude for the plots.
        plot_data = {}
        plot_data['xaxis']  = freqs/1e9
        plot_data['mags']   = np.transpose(   np.transpose(powerdat[0:powerind+1])/drive_amps_lin[0:powerind+1])
        plot_data['phases'] = phasedat[0:powerind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mag']   = powerdat[powerind]
        single_data['phase'] = phasedat[powerind]

        yaxis  = stab_points[0:powerind+1] - CAV_Attenuation
        labels = ['Freq (GHz)', 'Power (dBm)']
#        simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) #unnormalized plot
#        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=1) #normalized to drive level
        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=1, IQdata = False) #normalized to drive level
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

#        full_time = {}
#        full_time['xaxis']  = xaxis*1e6
#        full_time['mags']   = amps
#        full_time['phases'] = phases
        full_time = {}
        full_time['xaxis']  = xaxis*1e6
        full_time['Is']   = Is
        full_time['Qs'] = Qs

#        single_time = {}
#        single_time['xaxis'] = xaxis*1e6
#        single_time['mag']   = amp_full
#        single_time['phase'] = phase_full
        single_time = {}
        single_time['xaxis'] = xaxis*1e6
        single_time['I']   = I_full
        single_time['Q'] = Q_full

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Power: {}dBm'.format(stab_points[powerind]-CAV_Attenuation)
#        simplescan_plot(full_time, single_time, freqs/1e9, 
#                        'Raw_time_traces\n'+filename, 
#                        time_labels, 
#                        identifier, 
#                        fig_num=2)
#        simplescan_plot(full_time, single_time, freqs/1e9, 
#                        'Raw_time_traces\n'+filename, 
#                        time_labels, 
#                        identifier, 
#                        fig_num=2,
#                        IQdata = True)
        simplescan_plot(full_time, single_time, freqs/1e9, 
                        'Raw_time_traces\n'+filename, 
                        time_labels, 
                        identifier, 
                        fig_num=2,
                        IQdata = True)
        plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces.png'), dpi = 150)

        userfuncs.SaveFull(saveDir, filename, ['stab_points','freqs', 'powerdat', 'phasedat','xaxis','full_data', 'single_data', 'full_time', 'single_time'],
                             locals(), expsettings=settings, instruments=instruments, saveHWsettings=first_it)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
    
    return full_data