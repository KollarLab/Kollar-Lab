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
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process_2Readouts
from pdh_measurements.scheduler_pdh import scheduler_pdh

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'PulsedTrans'
    
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

def pulsed_2Readouts_trans(instruments, settings):
    ##Instruments used
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
    start_freq  = exp_settings['start_freq']
    stop_freq   = exp_settings['stop_freq']
    freq_points = exp_settings['freq_points']
    freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
    
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    start_power  = exp_settings['start_power'] + CAV_Attenuation
    stop_power   = exp_settings['stop_power'] + CAV_Attenuation
    power_points = exp_settings['power_points']
    powers = np.round(np.linspace(start_power,stop_power,power_points),2)
           
    # # # #
    amp_list = exp_settings['amp_list']
    phase_list = exp_settings['phase_list']
    mod_freq = exp_settings['mod_freq']
    # # # #    
    m_pulse = exp_globals['measurement_pulse']
    q_pulse = exp_globals['qubit_pulse']
    
    ## Generator settings
    cavitygen.freq   = freqs[0]
    cavitygen.power  = powers[0]
    
    cavitygen.enable_pulse()
    cavitygen.enable_IQ()
    cavitygen.output = 'On'

    LO.power  = 12
    LO.freq   = freqs[0] - exp_globals['IF']
    LO.output = 'On'
    
    ##Card settings
    configure_card(card, settings)

    ##HDAWG settings  
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    
    awg_sched = scheduler_pdh(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='Qubit_I')
    awg_sched.add_analog_channel(2, name='Qubit_Q')
    awg_sched.add_analog_channel(3, name='Cavity_I')
    awg_sched.add_analog_channel(4, name='Cavity_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(3, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
      
    cavity_I = awg_sched.analog_channels['Cavity_I']
    cavity_Q = awg_sched.analog_channels['Cavity_Q'] 
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    position = start_time
    cavity_I.add_pulse('pdh_I', position=position,
                           amp_list = amp_list, phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = window_time)

    cavity_Q.add_pulse('pdh_Q', position=position,
                           amp_list = amp_list, phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = window_time)
    cavity_marker.add_window(start_time, start_time+window_time)
       
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['blank1', 'blank2'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
    hdawg.AWGs[0].run_loop()
    # hdawg.AWGs[1].run_loop() #seeing if this is necessary
    time.sleep(0.1)
    
    ##create the digital down conversion filter if needed.
    if exp_globals['IF'] != 0:
        #create Chebychev type II digital filter
        filter_N = exp_globals['ddc_config']['order']
        filter_rs = exp_globals['ddc_config']['stop_atten']
        filter_cutoff = np.abs(exp_globals['ddc_config']['cutoff'])
        LPF = signal.cheby2(filter_N, filter_rs, filter_cutoff, btype='low', analog=False, output='sos', fs=card.sampleRate)
        settings['LPF'] = LPF
        
        xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
        
        # carrier
        digLO_sin_carr = np.sin(2*np.pi*exp_globals['IF']*xaxis)
        digLO_cos_carr = np.cos(2*np.pi*exp_globals['IF']*xaxis)
        settings['digLO_sin_carr'] = digLO_sin_carr 
        settings['digLO_cos_carr'] = digLO_cos_carr
        
        # lower sideband
        digLO_sin_sb = np.sin(2*np.pi*(-1)*(mod_freq-exp_globals['IF'])*xaxis )
        digLO_cos_sb = np.cos(2*np.pi*(-1)*(mod_freq-exp_globals['IF'])*xaxis )
        settings['digLO_sin_sb'] =  digLO_sin_sb 
        settings['digLO_cos_sb'] = digLO_cos_sb 
        

    ## Start main measurement loop 
    I_carr = np.zeros((len(powers), len(freqs)))
    Q_carr = np.zeros((len(powers), len(freqs)))
    powerdat_carr = np.zeros((len(powers), len(freqs)))
    phasedat_carr = np.zeros((len(powers), len(freqs)))
    
    I_sb = np.zeros((len(powers), len(freqs)))
    Q_sb = np.zeros((len(powers), len(freqs)))
    powerdat_sb = np.zeros((len(powers), len(freqs)))
    phasedat_sb = np.zeros((len(powers), len(freqs)))
    
    phasedat_diff = np.zeros((len(powers), len(freqs)))
    
    tstart = time.time()
    first_it = True

    drive_powers_lin = 10**(powers/10)
    drive_amps_lin = np.sqrt(drive_powers_lin)
    
    # ######hack test
    # cavitygen.disable_IQ() ###

    print('asdf')
    for powerind in range(len(powers)):
        cavitygen.power = powers[powerind]
        total_samples = card.samples 
        Is_carr  = np.zeros((len(freqs), total_samples))
        Qs_carr  = np.zeros((len(freqs), total_samples))
        
        Is_sb  = np.zeros((len(freqs), total_samples))
        Qs_sb  = np.zeros((len(freqs), total_samples))
        
        print('Current power:{}, max:{}'.format(powers[powerind] - CAV_Attenuation, powers[-1] - CAV_Attenuation))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]

            cavitygen.freq = freq
            LO.freq = freq - exp_globals['IF']
            LO.output = 'On'  
            cavitygen.phase = 0
            LO.phase = 0
    
            carr_data, sb_data, xaxis = read_and_process_2Readouts(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            
            I_final_carr = np.mean(carr_data['I_window']) #compute <I> in the data window
            Q_final_carr = np.mean(carr_data['Q_window']) #compute <Q> in the data window
            
            I_final_sb = np.mean(sb_data['I_window']) #compute <I> in the data window
            Q_final_sb = np.mean(sb_data['Q_window']) #compute <Q> in the data window
            
            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, len(powers)*len(freqs))
                first_it = False
                
            Is_carr[find,:] = carr_data['I_full']
            Qs_carr[find,:] = carr_data['Q_full']
            
            Is_sb[find,:] = sb_data['I_full']
            Qs_sb[find,:] = sb_data['Q_full']
            
            I_carr[powerind, find] = I_final_carr
            Q_carr[powerind, find] = Q_final_carr
            powerdat_carr[powerind, find] = np.sqrt(I_final_carr**2 + Q_final_carr**2)
            phasedat_carr[powerind, find] = np.arctan2(Q_final_carr, I_final_carr)*180/np.pi
            
            I_sb[powerind, find] = I_final_sb
            Q_sb[powerind, find] = Q_final_sb
            powerdat_sb[powerind, find] = np.sqrt(I_final_sb**2 + Q_final_sb**2)
            phasedat_sb[powerind, find] = np.arctan2(Q_final_sb, I_final_sb)*180/np.pi
            
            phasedat_diff[powerind, find] = phasedat_carr[powerind, find] - phasedat_sb[powerind, find]
        
        # save data
        IQdat = {}
        IQdat['I_carr'] = I_carr
        IQdat['Q_carr'] = Q_carr
        IQdat['I_sb'] = I_sb
        IQdat['Q_sb'] = Q_sb
                 
        full_data = {}
        full_data['xaxis']  = freqs/1e9
        full_data['carr_mags']   = powerdat_carr
        full_data['carr_phases'] = phasedat_carr
        full_data['sb_mags']   = powerdat_sb
        full_data['sb_phases'] = phasedat_sb
        full_data['phase_diff'] = phasedat_diff
        
        # Plot simple mag scan
        fig = plt.figure(1212, figsize=(13,8))
        plt.clf()
        ax = plt.subplot(1,2,1)
        plt.plot(full_data['xaxis'], full_data['carr_mags'][0])
        plt.xlabel('Freq (GHz)')
        plt.title('Carr mag')
        
        ax = plt.subplot(1,2,2)
        plt.plot(full_data['xaxis'], full_data['sb_mags'][0])
        plt.xlabel('Freq (GHz)')
        plt.title('sb mag')
        
        plt.suptitle('Filename: {}'.format(filename))
        plt.tight_layout()
        plt.show()
        plt.savefig(os.path.join(saveDir, filename+'_fullMagPlot.png'), dpi = 150)

        full_time = {}
        full_time['xaxis']  = xaxis*1e6
        full_time['Is_carr']   = Is_carr
        full_time['Qs_carr'] = Qs_carr
        full_time['Is_sb']   = Is_sb
        full_time['Qs_sb'] = Qs_sb

        identifier = 'Power: {}dBm'.format(powers[powerind]-CAV_Attenuation)
        new_filename = 'Raw_time_traces2\n'+filename
        
        # Plot raw time traces
        fig = plt.figure(1213, figsize=(13,8))
        plt.clf()
        ax = plt.subplot(1,2,1)
        plt.plot(full_time['xaxis'], carr_data['I_full'], label = 'I_full')
        plt.plot(full_time['xaxis'], carr_data['Q_full'], label = 'Q_full')
        plt.legend()
        plt.xlabel('Time (us)')
        plt.title('single carr trace')
        
        ax = plt.subplot(1,2,2)
        plt.plot(full_time['xaxis'], sb_data['I_full'], label = 'I_full')
        plt.plot(full_time['xaxis'], sb_data['Q_full'], label = 'Q_full')
        plt.legend()
        plt.xlabel('Time (us)')
        plt.title('single sb trace')
        
        plt.suptitle('Filename: {}'.format(new_filename))
        plt.tight_layout()
        plt.show()
        plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces2.png'), dpi = 150)

        # for normal plots similar to previous codes
        full_data_plot = {}
        full_data_plot['xaxis'] = freqs/1e9
        full_data_plot['mags'] = powerdat_carr
        full_data_plot['phases'] = phasedat_carr
        
        single_data_plot = {}
        single_data_plot['xaxis'] = freqs/1e9
        single_data_plot['mag'] = powerdat_carr[powerind, :]
        single_data_plot['phase'] = phasedat_carr[powerind, :]
        
        yaxis = powers - CAV_Attenuation
        labels = ['Freq (GHz)', 'Power (dBm)']
        simplescan_plot(full_data_plot, single_data_plot, yaxis, filename, labels, identifier='', fig_num=1)
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

        full_time_plot = {}
        full_time_plot['xaxis']  = xaxis*1e6
        full_time_plot['Is'] = Is_carr
        full_time_plot['Qs'] = Qs_carr
    
    
        single_time_plot = {}
        single_time_plot['xaxis'] = xaxis*1e6
        single_time_plot['I'] = carr_data['I_full']
        single_time_plot['Q'] = carr_data['Q_full']
        
        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Power: {}dBm'.format(powers[powerind]-CAV_Attenuation)
    
        simplescan_plot(full_time_plot, single_time_plot, freqs/1e9, 
                        'Raw_time_traces\n'+filename, 
                        time_labels, 
                        identifier, 
                        fig_num=2,
                        IQdata = True)
        plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)
        

        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'IQdat','full_data', 'full_time'],
                             locals(), expsettings=settings, instruments=instruments, saveHWsettings=first_it)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
    
    LO.output = 0
    cavitygen.output = 0
    
    return full_data