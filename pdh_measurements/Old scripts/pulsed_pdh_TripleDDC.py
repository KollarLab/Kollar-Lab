'''
8-25-21 AK modifying to normalize the amplitudes to the drive power.

8-25-21 AK modifying to return the raw data.

'''


import os
import time
import numpy as np 
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import simplescan_plot_TripleDDC
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process_TripleDDC
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
    
    settings['amp_list'] = [0,1,0]
    settings['phase_list'] = [1,0,0]
    settings['mod_freq'] = 20e6
    
    return settings

def pulsed_pdh_TripleDDC(instruments, settings):
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
    theta_lo = exp_settings['rotation_angle'] * np.pi / 180
    # # # #
    
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

    awg_sched.add_analog_channel(1, name='Cavity_I')
    awg_sched.add_analog_channel(2, name='Cavity_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=55e-9, HW_offset_off=0e-9)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
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
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['Qubit_enable', 'Cavity_enable'])
    
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
        settings['LPF'] = LPF
                
        xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
        
        
        # #lower sideband
        digLO_sin_lsb = np.sin(2*np.pi*(-1)*(mod_freq-exp_globals['IF'])*xaxis )
        digLO_cos_lsb = np.cos(2*np.pi*(-1)*(mod_freq-exp_globals['IF'])*xaxis )
        settings['digLO_sin_lsb'] =  digLO_sin_lsb 
        settings['digLO_cos_lsb'] = digLO_cos_lsb 
        
        # #carrier
        digLO_sin_carr = np.sin(2*np.pi*exp_globals['IF']*xaxis + theta_lo)
        digLO_cos_carr = np.cos(2*np.pi*exp_globals['IF']*xaxis + theta_lo)  
        settings['digLO_sin_carr'] = digLO_sin_carr 
        settings['digLO_cos_carr'] = digLO_cos_carr 
        
        # #upper sideband
        digLO_sin_usb = np.sin(2*np.pi*(mod_freq+exp_globals['IF'])*xaxis )
        digLO_cos_usb = np.cos(2*np.pi*(mod_freq+exp_globals['IF'])*xaxis )
        settings['digLO_sin_usb'] = digLO_sin_usb 
        settings['digLO_cos_usb'] = digLO_cos_usb 

        

    ## Start main measurement loop 
    powerdat_carr = np.zeros((len(powers), len(freqs)))
    phasedat_carr = np.zeros((len(powers), len(freqs)))
    
    powerdat_usb = np.zeros((len(powers), len(freqs)))
    phasedat_usb = np.zeros((len(powers), len(freqs)))
    
    powerdat_lsb = np.zeros((len(powers), len(freqs)))
    phasedat_lsb = np.zeros((len(powers), len(freqs)))
    
    # # # # 
    Idat = np.zeros((len(powers), len(freqs)))
    Qdat = np.zeros((len(powers), len(freqs)))
    # # # #
    
    tstart = time.time()
    first_it = True

    drive_powers_lin = 10**(powers/10)
    drive_amps_lin = np.sqrt(drive_powers_lin)

    print('asdf')
    for powerind in range(len(powers)):
        cavitygen.power = powers[powerind]
        total_samples = card.samples 
        Is_carr  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_carr  = np.zeros((len(freqs), total_samples))
        
        Is_usb  = np.zeros((len(freqs), total_samples)) 
        Qs_usb  = np.zeros((len(freqs), total_samples))
        
        Is_lsb  = np.zeros((len(freqs), total_samples)) 
        Qs_lsb  = np.zeros((len(freqs), total_samples))
        
        print('Current power:{}, max:{}'.format(powers[powerind] - CAV_Attenuation, powers[-1] - CAV_Attenuation))
    
        for find in range(0, len(freqs)):
            freq = freqs[find]

            cavitygen.freq = freq
            LO.freq = freq - exp_globals['IF']
            LO.output = 'On'       
            cavitygen.phase = 0
            LO.phase = 0

            carr_data, usb_data, lsb_data, xaxis = read_and_process_TripleDDC(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            # carr_data['I_window'], carr_data['Q_window'], carr_data['I_full'], carr_data['Q_full']                                                             
            I_final_carr = np.mean(carr_data['I_window'])
            Q_final_carr = np.mean(carr_data['Q_window'])
            
            I_final_usb = np.mean(usb_data['I_window'])
            Q_final_usb = np.mean(usb_data['Q_window'])
            
            I_final_lsb = np.mean(lsb_data['I_window'])
            Q_final_lsb = np.mean(lsb_data['Q_window'])
            
            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, len(powers)*len(freqs))
                first_it = False

            Is_carr[find,:] = carr_data['I_full']
            Qs_carr[find,:] = carr_data['Q_full']
            
            Is_usb[find,:] = usb_data['I_full']
            Qs_usb[find,:] = usb_data['Q_full']
            
            Is_lsb[find,:] = lsb_data['I_full']
            Qs_lsb[find,:] = lsb_data['Q_full']
            
            powerdat_carr[powerind, find] = np.sqrt(I_final_carr**2 + Q_final_carr**2)
            phasedat_carr[powerind, find] = np.arctan2(Q_final_carr, I_final_carr)*180/np.pi
            
            powerdat_usb[powerind, find] = np.sqrt(I_final_usb**2 + Q_final_usb**2)
            phasedat_usb[powerind, find] = np.arctan2(Q_final_usb, I_final_usb)*180/np.pi
            
            powerdat_lsb[powerind, find] = np.sqrt(I_final_lsb**2 + Q_final_lsb**2)
            phasedat_lsb[powerind, find] = np.arctan2(Q_final_lsb, I_final_lsb)*180/np.pi
            
            
            
            # # # # 
            Idat[powerind, find] = I_final_carr*I_final_usb - I_final_carr*I_final_lsb + Q_final_carr*Q_final_usb - Q_final_carr*Q_final_lsb
            Qdat[powerind, find] = -I_final_carr*Q_final_usb + Q_final_carr*I_final_usb - I_final_carr*Q_final_lsb + Q_final_carr*I_final_lsb
            # # # #
        
        powerdat = {}
        powerdat['carr'] = powerdat_carr
        powerdat['usb'] = powerdat_usb
        powerdat['lsb'] = powerdat_lsb
        
        phasedat = {}
        phasedat['carr'] = phasedat_carr
        phasedat['usb'] = phasedat_usb
        phasedat['lsb'] = phasedat_lsb
            
        full_data = {}
        full_data['xaxis']  = freqs/1e9
        full_data['mags_carr']   = powerdat_carr[0:powerind+1]
        full_data['phases_carr'] = phasedat_carr[0:powerind+1]
        full_data['mags_usb']   = powerdat_usb[0:powerind+1]
        full_data['phases_usb'] = phasedat_usb[0:powerind+1]
        full_data['mags_lsb']   = powerdat_lsb[0:powerind+1]
        full_data['phases_lsb'] = phasedat_lsb[0:powerind+1]
        
        #rescale to fractional amplitude for the plots.
        plot_data = {}
        plot_data['xaxis']  = freqs/1e9
        plot_data['mags_carr']   = np.transpose(np.transpose(powerdat_carr[0:powerind+1])/drive_amps_lin[0:powerind+1])
        plot_data['mags_usb']   = np.transpose(np.transpose(powerdat_usb[0:powerind+1])/drive_amps_lin[0:powerind+1])
        plot_data['mags_lsb']   = np.transpose(np.transpose(powerdat_lsb[0:powerind+1])/drive_amps_lin[0:powerind+1])
        plot_data['phases_carr'] = phasedat_carr[0:powerind+1]
        plot_data['phases_usb'] = phasedat_usb[0:powerind+1]
        plot_data['phases_lsb'] = phasedat_lsb[0:powerind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mags_carr']   = powerdat_carr[powerind]
        single_data['mags_usb']   = powerdat_usb[powerind]
        single_data['mags_lsb']   = powerdat_lsb[powerind]
        single_data['phases_carr'] = phasedat_carr[powerind]
        single_data['phases_usb'] = phasedat_usb[powerind]
        single_data['phases_lsb'] = phasedat_lsb[powerind]

        yaxis  = powers[0:powerind+1] - CAV_Attenuation
        labels = ['Freq (GHz)', 'Power (dBm)']
        ## CHECK AND REWRITE A CUSTOMIZED SIMPLESCAN PLOT FOR THE TRIPLE DDC DATA
        simplescan_plot_TripleDDC(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=11, IQdata = False) #normalized to drive level
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

        full_time = {}
        full_time['xaxis']  = xaxis*1e6
        full_time['I_carr']   = Is_carr 
        full_time['Q_carr'] = Qs_carr
        full_time['I_usb']   = Is_usb 
        full_time['Q_usb'] = Qs_usb
        full_time['I_lsb']   = Is_lsb 
        full_time['Q_lsb'] = Qs_lsb

        single_time = {}
        single_time['xaxis'] = xaxis*1e6
        single_time['I_carr']   = carr_data['I_full']
        single_time['Q_carr'] = carr_data['Q_full']
        single_time['I_usb']   = usb_data['I_full']
        single_time['Q_usb'] = usb_data['Q_full']
        single_time['I_lsb']   = lsb_data['I_full']
        single_time['Q_lsb'] = lsb_data['Q_full']

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = 'Power: {}dBm'.format(powers[powerind]-CAV_Attenuation)
        simplescan_plot_TripleDDC(full_time, single_time, freqs/1e9, 
                        'Raw_time_traces\n'+filename, 
                        time_labels, 
                        identifier, 
                        fig_num=12,
                        IQdata = True)
        plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces.png'), dpi = 150)

        userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','full_data', 'single_data', 'full_time', 'single_time'],
                             locals(), expsettings=settings, instruments=instruments, saveHWsettings=first_it)
        
        # # # # 
        fig = plt.figure(122)
        plt.clf()
        ax = plt.subplot(1,2,1)
        plt.plot(freqs/1e9, Idat[powerind], color = 'deepskyblue', label='I_full')
        plt.xlabel('Freq (GHz)')
        ax.legend()
        
        ax = plt.subplot(1,2,2)
        plt.plot(freqs/1e9, Qdat[powerind], color = 'deepskyblue', label='Q_full')
        plt.xlabel('Freq (GHz)')
        ax.legend()
        
        fig.suptitle(f'Final I and Q (TripleDDC)| Rot: {round(theta_lo *180/np.pi,2)} Deg')
        fig.tight_layout()
        plt.show()
        # # #
        
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
    # # # # # 
    # fig = plt.figure(121)
    # plt.clf()
    # ax = plt.subplot(1,2,1)
    # plt.plot(freqs/1e9, Idat[0:powerind+1], color = 'deepskyblue', label='I_full')
    # plt.xlabel('Freq (GHz)')
    # ax.legend()
    
    # ax = plt.subplot(1,2,2)
    # plt.plot(freqs/1e9, Qdat[0:powerind+1], color = 'deepskyblue', label='Q_full')
    # plt.xlabel('Freq (GHz)')
    # ax.legend()
    
    # fig.suptitle('Final I and Q')
    # fig.tight_layout()
    # plt.show()
    # # # #
    
    cavitygen.output = 'Off'
    LO.output = 'Off'
    
    return full_data