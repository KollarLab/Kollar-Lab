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
from utility.measurement_helpers import get_amp_comps, total_power

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'PDH_TripleDDC'
    
    #Sweep parameters
    settings['start_freq']   = 4.15*1e9  
    settings['stop_freq']    = 4.25*1e9 
    settings['freq_points']  = 50

    settings['upper_sideband_power']  = -90
    settings['carrier_power'] = -45
    settings['lower_sideband_power'] = -90

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
    qubitgen = instruments['qubitgen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    ## Readout Sweep settings
    start_freq  = exp_settings['start_freq']
    stop_freq   = exp_settings['stop_freq']
    freq_points = exp_settings['freq_points']
    freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
    cavitygen.freq = freqs[0]
    cavitygen.phase = 0
    cavitygen.enable_pulse()
    cavitygen.enable_IQ()

    LO.power = 12
    LO.freq = freqs[0] - exp_globals['IF']    
    LO.phase = 0

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

    # Set up card settings
    configure_card(card, settings)

    # carrier and sideband 
    phase_list = exp_settings['phase_list']
    mod_freq = exp_settings['mod_freq']
    theta_lo = exp_settings['rotation_angle'] * np.pi / 180
    power_points = exp_settings['power_points']
    CAV_Attenuation = exp_globals['CAV_Attenuation']

    usb_power_start = exp_settings['usb_power_start']
    usb_power_stop = exp_settings['usb_power_stop']
    usb_powers = np.round(np.linspace(usb_power_start, usb_power_stop,power_points), 2)

    lsb_power_start = exp_settings['lsb_power_start']
    lsb_power_stop = exp_settings['lsb_power_stop']
    lsb_powers = np.round(np.linspace(lsb_power_start, lsb_power_stop,power_points), 2)

    carrier_power_start = exp_settings['carrier_power_start']
    carrier_power_stop = exp_settings['carrier_power_stop']
    carrier_powers = np.round(np.linspace(carrier_power_start, carrier_power_stop,power_points), 2)

    ##HDAWG settings  
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    q_pulse = exp_globals['qubit_pulse']
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
    
    # data saving arrays
    powerdat_carr_g = np.zeros((power_points, len(freqs)))
    phasedat_carr_g = np.zeros((power_points, len(freqs)))
    powerdat_usb_g = np.zeros((power_points, len(freqs)))
    phasedat_usb_g = np.zeros((power_points, len(freqs)))
    powerdat_lsb_g = np.zeros((power_points, len(freqs)))
    phasedat_lsb_g = np.zeros((power_points, len(freqs)))

    powerdat_carr_e = np.zeros((power_points, len(freqs)))
    phasedat_carr_e = np.zeros((power_points, len(freqs)))
    powerdat_usb_e = np.zeros((power_points, len(freqs)))
    phasedat_usb_e = np.zeros((power_points, len(freqs)))
    powerdat_lsb_e = np.zeros((power_points, len(freqs)))
    phasedat_lsb_e = np.zeros((power_points, len(freqs)))
    
    # # # # 
    Idat_g = np.zeros((power_points, len(freqs)))
    Qdat_g = np.zeros((power_points, len(freqs)))
    Idat_e = np.zeros((power_points, len(freqs)))
    Qdat_e = np.zeros((power_points, len(freqs)))
    # # # #
    drive_powers_lin = 10**(power_points/10)
    drive_amps_lin_carr = np.sqrt(10**(carrier_powers/10))
    drive_amps_lin_usb = np.sqrt(10**(usb_powers/10))
    drive_amps_lin_lsb = np.sqrt(10**(lsb_powers/10))

    # drive_amps_lin = np.sqrt(drive_powers_lin)
    # For each power point for the 3 readout component, load the AWG for cavity readout
    for pind in range(power_points):
        usb_power = usb_powers[pind]
        lsb_power = lsb_powers[pind]
        carrier_power = carrier_powers[pind]
        power_list = [usb_power, carrier_power, lsb_power]

        # calculate the total power of the 3 components | and increase by 10 dBm
        total_comp_power = total_power(power_list) + 10
        # set gen power
        if total_comp_power > 0:
            gen_power = total_comp_power
        else:
            gen_power = 0
        # Use gen power to get amps of the components
        amp_list = get_amp_comps(power_list, gen_power)
        
        #update gen power to include cavity attenuation
        gen_power = gen_power + CAV_Attenuation

        if gen_power > 10:
            raise ValueError(f'SGS power of {gen_power} too high')
        
        # Load AWG for cavity and qubit
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
        
        qubit_I =  awg_sched.analog_channels['Qubit_I']
        qubit_marker  = awg_sched.digital_channels['Qubit_enable']
        delay = q_pulse['delay']
        sigma = q_pulse['sigma']
        num_sigma = q_pulse['num_sigma']
        
        position = start_time-delay-num_sigma*sigma-q_pulse['hold_time']
        qubit_I.add_pulse('gaussian_square', position=position, 
                                amplitude=q_pulse['piAmp'], length = q_pulse['hold_time'], 
                                ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
        
        qubit_marker.add_window(position-160e-9, position+2*160e-9+q_pulse['hold_time'])
        
        awg_sched.plot_waveforms()
        
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['blank1', 'blank2'])
        
        loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
        hdawg.AWGs[0].load_program(loadprog)
        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.2)
        
    
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

        
        tstart = time.time()
        first_it = True
    
        cavitygen.power = gen_power
        total_samples = card.samples 
        # For ground state 
        Is_carr_g  = np.zeros((len(freqs), total_samples)) 
        Qs_carr_g  = np.zeros((len(freqs), total_samples))
        Is_usb_g  = np.zeros((len(freqs), total_samples)) 
        Qs_usb_g  = np.zeros((len(freqs), total_samples))
        Is_lsb_g  = np.zeros((len(freqs), total_samples)) 
        Qs_lsb_g  = np.zeros((len(freqs), total_samples))
        # For excited state
        Is_carr_e  = np.zeros((len(freqs), total_samples)) 
        Qs_carr_e  = np.zeros((len(freqs), total_samples))
        Is_usb_e  = np.zeros((len(freqs), total_samples)) 
        Qs_usb_e  = np.zeros((len(freqs), total_samples))
        Is_lsb_e  = np.zeros((len(freqs), total_samples)) 
        Qs_lsb_e  = np.zeros((len(freqs), total_samples))
        
        print(f'Current Powers| carrier:{carrier_power}, sideband:{usb_power}, SGS:{gen_power-CAV_Attenuation}')
    
        for find in range(0, len(freqs)):

            if first_it:
                tstart = time.time()
            freq = freqs[find]
            cavitygen.freq = freq
            cavitygen.phase = 0
            cavitygen.output = 'On'
            LO.freq = freq - exp_globals['IF']
            LO.phase = 0
            LO.output = 'On'       
            
            # get ground state data
            qubitgen.output = 'Off'
            carr_data_g, usb_data_g, lsb_data_g, xaxis_g = read_and_process_TripleDDC(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            if exp_settings['subtract_background']:
                #Acquire background trace
                qubitgen.output='Off'
                carr_data_b, usb_data_b, lsb_data_b, xaxis_b = read_and_process_TripleDDC(card, settings, 
                                                                            plot=first_it, 
                                                                            IQstorage = True)
            else:
                carr_data_b = {key: np.zeros_like(value) for key, value in carr_data_g.items()}
                usb_data_b = {key: np.zeros_like(value) for key, value in usb_data_g.items()}
                lsb_data_b = {key: np.zeros_like(value) for key, value in lsb_data_g.items()}
            # get excited state data
            qubitgen.output = 'On'
            qubitgen.phase = 0
            time.sleep(0.1)
            carr_data_e, usb_data_e, lsb_data_e, xaxis_e = read_and_process_TripleDDC(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, power_points*len(freqs))
                first_it = False

            I_final_carr_g = np.mean(carr_data_g['I_window'])
            Q_final_carr_g = np.mean(carr_data_g['Q_window'])
            
            I_final_usb_g = np.mean(usb_data_g['I_window'])
            Q_final_usb_g = np.mean(usb_data_g['Q_window'])
            
            I_final_lsb_g = np.mean(lsb_data_g['I_window'])
            Q_final_lsb_g = np.mean(lsb_data_g['Q_window'])

            I_final_carr_e = np.mean(carr_data_e['I_window'])
            Q_final_carr_e = np.mean(carr_data_e['Q_window'])
            
            I_final_usb_e = np.mean(usb_data_e['I_window'])
            Q_final_usb_e = np.mean(usb_data_e['Q_window'])
            
            I_final_lsb_e = np.mean(lsb_data_e['I_window'])
            Q_final_lsb_e = np.mean(lsb_data_e['Q_window'])
            
            

            Is_carr_g[find,:] = carr_data_g['I_full']
            Qs_carr_g[find,:] = carr_data_g['Q_full']
            Is_usb_g[find,:] = usb_data_g['I_full']
            Qs_usb_g[find,:] = usb_data_g['Q_full']
            Is_lsb_g[find,:] = lsb_data_g['I_full']
            Qs_lsb_g[find,:] = lsb_data_g['Q_full']

            Is_carr_e[find,:] = carr_data_e['I_full']
            Qs_carr_e[find,:] = carr_data_e['Q_full']
            Is_usb_e[find,:] = usb_data_e['I_full']
            Qs_usb_e[find,:] = usb_data_e['Q_full']
            Is_lsb_e[find,:] = lsb_data_e['I_full']
            Qs_lsb_e[find,:] = lsb_data_e['Q_full']
            
            powerdat_carr_g[pind, find] = np.sqrt(I_final_carr_g**2 + Q_final_carr_g**2)
            phasedat_carr_g[pind, find] = np.arctan2(Q_final_carr_g, I_final_carr_g)*180/np.pi
            powerdat_usb_g[pind, find] = np.sqrt(I_final_usb_g**2 + Q_final_usb_g**2)
            phasedat_usb_g[pind, find] = np.arctan2(Q_final_usb_g, I_final_usb_g)*180/np.pi
            powerdat_lsb_g[pind, find] = np.sqrt(I_final_lsb_g**2 + Q_final_lsb_g**2)
            phasedat_lsb_g[pind, find] = np.arctan2(Q_final_lsb_g, I_final_lsb_g)*180/np.pi

            powerdat_carr_e[pind, find] = np.sqrt(I_final_carr_e**2 + Q_final_carr_e**2)
            phasedat_carr_e[pind, find] = np.arctan2(Q_final_carr_e, I_final_carr_e)*180/np.pi
            powerdat_usb_e[pind, find] = np.sqrt(I_final_usb_e**2 + Q_final_usb_e**2)
            phasedat_usb_e[pind, find] = np.arctan2(Q_final_usb_e, I_final_usb_e)*180/np.pi
            powerdat_lsb_e[pind, find] = np.sqrt(I_final_lsb_e**2 + Q_final_lsb_e**2)
            phasedat_lsb_e[pind, find] = np.arctan2(Q_final_lsb_e, I_final_lsb_e)*180/np.pi
            
            
            
            # # # # 
            Idat_g[pind, find] = I_final_carr_g*I_final_usb_g - I_final_carr_g*I_final_lsb_g + Q_final_carr_g*Q_final_usb_g - Q_final_carr_g*Q_final_lsb_g
            Qdat_g[pind, find] = -I_final_carr_g*Q_final_usb_g + Q_final_carr_g*I_final_usb_g - I_final_carr_g*Q_final_lsb_g + Q_final_carr_g*I_final_lsb_g

            Idat_e[pind, find] = I_final_carr_e*I_final_usb_e - I_final_carr_e*I_final_lsb_e + Q_final_carr_e*Q_final_usb_e - Q_final_carr_e*Q_final_lsb_e
            Qdat_e[pind, find] = -I_final_carr_e*Q_final_usb_e + Q_final_carr_e*I_final_usb_e - I_final_carr_e*Q_final_lsb_e + Q_final_carr_e*I_final_lsb_e
            # # # #
        IQdat = {}
        IQdat['Idat_g'] = Idat_g
        IQdat['Qdat_g'] = Qdat_g
        IQdat['Idat_e'] = Idat_e
        IQdat['Qdat_e'] = Qdat_e

        powerdat = {}
        powerdat['carr_g'] = powerdat_carr_g
        powerdat['usb_g'] = powerdat_usb_g
        powerdat['lsb_g'] = powerdat_lsb_g
        powerdat['carr_e'] = powerdat_carr_e
        powerdat['usb_e'] = powerdat_usb_e
        powerdat['lsb_e'] = powerdat_lsb_e
        
        phasedat = {}
        phasedat['carr_g'] = phasedat_carr_g
        phasedat['usb_g'] = phasedat_usb_g
        phasedat['lsb_g'] = phasedat_lsb_g
        phasedat['carr_e'] = phasedat_carr_e
        phasedat['usb_e'] = phasedat_usb_e
        phasedat['lsb_e'] = phasedat_lsb_e
            
        full_data = {}
        full_data['xaxis']  = freqs/1e9
        full_data['phasedat'] = phasedat
        full_data['powerdat'] = powerdat
        
        #rescale to fractional amplitude for the plots.
        plot_data = {}
        plot_data['xaxis']  = freqs/1e9
        plot_data['mags_carr']   = np.transpose(np.transpose(powerdat_carr_g[0:pind+1])/drive_amps_lin_carr[0:pind+1])
        plot_data['mags_usb']   = np.transpose(np.transpose(powerdat_usb_g[0:pind+1])/drive_amps_lin_usb[0:pind+1])
        plot_data['mags_lsb']   = np.transpose(np.transpose(powerdat_lsb_g[0:pind+1])/drive_amps_lin_lsb[0:pind+1])
        plot_data['phases_carr'] = phasedat_carr_g[0:pind+1]
        plot_data['phases_usb'] = phasedat_usb_g[0:pind+1]
        plot_data['phases_lsb'] = phasedat_lsb_g[0:pind+1]

        single_data = {}
        single_data['xaxis'] = freqs/1e9
        single_data['mags_carr']   = powerdat_carr_g[pind]
        single_data['mags_usb']   = powerdat_usb_g[pind]
        single_data['mags_lsb']   = powerdat_lsb_g[pind]
        single_data['phases_carr'] = phasedat_carr_g[pind]
        single_data['phases_usb'] = phasedat_usb_g[pind]
        single_data['phases_lsb'] = phasedat_lsb_g[pind]

        # yaxis  = power_points[0:pind+1] - CAV_Attenuation
        yaxis  = carrier_powers[0:pind+1] - CAV_Attenuation
        labels = ['Freq (GHz)', 'Power (dBm)']
        
        ## CHECK AND REWRITE A CUSTOMIZED SIMPLESCAN PLOT FOR THE TRIPLE DDC DATA
        # simplescan_plot_TripleDDC(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=11, IQdata = False) #normalized to drive level
        # plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

        full_time_g = {}
        full_time_g['xaxis']  = xaxis_g*1e6
        full_time_g['I_carr']   = Is_carr_g 
        full_time_g['Q_carr'] = Qs_carr_g
        full_time_g['I_usb']   = Is_usb_g 
        full_time_g['Q_usb'] = Qs_usb_g
        full_time_g['I_lsb']   = Is_lsb_g 
        full_time_g['Q_lsb'] = Qs_lsb_g

        single_time_g = {}
        single_time_g['xaxis'] = xaxis_g*1e6
        single_time_g['I_carr']   = carr_data_g['I_full']
        single_time_g['Q_carr'] = carr_data_g['Q_full']
        single_time_g['I_usb']   = usb_data_g['I_full']
        single_time_g['Q_usb'] = usb_data_g['Q_full']
        single_time_g['I_lsb']   = lsb_data_g['I_full']
        single_time_g['Q_lsb'] = lsb_data_g['Q_full']

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = f'Carrier:{carrier_power}dBm_Sideband{lsb_power}dBm'
        # simplescan_plot_TripleDDC(full_time_g, single_time_g, freqs/1e9, 
        #                 'Raw_time_traces\n'+filename, 
        #                 time_labels, 
        #                 identifier, 
        #                 fig_num=12,
        #                 IQdata = True)
        # plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces_g.png'), dpi = 150)

        full_time_e = {}
        full_time_e['xaxis']  = xaxis_e*1e6
        full_time_e['I_carr']   = Is_carr_e 
        full_time_e['Q_carr'] = Qs_carr_e
        full_time_e['I_usb']   = Is_usb_e 
        full_time_e['Q_usb'] = Qs_usb_e
        full_time_e['I_lsb']   = Is_lsb_e 
        full_time_e['Q_lsb'] = Qs_lsb_e

        single_time_e = {}
        single_time_e['xaxis'] = xaxis_e*1e6
        single_time_e['I_carr']   = carr_data_e['I_full']
        single_time_e['Q_carr'] = carr_data_e['Q_full']
        single_time_e['I_usb']   = usb_data_e['I_full']
        single_time_e['Q_usb'] = usb_data_e['Q_full']
        single_time_e['I_lsb']   = lsb_data_e['I_full']
        single_time_e['Q_lsb'] = lsb_data_e['Q_full']

        time_labels = ['Time (us)', 'Freq (GHz)']
        identifier = f'Carrier:{carrier_power}dBm_Sideband{lsb_power}dBm'
        # simplescan_plot_TripleDDC(full_time_e, single_time_e, freqs/1e9, 
        #                 'Raw_time_traces\n'+filename, 
        #                 time_labels, 
        #                 identifier, 
        #                 fig_num=121,
        #                 IQdata = True)
        # plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces_e.png'), dpi = 150)

        userfuncs.SaveFull(saveDir, filename, ['power_points','freqs','xaxis','full_data', 'IQdat', 'full_time_g', 'single_time_g', 'full_time_e', 'single_time_e'],
                             locals(), expsettings=settings, instruments=instruments, saveHWsettings=first_it)
        
        # # # # 
        fig = plt.figure(122)
        plt.clf()
        ax = plt.subplot(1,2,1)
        plt.plot(freqs/1e9, Idat_g[pind], label='I_full_g')
        plt.plot(freqs/1e9, Idat_e[pind], label='I_full_e')
        plt.xlabel('Freq (GHz)')
        ax.legend()
        
        ax = plt.subplot(1,2,2)
        plt.plot(freqs/1e9, Qdat_g[pind], label='Q_full_g')
        plt.plot(freqs/1e9, Qdat_e[pind], label='Q_full_e')
        plt.xlabel('Freq (GHz)')
        ax.legend()
        
        fig.suptitle(f'Final I and Q (TripleDDC)| Rot: {round(theta_lo *180/np.pi,2)} Deg \n Carrier: {carrier_power}dBm | Sideband: {lsb_power}dBm | ModFreq: {mod_freq/1e6}MHz')
        fig.tight_layout()
        plt.show()
        # # #
        
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
    
    cavitygen.output = 'Off'
    qubitgen.output = 'Off'
    LO.output = 'Off'
    
    return full_data