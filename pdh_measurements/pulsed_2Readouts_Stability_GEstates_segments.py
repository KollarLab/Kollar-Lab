'''
8-25-21 AK modifying to normalize the amplitudes to the drive power.

8-25-21 AK modifying to return the raw data.

'''


import os
import time
import numpy as np 
import matplotlib.pyplot as plt
import mplcursors
mplcursors.cursor(multiple=True)

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process_2Readouts_segments
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

def pulsed_2Readouts_Stability_GEstates_segments(instruments, settings):
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

    ##Sweep settings for qubit
    start_freq  = exp_settings['start_freq']
    stop_freq   = exp_settings['stop_freq']
    freq_points = exp_settings['freq_points']
    freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power  = exp_settings['Qbit_power'] + Qbit_Attenuation
    Qbit_freq = exp_settings['Qbit_freq']
    avgs = exp_settings['averages']
    
    
    # Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    gen_power  = exp_settings['gen_power'] + CAV_Attenuation
    CAV_freq = exp_settings['CAV_freq']

    
    # num_points = int(exp_settings['stability_points'])  # number of stability data points
    
           
    # # # #
    amp_list = exp_settings['amp_list']
    phase_list = exp_settings['phase_list']
    mod_freq = exp_settings['mod_freq']
    usb_power, carr_power, lsb_power = exp_settings['upper_sideband_power'], exp_settings['carrier_power'], exp_settings['lower_sideband_power']
    # # # #    
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = gen_power
    cavitygen.enable_pulse()
    cavitygen.enable_IQ()
    cavitygen.phase = 0
    cavitygen.Output = 'On'
    
    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.phase = 0
    LO.output = 'On'
    
    
    qubitgen.Freq = Qbit_freq
    qubitgen.Power = Qbit_power
    if exp_settings['Quasi_CW']:
        qubitgen.disable_pulse()
        qubitgen.disable_IQ()
    else:
        qubitgen.enable_pulse()
        qubitgen.enable_IQ()
    
    ##Card settings
    configure_card(card, settings)

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
    segments = card.segments
    I_carr_g, Q_carr_g = np.zeros((segments, len(freqs))), np.zeros((segments, len(freqs)))
    powerdat_carr_g, phasedat_carr_g = np.zeros((segments, len(freqs))), np.zeros((segments, len(freqs)))
    I_sb_g, Q_sb_g = np.zeros((segments, len(freqs))), np.zeros((segments, len(freqs)))
    powerdat_sb_g, phasedat_sb_g = np.zeros((segments, len(freqs))), np.zeros((segments, len(freqs)))
    phasedat_diff_g = np.zeros((segments, len(freqs)))
    
    I_carr_e, Q_carr_e = np.zeros((segments, len(freqs))), np.zeros((segments, len(freqs)))
    powerdat_carr_e, phasedat_carr_e = np.zeros((segments, len(freqs))), np.zeros((segments, len(freqs)))
    I_sb_e, Q_sb_e = np.zeros((segments, len(freqs))), np.zeros((segments, len(freqs)))
    powerdat_sb_e, phasedat_sb_e = np.zeros((segments, len(freqs))), np.zeros((segments, len(freqs)))
    phasedat_diff_e = np.zeros((segments, len(freqs)))
    
    
    # ts = np.zeros(num_points)
    

    tstart = time.time()
    first_it = True

    print('asdf')
    # for numind in range(num_points):
        
    qubitgen.power = Qbit_power
    qubitgen.freq = Qbit_freq
    cavitygen.power = gen_power
    total_samples = card.samples 
    Is_carr_g, Qs_carr_g  = np.zeros((len(freqs), total_samples)), np.zeros((len(freqs), total_samples))
    Is_sb_g, Qs_sb_g  = np.zeros((len(freqs), total_samples)), np.zeros((len(freqs), total_samples))
    
    Is_carr_e, Qs_carr_e  = np.zeros((len(freqs), total_samples)), np.zeros((len(freqs), total_samples))
    Is_sb_e, Qs_sb_e  = np.zeros((len(freqs), total_samples)), np.zeros((len(freqs), total_samples))
    
    
    # for background
    Is_carr_b, Qs_carr_b  = np.zeros((len(freqs), total_samples)), np.zeros((len(freqs), total_samples))
    Is_sb_b, Qs_sb_b  = np.zeros((len(freqs), total_samples)), np.zeros((len(freqs), total_samples))
    
    # print(f'SGS_power:{gen_power - CAV_Attenuation}, Qbit_power:{Qbit_power - Qbit_Attenuation}')
    print(f'SGS_power:{gen_power}, Qbit_power:{Qbit_power - Qbit_Attenuation}')

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
        qubitgen.output = 'Off'
        carr_data_g, sb_data_g, xaxis_g = read_and_process_2Readouts_segments(card, settings, 
                                                                        plot=first_it, 
                                                                        IQstorage = True)
        
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
            carr_data_b, sb_data_b, xaxis_b = read_and_process_2Readouts_segments(card, settings, 
                                                                            plot=first_it, 
                                                                            IQstorage = True)
        else:
            carr_data_b = {key: np.zeros_like(value) for key, value in carr_data_g.items()}
            sb_data_b = {key: np.zeros_like(value) for key, value in sb_data_g.items()}

        qubitgen.output = 'On'
        time.sleep(0.1)
        carr_data_e, sb_data_e, xaxis_e = read_and_process_2Readouts_segments(card, settings, 
                                                                        plot=first_it, 
                                                                        IQstorage = True)
        
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
            carr_data_b, sb_data_b, xaxis_b = read_and_process_2Readouts_segments(card, settings, 
                                                                            plot=first_it, 
                                                                            IQstorage = True)
        else:
            carr_data_b = {key: np.zeros_like(value) for key, value in carr_data_e.items()} 
            sb_data_b = {key: np.zeros_like(value) for key, value in sb_data_e.items()}
                        
        I_final_carr_g = np.mean(carr_data_g['I_window'], axis=1, keepdims=True)  - np.mean(carr_data_b['I_window'], axis=1, keepdims=True) #compute <I> in the data window
        Q_final_carr_g = np.mean(carr_data_g['Q_window'], axis=1, keepdims=True)  - np.mean(carr_data_b['Q_window'], axis=1, keepdims=True) #compute <Q> in the data window
        
        I_final_carr_e = np.mean(carr_data_e['I_window'], axis=1, keepdims=True)  - np.mean(carr_data_b['I_window'], axis=1, keepdims=True) #compute <I> in the data window
        Q_final_carr_e = np.mean(carr_data_e['Q_window'], axis=1, keepdims=True)  - np.mean(carr_data_b['Q_window'], axis=1, keepdims=True) #compute <Q> in the data window
        
        I_final_sb_g = np.mean(sb_data_g['I_window'], axis=1, keepdims=True) - np.mean(sb_data_b['I_window'], axis=1, keepdims=True) #compute <I> in the data window
        Q_final_sb_g = np.mean(sb_data_g['Q_window'], axis=1, keepdims=True) - np.mean(sb_data_b['Q_window'], axis=1, keepdims=True) #compute <Q> in the data window
        
        I_final_sb_e = np.mean(sb_data_e['I_window'], axis=1, keepdims=True) - np.mean(sb_data_b['I_window'], axis=1, keepdims=True) #compute <I> in the data window
        Q_final_sb_e = np.mean(sb_data_e['Q_window'], axis=1, keepdims=True) - np.mean(sb_data_b['Q_window'], axis=1, keepdims=True) #compute <Q> in the data window
        
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(freqs))
            first_it = False
        # get a segment of the whole data   
        Is_carr_g[find,:], Qs_carr_g[find,:] = carr_data_g['I_full'][0], carr_data_g['Q_full'][0]
        Is_sb_g[find,:], Qs_sb_g[find,:] = sb_data_g['I_full'][0], sb_data_g['Q_full'][0]
        
        Is_carr_e[find,:], Qs_carr_e[find,:] = carr_data_e['I_full'][0], carr_data_e['Q_full'][0]
        Is_sb_e[find,:], Qs_sb_e[find,:] = sb_data_e['I_full'][0], sb_data_e['Q_full'][0]
        
        Is_carr_b[find,:], Qs_carr_b[find,:] = carr_data_b['I_full'][0], carr_data_b['Q_full'][0]
        Is_sb_b[find,:], Qs_sb_b[find,:] = sb_data_b['I_full'][0], sb_data_b['Q_full'][0]
        
        I_final_carr_g, Q_final_carr_g = I_final_carr_g.flatten(), Q_final_carr_g.flatten()
        I_carr_g[:, find], Q_carr_g[:, find] = I_final_carr_g.flatten(), Q_final_carr_g.flatten()
        powerdat_carr_g[:, find], phasedat_carr_g[:, find] = np.sqrt(I_final_carr_g**2 + Q_final_carr_g**2), np.arctan2(Q_final_carr_g, I_final_carr_g)*180/np.pi
        
        I_final_carr_e, Q_final_carr_e = I_final_carr_e.flatten(), Q_final_carr_e.flatten()
        I_carr_e[:, find], Q_carr_e[:, find] = I_final_carr_e, Q_final_carr_e
        powerdat_carr_e[:, find], phasedat_carr_e[:, find] = np.sqrt(I_final_carr_e**2 + Q_final_carr_e**2), np.arctan2(Q_final_carr_e, I_final_carr_e)*180/np.pi
        
        I_final_sb_g, Q_final_sb_g = I_final_sb_g.flatten(), Q_final_sb_g.flatten()
        I_sb_g[:, find], Q_sb_g[:, find] = I_final_sb_g, Q_final_sb_g
        powerdat_sb_g[:, find], phasedat_sb_g[:, find] = np.sqrt(I_final_sb_g**2 + Q_final_sb_g**2), np.arctan2(Q_final_sb_g, I_final_sb_g)*180/np.pi
        
        I_final_sb_e, Q_final_sb_e = I_final_sb_e.flatten(), Q_final_sb_e.flatten()
        I_sb_e[:, find], Q_sb_e[:, find] = I_final_sb_e, Q_final_sb_e
        powerdat_sb_e[:, find], phasedat_sb_e[:, find] = np.sqrt(I_final_sb_e**2 + Q_final_sb_e**2), np.arctan2(Q_final_sb_e, I_final_sb_e)*180/np.pi
        
        
        phasedat_diff_g[:, find] = phasedat_carr_g[:, find] - phasedat_sb_g[:, find]
        phasedat_diff_e[:, find] = phasedat_carr_e[:, find] - phasedat_sb_e[:, find]
        
    tfinal = time.time() - tstart
    ts = np.linspace(0, tfinal, segments)

    # For saving
    IQdat = {}
    IQdat['I_carr_g'], IQdat['Q_carr_g'] = I_carr_g, Q_carr_g
    IQdat['I_sb_g'], IQdat['Q_sb_g'] = I_sb_g, Q_sb_g
    IQdat['I_carr_e'], IQdat['Q_carr_e'] = I_carr_e, Q_carr_e
    IQdat['I_sb_e'], IQdat['Q_sb_e'] = I_sb_e, Q_sb_e
                    
    full_data = {}
    full_data['xaxis']  = freqs/1e9
    full_data['carr_mags_g'], full_data['carr_phases_g']   = powerdat_carr_g, phasedat_carr_g
    full_data['sb_mags_g'], full_data['sb_phases_g']   = powerdat_sb_g, phasedat_sb_g
    full_data['carr_mags_e'], full_data['carr_phases_e']   = powerdat_carr_e, phasedat_carr_e
    full_data['sb_mags_e'], full_data['sb_phases_e']   = powerdat_sb_e, phasedat_sb_e
    
    full_data['phase_diff_g'], full_data['phase_diff_e'] = phasedat_diff_g, phasedat_diff_e
    
    full_time = {}
    full_time['xaxis']  = xaxis*1e6
    full_time['Is_carr_g'], full_time['Qs_carr_g']    = Is_carr_g, Qs_carr_g
    full_time['Is_sb_g'], full_time['Qs_sb_g']   = Is_sb_g, Qs_sb_g
    full_time['Is_carr_e'], full_time['Qs_carr_e']    = Is_carr_e, Qs_carr_e
    full_time['Is_sb_e'], full_time['Qs_sb_e']   = Is_sb_e, Qs_sb_e
    
    full_time['Is_carr_b']   = Is_carr_b
    full_time['Qs_carr_b'] = Qs_carr_b
    full_time['Is_sb_b']   = Is_sb_b
    full_time['Qs_sb_b'] = Qs_sb_b
    
    # PLOTS 
    if freq_points > 1:
        # Plot simple mag scan
        fig = plt.figure(1212, figsize=(13,8))
        plt.clf()
        ax = plt.subplot(1,2,1)
        plt.plot(full_data['xaxis'], full_data['carr_mags_g'][0], label = 'g')
        plt.plot(full_data['xaxis'], full_data['carr_mags_e'][0], label = 'e')
        plt.xlabel('Freq (GHz)')
        plt.legend()
        plt.title('Carr mag')
        
        ax = plt.subplot(1,2,2)
        plt.plot(full_data['xaxis'], full_data['sb_mags_g'][0], label = 'g')
        plt.plot(full_data['xaxis'], full_data['sb_mags_e'][0], label = 'e')
        plt.xlabel('Freq (GHz)')
        plt.legend()
        plt.title('sb mag')
        
        plt.suptitle('Filename: {}'.format(filename))
        plt.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_fullMagPlot.png'), dpi = 150)
        plt.show()
        
        # Plot simple phase scan
        fig = plt.figure(1217, figsize=(13,8))
        plt.clf()
        ax = plt.subplot(2,2,1)
        plt.plot(full_data['xaxis'], full_data['carr_phases_g'][0], label = 'g')
        plt.plot(full_data['xaxis'], full_data['carr_phases_e'][0], label = 'e')
        plt.xlabel('Freq (GHz)')
        plt.legend()
        plt.title('Carr phase')
        
        ax = plt.subplot(2,2,2)
        plt.plot(full_data['xaxis'], full_data['sb_phases_g'][0], label = 'g')
        plt.plot(full_data['xaxis'], full_data['sb_phases_e'][0], label = 'e')
        plt.xlabel('Freq (GHz)')
        plt.legend()
        plt.title('sb phases')
        
        ax = plt.subplot(2,2,3)
        plt.plot(full_data['xaxis'], full_data['phase_diff_g'][0], label = 'g')
        plt.plot(full_data['xaxis'], full_data['phase_diff_e'][0], label = 'e')
        plt.xlabel('Freq (GHz)')
        plt.legend()
        plt.title('diff phases')
        
        plt.suptitle('Filename: {}'.format(filename))
        plt.tight_layout()
        mplcursors.cursor(multiple=True)
        plt.savefig(os.path.join(saveDir, filename+'_fullPhasePlot.png'), dpi = 150)
        plt.show()

        identifier = 'Power: {}dBm'.format(gen_power-CAV_Attenuation)
        new_filename = 'Raw_time_traces\n'+filename
        
        # # Plot IQ for each freq
        # fig = plt.figure(1216, figsize=(13,8))
        # plt.clf()
        # ax = plt.subplot(2,2,1)
        # ax.scatter(I_carr_g, Q_carr_g, s=4, label = 'g')
        # ax.scatter(I_carr_e, Q_carr_e, s=4, label = 'e')
        # plt.xlabel('I')
        # plt.ylabel('Q')
        # plt.title('Carrier IQ')
        # plt.legend()
        # plt.axhline(y = 0, color ="black", linestyle ="-")
        # plt.axvline(x = 0, color ="black", linestyle ="-") 
        # y_max = np.max(np.abs(ax.get_ylim()))
        # x_max = np.max(np.abs(ax.get_xlim()))
        # list = np.array([x_max, y_max])
        # max = np.max(list)
        # ax.set_ylim(ymin=-max, ymax=max)
        # ax.set_xlim(xmin=-max, xmax=max)
        
        # ax = plt.subplot(2,2,2)
        # ax.scatter(I_sb_g, Q_sb_g, s=4, label = 'g')
        # ax.scatter(I_sb_e, Q_sb_e, s=4, label = 'e')
        # plt.xlabel('I')
        # plt.ylabel('Q')
        # plt.title('sideband IQ')
        # plt.legend()
        # plt.axhline(y = 0, color ="black", linestyle ="-")
        # plt.axvline(x = 0, color ="black", linestyle ="-") 
        # y_max = np.max(np.abs(ax.get_ylim()))
        # x_max = np.max(np.abs(ax.get_xlim()))
        # list = np.array([x_max, y_max])
        # max = np.max(list)
        # ax.set_ylim(ymin=-max, ymax=max)
        # ax.set_xlim(xmin=-max, xmax=max)
        
        # ax = plt.subplot(2,2,3)
        # ax.scatter(powerdat_carr_g[0,:]*np.cos(phasedat_carr_g[0,:]*np.pi/180), powerdat_carr_g[0,:]*np.sin(phasedat_carr_g[0,:]*np.pi/180), s=4, label = 'g')
        # ax.scatter(powerdat_carr_e[0,:]*np.cos(phasedat_carr_e[0,:]*np.pi/180), powerdat_carr_e[0,:]*np.sin(phasedat_carr_e[0,:]*np.pi/180), s=4, label = 'e')
        # plt.xlabel('I')
        # plt.ylabel('Q')
        # plt.legend()
        # plt.title('Carrier IQ reconstruct')
        # plt.axhline(y = 0, color ="black", linestyle ="-")
        # plt.axvline(x = 0, color ="black", linestyle ="-") 
        # y_max = np.max(np.abs(ax.get_ylim()))
        # x_max = np.max(np.abs(ax.get_xlim()))
        # list = np.array([x_max, y_max])
        # max = np.max(list)
        # ax.set_ylim(ymin=-max, ymax=max)
        # ax.set_xlim(xmin=-max, xmax=max)
        
        # ax = plt.subplot(2,2,4)
        # ax.scatter(powerdat_carr_g[0,:]*np.cos(phasedat_diff_g[0,:]*np.pi/180), powerdat_carr_g[0,:]*np.sin(phasedat_diff_g[0,:]*np.pi/180), s=4, label = 'g')
        # ax.scatter(powerdat_carr_e[0,:]*np.cos(phasedat_diff_e[0,:]*np.pi/180), powerdat_carr_e[0,:]*np.sin(phasedat_diff_e[0,:]*np.pi/180), s=4, label = 'e')
        # plt.xlabel('I')
        # plt.ylabel('Q')
        # plt.legend()
        # plt.title('Difference IQ')
        # plt.axhline(y = 0, color ="black", linestyle ="-")
        # plt.axvline(x = 0, color ="black", linestyle ="-") 
        # y_max = np.max(np.abs(ax.get_ylim()))
        # x_max = np.max(np.abs(ax.get_xlim()))
        # list = np.array([x_max, y_max])
        # max = np.max(list)
        # ax.set_ylim(ymin=-max, ymax=max)
        # ax.set_xlim(xmin=-max, xmax=max)
        
        # fig.suptitle('IQ freq sweep')
        # fig.tight_layout()
        # plt.show()
        # plt.savefig(os.path.join(saveDir, filename+'_IQFreqSweep.png'), dpi = 150)
        
        # Plot raw time traces
        fig = plt.figure(1215, figsize=(13,8))
        plt.clf()
        ax = plt.subplot(2,2,1)
        plt.plot(full_time['xaxis'], carr_data_g['I_full'][0], label = 'I_full_g')
        plt.plot(full_time['xaxis'], carr_data_g['Q_full'][0], label = 'Q_full_g')
        plt.legend()
        plt.xlabel('Time (us)')
        plt.title('single carr trace')
        
        ax = plt.subplot(2,2,2)
        plt.plot(full_time['xaxis'], sb_data_g['I_full'][0], label = 'I_full_g')
        plt.plot(full_time['xaxis'], sb_data_g['Q_full'][0], label = 'Q_full_g')
        plt.legend()
        plt.xlabel('Time (us)')
        plt.title('single sb trace')
        
        ax = plt.subplot(2,2,3)
        plt.plot(full_time['xaxis'], carr_data_e['I_full'][0], label = 'I_full_e')
        plt.plot(full_time['xaxis'], carr_data_e['Q_full'][0], label = 'Q_full_e')
        plt.legend()
        plt.xlabel('Time (us)')
        plt.title('single carr trace')
        
        ax = plt.subplot(2,2,4)
        plt.plot(full_time['xaxis'], sb_data_e['I_full'][0], label = 'I_full_e')
        plt.plot(full_time['xaxis'], sb_data_e['Q_full'][0], label = 'Q_full_e')
        plt.legend()
        plt.xlabel('Time (us)')
        plt.title('single sb trace')
        
        plt.suptitle('Filename: {}'.format(new_filename))
        plt.tight_layout()
        mplcursors.cursor(multiple=True)
        plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces.png'), dpi = 150)
        plt.show()

    else:
        fig = plt.figure(1213, figsize=(13,8))
        plt.clf()
        ax = plt.subplot(2,2,1)
        plt.plot(ts, phasedat_carr_g[:,0], label = 'g')
        plt.plot(ts, phasedat_carr_e[:,0], label = 'e')
        plt.xlabel('Time (s)')
        plt.legend()
        plt.title('Carrier phase')
        
        ax = plt.subplot(2,2,2)
        plt.plot(ts, phasedat_sb_g[:,0], label = 'g')
        plt.plot(ts, phasedat_sb_e[:,0], label = 'e')
        plt.xlabel('Time (s)')
        plt.legend()
        plt.title('Sideband phase')
        
        ax = plt.subplot(2,2,3)
        plt.plot(ts, phasedat_diff_g[:,0], label ='g')
        plt.plot(ts, phasedat_diff_e[:,0], label ='e')
        plt.xlabel('Time (s)')
        plt.title('Difference phase')
        
        fig.suptitle(f'{segments} points | avgs:{avgs} | Carrier: {carr_power}dBm | Sideband: {lsb_power}dBm | ModFreq: {mod_freq/1e6}MHz')
        fig.tight_layout()
        mplcursors.cursor(multiple=True)
        plt.savefig(os.path.join(saveDir, filename+'_PhaseStability.png'), dpi = 150)
        plt.show()
        
        fig2 = plt.figure(1214, figsize=(13,8))
        plt.clf()
        ax = plt.subplot(2,3,1)
        ax.scatter(I_carr_g, Q_carr_g, s=4, label = 'g')
        ax.scatter(I_carr_e, Q_carr_e, s=4, label = 'e')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('Carrier IQ')
        plt.legend()
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        y_max = np.max(np.abs(ax.get_ylim()))
        x_max = np.max(np.abs(ax.get_xlim()))
        list = np.array([x_max, y_max])
        max = np.max(list)
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        
        ax = plt.subplot(2,3,2)
        ax.scatter(I_sb_g, Q_sb_g, s=4, label = 'g')
        ax.scatter(I_sb_e, Q_sb_e, s=4, label = 'e')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('sideband IQ')
        plt.legend()
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        y_max = np.max(np.abs(ax.get_ylim()))
        x_max = np.max(np.abs(ax.get_xlim()))
        list = np.array([x_max, y_max])
        max = np.max(list)
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        
        ax = plt.subplot(2,3,3)
        ax.scatter(powerdat_carr_g[:,0]*np.cos(phasedat_carr_g[:,0]*np.pi/180), powerdat_carr_g[:,0]*np.sin(phasedat_carr_g[:,0]*np.pi/180), s=4, label = 'g')
        ax.scatter(powerdat_carr_e[:,0]*np.cos(phasedat_carr_e[:,0]*np.pi/180), powerdat_carr_e[:,0]*np.sin(phasedat_carr_e[:,0]*np.pi/180), s=4, label = 'e')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.legend()
        plt.title('Carrier IQ reconstruct')
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        y_max = np.max(np.abs(ax.get_ylim()))
        x_max = np.max(np.abs(ax.get_xlim()))
        list = np.array([x_max, y_max])
        max = np.max(list)
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        
        ax = plt.subplot(2,3,4)
        ax.scatter(powerdat_carr_g[:,0]*np.cos(phasedat_diff_g[:,0]*np.pi/180), powerdat_carr_g[:,0]*np.sin(phasedat_diff_g[:,0]*np.pi/180), s=4, label = 'g')
        ax.scatter(powerdat_carr_e[:,0]*np.cos(phasedat_diff_e[:,0]*np.pi/180), powerdat_carr_e[:,0]*np.sin(phasedat_diff_e[:,0]*np.pi/180), s=4, label = 'e')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.legend()
        plt.title('Difference IQ')
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        y_max = np.max(np.abs(ax.get_ylim()))
        x_max = np.max(np.abs(ax.get_xlim()))
        list = np.array([x_max, y_max])
        max = np.max(list)
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)

        ax = plt.subplot(2,3,5)
        ax.scatter(I_carr_g, Q_carr_g, s=4, label = 'g')
        # ax.scatter(I_carr_e, Q_carr_e, s=4, label = 'e')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('Carrier IQ (ground)')
        plt.legend()
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        y_max = np.max(np.abs(ax.get_ylim()))
        x_max = np.max(np.abs(ax.get_xlim()))
        list = np.array([x_max, y_max])
        max = np.max(list)
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)

        ax = plt.subplot(2,3,6)
        # ax.scatter(I_carr_g, Q_carr_g, s=4, label = 'g')
        ax.scatter(I_carr_e, Q_carr_e, s=4, label = 'e', color='orange')
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('Carrier IQ (excited)')
        plt.legend()
        plt.axhline(y = 0, color ="black", linestyle ="-")
        plt.axvline(x = 0, color ="black", linestyle ="-") 
        y_max = np.max(np.abs(ax.get_ylim()))
        x_max = np.max(np.abs(ax.get_xlim()))
        list = np.array([x_max, y_max])
        max = np.max(list)
        ax.set_ylim(ymin=-max, ymax=max)
        ax.set_xlim(xmin=-max, xmax=max)
        
        fig2.suptitle(f'{segments} points | avgs:{avgs} | Carrier: {carr_power}dBm | Sideband: {lsb_power}dBm | ModFreq: {mod_freq/1e6}MHz')
        fig2.tight_layout()
        mplcursors.cursor(multiple=True)
        plt.savefig(os.path.join(saveDir, filename+'_IQStability.png'), dpi = 150)
        plt.show()


    userfuncs.SaveFull(saveDir, filename, ['freqs', 'IQdat','full_data', 'full_time'],
                            locals(), expsettings=settings, instruments=instruments, saveHWsettings=first_it)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))

    LO.output = 0
    cavitygen.output = 0
    qubitgen.output = 0
    
    return full_data