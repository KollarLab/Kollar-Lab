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
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process
from pdh_measurements.scheduler_pdh import scheduler_pdh
from utility.measurement_helpers import get_amp_comps, total_power
from utility.userfits_v2 import fit_model

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'ac_stark'
    
    #Sweep parameters
    settings['cav_start_freq']   = 4.15*1e9  
    settings['cav_stop_freq']    = 4.25*1e9 
    settings['cav_freq_points']  = 50

    settings['Qbit_start_freq']   = 4.15*1e9  
    settings['Qbit_stop_freq']    = 4.25*1e9 
    settings['Qbit_freq_points']  = 50

    settings['cav_start_power']  = -90
    settings['cav_stop_power'] = -45
    settings['cav_power_points'] = 31

    #Card settings
    settings['segments'] = 1
    settings['reads']    = 1
    settings['averages'] = 5e3
    
    settings['mod_freq'] = 20e6
    
    return settings

def ac_stark(instruments, settings):
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

    # Cavity readout settings
    cav_start_freq  = exp_settings['cav_start_freq']
    cav_stop_freq   = exp_settings['cav_stop_freq']
    cav_freq_points = exp_settings['cav_freq_points']
    cav_freqs  = np.round(np.linspace(cav_start_freq,cav_stop_freq,cav_freq_points),-3)

    cav_start_power = exp_settings['cav_start_power']
    cav_stop_power = exp_settings['cav_stop_power']
    cav_power_points = exp_settings['cav_power_points']
    cav_powers = np.round(np.linspace(cav_start_power,cav_stop_power,cav_power_points),2)

    stimulation_pulse_length = exp_settings['stimulation_pulse_length']
    ring_down_time = exp_settings['ring_down_time']
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    cavitygen.freq = cav_freqs[0]
    cavitygen.phase = 0
    cavitygen.enable_pulse()
    cavitygen.enable_IQ()

    LO.power = 12
    LO.freq = cav_freqs[0] - exp_globals['IF']    
    LO.phase = 0

    # Qubit settings
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power  = exp_settings['Qbit_power'] + Qbit_Attenuation
    Qbit_start_freq = exp_settings['Qbit_start_freq']
    Qbit_stop_freq = exp_settings['Qbit_stop_freq']
    Qbit_freq_points = exp_settings['Qbit_freq_points']
    Qbit_freqs = np.round(np.linspace(Qbit_start_freq,Qbit_stop_freq,Qbit_freq_points),-3)

    qubitgen.freq = Qbit_freqs[0]
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
    phase_list = [0,0,0]
    mod_freq = exp_settings['mod_freq']

    ##HDAWG settings  
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    init_buffer  = m_pulse['init_buffer']
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
    cav_sweep_mag = np.zeros((len(cav_powers), len(cav_freqs)))
    cav_sweep_phase = np.zeros((len(cav_powers), len(cav_freqs)))
    cav_readout_freq = np.zeros(len(cav_powers))
    
    Qbit_sweep_mag = np.zeros((len(cav_powers), len(Qbit_freqs)))
    Qbit_sweep_phase = np.zeros((len(cav_powers), len(Qbit_freqs)))
    Qbit_max_freqs = np.zeros(len(cav_powers))

    # for time traces
    total_samples = int(card.samples) 
    Is_Qbit_sweep, Qs_Qbit_sweep = np.zeros((len(cav_powers), len(Qbit_freqs), total_samples)), np.zeros((len(cav_powers), len(Qbit_freqs), total_samples))
    Is_Qbit_sweep_b, Qs_Qbit_sweep_b = np.zeros((len(cav_powers), len(Qbit_freqs), total_samples)), np.zeros((len(cav_powers), len(Qbit_freqs), total_samples))

    drive_amps_lin = np.sqrt(10**(cav_powers/10))
    tstart = time.time()
    first_it = True

    readout_power = exp_settings['readout_power']
    hold = exp_settings['hold_time']

    #########################################

    # For each cav power point load the AWG for cavity readout and stimulation pulse
    for pind in range(cav_power_points):

        carr_power = cav_powers[pind]
        cav_power_list = [-100, carr_power, -100]

        # calculate the total power of the 3 components | and increase by 10 dBm
        total_comp_power = total_power(cav_power_list) + 10
        # set gen power
        if total_comp_power > 0:
            gen_power = total_comp_power
        else:
            gen_power = -5
            # gen_power = 5
        # Use gen power to get amps of the components
        cav_amp_list = get_amp_comps(cav_power_list, gen_power)
        readout_amp_list = get_amp_comps([-100, readout_power,-100], gen_power)
        
        #update gen power to include cavity attenuation
        gen_power = gen_power + CAV_Attenuation  # this is only compatible with 10db attenuation

        if gen_power > 12:
        # if gen_power > 20:
            raise ValueError(f'SGS power of {gen_power} too high')
        
        # Load AWG for cavity readout, stimulation pulse and qubit

        # For cavity readout / final measurement
        position = start_time
        cavity_I.add_pulse('pdh_I', position=position,
                           amp_list = readout_amp_list, phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = window_time)

        cavity_Q.add_pulse('pdh_Q', position=position,
                            amp_list = readout_amp_list, phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = window_time)
        # For stimulation pulse
        if exp_settings['with_stimulation'] == True:
            stim_position = start_time - ring_down_time - stimulation_pulse_length #+ window_time #- init_buffer
            if exp_settings['stimulation_on_resonance'] == True: 
                cavity_I.add_pulse('pdh_I', position=stim_position,
                           amp_list = cav_amp_list, phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = stimulation_pulse_length)
                cavity_Q.add_pulse('pdh_Q', position=stim_position,
                            amp_list = cav_amp_list, phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = stimulation_pulse_length)
            elif exp_settings['stimulation_on_resonance'] == False:
                cavity_I.add_pulse('pdh_I', position=stim_position,
                           amp_list = [cav_amp_list[1], cav_amp_list[0], cav_amp_list[2]], phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = stimulation_pulse_length)
                cavity_Q.add_pulse('pdh_Q', position=stim_position,
                            amp_list = cav_amp_list, phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = stimulation_pulse_length)
                
            cavity_marker.add_window(stim_position, start_time+window_time)

        cavity_marker.add_window(start_time, start_time+window_time)

        qubit_I =  awg_sched.analog_channels['Qubit_I']
        qubit_marker  = awg_sched.digital_channels['Qubit_enable']
        delay = q_pulse['delay']
        sigma = q_pulse['sigma']
        num_sigma = q_pulse['num_sigma']
        q_pulse['hold_time'] = hold 
        
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
            
            digLO_sin = np.sin(2*np.pi*exp_globals['IF']*xaxis)
            digLO_cos = np.cos(2*np.pi*exp_globals['IF']*xaxis)
            
            #store in settings so that the processing functions can get to them
            settings['digLO_sin'] = digLO_sin 
            settings['digLO_cos'] = digLO_cos
            settings['LPF'] = LPF


        cavitygen.power = gen_power
    
        # First, cavity sweep (make code faster by doing this only for the first power)
        print(f'Current Powers| Readout:{carr_power}, SGS:{gen_power}')
    
        for find in range(0, len(cav_freqs)):

            cav_freq = cav_freqs[find]
            cavitygen.freq = cav_freq
            cavitygen.phase = 0
            cavitygen.output = 'On'
            LO.freq = cav_freq - exp_globals['IF']
            LO.phase = 0
            LO.output = 'On'       
            qubitgen.output = 'Off'

            I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            
            I_final, Q_final = np.mean(I_window), np.mean(Q_window)

            cav_sweep_mag[pind, find] = np.sqrt(I_final**2 + Q_final**2)
            cav_sweep_phase[pind, find] = np.arctan2(Q_final, I_final)*180/np.pi

            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, len(cav_powers)*len(cav_freqs)*len(Qbit_freqs))
                first_it = False

        # plot the first cavity sweep
        if pind == 0:
            fig  = plt.figure(12321, figsize=(9,6))
            plt.clf()
            ax = plt.subplot(1,1,1)
            xaxis = cav_freqs/1e9
            yaxis = cav_sweep_mag[pind,:]
            params_cav = fit_model(xaxis*1e9, yaxis, 'lorenz', plot=False, ax=ax)
            linewidth = round(2*params_cav['sigma']/1e6,2)
            res_freq = round(params_cav['center']/1e9,5)
            plt.plot(xaxis, yaxis)
            plt.xlabel('Freq (GHz)')
            plt.title(f'Cavity Sweep | Res Freq:{res_freq}GHz | Linewidth:{linewidth}MHz')
            plt.tight_layout()
            mplcursors.cursor(multiple=True)
            plt.show()

        # get the cav resonance (min) for the above sweep
        # readout_freq_index = np.argmin(cav_sweep_mag[pind,:])
        # readout_freq = cav_freqs[readout_freq_index]
        readout_freq = res_freq*1e9
        cav_readout_freq[pind] = readout_freq

        # do qubit freq sweep for the readout frequency, readout power, and qubit freq
        print(f'Starting qubit scan for {cav_powers[pind]} readout power')
        first_it = True
        for find in range(0, len(Qbit_freqs)):
            cav_freq = exp_settings['cav_center'] #readout_freq - 10e6
            cavitygen.freq = cav_freq
            cavitygen.output = 'On'
            LO.freq = cav_freq - exp_globals['IF']
            LO.output = 'On'    
            qubitgen.freq = Qbit_freqs[find]   
            qubitgen.output = 'On'

            I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            if exp_settings['subtract_background']:
                #Acquire background trace
                qubitgen.output='Off'
#                time.sleep(0.1)
                I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process(card, settings, 
                                                                 plot=first_it, 
                                                                 IQstorage = True)
            else:
                I_window_b, Q_window_b, I_full_b, Q_full_b = 0,0,0,0
            
            if first_it:
                first_it=False

    
            I_sig, Q_sig   = [np.mean(I_window), np.mean(Q_window)] 
            I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)] 
           
            
            I_final = I_sig-I_back 
            Q_final = Q_sig-Q_back 
            I_full_net = I_full-I_full_b 
            Q_full_net = Q_full-Q_full_b 

            Is_Qbit_sweep[pind, find, :],  Qs_Qbit_sweep[pind, find, :] = I_full_net, Q_full_net
            Is_Qbit_sweep_b[pind, find, :],  Qs_Qbit_sweep_b[pind, find, :] = I_full_b, Q_full_b
            
            Qbit_sweep_mag[pind, find] = np.sqrt(I_final**2 + Q_final**2)
            Qbit_sweep_phase[pind, find] = np.arctan2(Q_final, I_final)*180/np.pi

        Qbit_max_index = np.argmax(Qbit_sweep_mag[pind,:])
        Qbit_max_freqs[pind] = Qbit_freqs[Qbit_max_index] 
        
    full_cav_data = {}
    full_cav_data['cav_freqs'] = cav_freqs
    full_cav_data['powers'] = cav_powers
    full_cav_data['cav_sweep_mag'] = cav_sweep_mag
    full_cav_data['cav_sweep_phase'] = cav_sweep_phase
    full_cav_data['cav_readout_freq'] = cav_readout_freq

    full_Qbit_data = {}
    full_Qbit_data['Qbit_freqs'] = Qbit_freqs
    full_Qbit_data['Qbit_sweep_mag'] = Qbit_sweep_mag
    full_Qbit_data['Qbit-sweep_phase'] = Qbit_sweep_phase
    full_Qbit_data['Qbit_max_freqs'] = Qbit_max_freqs
    full_Qbit_data['Is'], full_Qbit_data['Qs'] = Is_Qbit_sweep, Qs_Qbit_sweep
    full_Qbit_data['Is_back'], full_Qbit_data['Qs_back'] = Is_Qbit_sweep_b, Qs_Qbit_sweep_b

    # plot full cavity sweep for all powers
    plot_data_cav = {}
    plot_data_cav['xaxis']  = cav_freqs/1e9
    plot_data_cav['mags']   = np.transpose(np.transpose(cav_sweep_mag)/drive_amps_lin)
    plot_data_cav['phases'] = cav_sweep_phase

    single_data_cav = {}
    single_data_cav['xaxis'] = cav_freqs/1e9
    single_data_cav['mag']   = cav_sweep_mag[-1]
    single_data_cav['phase'] = cav_sweep_phase[-1]

    yaxis  = cav_powers
    labels = ['Freq (GHz)', 'Power (dBm)']
    filename_cav = filename+'_CavFreqSweep'
    simplescan_plot(plot_data_cav, single_data_cav, yaxis, filename_cav, labels, identifier='', fig_num=1020, IQdata = False) #normalized to drive level
    plt.savefig(os.path.join(saveDir, filename_cav+'_fullColorPlot.png'), dpi = 150)

    # plot full Qbit freq sweep
    plot_data_Qbit = {}
    plot_data_Qbit['xaxis']  = Qbit_freqs/1e9
    plot_data_Qbit['mags']   = np.transpose(np.transpose(Qbit_sweep_mag)/drive_amps_lin)
    plot_data_Qbit['phases'] = Qbit_sweep_phase

    single_data_Qbit = {}
    single_data_Qbit['xaxis'] = Qbit_freqs/1e9
    single_data_Qbit['mag']   = Qbit_sweep_mag[-1]
    single_data_Qbit['phase'] = Qbit_sweep_phase[-1]

    yaxis  = cav_powers
    labels = ['Freq (GHz)', 'Power (dBm)']
    filename_Qbit = filename+'_QbitFreqSweep'
    simplescan_plot(plot_data_Qbit, single_data_Qbit, yaxis, filename_Qbit, labels, identifier='', fig_num=1021, IQdata = False) #normalized to drive level
    plt.savefig(os.path.join(saveDir, filename_Qbit+'_fullColorPlot.png'), dpi = 150)

    # Quickly plot the first and last qubit spec scan
    fig = plt.figure(12322, figsize=(12,9))
    plt.clf()
    ax = plt.subplot(1,1,1)
    plt.plot(Qbit_freqs/1e9, Qbit_sweep_mag[0], label=f'{cav_powers[0]}')
    plt.plot(Qbit_freqs/1e9, Qbit_sweep_mag[-1], label=f'{cav_powers[-1]}')
    plt.xlabel('Spec Freq (GHz)')
    plt.ylabel('Amp (V)')
    plt.legend()
    plt.tight_layout()
    mplcursors.cursor(multiple=True)
    plt.show()


    userfuncs.SaveFull(saveDir, filename, ['full_cav_data', 'full_Qbit_data'],
                             locals(), expsettings=settings, instruments=instruments, saveHWsettings=first_it)
        
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
    
    cavitygen.output = 'Off'
    qubitgen.output = 'Off'
    LO.output = 'Off'
    
    return full_cav_data, full_Qbit_data