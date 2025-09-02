'''
8-25-21 AK modifying to normalize the amplitudes to the drive power. Undid that.

9-2-21 AK made it return the data

'''


import os
import time
import numpy as np
import matplotlib.pyplot as plt
import mplcursors
mplcursors.cursor(multiple=True)

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process, read_and_process_singleshot, read_and_process_2meas, configure_card_2meas
from utility.scheduler import scheduler
from pdh_measurements.scheduler_pdh import scheduler_pdh
from utility.measurement_helpers import get_amp_comps, total_power

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'QND_TwoReadout'
    
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
    settings['reverse']   = False
    settings['num_save'] = 1
    
    #background_subtraction (by taking reference trace with no qubit drive power)
    settings['subtract_background'] = False
    
    return settings


def QND_TwoReadout(instruments, settings):
    
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

    # Qubit settings 
    Qbit_Attenuation = exp_globals['Qbit_Attenuation'] 
    Qbit_freq = exp_settings['Qbit_freq']
    Qbit_power = exp_settings['Qbit_power'] + Qbit_Attenuation

    # Readout/Measurment settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = exp_settings['CAV_Power'] + CAV_Attenuation
    CAV_freq  = exp_settings['CAV_freq']

    # Stimulation pulse settings
    phase_list = [0,0,0]
    mod_freq = exp_settings['mod_freq']
    stim_start_power = exp_settings['stimulation_start_power'] + CAV_Attenuation
    stim_stop_power = exp_settings['stimulation_stop_power'] + CAV_Attenuation
    stim_power_points = exp_settings['stimulation_power_points']
    stim_powers = np.round(np.linspace(stim_start_power, stim_stop_power, stim_power_points),2)

    ## sweep settings (mainly for trans scan setup check)
    start_freq  = exp_settings['start_freq']
    stop_freq   = exp_settings['stop_freq']
    freq_points = exp_settings['freq_points']
    if exp_settings['reverse']:
        cav_freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points)[::-1],-3)
    else:
        cav_freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
    
    ## Generator settings
    cavitygen.freq   = CAV_freq
    cavitygen.power  = CAV_power
    cavitygen.phase = 0
    cavitygen.enable_pulse()
    cavitygen.enable_IQ()

    qubitgen.freq = Qbit_freq
    qubitgen.power = Qbit_power
    qubitgen.phase = 0
    if exp_settings['Quasi_CW']:
        qubitgen.disable_pulse()
        qubitgen.disable_IQ()
    else:
        qubitgen.enable_pulse()
        qubitgen.enable_IQ()
    
    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.phase = 0
    
    ## Card config
    configure_card_2meas(card, settings)

    ## HDAWG settings
    
    ## Sequencer program
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    init_buffer  = m_pulse['init_buffer']
    q_pulse = exp_globals['qubit_pulse']
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    stim_buffer = exp_settings['stimulation_time_buffer']
    
    awg_sched = scheduler_pdh(total_time=start_time+4*window_time+2*stim_buffer, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='Qubit_I')
    awg_sched.add_analog_channel(2, name='Qubit_Q')
    awg_sched.add_analog_channel(3, name='Cavity_I')
    awg_sched.add_analog_channel(4, name='Cavity_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(3, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)

    # data saving arrays
    first_it = True
    segments = int(card.segments)
    total_samples = int(card.samples)

    cav_sweep_mag1 = np.zeros((len(stim_powers), len(cav_freqs), segments))
    cav_sweep_phase1 = np.zeros((len(stim_powers), len(cav_freqs), segments))

    cav_sweep_mag2 = np.zeros((len(stim_powers), len(cav_freqs), segments))
    cav_sweep_phase2 = np.zeros((len(stim_powers), len(cav_freqs), segments))

    I1_all_powers = np.zeros((len(stim_powers), len(cav_freqs), segments))
    Q1_all_powers = np.zeros((len(stim_powers), len(cav_freqs), segments))
    I2_all_powers = np.zeros((len(stim_powers), len(cav_freqs), segments))
    Q2_all_powers = np.zeros((len(stim_powers), len(cav_freqs), segments))

    Is_all_powers, Qs_all_powers = np.zeros((len(stim_powers), len(cav_freqs), total_samples)), np.zeros((len(stim_powers), len(cav_freqs), total_samples))
    Is_all_powers_b, Qs_all_powers_b = np.zeros((len(stim_powers), len(cav_freqs), total_samples)), np.zeros((len(stim_powers), len(cav_freqs), total_samples))
    
    for pind in range(len(stim_powers)):
        # get amp for the corresponding readout and stimulation power
        readout_power_list = [-100, CAV_power, -100]
        stim_power = stim_powers[pind]
        stim_power_list = [-100, stim_power, -100]

        # Set SGS power
        gen_power = 10

        # get amps for the signals
        stim_amp_list = get_amp_comps(stim_power_list, gen_power)
        readout_amp_list = get_amp_comps(readout_power_list, gen_power)

        # Load HDAWG for cavity readout, stimulation pulse, and qubit
        cavity_I = awg_sched.analog_channels['Cavity_I']
        cavity_Q = awg_sched.analog_channels['Cavity_Q'] 
        cavity_marker = awg_sched.digital_channels['Cavity_enable']

        # get positions
        # final_meas_pos = start_time
        # stim_position = final_meas_pos - stim_buffer - window_time
        # first_meas_pos = stim_position - stim_buffer - window_time

        first_meas_pos = start_time
        stim_position = first_meas_pos + stim_buffer + window_time
        final_meas_pos = stim_position + stim_buffer + window_time

        # For first measurement
        position = first_meas_pos
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
            stim_position = stim_position
            if exp_settings['stimulation_on_resonance'] == True: 
                cavity_I.add_pulse('pdh_I', position=stim_position,
                           amp_list = stim_amp_list, phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = window_time)
                cavity_Q.add_pulse('pdh_Q', position=stim_position,
                            amp_list = stim_amp_list, phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = window_time)
            elif exp_settings['stimulation_on_resonance'] == False:
                cavity_I.add_pulse('pdh_I', position=stim_position,
                           amp_list = [stim_amp_list[1], stim_amp_list[0], stim_amp_list[2]], phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = window_time)
                cavity_Q.add_pulse('pdh_Q', position=stim_position,
                            amp_list = [stim_amp_list[1], stim_amp_list[0], stim_amp_list[2]], phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = window_time)
        # For second measurement
        position = final_meas_pos
        cavity_I.add_pulse('pdh_I', position=position,
                           amp_list = readout_amp_list, phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = window_time)

        cavity_Q.add_pulse('pdh_Q', position=position,
                            amp_list = readout_amp_list, phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = window_time)
        cavity_marker.add_window(first_meas_pos, final_meas_pos+window_time)

        # For qubit pulse
        if exp_settings['Qstate'] == 'Excited':
            qubit_I =  awg_sched.analog_channels['Qubit_I']
            qubit_marker  = awg_sched.digital_channels['Qubit_enable']
            delay = q_pulse['delay']
            sigma = q_pulse['sigma']
            num_sigma = q_pulse['num_sigma']
            position = first_meas_pos-delay-num_sigma*sigma-q_pulse['hold_time']
            qubit_I.add_pulse('gaussian_square', position=position, 
                                    amplitude=q_pulse['piAmp']/2, length = q_pulse['hold_time'], 
                                    ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
            
            qubit_marker.add_window(position-160e-9, position+2*160e-9+q_pulse['hold_time'])
        
        if first_it:
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
            
            xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
            digLO_sin = np.sin(2*np.pi*exp_globals['IF']*xaxis)
            digLO_cos = np.cos(2*np.pi*exp_globals['IF']*xaxis)
            
            #store in settings so that the processing functions can get to them
            settings['digLO_sin'] = digLO_sin 
            settings['digLO_cos'] = digLO_cos
            settings['LPF'] = LPF
         
        print('Current power:{}, max:{}'.format(stim_powers[pind]-Qbit_Attenuation, stim_powers[-1]-Qbit_Attenuation))
    
        for find in range(0, len(cav_freqs)):
            
            if first_it:
                tstart = time.time()
            
            ##Acquire signal
            if exp_settings['Qstate'] == 'Excited':
                qubitgen.freq = Qbit_freq
                qubitgen.power = Qbit_power
                qubitgen.phase = 0
                qubitgen.output = 'On'
            else:
                qubitgen.output = 'Off'
            
            cav_freq = cav_freqs[find]
            cavitygen.freq = cav_freq
            cavitygen.power = gen_power
            cavitygen.phase = 0
            cavitygen.output = 'On'
            
            LO.freq = cav_freq - exp_globals['IF']
            LO.phase = 0
            LO.output = 'On'

            I_window, Q_window, I_window2, Q_window2, I_full, Q_full, xaxis = read_and_process_2meas(card, settings, 
                                                                        plot=first_it, 
                                                                        IQstorage = True)
            if exp_settings['subtract_background']:
                #Acquire background trace
                qubitgen.output='Off'
                I_window_b, Q_window_b, I_window_b2, Q_window_b2, I_full_b, Q_full_b, xaxis_b = read_and_process_2meas(card, settings, 
                                                                    plot=first_it, 
                                                                    IQstorage = True)
            else:
                I_window_b, Q_window_b, I_full_b, Q_full_b = np.zeros(I_window.shape), np.zeros(Q_window.shape), np.zeros(I_full.shape), np.zeros(Q_full.shape)
                I_window_b2, Q_window_b2 = np.zeros(I_window2.shape), np.zeros(Q_window2.shape)
            
            if first_it:
                first_it=False

            ##Useful handles for variables
            I_final = np.mean(I_window, axis=1, keepdims=True) - np.mean(I_window_b, axis=1, keepdims=True)
            Q_final = np.mean(Q_window, axis=1, keepdims=True) - np.mean(Q_window_b, axis=1, keepdims=True)
            I_final2 = np.mean(I_window2, axis=1, keepdims=True) - np.mean(I_window_b2, axis=1, keepdims=True)
            Q_final2 = np.mean(Q_window2, axis=1, keepdims=True) - np.mean(Q_window_b2, axis=1, keepdims=True)

            I_final, Q_final = I_final.flatten(), Q_final.flatten()
            I_final2, Q_final2 = I_final2.flatten(), Q_final2.flatten()

            I1_all_powers[pind, find, :], Q1_all_powers[pind, find, :] = I_final, Q_final
            I2_all_powers[pind, find, :], Q2_all_powers[pind, find, :] = I_final2, Q_final2

            cav_sweep_mag1[pind, find, :] = np.sqrt(I_final**2 + Q_final**2)
            cav_sweep_phase1[pind, find, :] = np.arctan2(Q_final, I_final)*180/np.pi

            cav_sweep_mag2[pind, find, :] = np.sqrt(I_final2**2 + Q_final2**2)
            cav_sweep_phase2[pind, find, :] = np.arctan2(Q_final2, I_final2)*180/np.pi


            Is_all_powers[pind, find, :], Qs_all_powers[pind, find, :] = I_full[0], Q_full[0]
            Is_all_powers_b[pind, find, :], Qs_all_powers_b[pind, find, :] = I_full_b[0], Q_full_b[0]

        time.sleep(1)

    ## Plots
    ##########################
    if freq_points > 1:
        # fig = plt.figure(1113, figsize=(13,8))
        fig = plt.figure(1113)
        plt.clf()
        # ax = plt.subplot(1,1,1)
        xaxis = cav_freqs/1e9
        plt.plot(xaxis, cav_sweep_mag1[0,:,0], label='meas1')
        plt.plot(xaxis, cav_sweep_mag2[0,:,0], label='meas2')
        plt.xlabel('Freq (GHz)')
        plt.title('Cavity Sweep')
        plt.legend()
        plt.tight_layout()
        mplcursors.cursor(multiple=True)
        plt.savefig(os.path.join(saveDir, filename+'CavitySweep.png'), dpi=150)
        plt.show()

    else:
        if exp_settings['with_stimulation'] == True:
            stim_pulse = 'Y'
            if exp_settings['stimulation_on_resonance'] == True: 
                mod = 0
            else:
                mod = mod_freq
        else:
            stim_pulse = 'N'
            mod = mod_freq

    #################3########
        col = 3
        if len(stim_powers) > col:
            row = 2
        else:
            row = 1
        
        fig = plt.figure(1111, figsize=(16,10))
        plt.clf()
        for ind in range(len(stim_powers)):
            if ind <=5:
                ax = plt.subplot(row, col, ind+1)
                ax.scatter(I1_all_powers[ind, 0, :], Q1_all_powers[ind, 0, :], s=4, label = 'g')
                plt.xlabel('I')
                plt.ylabel('Q')
                plt.title(f'{stim_powers[ind]-CAV_Attenuation} dBm')
                plt.legend()
                plt.axhline(y = 0, color ="black", linestyle ="-")
                plt.axvline(x = 0, color ="black", linestyle ="-") 
                y_max = np.max(np.abs(ax.get_ylim()))
                x_max = np.max(np.abs(ax.get_xlim()))
                list = np.array([x_max, y_max])
                max = np.max(list)
                ax.set_ylim(ymin=-max, ymax=max)
                ax.set_xlim(xmin=-max, xmax=max)
                ax.set_aspect('equal')
        
        fig.suptitle(f'{segments} points | Readout:{round(CAV_freq/1e9,3)}GHz, {CAV_power - CAV_Attenuation}dBm | Stimulation:{stim_pulse} | $\Delta$ = {round((mod)/1e6,2)}MHz')
        fig.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_QNDPlot.png'), dpi = 150)
        plt.show()

        # plot above with imshow, cmap = hot
        bins = 100
        fig = plt.figure(1112, figsize=(16,10))
        plt.clf()
        for ind in range(len(stim_powers)):
            if ind <=5:
                ax = plt.subplot(row, col, ind+1)
                I_array, Q_array = I1_all_powers[ind, 0, :], Q1_all_powers[ind, 0, :]
                heatmap, xedges, yedges = np.histogram2d(I_array, Q_array, bins=bins)
                max_count = np.max(heatmap)
                heatmap_scaled = heatmap/max_count
                plt.imshow(heatmap_scaled.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        origin='lower', cmap='hot', vmin=0.05, vmax=0.5)
                plt.colorbar(label=fr'Counts $\times$ {int(max_count)}', fraction=0.04, pad=0.04)
                plt.xlabel('I')
                plt.ylabel('Q')
                plt.title(f'{stim_powers[ind]-CAV_Attenuation} dBm')
                plt.axhline(y = 0, color ="black", linestyle ="-")
                plt.axvline(x = 0, color ="black", linestyle ="-") 
                y_max = np.max(np.abs(ax.get_ylim()))
                x_max = np.max(np.abs(ax.get_xlim()))
                list = np.array([x_max, y_max])
                max = np.max(list)
                ax.set_ylim(ymin=-max, ymax=max)
                ax.set_xlim(xmin=-max, xmax=max)
                plt.gca().set_facecolor('black')
        fig.suptitle(f'{segments} points | Readout:{round(CAV_freq/1e9,3)}GHz, {CAV_power - CAV_Attenuation}dBm | Stimulation:{stim_pulse} | $\Delta$ = {round((mod)/1e6,2)}MHz')
        fig.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_QNDHotPlot.png'), dpi = 150)
        plt.show()
                

    # data to save
    full_data = {}
    full_data['I1_all_powers'] = I1_all_powers
    full_data['Q1_all_powers'] = Q1_all_powers
    full_data['I2_all_powers'] = I2_all_powers
    full_data['Q2_all_powers'] = Q2_all_powers
    full_data['Is_full_all_powers'] = Is_all_powers
    full_data['Qs_full_all_powers'] = Qs_all_powers
    full_data['Is_background'] = Is_all_powers_b
    full_data['Qs_background'] = Qs_all_powers_b

    userfuncs.SaveFull(saveDir, filename, ['stim_powers','cav_freqs', 'xaxis',
                                                       'full_data'], 
                                                     locals(), 
                                                     expsettings=settings, 
                                                     instruments=instruments, saveHWsettings=True)    


    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
       
    cavitygen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.output = 'Off'
        
    
    
    return full_data