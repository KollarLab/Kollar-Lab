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
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process, read_and_process_singleshot
from utility.scheduler import scheduler
from pdh_measurements.scheduler_pdh import scheduler_pdh

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'QND_normalReadout'
    
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


def QND_normalReadout(instruments, settings):
    
        ##Instruments used
    stimulationgen  = instruments['stimulationgen']
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
    
    ## stimulation settings
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
    cavitygen.Output = 'On'
    stimulationgen.Output  = 'On'
    
    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.disable_IQ()
    cavitygen.enable_pulse()
    cavitygen.output = 'On'
    
    ## Setting qubit generator to some safe starting point before we turn it on
    stimulationgen.Freq   = 4e9
    stimulationgen.Power  = -20
    
    if exp_settings['Quasi_CW']:
        stimulationgen.disable_pulse()
        stimulationgen.disable_IQ()
    else:
        stimulationgen.enable_pulse()
        stimulationgen.enable_IQ()
    
    ## Card config
    configure_card(card, settings)

    ## HDAWG settings
#    configure_hdawg(hdawg, settings)
    
    ## Sequencer program
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    q_pulse = exp_globals['qubit_pulse']
    
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    
    awg_sched = scheduler(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='stimulation_I')
    awg_sched.add_analog_channel(2, name='stimulation_Q')
    
    awg_sched.add_digital_channel(1, name='stimulation_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    stimulation_I       = awg_sched.analog_channels['stimulation_I']
    stimulation_marker  = awg_sched.digital_channels['stimulation_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    stimulation_pulse_length = exp_settings['stimulation_pulse_length']
    ring_down_time = exp_settings['ring_down_time']
    init_buffer = m_pulse['init_buffer']

    position = start_time-ring_down_time-delay-num_sigma*sigma-stimulation_pulse_length - init_buffer
    stimulation_I.add_pulse('gaussian_square', position=position, 
                              amplitude=q_pulse['piAmp'], length = stimulation_pulse_length, 
                              ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    
    total_pulse_length = stimulation_pulse_length + (num_sigma*sigma)
    marker_delay = 160e-9
    # stimulation_marker.add_window(position-160e-9, position+2*160e-9+stimulation_pulse_length+(num_sigma*sigma))
    stimulation_marker.add_window(position-marker_delay, position+marker_delay+total_pulse_length)

    ##
    cavity_marker.add_window(start_time, start_time+window_time)
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['stimulation_I', 'stimulation_Q'], ['stimulation_enable', 'Cavity_enable'])
    
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
    
    # powerdat = np.zeros((len(powers), len(freqs)))
    # phasedat = np.zeros((len(powers), len(freqs)))
    
    first_it = True
    segments = int(card.segments)
    total_samples = int(card.samples)

    I_all_powers = np.zeros((len(powers), len(freqs), segments))
    Q_all_powers = np.zeros((len(powers), len(freqs), segments))

    mags_all_powers = np.zeros((len(powers), len(freqs), segments))

    Is_full_all_powers, Qs_full_all_powers = np.zeros((len(powers), len(freqs), total_samples)), np.zeros((len(powers), len(freqs), total_samples))
    Is_back_all_powers, Qs_back_all_powers = np.zeros((len(powers), len(freqs), total_samples)), np.zeros((len(powers), len(freqs), total_samples))
    for powerind in range(len(powers)):
        stimulationgen.Power = powers[powerind]

        # I_carr, Q_carr = np.zeros((len(freqs), segments)), np.zeros((len(freqs), segments))

        # Is_full  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        # Qs_full  = np.zeros((len(freqs), total_samples))
        
        # Is_back = np.zeros(Is_full.shape)
        # Qs_back = np.zeros(Qs_full.shape)
        
        print('Current power:{}, max:{}'.format(powers[powerind]-Qbit_Attenuation, powers[-1]-Qbit_Attenuation))
    
        for find in range(0, len(freqs)):
            
            if first_it:
                tstart = time.time()
            
            ##Acquire signal
            freq = freqs[find]
            stimulationgen.Freq = freq
            stimulationgen.output='On'
            cavitygen.output = 'On'
            LO.output = 'On'
            stimulationgen.phase = 0
            cavitygen.phase = 0
            LO.phase = 0

            if exp_settings['sweep_cavity'] == True:
                cavitygen.freq = freq  # this is only needed to see if we can get trans data with this code, and show the code works as expected
                stimulationgen.output='Off'
                LO.freq = freq - exp_globals['IF']  
            else:
                cavitygen.freq = CAV_freq
                LO.freq = CAV_freq - exp_globals['IF']  
#            time.sleep(0.1)

            I_window, Q_window, I_full, Q_full, xaxis = read_and_process_singleshot(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            if exp_settings['subtract_background']:
                #Acquire background trace
                stimulationgen.output='Off'
#                time.sleep(0.1)
                I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process_singleshot(card, settings, 
                                                                 plot=first_it, 
                                                                 IQstorage = True)
            else:
                I_window_b, Q_window_b, I_full_b, Q_full_b = np.zeros(I_window.shape), np.zeros(Q_window.shape), np.zeros(I_full.shape), np.zeros(Q_full.shape)
            
            if first_it:
                first_it=False

            ##Useful handles for variables
            I_final = np.mean(I_window, axis=1, keepdims=True) - np.mean(I_window_b, axis=1, keepdims=True)
            Q_final = np.mean(Q_window, axis=1, keepdims=True) - np.mean(Q_window_b, axis=1, keepdims=True)

            I_final, Q_final = I_final.flatten(), Q_final.flatten()
            # I_carr[find, :] = I_final
            # Q_carr[find, :] = Q_final

            I_all_powers[powerind, find, :], Q_all_powers[powerind, find, :] = I_final, Q_final
            mags_all_powers[powerind, find, :] = np.sqrt(I_final**2 + Q_final**2)
            Is_full_all_powers[powerind, find, :], Qs_full_all_powers[powerind, find, :] = I_full[0], Q_full[0]
            Is_back_all_powers[powerind, find, :], Qs_back_all_powers[powerind, find, :] = I_full_b[0], Q_full_b[0]

            # Is_full[find, :], Qs_full[find, :] = I_full[0], Q_full[0]
            # Is_back[find, :], Qs_back[find, :] = I_full_b[0], Q_full_b[0]

        # I_all_powers[powerind, :, :], Q_all_powers[powerind, :, :] = I_carr.flatten(), Q_carr.flatten()
        # Is_all_powers[powerind, :], Qs_all_powers[powerind, :] = Is_full, Qs_full

    ## Plot IQ data points

    ##########################
    if freq_points > 1:
        fig = plt.figure(1113, figsize=(13,8))
        plt.clf()
        ax = plt.subplot(1,1,1)
        xaxis = freqs/1e9
        yaxis = np.sqrt(I_all_powers[0,:,0]**2 + Q_all_powers[0,:,0]**2)
        # plt.plot(xaxis, yaxis)
        plt.plot(xaxis, mags_all_powers[0,:,0])
        plt.xlabel('Freq (GHz)')
        plt.title('Cavity Sweep')
        plt.tight_layout()
        mplcursors.cursor(multiple=True)
        plt.savefig(os.path.join(saveDir, filename+'Sweep.png'), dpi=150)
        plt.show()


    else:
    #################3########
        col = 3
        if len(powers) > col:
            row = 2
        else:
            row = 1
        
        fig = plt.figure(1111, figsize=(16,10))
        plt.clf()
        for ind in range(len(powers)):
            if ind <=5:
                ax = plt.subplot(row, col, ind+1)
                ax.scatter(I_all_powers[ind, 0, :], Q_all_powers[ind, 0, :], s=4, label = 'g')
                plt.xlabel('I')
                plt.ylabel('Q')
                plt.title(f'{powers[ind]-Qbit_Attenuation} dBm')
                plt.legend()
                plt.axhline(y = 0, color ="black", linestyle ="-")
                plt.axvline(x = 0, color ="black", linestyle ="-") 
                y_max = np.max(np.abs(ax.get_ylim()))
                x_max = np.max(np.abs(ax.get_xlim()))
                list = np.array([x_max, y_max])
                max = np.max(list)
                ax.set_ylim(ymin=-max, ymax=max)
                ax.set_xlim(xmin=-max, xmax=max)
        
        fig.suptitle(f'{segments} points | Readout:{round(CAV_freq/1e9,3)}GHz, {CAV_power - CAV_Attenuation}dBm | Stimulation:{round(start_freq/1e9,3)}GHz | Ringdown:{round(ring_down_time/1e-6,2)}$\mu$s | $\Delta$ = {round((CAV_freq - start_freq)/1e6,2)}MHz')
        fig.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_QNDPlot.png'), dpi = 150)
        plt.show()

        # plot above with imshow, cmap = hot
        bins = 100
        fig = plt.figure(1112, figsize=(16,10))
        plt.clf()
        for ind in range(len(powers)):
            if ind <=5:
                ax = plt.subplot(row, col, ind+1)
                I_array, Q_array = I_all_powers[ind, 0, :], Q_all_powers[ind, 0, :]
                heatmap, xedges, yedges = np.histogram2d(I_array, Q_array, bins=bins)
                max_count = np.max(heatmap)
                heatmap_scaled = heatmap/max_count
                plt.imshow(heatmap_scaled.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        origin='lower', cmap='hot')
                plt.colorbar(label=fr'Counts $\times$ {int(max_count)}', fraction=0.04, pad=0.04)
                plt.xlabel('I')
                plt.ylabel('Q')
                plt.title(f'{powers[ind]-Qbit_Attenuation} dBm')
                plt.axhline(y = 0, color ="black", linestyle ="-")
                plt.axvline(x = 0, color ="black", linestyle ="-") 
                y_max = np.max(np.abs(ax.get_ylim()))
                x_max = np.max(np.abs(ax.get_xlim()))
                list = np.array([x_max, y_max])
                max = np.max(list)
                ax.set_ylim(ymin=-max, ymax=max)
                ax.set_xlim(xmin=-max, xmax=max)
                plt.gca().set_facecolor('black')
        fig.suptitle(f'{segments} points | Readout:{round(CAV_freq/1e9,3)}GHz, {CAV_power - CAV_Attenuation}dBm | Stimulation:{round(start_freq/1e9,3)}GHz | Ringdown:{round(ring_down_time/1e-6,2)}$\mu$s | $\Delta$= {round((CAV_freq - start_freq)/1e6,2)}MHz')
        fig.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_QNDHotPlot.png'), dpi = 150)
        plt.show()
                

    # data to save
    full_data = {}
    full_data['I_all_powers'] = I_all_powers
    full_data['Q_all_powers'] = Q_all_powers
    full_data['Is_full_all_powers'] = Is_full_all_powers
    full_data['Qs_full_all_powers'] = Qs_full_all_powers
    full_data['Is_background'] = Is_back_all_powers
    full_data['Qs_background'] = Qs_back_all_powers

    userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'xaxis',
                                                       'full_data'], 
                                                     locals(), 
                                                     expsettings=settings, 
                                                     instruments=instruments, saveHWsettings=True)    


    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
       
    cavitygen.Output = 'Off'
    stimulationgen.Output = 'Off'
    LO.output = 'Off'
        
            
    #     ##Packaging data for the plotting functions and saving 
    #     full_data = {}
    #     full_data['xaxis'] = freqs/1e9
    #     full_data['mags'] = powerdat[0:powerind+1]
    #     full_data['phases'] = phasedat[0:powerind+1]

    #     single_data = {}
    #     single_data['xaxis'] = freqs/1e9
    #     single_data['mag'] = powerdat[powerind, :]
    #     single_data['phase'] = phasedat[powerind, :]

    #     yaxis = powers[0:powerind+1] - Qbit_Attenuation
    #     labels = ['Freq (GHz)', 'Power (dBm)']
    #     simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1)
    #     if not powerind%exp_settings['num_save']:
    #         plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

    #     full_time = {}
    #     full_time['xaxis']  = xaxis*1e6
    #     full_time['Is'] = Is_full
    #     full_time['Qs'] = Qs_full
    #     full_time['Ib'] = Is_back
    #     full_time['Qb'] = Qs_back

    #     single_time = {}
    #     single_time['xaxis'] = xaxis*1e6
    #     single_time['I'] = I_full
    #     single_time['Q'] = Q_full

    #     time_labels = ['Time (us)', 'Freq (GHz)']
    #     identifier = 'Power: {}dBm'.format(powers[powerind]-Qbit_Attenuation)

    #     simplescan_plot(full_time, single_time, freqs/1e9, 
    #                     'Raw_time_traces\n'+filename, 
    #                     time_labels, 
    #                     identifier, 
    #                     fig_num=2,
    #                     IQdata = True)
    #     if not powerind%exp_settings['num_save']:
    #         plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)
        
    #     if not powerind%exp_settings['num_save']:
    #         userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'xaxis',
    #                                                'powerdat', 'phasedat',
    #                                                'full_data', 'single_data'], 
    #                                                #'full_time', 'single_time'],
    #                                              locals(), 
    #                                              expsettings=settings, 
    #                                              instruments=instruments, saveHWsettings=False)
    # userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'xaxis',
    #                                                    'powerdat', 'phasedat',
    #                                                    'full_data', 'single_data'], 
    #                                                    #'full_time', 'single_time'],
    #                                                  locals(), 
    #                                                  expsettings=settings, 
    #                                                  instruments=instruments, saveHWsettings=True)
    
    
    return full_data

def QND_normalReadout2(instruments, settings):
    
        ##Instruments used
    stimulationgen  = instruments['stimulationgen']
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
    
    ## stimulation settings
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
    cavitygen.Output = 'On'
    stimulationgen.Output  = 'On'
    
    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_IQ()
    cavitygen.enable_pulse()
    cavitygen.output = 'On'
    
    ## Setting qubit generator to some safe starting point before we turn it on
    stimulationgen.Freq   = 4e9
    stimulationgen.Power  = -20
    
    if exp_settings['Quasi_CW']:
        stimulationgen.disable_pulse()
        stimulationgen.disable_IQ()
    else:
        stimulationgen.enable_pulse()
        stimulationgen.enable_IQ()
    
    ## Card config
    configure_card(card, settings)

    ## HDAWG settings
#    configure_hdawg(hdawg, settings)
    
    ## Sequencer program
    # progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    q_pulse = exp_globals['qubit_pulse']
    
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']

    awg_sched = scheduler_pdh(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='stimulation_I')
    awg_sched.add_analog_channel(2, name='stimulation_Q')
    awg_sched.add_analog_channel(3, name='Cavity_I')
    awg_sched.add_analog_channel(4, name='Cavity_Q')
    
    awg_sched.add_digital_channel(1, name='stimulation_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(3, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)

    cavity_I = awg_sched.analog_channels['Cavity_I']
    cavity_Q = awg_sched.analog_channels['Cavity_Q'] 
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    position = start_time
    cavity_I.add_pulse('pdh_I', position=position,
                           amp_list = [0,1,0], phase_list = [0,0,0],
                           mod_freq = 0,
                           time = window_time)

    cavity_Q.add_pulse('pdh_Q', position=position,
                           amp_list = [0,1,0], phase_list = [0,0,0],
                           mod_freq = 0,
                           time = window_time)
    cavity_marker.add_window(start_time, start_time+window_time)
    
    stimulation_I       = awg_sched.analog_channels['stimulation_I']
    stimulation_marker  = awg_sched.digital_channels['stimulation_enable']
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    stimulation_pulse_length = exp_settings['stimulation_pulse_length']
    ring_down_time = exp_settings['ring_down_time']
    init_buffer = m_pulse['init_buffer']

    position = start_time-ring_down_time-delay-num_sigma*sigma-stimulation_pulse_length - init_buffer
    stimulation_I.add_pulse('gaussian_square', position=position, 
                              amplitude=q_pulse['piAmp'], length = stimulation_pulse_length, 
                              ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    
    total_pulse_length = stimulation_pulse_length + (num_sigma*sigma)
    marker_delay = 160e-9
    # stimulation_marker.add_window(position-160e-9, position+2*160e-9+stimulation_pulse_length+(num_sigma*sigma))
    stimulation_marker.add_window(position-marker_delay, position+marker_delay+total_pulse_length)

    ##
    cavity_marker.add_window(start_time, start_time+window_time)
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['stimulation_I', 'stimulation_Q'], ['stimulation_enable', 'Cavity_enable'])
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
        
        xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
        digLO_sin = np.sin(2*np.pi*exp_globals['IF']*xaxis)
        digLO_cos = np.cos(2*np.pi*exp_globals['IF']*xaxis)
        
        #store in settings so that the processing functions can get to them
        settings['digLO_sin'] = digLO_sin 
        settings['digLO_cos'] = digLO_cos
        settings['LPF'] = LPF
    
    # powerdat = np.zeros((len(powers), len(freqs)))
    # phasedat = np.zeros((len(powers), len(freqs)))
    
    first_it = True
    segments = int(card.segments)
    total_samples = int(card.samples)

    I_all_powers = np.zeros((len(powers), len(freqs), segments))
    Q_all_powers = np.zeros((len(powers), len(freqs), segments))

    mags_all_powers = np.zeros((len(powers), len(freqs), segments))

    Is_full_all_powers, Qs_full_all_powers = np.zeros((len(powers), len(freqs), total_samples)), np.zeros((len(powers), len(freqs), total_samples))
    Is_back_all_powers, Qs_back_all_powers = np.zeros((len(powers), len(freqs), total_samples)), np.zeros((len(powers), len(freqs), total_samples))
    for powerind in range(len(powers)):
        stimulationgen.Power = powers[powerind]

        # I_carr, Q_carr = np.zeros((len(freqs), segments)), np.zeros((len(freqs), segments))

        # Is_full  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        # Qs_full  = np.zeros((len(freqs), total_samples))
        
        # Is_back = np.zeros(Is_full.shape)
        # Qs_back = np.zeros(Qs_full.shape)
        
        print('Current power:{}, max:{}'.format(powers[powerind]-Qbit_Attenuation, powers[-1]-Qbit_Attenuation))
    
        for find in range(0, len(freqs)):
            
            if first_it:
                tstart = time.time()
            
            ##Acquire signal
            freq = freqs[find]
            stimulationgen.Freq = freq
            stimulationgen.output='On'
            cavitygen.output = 'On'
            LO.output = 'On'
            stimulationgen.phase = 0
            cavitygen.phase = 0
            LO.phase = 0

            if exp_settings['sweep_cavity'] == True:
                cavitygen.freq = freq  # this is only needed to see if we can get trans data with this code, and show the code works as expected
                stimulationgen.output='Off'
                LO.freq = freq - exp_globals['IF']  
            else:
                cavitygen.freq = CAV_freq
                LO.freq = CAV_freq - exp_globals['IF']  
#            time.sleep(0.1)

            I_window, Q_window, I_full, Q_full, xaxis = read_and_process_singleshot(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            if exp_settings['subtract_background']:
                #Acquire background trace
                stimulationgen.output='Off'
#                time.sleep(0.1)
                I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process_singleshot(card, settings, 
                                                                 plot=first_it, 
                                                                 IQstorage = True)
            else:
                I_window_b, Q_window_b, I_full_b, Q_full_b = np.zeros(I_window.shape), np.zeros(Q_window.shape), np.zeros(I_full.shape), np.zeros(Q_full.shape)
            
            if first_it:
                first_it=False

            ##Useful handles for variables
            I_final = np.mean(I_window, axis=1, keepdims=True) - np.mean(I_window_b, axis=1, keepdims=True)
            Q_final = np.mean(Q_window, axis=1, keepdims=True) - np.mean(Q_window_b, axis=1, keepdims=True)

            I_final, Q_final = I_final.flatten(), Q_final.flatten()
            # I_carr[find, :] = I_final
            # Q_carr[find, :] = Q_final

            I_all_powers[powerind, find, :], Q_all_powers[powerind, find, :] = I_final, Q_final
            mags_all_powers[powerind, find, :] = np.sqrt(I_final**2 + Q_final**2)
            Is_full_all_powers[powerind, find, :], Qs_full_all_powers[powerind, find, :] = I_full[0], Q_full[0]
            Is_back_all_powers[powerind, find, :], Qs_back_all_powers[powerind, find, :] = I_full_b[0], Q_full_b[0]

            # Is_full[find, :], Qs_full[find, :] = I_full[0], Q_full[0]
            # Is_back[find, :], Qs_back[find, :] = I_full_b[0], Q_full_b[0]

        # I_all_powers[powerind, :, :], Q_all_powers[powerind, :, :] = I_carr.flatten(), Q_carr.flatten()
        # Is_all_powers[powerind, :], Qs_all_powers[powerind, :] = Is_full, Qs_full

    ## Plot IQ data points

    ##########################
    if freq_points > 1:
        fig = plt.figure(1113, figsize=(13,8))
        plt.clf()
        ax = plt.subplot(1,1,1)
        xaxis = freqs/1e9
        yaxis = np.sqrt(I_all_powers[0,:,0]**2 + Q_all_powers[0,:,0]**2)
        # plt.plot(xaxis, yaxis)
        plt.plot(xaxis, mags_all_powers[0,:,0])
        plt.xlabel('Freq (GHz)')
        plt.title('Cavity Sweep')
        plt.tight_layout()
        mplcursors.cursor(multiple=True)
        plt.savefig(os.path.join(saveDir, filename+'Sweep.png'), dpi=150)
        plt.show()


    else:
    #################3########
        col = 3
        if len(powers) > col:
            row = 2
        else:
            row = 1
        
        fig = plt.figure(1111, figsize=(16,10))
        plt.clf()
        for ind in range(len(powers)):
            if ind <=5:
                ax = plt.subplot(row, col, ind+1)
                ax.scatter(I_all_powers[ind, 0, :], Q_all_powers[ind, 0, :], s=4, label = 'g')
                plt.xlabel('I')
                plt.ylabel('Q')
                plt.title(f'{powers[ind]-Qbit_Attenuation} dBm')
                plt.legend()
                plt.axhline(y = 0, color ="black", linestyle ="-")
                plt.axvline(x = 0, color ="black", linestyle ="-") 
                y_max = np.max(np.abs(ax.get_ylim()))
                x_max = np.max(np.abs(ax.get_xlim()))
                list = np.array([x_max, y_max])
                max = np.max(list)
                ax.set_ylim(ymin=-max, ymax=max)
                ax.set_xlim(xmin=-max, xmax=max)
        
        fig.suptitle(f'{segments} points | Readout:{round(CAV_freq/1e9,3)}GHz, {CAV_power - CAV_Attenuation}dBm | Stimulation:{round(start_freq/1e9,3)}GHz | Ringdown:{round(ring_down_time/1e-6,2)}$\mu$s | $\Delta$ = {round((CAV_freq - start_freq)/1e6,2)}MHz')
        fig.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_QNDPlot.png'), dpi = 150)
        plt.show()

        # plot above with imshow, cmap = hot
        bins = 100
        fig = plt.figure(1112, figsize=(16,10))
        plt.clf()
        for ind in range(len(powers)):
            if ind <=5:
                ax = plt.subplot(row, col, ind+1)
                I_array, Q_array = I_all_powers[ind, 0, :], Q_all_powers[ind, 0, :]
                heatmap, xedges, yedges = np.histogram2d(I_array, Q_array, bins=bins)
                max_count = np.max(heatmap)
                heatmap_scaled = heatmap/max_count
                plt.imshow(heatmap_scaled.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        origin='lower', cmap='hot')
                plt.colorbar(label=fr'Counts $\times$ {int(max_count)}', fraction=0.04, pad=0.04)
                plt.xlabel('I')
                plt.ylabel('Q')
                plt.title(f'{powers[ind]-Qbit_Attenuation} dBm')
                plt.axhline(y = 0, color ="black", linestyle ="-")
                plt.axvline(x = 0, color ="black", linestyle ="-") 
                y_max = np.max(np.abs(ax.get_ylim()))
                x_max = np.max(np.abs(ax.get_xlim()))
                list = np.array([x_max, y_max])
                max = np.max(list)
                ax.set_ylim(ymin=-max, ymax=max)
                ax.set_xlim(xmin=-max, xmax=max)
                plt.gca().set_facecolor('black')
        fig.suptitle(f'{segments} points | Readout:{round(CAV_freq/1e9,3)}GHz, {CAV_power - CAV_Attenuation}dBm | Stimulation:{round(start_freq/1e9,3)}GHz | Ringdown:{round(ring_down_time/1e-6,2)}$\mu$s | $\Delta$= {round((CAV_freq - start_freq)/1e6,2)}MHz')
        fig.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_QNDHotPlot.png'), dpi = 150)
        plt.show()
                

    # data to save
    full_data = {}
    full_data['I_all_powers'] = I_all_powers
    full_data['Q_all_powers'] = Q_all_powers
    full_data['Is_full_all_powers'] = Is_full_all_powers
    full_data['Qs_full_all_powers'] = Qs_full_all_powers
    full_data['Is_background'] = Is_back_all_powers
    full_data['Qs_background'] = Qs_back_all_powers

    userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'xaxis',
                                                       'full_data'], 
                                                     locals(), 
                                                     expsettings=settings, 
                                                     instruments=instruments, saveHWsettings=True)    


    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
       
    cavitygen.Output = 'Off'
    stimulationgen.Output = 'Off'
    LO.output = 'Off'
        
            
    #     ##Packaging data for the plotting functions and saving 
    #     full_data = {}
    #     full_data['xaxis'] = freqs/1e9
    #     full_data['mags'] = powerdat[0:powerind+1]
    #     full_data['phases'] = phasedat[0:powerind+1]

    #     single_data = {}
    #     single_data['xaxis'] = freqs/1e9
    #     single_data['mag'] = powerdat[powerind, :]
    #     single_data['phase'] = phasedat[powerind, :]

    #     yaxis = powers[0:powerind+1] - Qbit_Attenuation
    #     labels = ['Freq (GHz)', 'Power (dBm)']
    #     simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1)
    #     if not powerind%exp_settings['num_save']:
    #         plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

    #     full_time = {}
    #     full_time['xaxis']  = xaxis*1e6
    #     full_time['Is'] = Is_full
    #     full_time['Qs'] = Qs_full
    #     full_time['Ib'] = Is_back
    #     full_time['Qb'] = Qs_back

    #     single_time = {}
    #     single_time['xaxis'] = xaxis*1e6
    #     single_time['I'] = I_full
    #     single_time['Q'] = Q_full

    #     time_labels = ['Time (us)', 'Freq (GHz)']
    #     identifier = 'Power: {}dBm'.format(powers[powerind]-Qbit_Attenuation)

    #     simplescan_plot(full_time, single_time, freqs/1e9, 
    #                     'Raw_time_traces\n'+filename, 
    #                     time_labels, 
    #                     identifier, 
    #                     fig_num=2,
    #                     IQdata = True)
    #     if not powerind%exp_settings['num_save']:
    #         plt.savefig(os.path.join(saveDir, filename+'_Raw_time_traces.png'), dpi = 150)
        
    #     if not powerind%exp_settings['num_save']:
    #         userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'xaxis',
    #                                                'powerdat', 'phasedat',
    #                                                'full_data', 'single_data'], 
    #                                                #'full_time', 'single_time'],
    #                                              locals(), 
    #                                              expsettings=settings, 
    #                                              instruments=instruments, saveHWsettings=False)
    # userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'xaxis',
    #                                                    'powerdat', 'phasedat',
    #                                                    'full_data', 'single_data'], 
    #                                                    #'full_time', 'single_time'],
    #                                                  locals(), 
    #                                                  expsettings=settings, 
    #                                                  instruments=instruments, saveHWsettings=True)
    
    
    return full_data