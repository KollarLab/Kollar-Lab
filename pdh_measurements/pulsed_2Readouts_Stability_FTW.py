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
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process_2Readouts
from pdh_measurements.scheduler_pdh import scheduler_pdh

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'PhaseStability'
    
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
def generate_freq_array(start_freq, stop_freq, freq_step, num_points):
    if num_points == 1:
        freq_array = [start_freq]
    else:
        total_span = stop_freq - start_freq
        N_float = total_span/((num_points - 1)*freq_step)
        N_int = int(N_float)
        actual_step = N_int * freq_step
        freq_array = [start_freq + i*actual_step for i in range(num_points)]
    return np.array(freq_array)

def pulsed_2Readouts_Stability_FTW(instruments, settings):
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
    freq_tuning = exp_globals['freq_tuning'] 
    start_freq  = np.floor(exp_settings['start_freq']/freq_tuning) * freq_tuning
    stop_freq   = np.floor(exp_settings['stop_freq']/freq_tuning)*freq_tuning
    freq_points = exp_settings['freq_points']
    freqs = generate_freq_array(start_freq, stop_freq, freq_tuning, freq_points)


    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power  = exp_settings['Qbit_power'] + Qbit_Attenuation
    Qbit_freq = exp_settings['Qbit_freq']
    avgs = exp_settings['averages']
    
    
    # Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    gen_power  = exp_settings['gen_power'] + CAV_Attenuation
    CAV_freq = exp_settings['CAV_freq']

    
    num_points = int(exp_settings['stability_points'])  # number of stability data points
    stab_time_in_hr = exp_settings['stab_time_in_hr']
    wait_time_in_min = int(stab_time_in_hr * 60 / num_points)
    
           
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
    
    # ##create the digital down conversion filter if needed.
    # if exp_globals['IF'] != 0:
    #     #create Chebychev type II digital filter
    #     filter_N = exp_globals['ddc_config']['order']
    #     filter_rs = exp_globals['ddc_config']['stop_atten']
    #     filter_cutoff = np.abs(exp_globals['ddc_config']['cutoff'])
    #     LPF = signal.cheby2(filter_N, filter_rs, filter_cutoff, btype='low', analog=False, output='sos', fs=card.sampleRate)
    #     settings['LPF'] = LPF
        
    #     xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
        
    #     # carrier
    #     digLO_sin_carr = np.sin(2*np.pi*exp_globals['IF']*xaxis)
    #     digLO_cos_carr = np.cos(2*np.pi*exp_globals['IF']*xaxis)
    #     settings['digLO_sin_carr'] = digLO_sin_carr 
    #     settings['digLO_cos_carr'] = digLO_cos_carr
        
    #     # lower sideband
    #     digLO_sin_sb = np.sin(2*np.pi*(-1)*(mod_freq-exp_globals['IF'])*xaxis )
    #     digLO_cos_sb = np.cos(2*np.pi*(-1)*(mod_freq-exp_globals['IF'])*xaxis )
    #     settings['digLO_sin_sb'] =  digLO_sin_sb 
    #     settings['digLO_cos_sb'] = digLO_cos_sb 
        

    ## Start main measurement loop 
    I_carr = np.zeros((num_points, len(freqs)))
    Q_carr = np.zeros((num_points, len(freqs)))
    powerdat_carr = np.zeros((num_points, len(freqs)))
    phasedat_carr = np.zeros((num_points, len(freqs)))
    
    I_sb = np.zeros((num_points, len(freqs)))
    Q_sb = np.zeros((num_points, len(freqs)))
    powerdat_sb = np.zeros((num_points, len(freqs)))
    phasedat_sb = np.zeros((num_points, len(freqs)))
    
    phasedat_diff = np.zeros((num_points, len(freqs)))
    ts = np.zeros(num_points)
    

    tstart = time.time()
    first_it = True
    
    qubitgen.power = Qbit_power
    qubitgen.freq = Qbit_freq
    cavitygen.power = gen_power

    print('asdf')
    print(f'SGS_power:{gen_power - CAV_Attenuation}, Qbit_power:{Qbit_power - Qbit_Attenuation}')
    
    for numind in range(num_points):
        
        total_samples = card.samples 
        Is_carr  = np.zeros((len(freqs), total_samples))
        Qs_carr  = np.zeros((len(freqs), total_samples))
        
        Is_sb  = np.zeros((len(freqs), total_samples))
        Qs_sb  = np.zeros((len(freqs), total_samples))
        
        # for background
        Is_carr_b = np.zeros((len(freqs), total_samples))
        Qs_carr_b = np.zeros((len(freqs), total_samples))
        Is_sb_b = np.zeros((len(freqs), total_samples))
        Qs_sb_b = np.zeros((len(freqs), total_samples))
        
        # print(f'Scan run:{numind+1}, SGS_power:{gen_power - CAV_Attenuation}, Qbit_power:{Qbit_power - Qbit_Attenuation}')
        print(f'run {numind+1} of {num_points}')
    
        for find in range(0, len(freqs)):
                
            freq = freqs[find]
            
            cavitygen.freq = freq  
            cavitygen.phase = 0
            cavitygen.output = 'On'
            LO.freq = freq - exp_globals['IF']
            LO.phase = 0
            LO.output = 'On'
            if exp_settings['Qstate'] == 'Ground':
                qubitgen.output = 'Off'
            elif exp_settings['Qstate'] == 'Excited':
                qubitgen.output = 'On'
            else:
                print('Invalid Qubit state')

                ##create the digital down conversion filter if needed.
            if exp_globals['IF'] != 0:
                ftw_IF = freq - LO.freq
                #create Chebychev type II digital filter
                filter_N = exp_globals['ddc_config']['order']
                filter_rs = exp_globals['ddc_config']['stop_atten']
                filter_cutoff = np.abs(exp_globals['ddc_config']['cutoff'])
                LPF = signal.cheby2(filter_N, filter_rs, filter_cutoff, btype='low', analog=False, output='sos', fs=card.sampleRate)
                settings['LPF'] = LPF
                
                xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
                
                # carrier
                digLO_sin_carr = np.sin(2*np.pi*ftw_IF*xaxis)
                digLO_cos_carr = np.cos(2*np.pi*ftw_IF*xaxis)
                settings['digLO_sin_carr'] = digLO_sin_carr 
                settings['digLO_cos_carr'] = digLO_cos_carr
                
                # lower sideband
                digLO_sin_sb = np.sin(2*np.pi*(-1)*(mod_freq-ftw_IF)*xaxis )
                digLO_cos_sb = np.cos(2*np.pi*(-1)*(mod_freq-ftw_IF)*xaxis )
                settings['digLO_sin_sb'] =  digLO_sin_sb 
                settings['digLO_cos_sb'] = digLO_cos_sb 
    
            carr_data, sb_data, xaxis = read_and_process_2Readouts(card, settings, 
                                                                         plot=first_it, 
                                                                         IQstorage = True)
            
            if exp_settings['subtract_background']:
                #Acquire background trace
                qubitgen.output='Off'
                carr_data_b, sb_data_b, xaxis_b = read_and_process_2Readouts(card, settings, 
                                                                             plot=first_it, 
                                                                             IQstorage = True)
            else:
                carr_data_b = {key: 0 for key in carr_data}
                sb_data_b = {key: 0 for key in sb_data}
                           
            I_final_carr = np.mean(carr_data['I_window'])  - np.mean(carr_data_b['I_window']) #compute <I> in the data window
            Q_final_carr = np.mean(carr_data['Q_window'])  - np.mean(carr_data_b['Q_window']) #compute <Q> in the data window
            
            I_final_sb = np.mean(sb_data['I_window']) - np.mean(sb_data_b['I_window']) #compute <I> in the data window
            Q_final_sb = np.mean(sb_data['Q_window']) - np.mean(sb_data_b['Q_window']) #compute <Q> in the data window
            
            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, num_points*len(freqs))
                first_it = False
                
            Is_carr[find,:] = carr_data['I_full']
            Qs_carr[find,:] = carr_data['Q_full']
            
            Is_sb[find,:] = sb_data['I_full']
            Qs_sb[find,:] = sb_data['Q_full']
            
            Is_carr_b[find,:] = carr_data_b['I_full']
            Qs_carr_b[find,:] = carr_data_b['Q_full']
            
            Is_sb_b[find,:] = sb_data_b['I_full']
            Qs_sb_b[find,:] = sb_data_b['Q_full']
            
            I_carr[numind, find] = I_final_carr
            Q_carr[numind, find] = Q_final_carr
            powerdat_carr[numind, find] = np.sqrt(I_final_carr**2 + Q_final_carr**2)
            phasedat_carr[numind, find] = np.arctan2(Q_final_carr, I_final_carr)*180/np.pi
            
            I_sb[numind, find] = I_final_sb
            Q_sb[numind, find] = Q_final_sb
            powerdat_sb[numind, find] = np.sqrt(I_final_sb**2 + Q_final_sb**2)
            phasedat_sb[numind, find] = np.arctan2(Q_final_sb, I_final_sb)*180/np.pi
            
            phasedat_diff[numind, find] = phasedat_carr[numind, find] - phasedat_sb[numind, find]
            currT = time.time() - tstart
            ts[numind] = currT
        
        # wait for x minutes before the next point
        if num_points > 1:
            # time.sleep(5*60)
            time.sleep(wait_time_in_min * 60)
        # For saving
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
        full_data['ts'] = ts
        
        full_time = {}
        full_time['xaxis']  = xaxis*1e6
        full_time['Is_carr']   = Is_carr
        full_time['Qs_carr'] = Qs_carr
        full_time['Is_sb']   = Is_sb
        full_time['Qs_sb'] = Qs_sb
        full_time['Is_carr_b']   = Is_carr_b
        full_time['Qs_carr_b'] = Qs_carr_b
        full_time['Is_sb_b']   = Is_sb_b
        full_time['Qs_sb_b'] = Qs_sb_b
          
        Qstate = exp_settings['Qstate']
        if num_points > 1:
            fig = plt.figure(1213, figsize=(13,8))
            plt.clf()
            ax = plt.subplot(2,2,1)
            plt.plot(ts, phasedat_carr[:,0])
            plt.xlabel('Time (s)')
            #plt.ylim(-180,180)
            plt.title('Carrier phase')
            
            ax = plt.subplot(2,2,2)
            plt.plot(ts, phasedat_sb[:,0])
            plt.xlabel('Time (s)')
            #plt.ylim(-180,180)
            plt.title('Sideband phase')
            
            ax = plt.subplot(2,2,3)
            plt.plot(ts, phasedat_diff[:,0])
            plt.xlabel('Time (s)')
            #plt.ylim(-180,180)
            plt.title('Difference phase')
            
            # fig.suptitle(f'{num_points} points phase stability')
            fig.suptitle(f'{num_points} points | avgs:{avgs} | Carrier: {carr_power}dBm | Sideband: {lsb_power}dBm | ModFreq: {mod_freq/1e6}MHz | Qstate: {Qstate}')
            fig.tight_layout()
            mplcursors.cursor(multiple=True)
            plt.savefig(os.path.join(saveDir, filename+'_PhaseStability.png'), dpi = 150)
            plt.show()
            
            fig = plt.figure(1214, figsize=(13,8))
            plt.clf()
            ax = plt.subplot(2,2,1)
            ax.scatter(I_carr, Q_carr, s=4)
            plt.xlabel('I')
            plt.ylabel('Q')
            plt.title('Carrier IQ')
            plt.axhline(y = 0, color ="black", linestyle ="-")
            plt.axvline(x = 0, color ="black", linestyle ="-") 
            y_max = np.max(np.abs(ax.get_ylim()))
            x_max = np.max(np.abs(ax.get_xlim()))
            list = np.array([x_max, y_max])
            max = np.max(list)
            ax.set_ylim(ymin=-max, ymax=max)
            ax.set_xlim(xmin=-max, xmax=max)
            
            ax = plt.subplot(2,2,2)
            ax.scatter(I_sb, Q_sb, s=4)
            plt.xlabel('I')
            plt.ylabel('Q')
            plt.title('sideband IQ')
            plt.axhline(y = 0, color ="black", linestyle ="-")
            plt.axvline(x = 0, color ="black", linestyle ="-") 
            y_max = np.max(np.abs(ax.get_ylim()))
            x_max = np.max(np.abs(ax.get_xlim()))
            list = np.array([x_max, y_max])
            max = np.max(list)
            ax.set_ylim(ymin=-max, ymax=max)
            ax.set_xlim(xmin=-max, xmax=max)
            
            ax = plt.subplot(2,2,3)
            ax.scatter(powerdat_carr[:,0]*np.cos(phasedat_carr[:,0]*np.pi/180), powerdat_carr[:,0]*np.sin(phasedat_carr[:,0]*np.pi/180), s=4)
            plt.xlabel('I')
            plt.ylabel('Q')
            plt.title('Carrier IQ reconstruct')
            plt.axhline(y = 0, color ="black", linestyle ="-")
            plt.axvline(x = 0, color ="black", linestyle ="-") 
            y_max = np.max(np.abs(ax.get_ylim()))
            x_max = np.max(np.abs(ax.get_xlim()))
            list = np.array([x_max, y_max])
            max = np.max(list)
            ax.set_ylim(ymin=-max, ymax=max)
            ax.set_xlim(xmin=-max, xmax=max)
            
            ax = plt.subplot(2,2,4)
            ax.scatter(powerdat_carr[:,0]*np.cos(phasedat_diff[:,0]*np.pi/180), powerdat_carr[:,0]*np.sin(phasedat_diff[:,0]*np.pi/180), s=4)
            plt.xlabel('I')
            plt.ylabel('Q')
            plt.title('Difference IQ')
            plt.axhline(y = 0, color ="black", linestyle ="-")
            plt.axvline(x = 0, color ="black", linestyle ="-") 
            y_max = np.max(np.abs(ax.get_ylim()))
            x_max = np.max(np.abs(ax.get_xlim()))
            list = np.array([x_max, y_max])
            max = np.max(list)
            ax.set_ylim(ymin=-max, ymax=max)
            ax.set_xlim(xmin=-max, xmax=max)
            
            # fig.suptitle(f'{num_points} points IQ stability')
            fig.suptitle(f'{num_points} points | avgs:{avgs} | Carrier: {carr_power}dBm | Sideband: {lsb_power}dBm | ModFreq: {mod_freq/1e6}MHz | Qstate: {Qstate}')
            fig.tight_layout()
            mplcursors.cursor(multiple=True)
            plt.savefig(os.path.join(saveDir, filename+'_IQStability.png'), dpi = 150)
            plt.show()
            
            
        else:
        
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
            mplcursors.cursor(multiple=True)
            plt.savefig(os.path.join(saveDir, filename+'_fullMagPlot.png'), dpi = 150)
            plt.show()
    
            identifier = 'Power: {}dBm'.format(gen_power-CAV_Attenuation)
            new_filename = 'Raw_time_traces\n'+filename
            
            # Plot raw time traces
            fig = plt.figure(12155, figsize=(13,8))
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
            mplcursors.cursor(multiple=True)
            plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces.png'), dpi = 150)
            plt.show()
    
    
        userfuncs.SaveFull(saveDir, filename, ['freqs', 'IQdat','full_data', 'full_time'],
                             locals(), expsettings=settings, instruments=instruments, saveHWsettings=first_it)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))

    LO.output = 0
    cavitygen.output = 0
    qubitgen.output = 0
    
    return full_data