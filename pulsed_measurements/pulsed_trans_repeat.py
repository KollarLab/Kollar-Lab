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

def pulsed_trans_repeat(instruments, settings):
    
    ##Instruments used
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
    CAV_freqs = np.array([CAV_freq])

    
# =============================================================================
#     Qbit_Attenuation = exp_globals['Qbit_Attenuation']
#     Qbit_power  = exp_settings['Qbit_power'] + Qbit_Attenuation
#     power_points = exp_settings['power_points']
# =============================================================================
    num_points = int(exp_settings['num_points'])
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.Output = 'On'

    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_pulse()
    cavitygen.output = 'On'

    #######################################################################
    # First, do carrier transmission scan
    #######################################################################
    ##Card settings
    configure_card(card, settings)
    
    # The total measurement window ('Data Window Check' plot) is 2*meas_window + emp_delay + init_buffer + post_buffer
    # This might need to be reconfigured if the resonator Q is high enough that reflection is slow and make interference with input meas pulse
    # Or readout marker channel using power splitter might have a time delay.


    # progFile = open(r"C:\Users\kollarlab\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
    # oxford comp sequencer file location is different
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')

    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    
    awg_sched = scheduler(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='Qubit_I')
    awg_sched.add_analog_channel(2, name='Qubit_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=55e-9, HW_offset_off=0e-9)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
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

    ## Start main measurement loop 
    Idata = np.zeros(num_points)
    Qdata = np.zeros(num_points)
    ts = np.zeros(num_points)
    amp_data = np.zeros(num_points) 
    phase_data = np.zeros(num_points) 

    
    tstart = time.time()
    first_it = False#True

# =============================================================================
#     if first_it:
#         tstop = time.time()
# =============================================================================

    for avgind in range(num_points):
        cavitygen.power = CAV_power
        
        time.sleep(0.2)

        total_samples = card.samples 
        Is  = np.zeros((len(CAV_freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs  = np.zeros((len(CAV_freqs), total_samples))
    
        cavitygen.freq = CAV_freq
        LO.freq = CAV_freq - exp_globals['IF']
        LO.output = 'On'
        
        cavitygen.phase = 0
        LO.phase = 0
        time.sleep(0.2)

        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                        plot=first_it, 
                                                                        IQstorage = True)
        
        I_final = np.mean(I_window) #compute <I> in the data window
        Q_final = np.mean(Q_window) #compute <Q> in the data window
        
        if first_it:
            first_it = False
        
        # Store data
        Idata[avgind] = I_final
        Qdata[avgind] = Q_final
        amp_data[avgind] = np.sqrt(I_final**2 + Q_final**2)
        phase_data[avgind] = np.arctan2(Q_final, I_final)*180/np.pi
        currT = time.time() - tstart
        ts[avgind] = currT

    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
# =============================================================================
#     r_ctrans = find_pointsincircle(Iblob_ctrans, Qblob_ctrans, Iblob_ctrans_avg, Qblob_ctrans_avg, attempted_averages)
# =============================================================================
    Idata_avg = np.mean(Idata)   
    Qdata_avg = np.mean(Qdata)      
    trans_data = {}
    trans_data['xaxis'] = CAV_freq/1e9
    trans_data['amp'] = amp_data
    trans_data['phase'] = phase_data
# =============================================================================
#     trans_data['radius'] = r_ctrans
# =============================================================================

    if exp_settings['subtract_groundstate']:
        trans_data['Idata'] = Idata - Idata_avg
        trans_data['Qdata'] = Qdata - Qdata_avg
# =============================================================================
#         ctrans_data['center'] = [Iblob_ctrans_avg - Iblob_ctrans_avg, Qblob_ctrans_avg - Qblob_ctrans_avg]
# =============================================================================
    else:
        trans_data['Idata'] = Idata
        trans_data['Qdata'] = Qdata
# =============================================================================
#         ctrans_data['center'] = [Iblob_ctrans_avg, Qblob_ctrans_avg]
# =============================================================================

    # Plots
    # Need to plot the IQ points and the phase stability
    fig = plt.figure(101)
    plt.clf()
    ax = plt.subplot(2,2,1)
    plt.plot(ts, trans_data['amp'], label = "amp")
    plt.xlabel("Time (s)")
    plt.ylabel('v')
    plt.legend()
    
    ax = plt.subplot(2,2,2)
    plt.plot(ts, trans_data['phase'], label = "phase")
    plt.xlabel("Time (s)")
    plt.ylabel('Deg')
    plt.legend()
    
    ax = plt.subplot(2,2,3)
    ax.scatter(trans_data['Qdata'], trans_data['Idata'], s=4)
    plt.xlabel('Q')
    plt.ylabel('I')
    plt.axhline(y = 0, color ="black", linestyle ="-")
    plt.axvline(x = 0, color ="black", linestyle ="-") 
    y_max = np.max(np.abs(ax.get_ylim()))
    x_max = np.max(np.abs(ax.get_xlim()))
    list = np.array([x_max, y_max])
    max = np.max(list)
    ax.set_ylim(ymin=-max, ymax=max)
    ax.set_xlim(xmin=-max, xmax=max)
    
    fig.suptitle(f'{num_points} points')
    fig.tight_layout()
    plt.show()
    
    save_name = filename + 'TransRepeat'
    userfuncs.savefig(fig, save_name, saveDir, png=True)

    # Save data
    userfuncs.SaveFull(saveDir, filename, ['trans_data'],
                        locals(), 
                        expsettings=settings, 
                        instruments=instruments, saveHWsettings=True)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
       
    cavitygen.Output = 'Off'
    LO.output = 'Off'
    
    return trans_data

def pulsed_spec_repeat(instruments, settings):
    
    ##Instruments used
    cavitygen = instruments['cavitygen']
    qubitgen = instruments['qubitgen']
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
    CAV_freqs = np.array([CAV_freq])

    Qbit_freq  = exp_settings['Qbit_freq']
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power  = exp_settings['Qbit_power'] + Qbit_Attenuation
    power_points = exp_settings['power_points']
    
    
    qubitgen.Output  = 'On'
    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20
    
    if exp_settings['Quasi_CW']:
        qubitgen.disable_pulse()
        qubitgen.disable_IQ()
    else:
        qubitgen.enable_pulse()
        qubitgen.enable_IQ()
    
    
    num_points = int(exp_settings['num_points'])
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.Output = 'On'

    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_pulse()
    cavitygen.output = 'On'

    #######################################################################
    # First, do carrier transmission scan
    #######################################################################
    ##Card settings
    configure_card(card, settings)
    
    # The total measurement window ('Data Window Check' plot) is 2*meas_window + emp_delay + init_buffer + post_buffer
    # This might need to be reconfigured if the resonator Q is high enough that reflection is slow and make interference with input meas pulse
    # Or readout marker channel using power splitter might have a time delay.

    # progFile = open(r"C:\Users\kollarlab\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
    # oxford comp sequencer file location is different
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
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
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=55e-9, HW_offset_off=0e-9)
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
    Idata = np.zeros(num_points)
    Qdata = np.zeros(num_points)
    ts = np.zeros(num_points)
    amp_data = np.zeros(num_points) 
    phase_data = np.zeros(num_points) 

    
    tstart = time.time()
    first_it = False #True

# =============================================================================
#     if first_it:
#         tstop = time.time()
# =============================================================================

    for avgind in range(num_points):
        cavitygen.power = CAV_power
        qubitgen.power = Qbit_power
        time.sleep(0.2)

        total_samples = card.samples 
    
        cavitygen.freq = CAV_freq
        qubitgen.freq = Qbit_freq
        LO.freq = CAV_freq - exp_globals['IF']
        LO.output = 'On'
        
        cavitygen.phase = 0
        qubitgen.phase = 0
        LO.phase = 0
        time.sleep(0.2)

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
        
        I_final = np.mean(I_window) #compute <I> in the data window
        Q_final = np.mean(Q_window) #compute <Q> in the data window
        
        I_background = np.mean(I_window_b)
        Q_background = np.mean(Q_window_b)
        
        if first_it:
            first_it = False
        
        # Store data
        Idata[avgind] = I_final - I_background
        Qdata[avgind] = Q_final - Q_background
        amp_data[avgind] = np.sqrt(I_final**2 + Q_final**2)
        phase_data[avgind] = np.arctan2(Q_final, I_final)*180/np.pi
        currT = time.time() - tstart
        ts[avgind] = currT

    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
# =============================================================================
#     r_ctrans = find_pointsincircle(Iblob_ctrans, Qblob_ctrans, Iblob_ctrans_avg, Qblob_ctrans_avg, attempted_averages)
# =============================================================================
    Idata_avg = np.mean(Idata)   
    Qdata_avg = np.mean(Qdata)      
    trans_data = {}
    trans_data['xaxis'] = CAV_freq/1e9
    trans_data['amp'] = amp_data
    trans_data['phase'] = phase_data
# =============================================================================
#     trans_data['radius'] = r_ctrans
# =============================================================================

    if exp_settings['subtract_groundstate']:
        trans_data['Idata'] = Idata - Idata_avg
        trans_data['Qdata'] = Qdata - Qdata_avg
# =============================================================================
#         ctrans_data['center'] = [Iblob_ctrans_avg - Iblob_ctrans_avg, Qblob_ctrans_avg - Qblob_ctrans_avg]
# =============================================================================
    else:
        trans_data['Idata'] = Idata
        trans_data['Qdata'] = Qdata
# =============================================================================
#         ctrans_data['center'] = [Iblob_ctrans_avg, Qblob_ctrans_avg]
# =============================================================================

    # Plots
    # Need to plot the IQ points and the phase stability
    fig = plt.figure(101)
    plt.clf()
    ax = plt.subplot(2,2,1)
    plt.plot(ts, trans_data['amp'], label = "amp")
    plt.xlabel("Time (s)")
    plt.ylabel('v')
    plt.legend()
    
    ax = plt.subplot(2,2,2)
    plt.plot(ts, trans_data['phase'], label = "phase")
    plt.xlabel("Time (s)")
    plt.ylabel('Deg')
    plt.legend()
    
    ax = plt.subplot(2,2,3)
    ax.scatter(trans_data['Qdata'], trans_data['Idata'], s=4)
    plt.xlabel('Q')
    plt.ylabel('I')
    plt.axhline(y = 0, color ="black", linestyle ="-")
    plt.axvline(x = 0, color ="black", linestyle ="-") 
    y_max = np.max(np.abs(ax.get_ylim()))
    x_max = np.max(np.abs(ax.get_xlim()))
    list = np.array([x_max, y_max])
    max = np.max(list)
    ax.set_ylim(ymin=-max, ymax=max)
    ax.set_xlim(xmin=-max, xmax=max)
    
    fig.suptitle(f'{num_points} points')
    fig.tight_layout()
    plt.show()
    
    save_name = filename + 'SpecRepeat'
    userfuncs.savefig(fig, save_name, saveDir, png=True)

    # Save data
    userfuncs.SaveFull(saveDir, filename, ['trans_data'],
                        locals(), 
                        expsettings=settings, 
                        instruments=instruments, saveHWsettings=True)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
       
    cavitygen.Output = 'Off'
    LO.output = 'Off'
    
    return trans_data

def pulsed_spec_repeatOne(instruments, settings):
    
    ##Instruments used
    cavitygen = instruments['cavitygen']
    qubitgen = instruments['qubitgen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    
    ##Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = exp_settings['CAV_Power'] + CAV_Attenuation
    CAV_freq  = exp_settings['CAV_freq']

    Qbit_freq  = exp_settings['Qbit_freq']
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power  = exp_settings['Qbit_power'] + Qbit_Attenuation
    
    
    qubitgen.Output  = 'On'
    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20
    
    if exp_settings['Quasi_CW']:
        qubitgen.disable_pulse()
        qubitgen.disable_IQ()
    else:
        qubitgen.enable_pulse()
        qubitgen.enable_IQ()

    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.Output = 'On'

    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_pulse()
    cavitygen.output = 'On'

    ##Card settings
    configure_card(card, settings)
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
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
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=55e-9, HW_offset_off=0e-9)
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
    first_it = False #True
    
    cavitygen.power = CAV_power
    qubitgen.power = Qbit_power
    time.sleep(0.2)

    cavitygen.freq = CAV_freq
    qubitgen.freq = Qbit_freq
    LO.freq = CAV_freq - exp_globals['IF']
    LO.output = 'On'
    
    cavitygen.phase = 0
    qubitgen.phase = 0
    LO.phase = 0
    time.sleep(0.2)

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
    
    I_final = np.mean(I_window) #compute <I> in the data window
    Q_final = np.mean(Q_window) #compute <Q> in the data window
    
    I_background = np.mean(I_window_b)
    Q_background = np.mean(Q_window_b)
    
    if first_it:
        first_it = False
    
    # Store data
    Idata = I_final - I_background
    Qdata = Q_final - Q_background
    amp_data = np.sqrt(I_final**2 + Q_final**2)
    phase_data = np.arctan2(Q_final, I_final)*180/np.pi
     
    spec_data = {}
    spec_data['xaxis'] = CAV_freq/1e9
    spec_data['amp'] = amp_data
    spec_data['phase'] = phase_data
    spec_data['Idata'] = Idata
    spec_data['Qdata'] = Qdata
    
    cavitygen.output = 0
    qubitgen.output = 0
    LO.output = 0
    
    
    return spec_data

def pulsed_trans_repeatOne(instruments, settings):
    
    ##Instruments used
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    
    ##Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = exp_settings['CAV_Power'] + CAV_Attenuation
    CAV_freq  = exp_settings['CAV_freq']
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.Output = 'On'

    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_pulse()
    cavitygen.output = 'On'
    ##Card settings
    configure_card(card, settings)
    progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')

    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    
    m_pulse = exp_globals['measurement_pulse']
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    
    awg_sched = scheduler(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='Qubit_I')
    awg_sched.add_analog_channel(2, name='Qubit_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=55e-9, HW_offset_off=0e-9)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
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

    ## Start main measurement loop   
    first_it = False
    #for avgind in range(num_points):
    cavitygen.power = CAV_power
    
    time.sleep(0.2)
    cavitygen.freq = CAV_freq
    LO.freq = CAV_freq - exp_globals['IF']
    LO.output = 'On'
    
    cavitygen.phase = 0
    LO.phase = 0
    time.sleep(0.2)

    I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                    plot=first_it, 
                                                                    IQstorage = True)
    
    I_final = np.mean(I_window) #compute <I> in the data window
    Q_final = np.mean(Q_window) #compute <Q> in the data window
    
    if first_it:
        first_it = False
    
    # Store data
    Idata = I_final
    Qdata = Q_final
    amp_data = np.sqrt(I_final**2 + Q_final**2)
    phase_data = np.arctan2(Q_final, I_final)*180/np.pi
     
    trans_data = {}
    trans_data['xaxis'] = CAV_freq/1e9
    trans_data['amp'] = amp_data
    trans_data['phase'] = phase_data
    trans_data['Idata'] = Idata
    trans_data['Qdata'] = Qdata
    
    cavitygen.output = 0
    LO.output = 0

    return trans_data


def gestate(instruments, settings):
    
    cavitygen = instruments['cavitygen']
    qubitgen = instruments['qubitgen']
    card      = instruments['card']
    LO        = instruments['LO']
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    ##Data saving and naming
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp
    
    num_points = int(exp_settings['num_points'])
    
    Idata_trans = np.zeros(num_points)
    Qdata_trans = np.zeros(num_points)
    amp_data_trans = np.zeros(num_points)
    phase_data_trans = np.zeros(num_points)
    ts_trans = np.zeros(num_points)
    
    Idata_spec = np.zeros(num_points)
    Qdata_spec = np.zeros(num_points)
    amp_data_spec = np.zeros(num_points)
    phase_data_spec = np.zeros(num_points)
    ts_spec = np.zeros(num_points)
    
    tstart = time.time()
    for ind in range(num_points):
        
        trans_data = pulsed_trans_repeatOne(instruments, settings)
        currT_trans = time.time() - tstart
        time.sleep(1)
        
        spec_data = pulsed_spec_repeatOne(instruments, settings)
        currT_spec = time.time() - tstart
        time.sleep(1)
        
        
        Idata_trans[ind] = trans_data['Idata']
        Qdata_trans[ind] = trans_data['Qdata']
        amp_data_trans[ind] = trans_data['amp']
        phase_data_trans[ind] = trans_data['phase']
        ts_trans[ind] = currT_trans
        
        Idata_spec[ind] = spec_data['Idata']
        Qdata_spec[ind] = spec_data['Qdata']
        amp_data_spec[ind] = spec_data['amp']
        phase_data_spec[ind] = spec_data['phase']
        ts_spec[ind] = currT_spec
    

    if exp_settings['subtract_groundstate']:
        Idata_trans = Idata_trans - np.mean(Idata_trans)
        Qdata_trans = Qdata_trans - np.mean(Qdata_trans)
        
        Idata_spec = Idata_spec - np.mean(Idata_spec)
        Qdata_spec = Qdata_spec - np.mean(Qdata_spec)
    
    else:
        Idata_trans = Idata_trans 
        Qdata_trans = Qdata_trans 
        
        Idata_spec = Idata_spec 
        Qdata_spec = Qdata_spec 
        
    # PLOT
    fig = plt.figure(101)
    plt.clf()
    ax = plt.subplot(2,2,1)
    plt.plot(ts_trans, amp_data_trans, label = "amp_trans")
    plt.plot(ts_spec, amp_data_spec, label = "amp_spec")
    plt.xlabel("Time (s)")
    plt.ylabel('v')
    plt.legend()
    
    ax = plt.subplot(2,2,2)
    plt.plot(ts_trans, phase_data_trans, label = "phase_trans")
    plt.plot(ts_spec, phase_data_spec, label = "phase_spec")
    plt.xlabel("Time (s)")
    plt.ylabel('Deg')
    plt.legend()
    
    ax = plt.subplot(2,2,3)
    ax.scatter(Qdata_trans, Idata_trans, s=4, label = 'trans')
    ax.scatter(Qdata_spec, Idata_spec, s=4, label = 'spec')
    plt.xlabel('Q')
    plt.ylabel('I')
    plt.legend()
    plt.axhline(y = 0, color ="black", linestyle ="-")
    plt.axvline(x = 0, color ="black", linestyle ="-") 
    y_max = np.max(np.abs(ax.get_ylim()))
    x_max = np.max(np.abs(ax.get_xlim()))
    list = np.array([x_max, y_max])
    max = np.max(list)
    ax.set_ylim(ymin=-max, ymax=max)
    ax.set_xlim(xmin=-max, xmax=max)
    
    fig.suptitle(f'{num_points} points_{card.averages} Avgs')
    fig.tight_layout()
    plt.show()
    
    save_name = filename + 'TransSpecRepeat'
    userfuncs.savefig(fig, save_name, saveDir, png=True)

    # Save data
    userfuncs.SaveFull(saveDir, filename, ['trans_data'],
                        locals(), 
                        expsettings=settings, 
                        instruments=instruments, saveHWsettings=True)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
       
       
    trans_data_final = {}
    trans_data_final['Idata'] = Idata_trans
    trans_data_final['Qdata'] = Qdata_trans
    trans_data_final['amp_data'] = amp_data_trans
    trans_data_final['phase_data'] = phase_data_trans
    trans_data_final['time'] = ts_trans
    
    spec_data_final = {}
    spec_data_final['Idata'] = Idata_spec
    spec_data_final['Qdata'] = Qdata_spec
    spec_data_final['amp_data'] = amp_data_spec
    spec_data_final['phase_data'] = phase_data_spec
    spec_data_final['time'] = ts_spec
    
    full_data = {}
    full_data['trans_data_final'] = trans_data_final
    full_data['spec_data_final'] = spec_data_final
    
    cavitygen.output = 0
    qubitgen.output = 0
    LO.ouput = 0
    
    return full_data

        
        
    
    