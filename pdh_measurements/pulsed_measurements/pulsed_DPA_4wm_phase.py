import os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import userfuncs
from utility.plotting_tools import simplescan_plot, load_fig
from pdh_measurements.utility.measurement_helpers_pdh import configure_card, estimate_time, read_and_process_JPA, get_amp_comps, V_to_dBm
from pdh_measurements.utility.scheduler_pdh import scheduler_pdh

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'pulsed_DPA_4wm_phase'

    # vbias parameters
    settings['voltage']          = 0.4

    # cavity parameters
    settings['carr_freq']         = 5e9
    settings['del_carr']          = 20e6
    settings['total_power']       = -30

    # Three tone parameters, assume cav power is const
    settings['start_power']      = [-40, -45, -40]
    settings['stop_power']       = [-35, -45, -35]
    settings['power_points']     = 10

    settings['del_cav']          = 2e6
    settings['del_sb']           = 40e6
    
    settings['start_phase']      = [1,0,0]
    settings['stop_phase']       = [0,0,0]
    settings['phase_points']     = 10

    settings['pump_window']      = 10e-6

    settings['plot_phase']       = 'scissors'
    
    # card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 1e3
    
    # measurement settings
    settings['Quasi_CW']            = False
    settings['reverse']             = False
    settings['num_save']            = 1
    settings['subtract_background'] = False
    
    return settings



def pulsed_DPA_4wm_phase(instruments, settings):
    cavitygen = instruments['cavitygen']
    # qubitgen  = instruments['qubitgen']
    LO        = instruments['LO']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    pumpbias  = instruments['pumpbias']
    # qubitbias = instruments['qubitbias']
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    # Data saving and naming
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename_template = exp_settings['scanname'] + '{}V_C{}GHz_\u0394{}MHz_\u03B4{}MHz' + stamp

    # vbias settings
    voltage = exp_settings['voltage']
    
    # Cavity settings
    carr_freq = exp_settings['carr_freq']
    del_carr = exp_settings['del_carr']

    # PDH settings
    total_power = exp_settings['total_power'] + exp_globals['CAV_Attenuation']
    start_power = np.array(exp_settings['start_power']) + exp_globals['CAV_Attenuation']
    stop_power = np.array(exp_settings['stop_power']) + exp_globals['CAV_Attenuation']
    power_points = exp_settings['power_points']
    power_list = np.round(np.linspace(start_power, stop_power, power_points),2)

    amp_off_list = get_amp_comps([-150, start_power[1], -150], total_power)
    amp_list = np.zeros((power_points, 3))
    print(amp_off_list)

    for i in range(power_points):
        amp_list[i] = get_amp_comps(power_list[i], total_power)

    del_cav = exp_settings['del_cav']
    del_sb = exp_settings['del_sb']

    if del_sb <= 0 or del_carr <= 0:
        raise ValueError("Detunings SB and carrier should be positive.")

    freq_list = np.array([del_carr - del_sb, del_carr + del_cav, del_carr + del_sb])

    start_phase = np.array(exp_settings['start_phase'])
    stop_phase  = np.array(exp_settings['stop_phase'])
    phase_points = exp_settings['phase_points']
    phase_list = np.linspace(start_phase, stop_phase, phase_points)
    # scissors_phase = 2*phase_list[:,1] - (phase_list[:,0] + phase_list[:,2])


    # calculate phase vector
    # 1. get from inverse matrix
    # B = np.array([
    #     [ 1, -1, -1],
    #     [ 1,  0,  2],
    #     [ 1,  1, -1]
    # ], dtype=float)

    # coeff = np.linalg.solve(B, phase_list.T).T

    # alpha = coeff[:, 0]
    # beta  = coeff[:, 1]
    # gamma = coeff[:, 2]

    # 2. Use the known inverse matrix
    phi_m = phase_list[:, 0]
    phi_0 = phase_list[:, 1]
    phi_p = phase_list[:, 2]

    alpha = (phi_m + phi_0 + phi_p) / 3
    beta  = (phi_p - phi_m) / 2
    gamma = (2*phi_0 - phi_m - phi_p) / 6

    # what gain is dependent on is actual phase diff, not normalized vector
    alpha = 3*alpha
    beta = 2*beta
    gamma = 6*gamma


    # voltage bias settings
    pumpbias.voltage_ramp(0)
    pumpbias.Output = 'On'
    pumpbias.voltage_ramp(voltage)
    
    # Generator settings
    cavitygen.Freq   = carr_freq
    cavitygen.Power  = total_power
    cavitygen.Phase  = 0
    cavitygen.enable_IQ()
    cavitygen.enable_pulse()
    cavitygen.Output = 'On'

    LO.Power = 12
    LO.Freq = carr_freq + del_carr + del_cav +  exp_globals['IF']
    LO.Phase = 0
    LO.Output = 'On'
    LO.disable_pulse()

    # exp_globals measurement settings
    m_pulse = exp_globals['measurement_pulse']
    window_time = m_pulse['meas_window']
    pump_window = exp_settings['pump_window']
    cav_position  = m_pulse['meas_pos']
    pump_position = cav_position - pump_window/2
    cycle_buffer = 20e-6 #padding time to give the AWG a break between cycles
    cycle_total_time = 1/exp_globals['trigger_rate'] - cycle_buffer

    q_pulse = exp_globals['qubit_pulse']
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    qbit_position = cav_position-delay-num_sigma*sigma-q_pulse['hold_time']


    
    ###################################################
    # phase scan
    ###################################################
    tstart = time.time()
    filename = filename_template.format(voltage, (carr_freq+del_carr)/1e9, del_sb/1e6, del_cav/1e6)
    first_it = True

    I_full = np.zeros((power_points, phase_points))
    Q_full = np.zeros((power_points, phase_points))
    powerdat_full = np.zeros((power_points, phase_points))
    phasedat_full = np.zeros((power_points, phase_points))
    diff_full = np.zeros((power_points, phase_points))

    I_off = np.zeros(1)
    Q_off = np.zeros(1)
    powerdat_off = np.zeros(1)
    phasedat_off = np.zeros(1)



    ###################################################
    # pump off meas
    ###################################################
    configure_card(card, settings)

    # Digital Filter settings
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


    progFile = open(r"C:\Users\BF2-meas\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()

    awg_sched = scheduler_pdh(total_time=cycle_total_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='blank1')
    awg_sched.add_analog_channel(2, name='blank2')
    awg_sched.add_analog_channel(3, name='Cavity_I')
    awg_sched.add_analog_channel(4, name='Cavity_Q')
    
    awg_sched.add_digital_channel(1, name='blank3', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(2, name='blank4', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(3, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(4, name='blank5', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    cavity_I      = awg_sched.analog_channels['Cavity_I']
    cavity_Q      = awg_sched.analog_channels['Cavity_Q']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']

    # cav pulse
    cavity_I.add_pulse(
        type='DPA_I', position=pump_position,
        amp_list = amp_off_list,
        freq_list = freq_list, 
        phase_list = phase_list[0],
        pump_win = pump_window,
        cav_win = window_time
    )
    cavity_Q.add_pulse(
        type='DPA_Q', position=pump_position,
        amp_list = amp_off_list,
        freq_list = freq_list, 
        phase_list = phase_list[0],
        pump_win = pump_window,
        cav_win = window_time
    )

    # marker settings
    cavity_marker.add_window(pump_position, pump_position+window_time+pump_window)

    if first_it:
        awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['blank1', 'blank2'], ['blank3', 'blank4'])
    [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['Cavity_enable', 'blank5'])

    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    
    hdawg.AWGs[0].stop()
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
    hdawg.AWGs[0].run_loop()
    hdawg.AWGs[1].run_loop()
    time.sleep(0.1)



    # substract pump off data for each vbias
    cavitygen.Freq = carr_freq
    LO.Freq = carr_freq + del_carr + del_cav +  exp_globals['IF']
    cavitygen.Phase = 0
    LO.Phase = 0

    I_window, Q_window, _, _, _ = read_and_process_JPA(card, settings, plot=first_it, IQstorage = True)

    I_sig, Q_sig = [np.mean(I_window), np.mean(Q_window)]

    I_off[0] = I_sig
    Q_off[0] = Q_sig
    powerdat_off[0] = V_to_dBm(np.sqrt(I_sig**2 + Q_sig**2))
    # powerdat_off[0] = np.sqrt(I_sig**2 + Q_sig**2)
    phasedat_off[0] = np.arctan2(Q_sig,I_sig)*180/np.pi

    full_data_off = {}
    full_data_off['mags'] = powerdat_off
    full_data_off['phases'] = phasedat_off
    full_data_off['I'] = I_off
    full_data_off['Q'] = Q_off


    
    ###################################################
    # pump on meas
    ###################################################
    for pind in range(power_points):
        print(f'Current power: {power_list[pind, 0] - exp_globals["CAV_Attenuation"]}.')

        for phind in range(phase_points):
            progFile = open(r"C:\Users\BF2-meas\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
            rawprog  = progFile.read()
            loadprog = rawprog
            progFile.close()
        
            awg_sched = scheduler_pdh(total_time=cycle_total_time, sample_rate=2.4e9)

            awg_sched.add_analog_channel(1, name='blank1')
            awg_sched.add_analog_channel(2, name='blank2')
            awg_sched.add_analog_channel(3, name='Cavity_I')
            awg_sched.add_analog_channel(4, name='Cavity_Q')
            
            awg_sched.add_digital_channel(1, name='blank3', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(2, name='blank4', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(3, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(4, name='blank5', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            
            cavity_I      = awg_sched.analog_channels['Cavity_I']
            cavity_Q      = awg_sched.analog_channels['Cavity_Q']
            cavity_marker = awg_sched.digital_channels['Cavity_enable']

            # cav pulse
            cavity_I.add_pulse(
                type='DPA_I', position=pump_position,
                amp_list = amp_list[pind],
                freq_list = freq_list,
                phase_list = phase_list[phind],
                pump_win = pump_window,
                cav_win = window_time
            )
            cavity_Q.add_pulse(
                type='DPA_Q', position=pump_position,
                amp_list = amp_list[pind],
                freq_list = freq_list,
                phase_list = phase_list[phind],
                pump_win = pump_window,
                cav_win = window_time
            )

            # marker settings
            cavity_marker.add_window(pump_position, pump_position+window_time+pump_window)

            if first_it:
                awg_sched.plot_waveforms()
            
            [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['blank1', 'blank2'], ['blank3', 'blank4'])
            [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['Cavity_enable', 'blank5'])

            loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
            hdawg.AWGs[0].load_program(loadprog)
            
            hdawg.AWGs[0].stop()
            hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
            hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
            hdawg.AWGs[0].run_loop()
            hdawg.AWGs[1].run_loop()
            time.sleep(0.1)



            # cavitygen.Freq = cav_freq
            # LO.Freq = cav_freq + exp_globals['IF']
            cavitygen.Phase = 0
            LO.Phase = 0

            I_window, Q_window, _, _, _ = read_and_process_JPA(card, settings, plot=first_it, IQstorage = True)

            I_sig, Q_sig = [np.mean(I_window), np.mean(Q_window)]

            I_full[pind, phind] = I_sig
            Q_full[pind, phind] = Q_sig
            powerdat_full[pind, phind] = V_to_dBm(np.sqrt(I_sig**2 + Q_sig**2))
            # powerdat_full[pind, phind] = np.sqrt(I_sig**2 + Q_sig**2)
            phasedat_full[pind, phind] = np.arctan2(Q_sig, I_sig)*180/np.pi

            if first_it:
                first_it = False
            
            diff_full[pind] = powerdat_full[pind] - powerdat_off



        ## Save Full data
        if exp_settings['plot_phase'] == 'common':
            xaxis = alpha
            labels = ['Common Mode', 'Pump Power Sweep (dBm)']
        elif exp_settings['plot_phase'] == 'differential':
            xaxis = beta
            labels = ['Differential Mode', 'Pump Power Sweep (dBm)']
        else:
            xaxis = gamma
            labels = ['Scissors Phase', 'Pump Power Sweep (dBm)']
        yaxis = power_list[0:pind+1, 0] - exp_globals['CAV_Attenuation']

        full_data = {}
        full_data['xaxis'] = xaxis
        full_data['yaxis'] = yaxis
        full_data['mags'] = diff_full[0:pind+1] #powerdat_full[0:pfind+1]
        full_data['phases'] = phasedat_full[0:pind+1]
        full_data['abs_power'] = powerdat_full[0:pind+1]
        full_data['I'] = I_full[0:pind+1]
        full_data['Q'] = Q_full[0:pind+1]

        single_data = {}
        single_data['xaxis'] = xaxis
        single_data['mag'] = powerdat_full[pind, :]
        single_data['phase'] = phasedat_full[pind, :]

        
        simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1)
        plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
    
    userfuncs.SaveFull(saveDir, filename, ['full_data', 'full_data_off', 'voltage', 'amp_list', 'freq_list', 'phase_list'], 
                        locals(), expsettings=settings, instruments=instruments, saveHWsettings=True)


    tstop = time.time()
    print('elapsed time = ' + str(tstop-tstart))
       
    cavitygen.Output = 'Off'
    LO.Output = 'Off'
    
    return full_data, full_data_off