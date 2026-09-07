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
    settings['meas_type'] = 'pulsed_DPA_4wm_qubit'

    # JPA vbias parameters
    settings['JPA_voltage']      = 0.4

    # carrier parameters
    settings['start_sig_freq']  = 5e9
    settings['stop_sig_freq']   = 6e9
    settings['sig_freq_points'] = 11
    settings['del_carr']         = 20e6
    settings['total_power']      = -30

    # Three tone parameters, assume cav power is const
    settings['pdh_power']        = [-40, -60, -40]

    settings['del_cav']          = 2e6
    settings['del_sb']           = 40e6

    settings['start_phase']      = [1,0,0]
    settings['stop_phase']       = [0,0,0]
    settings['phase_points']     = 10

    settings['pump_window']      = 10e-6

    # qubit parameters
    settings['qubit_meas']       = False
    settings['qubit_voltage']    = 0.1
    settings['qubit_freq']       = 4.5e9
    settings['qubit_power']      = -10
    
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



def pulsed_DPA_4wm_qubit(instruments, settings):
    cavitygen = instruments['cavitygen']
    qubitgen  = instruments['qubitgen']
    LO        = instruments['LO']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    pumpbias  = instruments['pumpbias']
    qubitbias = instruments['qubitbias']
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    # Data saving and naming
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename_template = exp_settings['scanname'] + 'qubit_{}' + stamp

    # JPA vbias settings
    JPA_voltage = exp_settings['JPA_voltage']
    
    # Carrier / sig settings
    sig_freq_points = exp_settings['sig_freq_points']
    del_carr = exp_settings['del_carr']
    carr_freqs = np.round(np.linspace(exp_settings['start_sig_freq'],exp_settings['stop_sig_freq'],sig_freq_points),-3) - del_carr
    total_power = exp_settings['total_power'] + exp_globals['CAV_Attenuation']

    # PDH settings
    pdh_power = np.array(exp_settings['pdh_power']) + exp_globals['CAV_Attenuation']

    amp_off_list = get_amp_comps([-150, pdh_power[1], -150], total_power)
    amp_list = get_amp_comps(pdh_power, total_power)

    del_cav = exp_settings['del_cav']
    del_sb = exp_settings['del_sb']

    if del_sb <= 0 or del_carr <= 0:
        raise ValueError("Detunings SB and carrier should be positive.")

    freq_list = np.array([del_carr - del_sb, del_carr + del_cav, del_carr + del_sb])

    phase_points = exp_settings['phase_points']
    phase_list = np.linspace(np.array(exp_settings['start_phase']), np.array(exp_settings['stop_phase']), phase_points)

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
    alpha = np.round(3*alpha,3)
    beta = np.round(2*beta,3)
    gamma = np.round(6*gamma,3)

    # qubit settings
    if exp_settings['qubit_meas']:
        qubit_voltage = exp_settings['qubit_voltage']
        qubit_freq = exp_settings['qubit_freq']
        qubit_power = exp_settings['qubit_power'] + exp_globals['Qbit_Attenuation']


    # voltage bias settings
    pumpbias.voltage_ramp(0)
    pumpbias.Output = 'On'
    pumpbias.voltage_ramp(JPA_voltage)

    qubitbias.voltage_ramp(0)
    qubitbias.Output = 'On'
    if exp_settings['qubit_meas']:
        qubitbias.voltage_ramp(qubit_voltage)
    
    # Generator settings
    cavitygen.Freq   = carr_freqs[0]
    cavitygen.Power  = total_power
    cavitygen.Phase  = 0
    cavitygen.enable_IQ()
    cavitygen.enable_pulse()
    cavitygen.Output = 'On'

    qubitgen.Output = 'Off'
    if exp_settings['qubit_meas']:
        qubitgen.Freq   = qubit_freq
        qubitgen.Power  = qubit_power
        qubitgen.Phase  = 0
        qubitgen.enable_IQ()
        qubitgen.enable_pulse()
        qubitgen.Output = 'Off'

    LO.Power = 12
    LO.Freq = carr_freqs[0] + del_carr + del_cav +  exp_globals['IF']
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
    # 1. for loop phase
    # 2. for loop carr freq
    ###################################################
    tstart = time.time()
    first_it = True

    I_full_g = np.zeros((phase_points, sig_freq_points))
    Q_full_g = np.zeros((phase_points, sig_freq_points))
    powerdat_full_g = np.zeros((phase_points, sig_freq_points))
    phasedat_full_g = np.zeros((phase_points, sig_freq_points))
    diff_full_g = np.zeros((phase_points, sig_freq_points))

    I_full_e = np.zeros((phase_points, sig_freq_points))
    Q_full_e = np.zeros((phase_points, sig_freq_points))
    powerdat_full_e = np.zeros((phase_points, sig_freq_points))
    phasedat_full_e = np.zeros((phase_points, sig_freq_points))
    diff_full_e = np.zeros((phase_points, sig_freq_points))

    I_off = np.zeros(sig_freq_points)
    Q_off = np.zeros(sig_freq_points)
    powerdat_off = np.zeros(sig_freq_points)
    phasedat_off = np.zeros(sig_freq_points)



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



    ###################################################
    # pump off meas
    ###################################################
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

    for i in range(sig_freq_points):
        cavitygen.Freq = carr_freqs[i]
        LO.Freq = carr_freqs[i] + del_carr + del_cav +  exp_globals['IF']
        cavitygen.Phase = 0
        LO.Phase = 0

        I_window, Q_window, _, _, _ = read_and_process_JPA(card, settings, plot=first_it, IQstorage = True)

        I_sig, Q_sig = [np.mean(I_window), np.mean(Q_window)]

        I_off[i] = I_sig
        Q_off[i] = Q_sig
        powerdat_off[i] = V_to_dBm(np.sqrt(I_sig**2 + Q_sig**2))
        phasedat_off[i] = np.arctan2(Q_sig,I_sig)*180/np.pi

        if first_it:
            first_it = False

    first_it = True

    labels = ['signal freq (GHz)', 'Pump off']
    xaxis = carr_freqs + del_carr
    yaxis = [0]

    full_data_off = {}
    full_data_off['xaxis'] = xaxis
    full_data_off['mags'] = np.reshape(powerdat_off, (1, sig_freq_points))
    full_data_off['phases'] = np.reshape(phasedat_off, (1, sig_freq_points))
    full_data_off['I'] = I_off
    full_data_off['Q'] = Q_off

    single_data_off = {}
    single_data_off['xaxis'] = xaxis
    single_data_off['mag'] = powerdat_off
    single_data_off['phase'] = phasedat_off

    simplescan_plot(full_data_off, single_data_off, yaxis, filename_template.format('pumpoff'), labels, identifier='', fig_num=1)
    plt.savefig(os.path.join(saveDir, filename_template.format('pumpoff')+'.png'), dpi = 150)


    
    ###################################################
    # qubit g state meas = cav meas
    ###################################################
    for phind in range(phase_points):
        print(f'Current scissors phase: {gamma[phind]}.')

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
            amp_list = amp_list,
            freq_list = freq_list,
            phase_list = phase_list[phind],
            pump_win = pump_window,
            cav_win = window_time
        )
        cavity_Q.add_pulse(
            type='DPA_Q', position=pump_position,
            amp_list = amp_list,
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


        for find in range(sig_freq_points):
            cavitygen.Freq = carr_freqs[find]
            LO.Freq = carr_freqs[find] + del_carr + del_cav +  exp_globals['IF']
            cavitygen.Phase = 0
            LO.Phase = 0

            I_window, Q_window, _, _, _ = read_and_process_JPA(card, settings, plot=first_it, IQstorage = True)

            I_sig, Q_sig = [np.mean(I_window), np.mean(Q_window)]

            I_full_g[phind, find] = I_sig
            Q_full_g[phind, find] = Q_sig
            powerdat_full_g[phind, find] = V_to_dBm(np.sqrt(I_sig**2 + Q_sig**2))
            phasedat_full_g[phind, find] = np.arctan2(Q_sig, I_sig)*180/np.pi

            diff_full_g[phind] = powerdat_full_g[phind] - powerdat_off

            if first_it:
                first_it = False



        ## Save Full data
        labels = ['signal freq (GHz)', 'Scissors Phase']
        xaxis = carr_freqs + del_carr
        yaxis = gamma[0:phind+1]

        full_data_g = {}
        full_data_g['xaxis'] = xaxis
        full_data_g['yaxis'] = yaxis
        full_data_g['mags'] = diff_full_g[0:phind+1] #powerdat_full[0:pfind+1]
        full_data_g['phases'] = phasedat_full_g[0:phind+1]
        full_data_g['abs_power'] = powerdat_full_g[0:phind+1]
        full_data_g['I'] = I_full_g[0:phind+1]
        full_data_g['Q'] = Q_full_g[0:phind+1]

        single_data = {}
        single_data['xaxis'] = xaxis
        single_data['mag'] = diff_full_g[phind, :]
        single_data['phase'] = phasedat_full_g[phind, :]

        
        simplescan_plot(full_data_g, single_data, yaxis, filename_template.format('g'), labels, identifier='', fig_num=2)
        plt.savefig(os.path.join(saveDir, filename_template.format('g')+'.png'), dpi = 150)
    
    userfuncs.SaveFull(saveDir, filename_template.format('g'), ['full_data_g', 'full_data_off', 'JPA_voltage'], 
                        locals(), expsettings=settings, instruments=instruments, saveHWsettings=True)
    
    first_it = True

    ###################################################
    # qubit e state meas = qubit meas
    ###################################################
    if exp_settings['qubit_meas']:
        for phind in range(phase_points):
            print(f'Current scissors phase: {gamma[phind]}.')

            progFile = open(r"C:\Users\BF2-meas\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
            rawprog  = progFile.read()
            loadprog = rawprog
            progFile.close()
        
            awg_sched = scheduler_pdh(total_time=cycle_total_time, sample_rate=2.4e9)

            awg_sched.add_analog_channel(1, name='Qubit_I')
            awg_sched.add_analog_channel(2, name='Qubit_Q')
            awg_sched.add_analog_channel(3, name='Cavity_I')
            awg_sched.add_analog_channel(4, name='Cavity_Q')
            
            awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(2, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(3, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            
            qubit_I       = awg_sched.analog_channels['Qubit_I']
            cavity_I      = awg_sched.analog_channels['Cavity_I']
            cavity_Q      = awg_sched.analog_channels['Cavity_Q']
            qubit_marker  = awg_sched.digital_channels['Qubit_enable']
            cavity_marker = awg_sched.digital_channels['Cavity_enable']

            # cav pulse
            cavity_I.add_pulse(
                type='DPA_I', position=pump_position,
                amp_list = amp_list,
                freq_list = freq_list,
                phase_list = phase_list[phind],
                pump_win = pump_window,
                cav_win = window_time
            )
            cavity_Q.add_pulse(
                type='DPA_Q', position=pump_position,
                amp_list = amp_list,
                freq_list = freq_list,
                phase_list = phase_list[phind],
                pump_win = pump_window,
                cav_win = window_time
            )
            qubit_I.add_pulse(
                type='gaussian_square', 
                position=qbit_position, 
                amplitude=q_pulse['piAmp'],
                length = q_pulse['hold_time'], 
                ramp_sigma=q_pulse['sigma'], 
                num_sigma=q_pulse['num_sigma']
            )
            # marker settings
            qubit_marker.add_window(qbit_position-160e-9, qbit_position+2*160e-9+q_pulse['hold_time'])
            cavity_marker.add_window(pump_position, pump_position+window_time+pump_window)

            if first_it:
                awg_sched.plot_waveforms()
            
            [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'blank1'])
            [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['Cavity_enable', 'blank2'])

            loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
            hdawg.AWGs[0].load_program(loadprog)
            
            hdawg.AWGs[0].stop()
            hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
            hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
            hdawg.AWGs[0].run_loop()
            hdawg.AWGs[1].run_loop()
            time.sleep(0.1)

            qubitgen.Output = 'On'
            for find in range(sig_freq_points):
                cavitygen.Freq = carr_freqs[find]
                LO.Freq = carr_freqs[find] + del_carr + del_cav +  exp_globals['IF']
                cavitygen.Phase = 0
                LO.Phase = 0

                I_window, Q_window, _, _, _ = read_and_process_JPA(card, settings, plot=first_it, IQstorage = True)

                I_sig, Q_sig = [np.mean(I_window), np.mean(Q_window)]

                I_full_e[phind, find] = I_sig
                Q_full_e[phind, find] = Q_sig
                powerdat_full_e[phind, find] = V_to_dBm(np.sqrt(I_sig**2 + Q_sig**2))
                phasedat_full_e[phind, find] = np.arctan2(Q_sig, I_sig)*180/np.pi

                diff_full_e[phind] = powerdat_full_e[phind] - powerdat_off

                if first_it:
                    first_it = False



            ## Save Full data
            labels = ['signal freq (GHz)', 'Scissors Phase']
            xaxis = carr_freqs + del_carr
            yaxis = gamma[0:phind+1]

            full_data_e = {}
            full_data_e['xaxis'] = xaxis
            full_data_e['yaxis'] = yaxis
            full_data_e['mags'] = diff_full_e[0:phind+1] #powerdat_full[0:pfind+1]
            full_data_e['phases'] = phasedat_full_e[0:phind+1]
            full_data_e['abs_power'] = powerdat_full_e[0:phind+1]
            full_data_e['I'] = I_full_e[0:phind+1]
            full_data_e['Q'] = Q_full_e[0:phind+1]

            single_data = {}
            single_data['xaxis'] = xaxis
            single_data['mag'] = diff_full_e[phind, :]
            single_data['phase'] = phasedat_full_e[phind, :]

            
            simplescan_plot(full_data_e, single_data, yaxis, filename_template.format('e'), labels, identifier='', fig_num=3)
            plt.savefig(os.path.join(saveDir, filename_template.format('e')+'.png'), dpi = 150)
        
        userfuncs.SaveFull(saveDir, filename_template.format('e'), ['full_data_e', 'full_data_off', 'JPA_voltage', 'qubit_voltage'], 
                            locals(), expsettings=settings, instruments=instruments, saveHWsettings=True)

    else:
        full_data_e = {}


    tstop = time.time()
    print('elapsed time = ' + str(tstop-tstart))
       
    cavitygen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.Output = 'Off'
    
    return full_data_g, full_data_e, full_data_off