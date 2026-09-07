import os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import userfuncs
from utility.plotting_tools import simplescan_plot
from pdh_measurements.utility.measurement_helpers_pdh import configure_card, estimate_time, read_and_process_JPA, get_amp_comps, total_power, V_to_dBm
from pdh_measurements.utility.scheduler_pdh import scheduler_pdh

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'pulsed_NDPA_doublepump_powercal'

    # vbias parameters
    settings['start_voltage']    = 0.4
    settings['stop_voltage']     = 0.9
    settings['voltage_points']   = 15

    # pump parameters
    settings['total_pump_power'] = -30
    settings['start_pump_power'] = -50
    settings['stop_pump_power']  = -48
    settings['pump_power_points']= 10

    settings['pump_detuning']  = 10e6
    
    # cavity parameters
    settings['carr_power']         = -60
    settings['carr_freq']          = 7e9
    settings['sig_detuning']       = 5e6 #delta

    # special flags
    settings['upper_pump'] = True
    settings['lower_pump'] = True

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



def pulsed_NDPA_doublepump_powercal(instruments, settings):
    cavitygen = instruments['cavitygen']
    pumpgen   = instruments['pumpgen']
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
    filename_template = exp_settings['scanname'] + '{}V' + stamp

    # vbias settings
    voltage_points = exp_settings['voltage_points']
    voltages = np.round(np.linspace(exp_settings['start_voltage'],exp_settings['stop_voltage'],voltage_points),5)
    
    # Cavity settings
    carr_freq = exp_settings['carr_freq']
    sig_detuning = exp_settings['sig_detuning']
    cav_freq = carr_freq + sig_detuning

    cav_power = exp_settings['carr_power'] + exp_globals['CAV_Attenuation']
    cav_amp_list = get_amp_comps([-150, cav_power, -150], cav_power+7)
    
    # Pump settings
    pump_detuning = exp_settings['pump_detuning']

    pump_power_start    = exp_settings['start_pump_power'] + exp_globals['Qbit_Attenuation']
    pump_power_stop     = exp_settings['stop_pump_power'] + exp_globals['Qbit_Attenuation']
    pump_power_points   = exp_settings['pump_power_points']
    pump_powers         = np.round(np.linspace(pump_power_start,pump_power_stop, pump_power_points),2)
    total_pump_power    = exp_settings['total_pump_power'] + exp_globals['Qbit_Attenuation']

    pump_amp_list  = np.zeros((pump_power_points, pump_power_points, 3))

    # usb = yaxis
    # lsb = xaxis
    for i in range(pump_power_points):
        for j in range(pump_power_points):
            if exp_settings['upper_pump']:
                usb_power = pump_powers[i]
            else:
                usb_power = -150

            if exp_settings['lower_pump']:
                lsb_power = pump_powers[j]
            else:
                lsb_power = -150  
            pump_amp_list[i,j] = get_amp_comps([lsb_power, -150, usb_power], total_pump_power)



    # voltage bias settings
    pumpbias.voltage_ramp(0)
    pumpbias.Output = 'On'
    
    # Generator settings
    cavitygen.Freq   = cav_freq
    cavitygen.Power  = cav_power+7
    cavitygen.Phase  = 0
    cavitygen.enable_IQ()
    cavitygen.enable_pulse()
    cavitygen.Output = 'On'
    
    pumpgen.Freq  = carr_freq
    pumpgen.Power = total_pump_power
    pumpgen.Phase = 0
    pumpgen.enable_pulse()
    pumpgen.enable_IQ()
    pumpgen.Output = 'Off'

    LO.Power = 12
    LO.Freq = cav_freq + exp_globals['IF']
    LO.Phase = 0
    LO.Output = 'On'
    LO.disable_pulse()

    # exp_globals measurement settings
    m_pulse = exp_globals['measurement_pulse']
    cav_position  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    cycle_buffer = 20e-6 #padding time to give the AWG a break between cycles
    cycle_total_time = 1/exp_globals['trigger_rate'] - cycle_buffer

    q_pulse = exp_globals['qubit_pulse']
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    qbit_position = cav_position-delay-num_sigma*sigma-q_pulse['hold_time']
    



    ###################################################
    # 1. for loop vbias
    # 2. for loop pump power
    # 3. for loop pump detuning
    ###################################################
    tstart = time.time()
    for vind in range(voltage_points):
        current_voltage = voltages[vind]
        pumpbias.voltage_ramp(voltages[vind])
        print(f'Current voltage: {current_voltage}.')

        filename = filename_template.format(current_voltage)
        first_it = True

        I_full = np.zeros((pump_power_points, pump_power_points))
        Q_full = np.zeros((pump_power_points, pump_power_points))
        powerdat_full = np.zeros((pump_power_points, pump_power_points))
        phasedat_full = np.zeros((pump_power_points, pump_power_points))
        diff_full = np.zeros((pump_power_points, pump_power_points))

        I_off = np.zeros(1)
        Q_off = np.zeros(1)
        powerdat_off = np.zeros(1)
        phasedat_off = np.zeros(1)

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



        # pump off
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
            type='pdh_I', position=cav_position,
            amp_list = cav_amp_list, phase_list = [0,0,0],
            mod_freq = 10e6,
            time = window_time
        )
        cavity_Q.add_pulse(
            type='pdh_Q', position=cav_position,
            amp_list = cav_amp_list, phase_list = [0,0,0],
            mod_freq = 10e6,
            time = window_time
        )

        # marker settings
        cavity_marker.add_window(cav_position, cav_position+window_time)

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

        # pump off measurement
        pumpgen.Output='Off'
        time.sleep(0.1)
        cavitygen.Freq = cav_freq
        LO.Freq = cav_freq + exp_globals['IF']
        cavitygen.Phase = 0
        pumpgen.Phase = 0
        LO.Phase = 0

        I_window, Q_window, _, _, _ = read_and_process_JPA(card, settings, plot=first_it, IQstorage = True)

        I_sig, Q_sig = [np.mean(I_window), np.mean(Q_window)]

        I_off[0] = I_sig
        Q_off[0] = Q_sig
        powerdat_off[0] = V_to_dBm(np.sqrt(I_sig**2 + Q_sig**2))
        phasedat_off[0] = np.arctan2(Q_sig,I_sig)*180/np.pi

        full_data_off = {}
        full_data_off['mags'] = powerdat_off[0]
        full_data_off['phases'] = phasedat_off[0]
        full_data_off['I'] = I_off
        full_data_off['Q'] = Q_off



        for usbind in range(pump_power_points):
            current_pump_power = pump_powers[usbind] - exp_globals['Qbit_Attenuation']
            print(f'Current pump power: {current_pump_power}.')
            
            for lsbind in range(pump_power_points):
                progFile = open(r"C:\Users\BF2-meas\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
                rawprog  = progFile.read()
                loadprog = rawprog
                progFile.close()
            
                awg_sched = scheduler_pdh(total_time=cycle_total_time, sample_rate=2.4e9)

                awg_sched.add_analog_channel(1, name='pump_I')
                awg_sched.add_analog_channel(2, name='pump_Q')
                awg_sched.add_analog_channel(3, name='Cavity_I')
                awg_sched.add_analog_channel(4, name='Cavity_Q')
                
                awg_sched.add_digital_channel(1, name='pump_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
                awg_sched.add_digital_channel(2, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
                awg_sched.add_digital_channel(3, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
                awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
                
                cavity_I      = awg_sched.analog_channels['Cavity_I']
                cavity_Q      = awg_sched.analog_channels['Cavity_Q']
                cavity_marker = awg_sched.digital_channels['Cavity_enable']
                pump_I        = awg_sched.analog_channels['pump_I']
                pump_Q        = awg_sched.analog_channels['pump_Q']
                pump_marker   = awg_sched.digital_channels['pump_enable']

                # cav pulse
                cavity_I.add_pulse(
                    type='pdh_I', position=cav_position,
                    amp_list = cav_amp_list, phase_list = [0,0,0],
                    mod_freq = 10e6,
                    time = window_time
                )
                cavity_Q.add_pulse(
                    type='pdh_Q', position=cav_position,
                    amp_list = cav_amp_list, phase_list = [0,0,0],
                    mod_freq = 10e6,
                    time = window_time
                )
                
                # pump pulse
                pump_I.add_pulse(
                    type='pdh_I', position=0,
                    amp_list = pump_amp_list[usbind, lsbind], phase_list = [0,0,0],
                    mod_freq = pump_detuning,
                    time = cycle_total_time
                )
                pump_Q.add_pulse(
                    type='pdh_Q', position=0,
                    amp_list = pump_amp_list[usbind, lsbind], phase_list = [0,0,0],
                    mod_freq = pump_detuning,
                    time = cycle_total_time
                )

                # marker settings
                cavity_marker.add_window(cav_position, cav_position+window_time)
                pump_marker.add_window(0, cycle_total_time)

                if first_it:
                    awg_sched.plot_waveforms()
                
                [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['pump_I', 'pump_Q'], ['pump_enable', 'blank1'])
                [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['Cavity_enable', 'blank2'])

                loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
                hdawg.AWGs[0].load_program(loadprog)
                
                hdawg.AWGs[0].stop()
                hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
                hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
                hdawg.AWGs[0].run_loop()
                hdawg.AWGs[1].run_loop()
                time.sleep(0.1)

            
            
                # Ready to do actual scans
                pumpgen.Output = 'On'
                time.sleep(0.1)

                cavitygen.Freq = cav_freq
                pumpgen.Freq = carr_freq
                LO.Freq = cav_freq + exp_globals['IF']
                cavitygen.Phase = 0
                pumpgen.Phase = 0
                LO.Phase = 0

                I_window, Q_window, _, _, _ = read_and_process_JPA(card, settings, plot=first_it, IQstorage = True)

                I_sig, Q_sig = [np.mean(I_window), np.mean(Q_window)]

                I_full[usbind, lsbind] = I_sig
                Q_full[usbind, lsbind] = Q_sig
                powerdat_full[usbind, lsbind] = V_to_dBm(np.sqrt(I_sig**2 + Q_sig**2))
                phasedat_full[usbind, lsbind] = np.arctan2(Q_sig, I_sig)*180/np.pi

                if first_it:
                    first_it = False
                    
            diff_full[usbind] = powerdat_full[usbind] - powerdat_off[0]
        

            ## Save Full data
            full_data = {}
            full_data['xaxis'] = pump_powers - exp_globals['Qbit_Attenuation']
            full_data['yaxis'] = pump_powers[0:usbind+1] - exp_globals['Qbit_Attenuation']
            full_data['mags'] = diff_full[0:usbind+1] #powerdat_full[0:pind+1]
            full_data['phases'] = phasedat_full[0:usbind+1]
            full_data['abs_power'] = powerdat_full[0:usbind+1]
            full_data['I'] = I_full[0:usbind+1]
            full_data['Q'] = Q_full[0:usbind+1]

            single_data = {}
            single_data['xaxis'] = pump_powers - exp_globals['Qbit_Attenuation']
            single_data['mag'] = powerdat_full[usbind, :]
            single_data['phase'] = phasedat_full[usbind, :]

            yaxis = full_data['yaxis']
            labels = ['lsb Power (dBm)', 'usb Power (dBm)']
            simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier='', fig_num=1)
            plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
        userfuncs.SaveFull(saveDir, filename, ['full_data', 'full_data_off', 'current_voltage'], 
                            locals(), expsettings=settings, instruments=instruments, saveHWsettings=True)


    tstop = time.time()
    print('elapsed time = ' + str(tstop-tstart))

    cavitygen.Output = 'Off'
    pumpgen.Output = 'Off'
    LO.Output = 'Off'
    
    return full_data, full_data_off