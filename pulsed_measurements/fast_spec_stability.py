'''
8-25-21 AK modifying to normalize the amplitudes to the drive power. Undid that.

9-2-21 AK made it return the data

'''


import os
import time
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, generate_filter, read_and_process
from utility.scheduler import scheduler

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'pulsed_spec_stability'
    
    #Qubit parameters
    settings['span']   = 50e6 
    settings['freq_points'] = 51 
    settings['spec_power'] = -15
    
    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 5e3
    
    #Measurement settings
    settings['Quasi_CW']    = False
    settings['reverse']   = False
    settings['num_save'] = 1
    settings['repeats'] = 20
    settings['hold_time']      = 30e-6
    #background_subtraction (by taking reference trace with no qubit drive power)
    settings['subtract_background'] = False
    
    return settings

def pulsed_spec_stability(instruments, settings):
    
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    calibration  = settings['calibration']
    
    ##Data saving and naming
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp
    
    ##Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = calibration['CAV_Power'] + CAV_Attenuation
    CAV_freq  = calibration['CAV_Freq']
    
    ##Qubit settings
    start_freq  = calibration['Q_Freq']-exp_settings['span']/2
    stop_freq   = calibration['Q_Freq']+exp_settings['span']/2
    freq_points = exp_settings['freq_points']

    freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
    
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Q_power = exp_settings['spec_power'] + Qbit_Attenuation
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.Output = 'On'
    qubitgen.Output  = 'On'
    
    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_pulse()
    cavitygen.output = 'On'
    
    ## Setting qubit generator to some safe starting point before we turn it on
    qubitgen.Freq   = 4e9
    qubitgen.Power  = Q_power

    qubitgen.enable_pulse()
    qubitgen.enable_IQ()
    
    ## Card config
    configure_card(card, settings)
    
    ## Sequencer program
    progFile = open(r"C:\Users\kollarlab\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
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
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=50e-9, HW_offset_off=50e-9)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    
    qubit_I       = awg_sched.analog_channels['Qubit_I']
    qubit_marker  = awg_sched.digital_channels['Qubit_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']

    position = start_time-delay-num_sigma*sigma-exp_settings['hold_time']
    qubit_I.add_pulse('gaussian_square', position=position, 
                              amplitude=q_pulse['piAmp'], length = exp_settings['hold_time'], 
                              ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    
    qubit_marker.add_window(position, position+exp_settings['hold_time'])
    cavity_marker.add_window(start_time, start_time+window_time)
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)
    
    
    ##create the digital down conversion filter if needed.
    generate_filter(card, settings)
    
    mags = np.zeros((exp_settings['repeats'], exp_settings['freq_points']))
    mag_temp = np.zeros(exp_settings['freq_points'])
    t1 = time.time()
    for ind in range(exp_settings['repeats']):
        print(ind)
        for find in range(0, len(freqs)):
            qubitgen.freq = freqs[find]
            
            I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                         plot=False, 
                                                                         IQstorage = True)
            ##Useful handles for variables
            I_sig, Q_sig   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
            mag_temp[find] = np.sqrt(I_sig**2+Q_sig**2)
        mags[ind] = mag_temp
    
    t2 = time.time()
    f_step = (freqs[1]-freqs[0])/2
    elapsed_time = t2-t1
    plt.figure()
    plt.imshow(mags, origin='lower', aspect='auto',extent=[(freqs[0]-f_step)/1e9, 
                                                           (freqs[-1]+f_step)/1e9, 0, elapsed_time])
    plt.xlabel('Freq (GHz)')
    plt.ylabel('Elapsed time (s)')
    plt.title('Stability started at: '+stamp)
    plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
    userfuncs.SaveFull(saveDir, filename, ['mags','freqs'],
                                                     locals(), 
                                                     expsettings=settings, 
                                                     instruments=instruments, saveHWsettings=True)

    print('elapsed time = ' + str(t2-t1))
       
    cavitygen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.output = 'Off'
    full_data = {
        'mag_array':mags,
        'freqs':freqs
        }
    return full_data