'''
8-25-21 AK modifying to normalize the amplitudes to the drive power. Undid that.

9-2-21 AK made it return the data

'''


import os
import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import userfuncs
from utility.measurement_helpers import configure_card, generate_filter, read_and_process
from utility.scheduler import scheduler

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'coherent_state_cal'
    
    #Qubit parameters
    settings['span']   = 50e6 
    settings['freq_points'] = 51 
    settings['spec_power'] = -15
    
    #Coherent state drive
    settings['start_delay'] = 0
    settings['stop_delay'] = 21e-6
    settings['delay_points'] = 11
    
    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 5e3
    
    #Measurement settings
    settings['num_save']  = 1
    settings['hold_time'] = 100e-6
    settings['drive_time'] = 50e-6
    settings['delay'] = 1e-6
    
    settings['subtract_background'] = True
    
    return settings

def upload_schedule(hdawg, settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    
    ## Sequencer program
    root_dir = Path(__file__).parent
    progFile = open(os.path.join(root_dir, r"HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp"),'r')
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
    awg_sched.add_analog_channel(3, name='blank')
    awg_sched.add_analog_channel(4, name='blank2')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=50e-9, HW_offset_off=50e-9)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(3, name='Boost_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(4, name='digital_blank', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    qubit_I       = awg_sched.analog_channels['Qubit_I']
    qubit_marker  = awg_sched.digital_channels['Qubit_enable']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    boost_marker  = awg_sched.digital_channels['Boost_enable']
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']

    position = start_time-delay-num_sigma*sigma-exp_settings['hold_time']
    qubit_I.add_pulse('gaussian_square', position=position, 
                              amplitude=q_pulse['piAmp'], length = exp_settings['hold_time'], 
                              ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    
    qubit_marker.add_window(position, position+exp_settings['hold_time'])
    cavity_marker.add_window(start_time, start_time+window_time)
    boost_tone_pos = position+exp_settings['hold_time']+num_sigma*sigma
    boost_marker.add_window(boost_tone_pos-exp_settings['drive_time']-exp_settings['delay'], 
                            boost_tone_pos-exp_settings['delay'])
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['blank', 'blank2'], ['Boost_enable', 'digital_blank'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)
    
def coherent_state_evolution(instruments, settings):
    
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    card      = instruments['card']
    hdawg     = instruments['AWG']
    LO        = instruments['LO']
    boostgen  = instruments['boostgen']
    
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
    
    ##Boost settings
    Boost_Attenuation = exp_globals['Boost_Attenuation']
    Boost_freq = calibration['Boost_Freq']
    
    ##Qubit settings
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Q_power = exp_settings['spec_power'] + Qbit_Attenuation
    
    ##Qbit freq sweep
    start_freq  = calibration['Q_Freq']-exp_settings['span']/2
    stop_freq   = calibration['Q_Freq']+exp_settings['span']/2
    freq_points = exp_settings['freq_points']
    
    ##Boost power sweep
    start_delay  = exp_settings['start_delay']
    stop_delay   = exp_settings['stop_delay']
    delay_points = exp_settings['delay_points']
    
    #Configure sweeps 
    freqs  = np.round(np.linspace(start_freq,stop_freq,freq_points),-3)
    delays = np.round(np.linspace(start_delay, stop_delay, delay_points),9)
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.Output = 'On'
    
    boostgen.Freq = Boost_freq
    boostgen.power = calibration['Boost_Power']+Boost_Attenuation
    boostgen.Output = 'On'
    
    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_pulse()
    boostgen.enable_pulse()
    
    ## Setting qubit generator to some safe starting point before we turn it on
    qubitgen.Output  = 'On'
    qubitgen.Freq   = 4e9
    qubitgen.Power  = Q_power

    qubitgen.enable_pulse()
    qubitgen.enable_IQ()
    
    ## Card config
    configure_card(card, settings)
    
    ## Upload schedule
    upload_schedule(hdawg, settings)
    
    ##create the digital down conversion filter if needed.
    generate_filter(card, settings)
    
    mags = np.zeros((exp_settings['delay_points'], exp_settings['freq_points']))
    mag_temp = np.zeros(exp_settings['freq_points'])
    t1 = time.time()
    for dind in range(len(delays)):
        print('Delay:{}us'.format(delays[dind]*1e6))
        exp_settings['delay'] = delays[dind]
        upload_schedule(hdawg, settings)
        
        for find in range(0, len(freqs)):
            qubitgen.freq = freqs[find]
            
            I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                         plot=False, 
                                                                         IQstorage = True)
            ##Useful handles for variables
            I_sig, Q_sig   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
            mag_temp[find] = np.sqrt(I_sig**2+Q_sig**2)
        mags[dind] = mag_temp
        
        f_step = (freqs[1]-freqs[0])/2
        try:
            d_step = ((delays[1]-delays[0])/2)*1e6
        except:
            d_step = 1
        
        fig = plt.figure(151)
        plt.clf()
        plt.imshow(mags[:dind+1], origin='lower', aspect='auto',extent=[(freqs[0]-f_step)/1e9, 
                                                               (freqs[-1]+f_step)/1e9, 
                                                               delays[0]*1e6-d_step, 
                                                               delays[dind]*1e6+d_step])
        plt.xlabel('Freq (GHz)')
        plt.ylabel('Delay (us)')
        plt.title('Coherent state evolution, started at: '+stamp+'\n'+saveDir)
        fig.canvas.draw()
        fig.canvas.flush_events()
        
        plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
        userfuncs.SaveFull(saveDir, filename, ['mags','freqs','delays'],
                                                         locals(), 
                                                         expsettings=settings, 
                                                         instruments=instruments, saveHWsettings=True)
    t2 = time.time()
    elapsed_time = t2-t1
    print('elapsed time = ' + str(elapsed_time))
    
    cavitygen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.output = 'Off'
    boostgen.output = 'Off'
    
    full_data = {
        'mag_array':mags,
        'freqs':freqs,
        'delays':delays
        }
    return full_data