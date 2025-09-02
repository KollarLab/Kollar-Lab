import os
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process
from utility.scheduler import scheduler

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'pulsed_gefstate'
    
    #Cavity parameters
    settings['CAV_Power']       = -60
    settings['CAV_freq']        = 7e9
    settings['sideband_power']  = -60
    settings['mod_freq']        = 40e6
    
    #Qubit parameters
    settings['Qbit_freq']       = 4.15e9  
    settings['freq_points']     = 1

    settings['Qbit_power']     = -30
    settings['power_points']    = 1

    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 1
    settings['attempted_averages'] = int(2e3)
    
    #Measurement settings
    settings['Quasi_CW']    = False
    settings['reverse']   = False
    settings['num_save'] = 1
    
    #background_subtraction (by taking reference trace with no qubit drive power)
    settings['subtract_background'] = False
    settings['subtract_groundstate'] = True
    
    return settings

def find_pointsincircle(x, y, centerx, centery, samples):
    rxmax = np.max(x-centerx)
    rymax = np.max(y-centery)
    rmax = np.sqrt(rxmax**2 + rymax**2)
    rlist = np.linspace(0, rmax, samples)

    for i in range(samples):
        point_count = 0
        for px, py in zip(x, y):
            if (px-centerx)**2 + (py-centery)**2 <= rlist[i]**2:
                point_count += 1
        r = rlist[i]
        if point_count > samples*0.7609:
            break

    return r
            


def pulsed_gefstate(instruments, settings):
    
    ##Instruments used
    qubitgen  = instruments['qubitgen']
    cavitygen = instruments['cavitygen']
    sidebandgen = instruments['sidebandgen']
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
    sideband_power = exp_settings['sideband_power'] + CAV_Attenuation
    sideband_freq = exp_settings['CAV_freq'] - exp_settings['mod_freq']
    sideband_freqs = np.array([sideband_freq])
    
    ##Qubit settings
    Qbit_freq  = exp_settings['Qbit_freq']
    freq_points = exp_settings['freq_points']
    freqs = np.array([Qbit_freq])
    
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power  = exp_settings['Qbit_power'] + Qbit_Attenuation
    power_points = exp_settings['power_points']
    powers = np.array([Qbit_power])
    attempted_averages = int(exp_settings['attempted_averages'])
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = CAV_power
    cavitygen.Output = 'On'

    sidebandgen.Freq = sideband_freq
    sidebandgen.Power = sideband_power
    sidebandgen.Output = 'On'

    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_pulse()
    cavitygen.output = 'On'
    sidebandgen.enable_pulse()
    sidebandgen.output = 'On'






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
    
    # Make sure sideband signal is off for this carrier only transmission scan
    qubitgen.Output = 'Off'
    sidebandgen.output = 'Off'

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
    Iblob_ctrans = np.zeros(attempted_averages)
    Qblob_ctrans = np.zeros(attempted_averages)
    powerdat_ctrans = np.zeros(attempted_averages)
    phasedat_ctrans = np.zeros(attempted_averages)
    
    tstart = time.time()
    first_it = True

    if first_it:
        tstop = time.time()

    for avgind in range(attempted_averages):
        cavitygen.power = CAV_power
        
        time.sleep(0.2)

        total_samples = card.samples 
        Is  = np.zeros((len(CAV_freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs  = np.zeros((len(CAV_freqs), total_samples))
    
        cavitygen.freq = CAV_freq

        #LO.freq   = freq-exp_globals['IF']
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
        

        
        # Store data
        Iblob_ctrans[avgind] = I_final
        Qblob_ctrans[avgind] = Q_final
        Iblob_ctrans_avg = np.mean(Iblob_ctrans)
        Qblob_ctrans_avg = np.mean(Qblob_ctrans)
        powerdat_ctrans[avgind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_ctrans[avgind] = np.arctan2(Q_final, I_final)*180/np.pi

    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
    r_ctrans = find_pointsincircle(Iblob_ctrans, Qblob_ctrans, Iblob_ctrans_avg, Qblob_ctrans_avg, attempted_averages)
            
    ctrans_data = {}
    ctrans_data['xaxis'] = CAV_freq/1e9
    ctrans_data['mags'] = powerdat_ctrans
    ctrans_data['phases'] = phasedat_ctrans
    ctrans_data['radius'] = r_ctrans

    if exp_settings['subtract_groundstate']:
        ctrans_data['Iblob'] = Iblob_ctrans - Iblob_ctrans_avg
        ctrans_data['Qblob'] = Qblob_ctrans - Qblob_ctrans_avg
        ctrans_data['center'] = [Iblob_ctrans_avg - Iblob_ctrans_avg, Qblob_ctrans_avg - Qblob_ctrans_avg]
    else:
        ctrans_data['Iblob'] = Iblob_ctrans
        ctrans_data['Qblob'] = Qblob_ctrans
        ctrans_data['center'] = [Iblob_ctrans_avg, Qblob_ctrans_avg]





    #######################################################################
    # Second, do carrier spec scan
    #######################################################################
    ## Setting qubit generator to some safe starting point before we turn it on
    qubitgen.Output  = 'On'
    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20
    
    if exp_settings['Quasi_CW']:
        qubitgen.disable_pulse()
        qubitgen.disable_IQ()
    else:
        qubitgen.enable_pulse()
        qubitgen.enable_IQ()
    
    ## Card config
    configure_card(card, settings)

    ## HDAWG settings
#    configure_hdawg(hdawg, settings)
    
    ## Sequencer program
    # progFile = open(r"C:\Users\kollarlab\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
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
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
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
#    qubitgen.disable_IQ()
#    qubit_marker.add_window(0,start_time+2e-6)
    ##
    cavity_marker.add_window(start_time, start_time+window_time)
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)

    # Make sure sideband signal is off for this carrier only spec scan
    sidebandgen.output = 'Off'
    
    
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
    
    #powerdat = np.zeros((len(powers), len(freqs)))
    #phasedat = np.zeros((len(powers), len(freqs)))
    
    first_it = True
    Iblob_cspec = np.zeros(attempted_averages)
    Qblob_cspec = np.zeros(attempted_averages)
    powerdat_cspec = np.zeros(attempted_averages)
    phasedat_cspec = np.zeros(attempted_averages)
    if first_it:
        tstart = time.time()
        
    for avgind in range(attempted_averages):
        qubitgen.Power = Qbit_power
        
        total_samples = card.samples

        Is_full  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_full  = np.zeros((len(freqs), total_samples))
        
        Is_back = np.zeros(Is_full.shape)
        Qs_back = np.zeros(Qs_full.shape)

    
        # start from here         

        
        ##Acquire signal
        qubitgen.Freq = Qbit_freq
        qubitgen.output='On'
#            time.sleep(0.1)

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
        ##Useful handles for variables
        I_sig, Q_sig   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
        I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)] #<I>, <Q> for background trace
        theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
        theta_back = np.arctan2(Q_back, I_back)*180/np.pi #angle relative to x axis in IQ plane 
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back #compute <Q_net> in the data window
        I_full_net = I_full-I_full_b #full I data with background subtracted
        Q_full_net = Q_full-Q_full_b #full Q data with background subtracted
        
        # Store data
        Iblob_cspec[avgind] = I_final
        Qblob_cspec[avgind] = Q_final
        Iblob_cspec_avg = np.mean(Iblob_cspec)
        Qblob_cspec_avg = np.mean(Qblob_cspec)
        powerdat_cspec[avgind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_cspec[avgind] = np.arctan2(Q_final, I_final)*180/np.pi

    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
    r_cspec = find_pointsincircle(Iblob_cspec, Qblob_cspec, Iblob_cspec_avg, Qblob_cspec_avg, attempted_averages)
        
    ##Packaging data for the plotting functions and saving 
    cspec_data = {}
    cspec_data['xaxis'] = Qbit_freq/1e9
    cspec_data['mags'] = powerdat_cspec
    cspec_data['phases'] = phasedat_cspec
    cspec_data['radius'] = r_cspec

    if exp_settings['subtract_groundstate']:
        cspec_data['Iblob'] = Iblob_cspec - Iblob_ctrans_avg
        cspec_data['Qblob'] = Qblob_cspec - Qblob_ctrans_avg
        cspec_data['center'] = [Iblob_cspec_avg - Iblob_ctrans_avg, Qblob_cspec_avg - Qblob_ctrans_avg]
    else:
        cspec_data['Iblob'] = Iblob_cspec
        cspec_data['Qblob'] = Qblob_cspec
        cspec_data['center'] = [Iblob_cspec_avg, Qblob_cspec_avg]





    #######################################################################
    # Third, do carrier+sideband transmission scan
    #######################################################################

    ##Card settings
    configure_card(card, settings)
    
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

    # Make sure sideband signal is on for this carrier+sideband transmission scan
    sidebandgen.output = 'On'
    qubitgen.Output = 'Off'
    
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
    Iblob_cstrans = np.zeros(attempted_averages)
    Qblob_cstrans = np.zeros(attempted_averages)
    powerdat_cstrans = np.zeros(attempted_averages)
    phasedat_cstrans = np.zeros(attempted_averages)
    
    tstart = time.time()
    first_it = True

    if first_it:
        tstop = time.time()

    for avgind in range(attempted_averages):
        cavitygen.power = CAV_power
        sidebandgen.power = sideband_power
        
        time.sleep(0.2)

        total_samples = card.samples 
        Is  = np.zeros((len(CAV_freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs  = np.zeros((len(CAV_freqs), total_samples))
    
        cavitygen.freq = CAV_freq
        sidebandgen.freq = sideband_freq

        #LO.freq   = freq-exp_globals['IF']
        LO.freq = CAV_freq - exp_globals['IF']
        LO.output = 'On'
        
        cavitygen.phase = 0
        sidebandgen.phase = 0
        LO.phase = 0
        time.sleep(0.2)

        I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
                                                                        plot=first_it, 
                                                                        IQstorage = True)
        
        I_final = np.mean(I_window) #compute <I> in the data window
        Q_final = np.mean(Q_window) #compute <Q> in the data window
        

        
        # Store data
        Iblob_cstrans[avgind] = I_final
        Qblob_cstrans[avgind] = Q_final
        Iblob_cstrans_avg = np.mean(Iblob_cstrans)
        Qblob_cstrans_avg = np.mean(Qblob_cstrans)
        powerdat_cstrans[avgind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_cstrans[avgind] = np.arctan2(Q_final, I_final)*180/np.pi

    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
    r_cstrans = find_pointsincircle(Iblob_cstrans, Qblob_cstrans, Iblob_cstrans_avg, Qblob_cstrans_avg, attempted_averages)
            
    cstrans_data = {}
    cstrans_data['xaxis'] = CAV_freq/1e9
    cstrans_data['mags'] = powerdat_cstrans
    cstrans_data['phases'] = phasedat_cstrans
    cstrans_data['radius'] = r_cstrans

    if exp_settings['subtract_groundstate']:
        cstrans_data['Iblob'] = Iblob_cstrans - Iblob_ctrans_avg
        cstrans_data['Qblob'] = Qblob_cstrans - Qblob_ctrans_avg
        cstrans_data['center'] = [Iblob_cstrans_avg - Iblob_ctrans_avg, Qblob_cstrans_avg - Qblob_ctrans_avg]
    else:
        cstrans_data['Iblob'] = Iblob_cstrans
        cstrans_data['Qblob'] = Qblob_cstrans
        cstrans_data['center'] = [Iblob_cstrans_avg, Qblob_cstrans_avg]




    #######################################################################
    # Fourth, do carrier+sideband spec scan
    #######################################################################
    ## Setting qubit generator to some safe starting point before we turn it on
    qubitgen.Output  = 'On'
    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20
    
    if exp_settings['Quasi_CW']:
        qubitgen.disable_pulse()
        qubitgen.disable_IQ()
    else:
        qubitgen.enable_pulse()
        qubitgen.enable_IQ()
    
    ## Card config
    configure_card(card, settings)

    ## HDAWG settings
#    configure_hdawg(hdawg, settings)
    
    ## Sequencer program
    # progFile = open(r"C:\Users\kollarlab\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
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
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
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
#    qubitgen.disable_IQ()
#    qubit_marker.add_window(0,start_time+2e-6)
    ##
    cavity_marker.add_window(start_time, start_time+window_time)
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)
    
    # Make sure sideband signal is on for this carrier+sideband transmission scan
    sidebandgen.output = 'On'
    
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
    
    #powerdat = np.zeros((len(powers), len(freqs)))
    #phasedat = np.zeros((len(powers), len(freqs)))
    
    first_it = True
    Iblob_csspec = np.zeros(attempted_averages)
    Qblob_csspec = np.zeros(attempted_averages)
    powerdat_csspec = np.zeros(attempted_averages)
    phasedat_csspec = np.zeros(attempted_averages)
    if first_it:
        tstart = time.time()
        
    for avgind in range(attempted_averages):
        qubitgen.Power = Qbit_power
        
        total_samples = card.samples

        Is_full  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_full  = np.zeros((len(freqs), total_samples))
        
        Is_back = np.zeros(Is_full.shape)
        Qs_back = np.zeros(Qs_full.shape)

    
        # start from here         

        
        ##Acquire signal
        qubitgen.Freq = Qbit_freq
        qubitgen.output='On'
#            time.sleep(0.1)

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
        ##Useful handles for variables
        I_sig, Q_sig   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
        I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)] #<I>, <Q> for background trace
        theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
        theta_back = np.arctan2(Q_back, I_back)*180/np.pi #angle relative to x axis in IQ plane 
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back #compute <Q_net> in the data window
        I_full_net = I_full-I_full_b #full I data with background subtracted
        Q_full_net = Q_full-Q_full_b #full Q data with background subtracted
        
        # Store data
        Iblob_csspec[avgind] = I_final
        Qblob_csspec[avgind] = Q_final
        Iblob_csspec_avg = np.mean(Iblob_csspec)
        Qblob_csspec_avg = np.mean(Qblob_csspec)
        powerdat_csspec[avgind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_csspec[avgind] = np.arctan2(Q_final, I_final)*180/np.pi

    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
    r_csspec = find_pointsincircle(Iblob_csspec, Qblob_csspec, Iblob_csspec_avg, Qblob_csspec_avg, attempted_averages)
        
    ##Packaging data for the plotting functions and saving 
    csspec_data = {}
    csspec_data['xaxis'] = Qbit_freq/1e9
    csspec_data['mags'] = powerdat_csspec
    csspec_data['phases'] = phasedat_csspec
    csspec_data['radius'] = r_csspec

    if exp_settings['subtract_groundstate']:
        csspec_data['Iblob'] = Iblob_csspec - Iblob_ctrans_avg
        csspec_data['Qblob'] = Qblob_csspec - Qblob_ctrans_avg
        csspec_data['center'] = [Iblob_csspec_avg - Iblob_ctrans_avg, Qblob_csspec_avg - Qblob_ctrans_avg]
    else:
        csspec_data['Iblob'] = Iblob_csspec
        csspec_data['Qblob'] = Qblob_csspec
        csspec_data['center'] = [Iblob_csspec_avg, Qblob_csspec_avg]


    full_data = {}
    full_data['carrier_trans'] = ctrans_data
    full_data['carrier_spec'] = cspec_data
    full_data['carrier_and_sideband_trans'] = cstrans_data
    full_data['carrier_and_sideband_spec'] = csspec_data

    # Plot
    fig = plt.figure(figsize=(8,8))
    plt.clf()
    ax = plt.subplot(1,1,1)
    ax.scatter(ctrans_data['Iblob'], ctrans_data['Qblob'], s=4, c='k', label='carrier_trans')
    ax.scatter(cspec_data['Iblob'], cspec_data['Qblob'], s=4, c='g', label='carrier_spec')
    ax.scatter(cstrans_data['Iblob'], cstrans_data['Qblob'], s=4, c='b', label='carrier and sideband_trans')
    ax.scatter(csspec_data['Iblob'], csspec_data['Qblob'], s=4, c='r', label='carrier and sideband_spec')
    
    # circle_trans = Circle((Iblob_trans_avg, Qblob_trans_avg), radius=r_trans, color ='black', fill =False)
    # plt.gca().add_patch(circle_trans)
    # circle_trans_path = circle_trans.get_path()

    # circle_spec = Circle((Iblob_spec_avg, Qblob_spec_avg), radius=r_spec, color ='red', fill =False)
    # plt.gca().add_patch(circle_spec)
    # circle_path = circle_spec.get_path()

    plt.title('Modulation Freq(MHz): {}, Number of attempted shots: {}, meas time: {}'.format(exp_settings['mod_freq']/1e6, attempted_averages, window_time))
    plt.xlabel("I")
    plt.ylabel("Q")
    plt.axhline(y = 0, color ="black", linestyle ="-")
    plt.axvline(x = 0, color ="black", linestyle ="-") 
    y_max = np.max(np.abs(ax.get_ylim()))
    x_max = np.max(np.abs(ax.get_xlim()))
    list = np.array([x_max, y_max])
    max = np.max(list)
    ax.set_ylim(ymin=-max, ymax=max)
    ax.set_xlim(xmin=-max, xmax=max)
    plt.legend()
    plt.show()

    save_name = filename + '.png'
    userfuncs.savefig(fig, save_name, saveDir, png=True)

    # print averages result
    # if r_trans+r_spec > np.sqrt((np.abs(trans_data['center'][0] - spec_data['center'][0]))**2 + (np.abs(trans_data['center'][1] - spec_data['center'][1]))**2):
    #     print('Single shots should be higher. I,Q blobs overlap.')
    # else:
    #     print('Single Shots can be lower. I,Q blobs do not overlap.')

    # Save data
    userfuncs.SaveFull(saveDir, filename, ['powers', 'freqs', 'full_data'],
                        locals(), 
                        expsettings=settings, 
                        instruments=instruments, saveHWsettings=True)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart))
       
    cavitygen.Output = 'Off'
    sidebandgen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.output = 'Off'
    
    return full_data