import os
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.optimize import fsolve

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process
from pdh_measurements.scheduler_pdh import scheduler_pdh

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'pulsed_gefstate_SGS'
    
    #Cavity parameters
    settings['CAV_power']  = -53
    settings['SB_start_power'] = -50
    settings['SB_stop_power']  = -40
    settings['power_points']    = 5
    settings['SGS_power'] = -30
    
    settings['CAV_freq']        = 7e9
    settings['CAV_estate_freq'] = 6.99e9
    settings['mod_freq']        = 40e6

    settings['phase_list']      = [0,0,0]
    
    #Qubit parameters
    settings['Qbit_freq']       = 4.15e9
    settings['Qbit_power']      = -10  
    settings['freq_points']     = 1


    #Card settings
    settings['segments']         = 1
    settings['reads']            = 1
    settings['averages']         = 1
    settings['single_shots'] = int(2e3)
    
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
            

def total_power(power_list_dBm):
    # calculate the power in mW for each power in power_list_dBm
    mW_list = [10**(ind/10) for ind in power_list_dBm]
    # get total power in mW
    mW_sum = np.sum(mW_list)
    # convert mW power to dBm
    dBm = 10*np.log10(mW_sum)

    return dBm

def get_amp_comps(power_list_dBm, gen_dBm):

    # confirm if the sum of the components is less than the generator's
    gen_dBm = gen_dBm + 10
    total_comp_power = total_power(power_list_dBm)
    if total_comp_power > gen_dBm:
        raise ValueError('Sum of the components is greater than the available power. Use different combinations.')
    else:
        # calculate the corresponding power factor for each component / getting ratio in units of power
        power_factor = [10**(0.1*(pind - gen_dBm)) for pind in power_list_dBm]
        # now calculate the amplitudes list
        amp_list = [np.round(np.sqrt(ind),3) for ind in power_factor]

    return amp_list



def pulsed_gefstate_SGS_powercali(instruments, settings):
    
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

    
    ##Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    SGS_power       = exp_settings['SGS_power'] + CAV_Attenuation
    CAV_power       = exp_settings['CAV_power'] + CAV_Attenuation
    SB_start_power  = exp_settings['SB_start_power'] + CAV_Attenuation
    SB_stop_power   = exp_settings['SB_stop_power']  + CAV_Attenuation
    power_points    = exp_settings['power_points']
    SB_power        = np.round(np.linspace(SB_start_power,SB_stop_power,power_points),2)

    CAV_freq        = exp_settings['CAV_freq']
    CAV_estate_freq = exp_settings['CAV_estate_freq']
    CAV_freqs       = np.array([CAV_freq])
    mod_freq        = exp_settings['mod_freq']


    ##Carrier and Sideband AWG settings    
    carrier_amp_list = get_amp_comps([-100,CAV_power,-100], SGS_power)
    SB_amp_list      = np.zeros((power_points, 3))
    SB_power_list    = np.zeros((power_points, 3))
    phase_list       = exp_settings['phase_list']

    for i in range(power_points):
        SB_power_list[i][0] = SB_power[i]
        SB_power_list[i][1] = CAV_power
        SB_power_list[i][2] = SB_power[i]

        SB_amp_list[i] = get_amp_comps(SB_power_list[i], SGS_power)


    ##Qubit settings
    Qbit_freq  = exp_settings['Qbit_freq']
    freq_points = exp_settings['freq_points']
    freqs = np.array([Qbit_freq])

    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power = exp_settings['Qbit_power'] + Qbit_Attenuation
    

    single_shots = int(exp_settings['single_shots'])
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = SGS_power

    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_pulse()
    cavitygen.output = 'On'





    #######################################################################
    # First, do carrier spec scan with e state cav freq
    #######################################################################
    ## Setting qubit generator to some safe starting point before we turn it on
    print('Carrier with qubit at e state')
    qubitgen.Output  = 'On'
    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20
    
    # Park cavity at e state
    cavitygen.Freq   = CAV_estate_freq
    cavitygen.Power  = SGS_power
    cavitygen.Output = 'On'

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
    
    awg_sched = scheduler_pdh(total_time=start_time+2*window_time, sample_rate=2.4e9)

    awg_sched.add_analog_channel(1, name='Qubit_I')
    awg_sched.add_analog_channel(2, name='Qubit_Q')
    awg_sched.add_analog_channel(3, name='Cavity_I')
    awg_sched.add_analog_channel(4, name='Cavity_Q')
    
    awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0e-9, HW_offset_off=0e-9)
    awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(3, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
    
    cavity_I      = awg_sched.analog_channels['Cavity_I']
    cavity_Q      = awg_sched.analog_channels['Cavity_Q']
    cavity_marker = awg_sched.digital_channels['Cavity_enable']
    qubit_I       = awg_sched.analog_channels['Qubit_I']
    qubit_marker  = awg_sched.digital_channels['Qubit_enable']

    
    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    cav_position = start_time
    qbit_position = start_time-delay-num_sigma*sigma-q_pulse['hold_time']
    cavity_I.add_pulse('pdh_I', position=cav_position,
                           amp_list = carrier_amp_list, phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = window_time)

    cavity_Q.add_pulse('pdh_Q', position=cav_position,
                           amp_list = carrier_amp_list, phase_list = phase_list,
                           mod_freq = mod_freq,
                           time = window_time)
    
    qubit_I.add_pulse('gaussian_square', position=qbit_position, 
                              amplitude=q_pulse['piAmp'], length = q_pulse['hold_time'], 
                              ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
    
    cavity_marker.add_window(start_time, start_time+window_time)
    qubit_marker.add_window(qbit_position-160e-9, qbit_position+2*160e-9+q_pulse['hold_time'])
    
    
    awg_sched.plot_waveforms()
    
    [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
    [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['blank1', 'blank2'])
    
    loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].stop()
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
    
    #powerdat = np.zeros((len(powers), len(freqs)))
    #phasedat = np.zeros((len(powers), len(freqs)))
    
    first_it = True
    Iblob_cspec_e = np.zeros(single_shots)
    Qblob_cspec_e = np.zeros(single_shots)
    powerdat_cspec_e = np.zeros(single_shots)
    phasedat_cspec_e = np.zeros(single_shots)
    if first_it:
        tstart0 = time.time()
        


    for avgind in range(single_shots):
        ##Acquire signal
        qubitgen.Power = Qbit_power
        qubitgen.Freq = Qbit_freq
        qubitgen.output='On'

        total_samples = card.samples

        Is_full  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_full  = np.zeros((len(freqs), total_samples))
        
        Is_back = np.zeros(Is_full.shape)
        Qs_back = np.zeros(Qs_full.shape)

        cavitygen.phase = 0
        LO.phase = 0

    
        # start from here         
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
        Iblob_cspec_e[avgind] = I_final
        Qblob_cspec_e[avgind] = Q_final
        powerdat_cspec_e[avgind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_cspec_e[avgind] = np.arctan2(Q_final, I_final)*180/np.pi

    
    # Park cavity to the normal freq
    cavitygen.Freq   = CAV_freq
    cavitygen.Output = 'On'
    
    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
    Iblob_cspec_e_avg = np.mean(Iblob_cspec_e)
    Qblob_cspec_e_avg = np.mean(Qblob_cspec_e)
    r_cspec_e = find_pointsincircle(Iblob_cspec_e, Qblob_cspec_e, Iblob_cspec_e_avg, Qblob_cspec_e_avg, single_shots)
        
    ##Packaging data for the plotting functions and saving 
    cspec_e_data = {}
    cspec_e_data['xaxis'] = Qbit_freq/1e9
    cspec_e_data['mags'] = powerdat_cspec_e
    cspec_e_data['phases'] = phasedat_cspec_e
    cspec_e_data['radius'] = r_cspec_e


    if exp_settings['subtract_groundstate']:
        cspec_e_data['Iblob'] = Iblob_cspec_e - Iblob_cspec_e_avg
        cspec_e_data['Qblob'] = Qblob_cspec_e - Qblob_cspec_e_avg
        cspec_e_data['center'] = [Iblob_cspec_e_avg - Iblob_cspec_e_avg, Qblob_cspec_e_avg - Qblob_cspec_e_avg]
    else:
        cspec_e_data['Iblob'] = Iblob_cspec_e
        cspec_e_data['Qblob'] = Qblob_cspec_e
        cspec_e_data['center'] = [Iblob_cspec_e_avg, Qblob_cspec_e_avg]








    #######################################################################
    # Second, do carrier+sideband spec scan with e state cav freq
    #######################################################################
    ## Setting qubit generator to some safe starting point before we turn it on
    print('Carrier + Sideband with qubit at e state')
    qubitgen.Output  = 'On'
    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20

    # Park cavity at e state
    cavitygen.Power = SGS_power
    cavitygen.Freq   = CAV_estate_freq
    cavitygen.Output = 'On'
    
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
    Iblob_csspec_e = np.zeros((power_points, single_shots))
    Qblob_csspec_e = np.zeros((power_points, single_shots))
    powerdat_csspec_e = np.zeros((power_points, single_shots))
    phasedat_csspec_e = np.zeros((power_points, single_shots))

    if first_it:
        tstart1 = time.time()
        
    for pow_ind in range(power_points):
        print('Current CAV power:{}, Current SB power:{}, max:{}'.format(CAV_power-CAV_Attenuation, SB_power[pow_ind]-CAV_Attenuation, SB_power[-1]-CAV_Attenuation))

        awg_sched = scheduler_pdh(total_time=start_time+2*window_time, sample_rate=2.4e9)

        awg_sched.add_analog_channel(1, name='Qubit_I')
        awg_sched.add_analog_channel(2, name='Qubit_Q')
        awg_sched.add_analog_channel(3, name='Cavity_I')
        awg_sched.add_analog_channel(4, name='Cavity_Q')
        
        awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0e-9, HW_offset_off=0e-9)
        awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        awg_sched.add_digital_channel(3, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        
        cavity_I      = awg_sched.analog_channels['Cavity_I']
        cavity_Q      = awg_sched.analog_channels['Cavity_Q']
        cavity_marker = awg_sched.digital_channels['Cavity_enable']
        qubit_I       = awg_sched.analog_channels['Qubit_I']
        qubit_marker  = awg_sched.digital_channels['Qubit_enable']

        
        delay = q_pulse['delay']
        sigma = q_pulse['sigma']
        num_sigma = q_pulse['num_sigma']
        cav_position = start_time
        qbit_position = start_time-delay-num_sigma*sigma-q_pulse['hold_time']
        cavity_I.add_pulse('pdh_I', position=cav_position,
                            amp_list = SB_amp_list[pow_ind], phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = window_time)

        cavity_Q.add_pulse('pdh_Q', position=cav_position,
                            amp_list = SB_amp_list[pow_ind], phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = window_time)
        
        qubit_I.add_pulse('gaussian_square', position=qbit_position, 
                                amplitude=q_pulse['piAmp'], length = q_pulse['hold_time'], 
                                ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])
        
        cavity_marker.add_window(start_time, start_time+window_time)
        qubit_marker.add_window(qbit_position-160e-9, qbit_position+2*160e-9+q_pulse['hold_time'])
        
        if first_it:
            awg_sched.plot_waveforms()

        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['blank1', 'blank2'])
        
        loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
        hdawg.AWGs[0].load_program(loadprog)
        hdawg.AWGs[0].stop()
        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)

        for avgind in range(single_shots):
            ##Acquire signal
            qubitgen.Power = Qbit_power
            qubitgen.Freq = Qbit_freq
            qubitgen.output='On'
            
            total_samples = card.samples

            Is_full  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
            Qs_full  = np.zeros((len(freqs), total_samples))
            
            Is_back = np.zeros(Is_full.shape)
            Qs_back = np.zeros(Qs_full.shape)

            cavitygen.Phase = 0
            qubitgen.Phase = 0
        
            # start from here         
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
            Iblob_csspec_e[pow_ind, avgind] = I_final
            Qblob_csspec_e[pow_ind, avgind] = Q_final
            powerdat_csspec_e[pow_ind, avgind] = np.sqrt(I_final**2 + Q_final**2)
            phasedat_csspec_e[pow_ind, avgind] = np.arctan2(Q_final, I_final)*180/np.pi

    # Park cavity to the normal freq
    cavitygen.Freq   = CAV_freq
    cavitygen.Output = 'On'

    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
    Iblob_csspec_e_avg = np.mean(Iblob_csspec_e, axis=1)
    Qblob_csspec_e_avg = np.mean(Qblob_csspec_e, axis=1)
    r_csspec_e = np.zeros((power_points,1))
    for i in range(power_points):
        r_csspec_e[i][0] = find_pointsincircle(Iblob_csspec_e[i], Qblob_csspec_e[i], Iblob_csspec_e_avg[i], Qblob_csspec_e_avg[i], single_shots)
        
    ##Packaging data for the plotting functions and saving 
    csspec_e_data = {}
    csspec_e_data['xaxis'] = Qbit_freq/1e9
    csspec_e_data['mags'] = powerdat_csspec_e
    csspec_e_data['phases'] = phasedat_csspec_e
    csspec_e_data['radius'] = r_csspec_e

    if exp_settings['subtract_groundstate']:
        csspec_e_data['Iblob'] = Iblob_csspec_e - Iblob_cspec_e_avg
        csspec_e_data['Qblob'] = Qblob_csspec_e - Qblob_cspec_e_avg
        csspec_e_data['center'] = [Iblob_csspec_e_avg - Iblob_cspec_e_avg, Qblob_csspec_e_avg - Qblob_cspec_e_avg]
    else:
        csspec_e_data['Iblob'] = Iblob_csspec_e
        csspec_e_data['Qblob'] = Qblob_csspec_e
        csspec_e_data['center'] = [Iblob_csspec_e_avg, Qblob_csspec_e_avg]



    full_data = {}
    full_data['carrier_spec_estate'] = cspec_e_data
    full_data['carrier_and_sideband_spec_estate'] = csspec_e_data



    ############################################
    # Plot total : plot all does not make a good visualization
    ############################################
    
    # fig = plt.figure(num=31, figsize=(8,8))
    # plt.clf()
    # ax = plt.subplot(1,1,1)
    # ax.scatter(Iblob_cspec_e, Qblob_cspec_e, s=4, c='k', label='carrier_spec_estate')
    # ax.scatter(Iblob_csspec_e[0], Qblob_csspec_e[0], s=4, c='r', label='carrier and sideband_spec_estate')
    
    # # circle_trans = Circle((Iblob_trans_avg, Qblob_trans_avg), radius=r_trans, color ='black', fill =False)
    # # plt.gca().add_patch(circle_trans)
    # # circle_trans_path = circle_trans.get_path()

    # # circle_spec = Circle((Iblob_spec_avg, Qblob_spec_avg), radius=r_spec, color ='red', fill =False)
    # # plt.gca().add_patch(circle_spec)
    # # circle_path = circle_spec.get_path()

    # plt.title('Modulation Freq(MHz): {}, Number of single shots: {}, meas time: {}'.format(exp_settings['mod_freq']/1e6, single_shots, window_time))
    # plt.xlabel("I")
    # plt.ylabel("Q")
    # plt.axhline(y = 0, color ="black", linestyle ="-")
    # plt.axvline(x = 0, color ="black", linestyle ="-") 
    # y_max = np.max(np.abs(ax.get_ylim()))
    # x_max = np.max(np.abs(ax.get_xlim()))
    # list = np.array([x_max, y_max])
    # max = np.max(list)
    # ax.set_ylim(ymin=-max, ymax=max)
    # ax.set_xlim(xmin=-max, xmax=max)
    # plt.legend()
    # plt.show()

    # save_name = filename + '_total'
    # userfuncs.savefig(fig, save_name, saveDir, png=True)


    # Plot fstate plot
    arb = 100
    for i in range(power_points):
        fig = plt.figure(num=arb+i, figsize=(8,8))
        plt.clf()
        ax = plt.subplot(1,1,1)
        ax.scatter(Iblob_cspec_e, Qblob_cspec_e, s=4, c='k', label='carrier_spec_estate')
        ax.scatter(Iblob_csspec_e[i], Qblob_csspec_e[i], s=4, c='r', label='carrier and sideband_spec_estate')
        
        # circle_trans = Circle((Iblob_trans_avg, Qblob_trans_avg), radius=r_trans, color ='black', fill =False)
        # plt.gca().add_patch(circle_trans)
        # circle_trans_path = circle_trans.get_path()

        # circle_spec = Circle((Iblob_spec_avg, Qblob_spec_avg), radius=r_spec, color ='red', fill =False)
        # plt.gca().add_patch(circle_spec)
        # circle_path = circle_spec.get_path()

        plt.title('Modulation Freq(MHz): {}, Single_shots: {}, Sideband drive power: {}'.format(exp_settings['mod_freq']/1e6, single_shots, SB_power[i]))
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

        save_name = filename + '_fstate_' + 'SBpower: {}'.format(SB_power[i]-CAV_Attenuation)
        userfuncs.savefig(fig, save_name, saveDir, png=True)




    # print averages result
    # if r_trans+r_spec > np.sqrt((np.abs(trans_data['center'][0] - spec_data['center'][0]))**2 + (np.abs(trans_data['center'][1] - spec_data['center'][1]))**2):
    #     print('Single shots should be higher. I,Q blobs overlap.')
    # else:
    #     print('Single Shots can be lower. I,Q blobs do not overlap.')

    # Save data
    userfuncs.SaveFull(saveDir, filename, ['SGS_power', 'CAV_power', 'SB_power', 'full_data'],
                        locals(), 
                        expsettings=settings, 
                        instruments=instruments, saveHWsettings=True)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart0))
       
    cavitygen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.output = 'Off'
    
    return full_data