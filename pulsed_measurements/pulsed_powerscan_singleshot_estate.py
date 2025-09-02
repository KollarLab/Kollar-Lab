import os
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.optimize import fsolve

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card, configure_hdawg, estimate_time, read_and_process, read_and_process_singleshot
from pdh_measurements.scheduler_pdh import scheduler_pdh
from utility.scheduler import scheduler

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'pulsed_powerscan'
    
    #Cavity parameters
    settings['CAV_start_power']  = -53
    settings['CAV_stop_power']   = -50
    settings['CAV_power_points'] = 5
    settings['SB_start_power']   = -50
    settings['SB_stop_power']    = -40
    settings['SB_power_points']  = 5
    settings['SGS_power']        = -30
    settings['sidebandonoff']    = 'Off' 
    
    settings['CAV_freq']        = 7e9
    settings['chi_shift']       = 1e6
    settings['mod_freq']        = 40e6

    settings['phase_list']      = [0,0,0]
    
    #Qubit parameters
    settings['Qbit_freq']       = 4.15e9
    settings['Qbit_power']      = -10  
    settings['freq_points']     = 1


    #Card settings
    settings['segments']         = 2e3
    settings['reads']            = 1
    settings['averages']         = 1
    
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
    total_comp_power = total_power(power_list_dBm)
    if total_comp_power > gen_dBm:
        raise ValueError('Sum of the components is greater than the available power. Use different combinations.')
    else:
        # calculate the corresponding power factor for each component / getting ratio in units of power
        power_factor = [10**(0.1*(pind - gen_dBm)) for pind in power_list_dBm]
        # now calculate the amplitudes list
        amp_list = [np.round(np.sqrt(ind),3) for ind in power_factor]

    return amp_list



def pulsed_powerscan_singleshot_estate(instruments, settings):
    
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
    CAV_Attenuation  = exp_globals['CAV_Attenuation']
    EnableIQ_loss    = exp_globals['EnableIQ_loss']
    SGS_power        = exp_settings['SGS_power'] + EnableIQ_loss + CAV_Attenuation
    CAV_start_power  = exp_settings['CAV_start_power'] 
    CAV_stop_power   = exp_settings['CAV_stop_power'] 
    SB_start_power   = exp_settings['SB_start_power'] 
    SB_stop_power    = exp_settings['SB_stop_power']
    CAV_power_points = exp_settings['CAV_power_points']
    SB_power_points  = exp_settings['SB_power_points']
    CAV_power        = np.round(np.linspace(CAV_start_power,CAV_stop_power,CAV_power_points),2)
    SB_power         = np.round(np.linspace(SB_start_power,SB_stop_power,SB_power_points),2)

    CAV_freq         = exp_settings['CAV_freq']
    #CAV_freq_list    = exp_settings['CAV_freq_list'] 
    chi              = exp_settings['chi_shift']
    CAV_e_freq       = CAV_freq - chi
    #CAV_e_freq_list  = CAV_freq_list - chi
    mod_freq         = exp_settings['mod_freq']


    ##Carrier and Sideband AWG settings    
    carrier_amp_list = np.zeros((CAV_power_points,3))
    SB_amp_list      = np.zeros((SB_power_points, 3))
    phase_list       = exp_settings['phase_list']
    
    for i in range(CAV_power_points):
        carrier_amp_list[i] = get_amp_comps([-100,CAV_power[i],-100], exp_settings['SGS_power'])
    
    for i in range(SB_power_points):
        SB_amp_list[i] = get_amp_comps([SB_power[i], CAV_power[0], SB_power[i]], exp_settings['SGS_power'])


    ## carrier / sideband settings 
    if exp_settings['sidebandonoff'] == 'On':
        powers = SB_power
        power_points = SB_power_points
        
    if exp_settings['sidebandonoff'] == 'Off':
        raise ValueError('Turn on the sidebands')

            

    ##Qubit settings
    Qbit_freq  = exp_settings['Qbit_freq']
    freq_points = exp_settings['freq_points']
    freqs = np.array([Qbit_freq])

    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power = exp_settings['Qbit_power'] + Qbit_Attenuation
    
    ##Defining the number of single shots
    segments = int(exp_settings['segments'])
    
    ## Generator settings
    cavitygen.Freq   = CAV_e_freq
    cavitygen.Power  = SGS_power
    cavitygen.enable_IQ()
    cavitygen.enable_pulse()
    cavitygen.output = 'On'

    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'


    if exp_settings['Quasi_CW']:
        qubitgen.disable_pulse()
        qubitgen.disable_IQ()
    else:
        qubitgen.enable_pulse()
        qubitgen.enable_IQ()


    ## exp_globals measurement settings
    m_pulse = exp_globals['measurement_pulse']
    q_pulse = exp_globals['qubit_pulse']
    
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']

    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    cav_position = start_time
    qbit_position = start_time-delay-num_sigma*sigma-q_pulse['hold_time']



    ## Digital Filter settings
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


    ## Set full_data
    Iblob_trans = np.zeros((power_points, segments))
    Qblob_trans = np.zeros((power_points, segments))
    powerdat_trans = np.zeros((power_points, segments))
    phasedat_trans = np.zeros((power_points, segments))

    Iblob_spec = np.zeros((power_points, segments))
    Qblob_spec = np.zeros((power_points, segments))
    powerdat_spec = np.zeros((power_points, segments))
    phasedat_spec = np.zeros((power_points, segments))





    for pind in range(power_points):
        #######################################################################
        # First, spec at e state with carrier
        #######################################################################
        ## Setting generators to some safe starting point before we turn it on
        print('Current power:{}'.format(powers[pind]))

        qubitgen.Output  = 'On'
        qubitgen.Freq   = 4e9
        qubitgen.Power  = -20
        
        cavitygen.Freq = CAV_e_freq
        LO.freq = CAV_e_freq - exp_globals['IF']
        cavitygen.Power  = SGS_power
        cavitygen.enable_pulse()
        cavitygen.enable_IQ()
        cavitygen.Output = 'On'
        
        
        ## Card config
        configure_card(card, settings)
        
        ## Sequencer program
        # progFile = open(r"C:\Users\kollarlab\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
        # oxford comp sequencer file location is different
        progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
        rawprog  = progFile.read()
        loadprog = rawprog
        progFile.close()
        
        
        first_it = True
        if first_it:
            tstart0 = time.time()

        
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

    
        cavity_I.add_pulse('pdh_I', position=cav_position,
                            amp_list = carrier_amp_list[0], phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = window_time)

        cavity_Q.add_pulse('pdh_Q', position=cav_position,
                            amp_list = carrier_amp_list[0], phase_list = phase_list,
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

    
        
        
        ##Acquire signal
        qubitgen.Power  = Qbit_power
        qubitgen.Freq   = Qbit_freq
        qubitgen.output ='On'


        total_samples = card.samples

        Is_full  = np.zeros((segments, total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_full  = np.zeros((segments, total_samples))
        
        Is_back = np.zeros(Is_full.shape)
        Qs_back = np.zeros(Qs_full.shape)

        cavitygen.Phase = 0
        qubitgen.Phase = 0
        LO.Phase = 0

    
        # start from here         
        I_window, Q_window, I_full, Q_full, xaxis = read_and_process_singleshot(card, settings, 
                                                                        plot=first_it, 
                                                                        IQstorage = True)
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
#                time.sleep(0.1)
            I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process_singleshot(card, settings, 
                                                                plot=first_it, 
                                                                IQstorage = True)
        else:
            I_window_b, Q_window_b, I_full_b, Q_full_b = 0,0,0,0
        
        ##Useful handles for variables
        I_sig, Q_sig   = [np.mean(I_window, axis=1), np.mean(Q_window, axis=1)] #<I>, <Q> for signal trace
        I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)] #<I>, <Q> for background trace
        theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
        theta_back = np.arctan2(Q_back, I_back)*180/np.pi #angle relative to x axis in IQ plane 
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back #compute <Q_net> in the data window
        I_full_net = I_full-I_full_b #full I data with background subtracted
        Q_full_net = Q_full-Q_full_b #full Q data with background subtracted
    
        # Store data
        Iblob_trans[pind] = I_final
        Qblob_trans[pind] = Q_final
        powerdat_trans[pind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_trans[pind] = np.arctan2(Q_final, I_final)*180/np.pi

        
        # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
        Iblob_trans_avg = np.mean(Iblob_trans, axis=1, keepdims=True)
        Qblob_trans_avg = np.mean(Qblob_trans, axis=1, keepdims=True)
        r_trans = np.zeros((power_points,1))
        for i in range(power_points):
            r_trans[i][0] = find_pointsincircle(Iblob_trans[i], Qblob_trans[i], Iblob_trans_avg[i], Qblob_trans_avg[i], segments)
            
        ##Packaging data for the plotting functions and saving 
        cspec_data = {}
        cspec_data['xaxis'] = powers
        cspec_data['sidebandonoff'] = exp_settings['sidebandonoff']
        cspec_data['mags'] = powerdat_trans
        cspec_data['phases'] = phasedat_trans
        cspec_data['radius'] = r_trans


        if exp_settings['subtract_groundstate']:
            cspec_data['Iblob'] = Iblob_trans - Iblob_trans_avg
            cspec_data['Qblob'] = Qblob_trans - Qblob_trans_avg
            cspec_data['center'] = [Iblob_trans_avg - Iblob_trans_avg, Qblob_trans_avg - Qblob_trans_avg]
        else:
            cspec_data['Iblob'] = Iblob_trans
            cspec_data['Qblob'] = Qblob_trans
            cspec_data['center'] = [Iblob_trans_avg, Qblob_trans_avg]





        #######################################################################
        # Second, spec at e state with carrier+sideband
        #######################################################################
        ## Setting qubit generator to some safe starting point before we turn it on
        qubitgen.Output  = 'On'
        qubitgen.Freq   = 4e9
        qubitgen.Power  = -20
        

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

    
        cavity_I.add_pulse('pdh_I', position=cav_position,
                            amp_list = SB_amp_list[pind], phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = window_time)

        cavity_Q.add_pulse('pdh_Q', position=cav_position,
                            amp_list = SB_amp_list[pind], phase_list = phase_list,
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



        ##Acquire signal
        qubitgen.Power  = Qbit_power
        qubitgen.Freq   = Qbit_freq
        qubitgen.output ='On'
        
        total_samples = card.samples

        Is_full  = np.zeros((segments, total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_full  = np.zeros((segments, total_samples))
        
        Is_back = np.zeros(Is_full.shape)
        Qs_back = np.zeros(Qs_full.shape)

        cavitygen.Phase = 0
        qubitgen.Phase = 0
        LO.Phase = 0
    
        # start from here         
        I_window, Q_window, I_full, Q_full, xaxis = read_and_process_singleshot(card, settings, 
                                                                        plot=first_it, 
                                                                        IQstorage = True)
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
#                time.sleep(0.1)
            I_window_b, Q_window_b, I_full_b, Q_full_b, xaxis_b = read_and_process_singleshot(card, settings, 
                                                                plot=first_it, 
                                                                IQstorage = True)
        else:
            I_window_b, Q_window_b, I_full_b, Q_full_b = 0,0,0,0
        
        if first_it:
            first_it=False
        ##Useful handles for variables
        I_sig, Q_sig   = [np.mean(I_window, axis=1), np.mean(Q_window, axis=1)] #<I>, <Q> for signal trace
        I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)] #<I>, <Q> for background trace
        theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
        theta_back = np.arctan2(Q_back, I_back)*180/np.pi #angle relative to x axis in IQ plane 
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back #compute <Q_net> in the data window
        I_full_net = I_full-I_full_b #full I data with background subtracted
        Q_full_net = Q_full-Q_full_b #full Q data with background subtracted
        
        # Store data
        Iblob_spec[pind] = I_final
        Qblob_spec[pind] = Q_final
        powerdat_spec[pind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_spec[pind] = np.arctan2(Q_final, I_final)*180/np.pi


        # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
        Iblob_spec_avg = np.mean(Iblob_spec, axis=1, keepdims=True)
        Qblob_spec_avg = np.mean(Qblob_spec, axis=1, keepdims=True)
        r_spec = np.zeros((power_points,1))
        for i in range(power_points):
            r_spec[i][0] = find_pointsincircle(Iblob_spec[i], Qblob_spec[i], Iblob_spec_avg[i], Qblob_spec_avg[i], segments)
            
        ##Packaging data for the plotting functions and saving 
        csspec_data = {}
        csspec_data['xaxis'] = powers
        csspec_data['sidebandonoff'] = exp_settings['sidebandonoff']
        csspec_data['mags'] = powerdat_spec
        csspec_data['phases'] = phasedat_spec
        csspec_data['radius'] = r_spec

        if exp_settings['subtract_groundstate']:
            csspec_data['Iblob'] = Iblob_spec - Iblob_trans_avg
            csspec_data['Qblob'] = Qblob_spec - Qblob_trans_avg
            csspec_data['center'] = [Iblob_spec_avg - Iblob_trans_avg, Qblob_spec_avg - Qblob_trans_avg]
        else:
            csspec_data['Iblob'] = Iblob_spec
            csspec_data['Qblob'] = Qblob_spec
            csspec_data['center'] = [Iblob_spec_avg, Qblob_spec_avg]



    ## Save full data
    full_data = {}
    full_data['cspec'] = cspec_data
    full_data['csspec'] = csspec_data



    ############################################
    # Plot total
    ############################################
    
    # fig = plt.figure(num=31, figsize=(8,8))
    # plt.clf()
    # ax = plt.subplot(1,1,1)
    # ax.scatter(Iblob_trans[0], Qblob_trans[0], s=4, c='k', label='carrier_spec_estate')
    # ax.scatter(Iblob_spec[0], Qblob_spec_e[0], s=4, c='r', label='carrier and sideband_spec_estate')
    
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


    # Plot
    arb = 100
    
    for i in range(power_points):
        fig = plt.figure(num=arb+i, figsize=(8,8))
        plt.clf()
        ax = plt.subplot(1,1,1)
        ax.scatter(Iblob_trans[i], Qblob_trans[i], s=4, c='k', label='carrier')
        ax.scatter(Iblob_spec[i], Qblob_spec[i], s=4, c='r', label='carrier+sideband')
        
        # circle_trans = Circle((Iblob_trans_avg, Qblob_trans_avg), radius=r_trans, color ='black', fill =False)
        # plt.gca().add_patch(circle_trans)
        # circle_trans_path = circle_trans.get_path()

        # circle_spec = Circle((Iblob_spec_avg, Qblob_spec_avg), radius=r_spec, color ='red', fill =False)
        # plt.gca().add_patch(circle_spec)
        # circle_path = circle_spec.get_path()

        
        plt.title('Modulation Freq(MHz): {}, Single_shots: {}, Scanning SB power: {}'.format(exp_settings['mod_freq']/1e6, segments, powers[i]))
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

        save_name = filename + '_power scan_estate_' + str(powers[i])
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