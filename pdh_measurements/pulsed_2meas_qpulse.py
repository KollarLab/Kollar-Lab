'''
8-25-21 AK modifying to normalize the amplitudes to the drive power. Undid that.

9-2-21 AK made it return the data

'''


import os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import configure_card_2meas, configure_hdawg, estimate_time, read_and_process_singleshot, read_and_process_2meas
from pdh_measurements.scheduler_pdh import scheduler_pdh

import scipy.signal as signal

def get_default_settings():
    settings = {}
    
    settings['scanname'] = 'scanname'
    settings['meas_type'] = 'pulsed_2meas'
    
    #Cavity parameters
    settings['SGS_power']        = -30
    settings['CAV_Power']        = -60
    settings['CAV_freq']         = 7e9

    #Sideband parameters
    settings['SB_start_power']   = -60
    settings['SB_stop_power']    = -50
    settings['SB_power_points']  = 11
    settings['Mod_freq']         = 6.98e9
    settings['stimulation_time'] = 1e-6

    settings['qubit_pulse_enable'] = 'Off'
    
    settings['Qbit_power']      = -20
    settings['Qbit_freq']       = 5.6e9
    
    #Card settings
    settings['segments']         = 1e3
    settings['reads']            = 1
    settings['averages']         = 1
    
    #Measurement settings
    settings['Quasi_CW']         = False
    settings['reverse']          = False
    settings['num_save']         = 1
    
    #background_subtraction (by taking reference trace with no qubit drive power)
    settings['subtract_background']  = False
    settings['subtract_groundstate'] = True
    
    return settings



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

def find_pointsincircle(x, y, centerx, centery, samples):
    rxmax = np.max(x-centerx)
    rymax = np.max(y-centery)
    rmax = np.sqrt(rxmax**2 + rymax**2)
    rlist = np.linspace(0, rmax, num=samples)

    for i in range(samples):
        point_count = 0
        for px, py in zip(x, y):
            if (px-centerx)**2 + (py-centery)**2 <= rlist[i]**2:
                point_count += 1
        r = rlist[i]
        if point_count > samples*0.7609:
            break

    return r



def pulsed_2meas(instruments, settings):
    
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

    CAV_freq  = exp_settings['CAV_freq']
    CAV_power = exp_settings['CAV_Power']
    
    ##SB settings
    mod_freq        = exp_settings['Mod_freq']
    SB_start_power  = exp_settings['SB_start_power']
    SB_stop_power   = exp_settings['SB_stop_power']
    SB_power_points = exp_settings['SB_power_points']
    SB_powers = np.round(np.linspace(SB_start_power, SB_stop_power, SB_power_points),2)

    CAV_amp_list = np.zeros((1,3))
    SB_amp_list      = np.zeros((SB_power_points, 3))
    for i in range(1):
        CAV_amp_list[i] = get_amp_comps([-100, CAV_power, -100], exp_settings['SGS_power'])
    
    for i in range(SB_power_points):
        SB_amp_list[i] = get_amp_comps([SB_powers[i], -100, -100], exp_settings['SGS_power'])

    


    ##Qubit settings
    Qbit_freq  = exp_settings['Qbit_freq']
    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power = exp_settings['Qbit_power'] + Qbit_Attenuation
    
    ##Defining the number of single shots
    segments = int(exp_settings['segments'])


    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = SGS_power
    cavitygen.enable_IQ()
    cavitygen.enable_pulse()
    cavitygen.Output = 'On'
    
    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'

    qubitgen.freq = Qbit_freq
    qubitgen.power = Qbit_power
    qubitgen.enable_pulse()
    qubitgen.enable_IQ()
    qubitgen.Output = 'Off'


    ## exp_globals measurement settings
    m_pulse = exp_globals['measurement_pulse']
    q_pulse = exp_globals['qubit_pulse']
    
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    SB_delay    = m_pulse['SB_delay'] 
    cav2_delay  = m_pulse['meas2_delay'] 
    stm_time    = exp_settings['stimulation_time']

    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    cav_position = start_time
    SB_position = start_time + SB_delay
    cav2_position = start_time + cav2_delay
    qbit_position = start_time-delay-num_sigma*sigma-q_pulse['hold_time']





    ###################################################
    # No sideband measurement config
    ###################################################
    ## Card config
    configure_card_2meas(card, settings)

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

    
    if exp_settings['qubit_pulse_enable'] == 'Off' or exp_settings['qubit_pulse_enable'] == 'off':
        awg_sched = scheduler_pdh(total_time=start_time+3*window_time+cav2_delay, sample_rate=2.4e9)

        awg_sched.add_analog_channel(1, name='blank1')
        awg_sched.add_analog_channel(2, name='blank2')
        awg_sched.add_analog_channel(3, name='Cavity_I')
        awg_sched.add_analog_channel(4, name='Cavity_Q')
        
        awg_sched.add_digital_channel(1, name='blank3', polarity='Pos', HW_offset_on=0e-9, HW_offset_off=0e-9)
        awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        awg_sched.add_digital_channel(3, name='blank4', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        awg_sched.add_digital_channel(4, name='blank5', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        
        cavity_I      = awg_sched.analog_channels['Cavity_I']
        #cavity_Q      = awg_sched.analog_channels['Cavity_Q']
        cavity_marker = awg_sched.digital_channels['Cavity_enable']


        # first measurement pulse
        cavity_I.add_pulse('pdh_I', position=cav_position,
                            amp_list = CAV_amp_list[0], phase_list = [0,0,0],
                            mod_freq = 0,
                            time = window_time)
        
        # second measurement pulse
        cavity_I.add_pulse('pdh_I', position=cav2_position,
                    amp_list = CAV_amp_list[0], phase_list = [0,0,0],
                    mod_freq = 0,
                    time = window_time)

        # cavity marker on for whole time
        cavity_marker.add_window(start_time, start_time+cav2_delay+2*window_time)


        if first_it:
            awg_sched.plot_waveforms()
        
        [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['blank1', 'blank2'], ['blank3', 'Cavity_enable'])
        [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['blank4', 'blank5'])
        
        loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
        hdawg.AWGs[0].load_program(loadprog)
        hdawg.AWGs[0].stop()
        hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
        hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)
    
    else:
        awg_sched = scheduler_pdh(total_time=start_time+3*window_time+cav2_delay, sample_rate=2.4e9)

        awg_sched.add_analog_channel(1, name='Qubit_I')
        awg_sched.add_analog_channel(2, name='Qubit_Q')
        awg_sched.add_analog_channel(3, name='Cavity_I')
        awg_sched.add_analog_channel(4, name='Cavity_Q')
        
        awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0e-9, HW_offset_off=0e-9)
        awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        awg_sched.add_digital_channel(3, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        
        cavity_I      = awg_sched.analog_channels['Cavity_I']
        #cavity_Q      = awg_sched.analog_channels['Cavity_Q']
        cavity_marker = awg_sched.digital_channels['Cavity_enable']
        qubit_I       = awg_sched.analog_channels['Qubit_I']
        qubit_marker  = awg_sched.digital_channels['Qubit_enable']


        # Qubit pulse pi/2
        qubit_I.add_pulse('gaussian_square', position=qbit_position, 
                            amplitude=q_pulse['piAmp']/2, length = q_pulse['hold_time'], 
                            ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])

        # first measurement pulse
        cavity_I.add_pulse('pdh_I', position=cav_position,
                            amp_list = CAV_amp_list[0], phase_list = [0,0,0],
                            mod_freq = 0,
                            time = window_time)
        
        # second measurement pulse
        cavity_I.add_pulse('pdh_I', position=cav2_position,
                    amp_list = CAV_amp_list[0], phase_list = [0,0,0],
                    mod_freq = 0,
                    time = window_time)

        # cavity marker on for whole time
        cavity_marker.add_window(start_time, start_time+cav2_delay+2*window_time)
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
    # trans: SB off, spec: SB on
    Iblob_trans = np.zeros((SB_power_points, segments))
    Qblob_trans = np.zeros((SB_power_points, segments))
    powerdat_trans = np.zeros((SB_power_points, segments))
    phasedat_trans = np.zeros((SB_power_points, segments))
    Iblob_trans2 = np.zeros((SB_power_points, segments))
    Qblob_trans2 = np.zeros((SB_power_points, segments))
    powerdat_trans2 = np.zeros((SB_power_points, segments))
    phasedat_trans2 = np.zeros((SB_power_points, segments))

    Iblob_spec = np.zeros((SB_power_points, segments))
    Qblob_spec = np.zeros((SB_power_points, segments))
    powerdat_spec = np.zeros((SB_power_points, segments))
    phasedat_spec = np.zeros((SB_power_points, segments))
    Iblob_spec2 = np.zeros((SB_power_points, segments))
    Qblob_spec2 = np.zeros((SB_power_points, segments))
    powerdat_spec2 = np.zeros((SB_power_points, segments))
    phasedat_spec2 = np.zeros((SB_power_points, segments))


    
    for pind in range(SB_power_points):
        #######################################################################
        # First, do SB off
        #######################################################################
        print('SB off,  Current power:{}'.format(SB_powers[pind]))
        if exp_settings['qubit_pulse_enable'] == 'Off' or exp_settings['qubit_pulse_enable'] == 'off':
            qubitgen.Output  = 'Off'
        else:
            qubitgen.Output = 'On'
        
        
        ##Acquire signal
        total_samples = card.samples

        Is_full  = np.zeros((segments, total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_full  = np.zeros((segments, total_samples))
        
        Is_back = np.zeros(Is_full.shape)
        Qs_back = np.zeros(Qs_full.shape)

        cavitygen.phase = 0
        qubitgen.phase = 0
        LO.phase = 0

    
        # read and process
        # I,Q : dim [segment, time trace]
        I_window, Q_window, I_window2, Q_window2, I_full, Q_full, xaxis = read_and_process_2meas(card, settings, 
                                                                        plot=first_it, 
                                                                        IQstorage = True)
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
#                time.sleep(0.1)
            I_window_b, Q_window_b, I_window_b2, Q_window_b2, I_full_b, Q_full_b, xaxis_b = read_and_process_2meas(card, settings, 
                                                                plot=first_it, 
                                                                IQstorage = True)
        else:
            I_window_b, Q_window_b, I_window_b2, Q_window_b2, I_full_b, Q_full_b = 0,0,0,0,0,0
        
        ##Useful handles for variables
        #<I>, <Q> for signal time trace and get all segments
        I_sig, Q_sig   = [np.mean(I_window, axis=1), np.mean(Q_window, axis=1)] 
        I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)]
        I_sig2, Q_sig2   = [np.mean(I_window2, axis=1), np.mean(Q_window2, axis=1)] 
        I_back2, Q_back2 = [np.mean(I_window_b2), np.mean(Q_window_b2)]

        # theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
        # theta_back = np.arctan2(Q_back, I_back)*180/np.pi #angle relative to x axis in IQ plane 
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back #compute <Q_net> in the data window
        I_final2 = I_sig2-I_back2 #compute <I_net> in the data window
        Q_final2 = Q_sig2-Q_back2 #compute <Q_net> in the data window
        I_full_net = I_full-I_full_b #full I data with background subtracted
        Q_full_net = Q_full-Q_full_b #full Q data with background subtracted

    
        # Store data
        Iblob_trans[pind] = I_final
        Qblob_trans[pind] = Q_final
        powerdat_trans[pind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_trans[pind] = np.arctan2(Q_final, I_final)*180/np.pi

        Iblob_trans2[pind] = I_final2
        Qblob_trans2[pind] = Q_final2
        powerdat_trans2[pind] = np.sqrt(I_final2**2 + Q_final2**2)
        phasedat_trans2[pind] = np.arctan2(Q_final2, I_final2)*180/np.pi

        
        # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
        Iblob_trans_avg = np.mean(Iblob_trans, axis=1, keepdims=True)
        Qblob_trans_avg = np.mean(Qblob_trans, axis=1, keepdims=True)
        Iblob_trans_avg2 = np.mean(Iblob_trans2, axis=1, keepdims=True)
        Qblob_trans_avg2 = np.mean(Qblob_trans2, axis=1, keepdims=True)

        # r_trans = np.zeros((power_points,1))        
        # for i in range(power_points):
        #     r_trans[i][0] = find_pointsincircle(Iblob_trans[i], Qblob_trans[i], Iblob_trans_avg[i], Qblob_trans_avg[i], segments)
            
        ##Packaging data for the plotting functions and saving 
        trans_data = {}
        trans_data['xaxis'] = SB_powers
        trans_data['mags'] = powerdat_trans
        trans_data['phases'] = phasedat_trans
        trans_data['mags2'] = powerdat_trans2
        trans_data['phases2'] = phasedat_trans2
        #trans_data['radius'] = r_trans


        if exp_settings['subtract_groundstate']:
            trans_data['Iblob'] = Iblob_trans - Iblob_trans_avg
            trans_data['Qblob'] = Qblob_trans - Qblob_trans_avg
            trans_data['center'] = [Iblob_trans_avg - Iblob_trans_avg, Qblob_trans_avg - Qblob_trans_avg]
            trans_data['Iblob2'] = Iblob_trans2 - Iblob_trans_avg2
            trans_data['Qblob2'] = Qblob_trans2 - Qblob_trans_avg2
            trans_data['center2'] = [Iblob_trans_avg2 - Iblob_trans_avg2, Qblob_trans_avg2 - Qblob_trans_avg2]
        else:
            trans_data['Iblob'] = Iblob_trans
            trans_data['Qblob'] = Qblob_trans
            trans_data['center'] = [Iblob_trans_avg, Qblob_trans_avg]
            trans_data['Iblob2'] = Iblob_trans2
            trans_data['Qblob2'] = Qblob_trans2
            trans_data['center2'] = [Iblob_trans_avg2, Qblob_trans_avg2]
        
        if first_it:
            first_it = False
    
    time.sleep(1)





    ###################################################
    # Sideband measurement config; need to move configure_card_2meas and pdh scheduler inside of for loop to change the SB amp list
    ###################################################



    for pind in range(SB_power_points):
        #######################################################################
        # Second, do SB on using qubitgen
        #######################################################################
        ## Setting qubit generator to some safe starting point before we turn it on
        print('SB on,  Current power:{}'.format(SB_powers[pind]))
        if exp_settings['qubit_pulse_enable'] == 'Off' or exp_settings['qubit_pulse_enable'] == 'off':
            qubitgen.Output  = 'Off'
        else:
            qubitgen.Output = 'On'
        


        ## Card config
        configure_card_2meas(card, settings)

        ## Sequencer program
        # progFile = open(r"C:\Users\kollarlab\Documents\GitHub\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder.cpp",'r')
        # oxford comp sequencer file location is different
        progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
        rawprog  = progFile.read()
        loadprog = rawprog
        progFile.close()
        
        first_it = True
        
        if exp_settings['qubit_pulse_enable'] == 'Off' or exp_settings['qubit_pulse_enable'] == 'off':
            awg_sched = scheduler_pdh(total_time=start_time+3*window_time+cav2_delay, sample_rate=2.4e9)

            awg_sched.add_analog_channel(1, name='blank1')
            awg_sched.add_analog_channel(2, name='blank2')
            awg_sched.add_analog_channel(3, name='Cavity_I')
            awg_sched.add_analog_channel(4, name='Cavity_Q')
            
            awg_sched.add_digital_channel(1, name='blank3', polarity='Pos', HW_offset_on=0e-9, HW_offset_off=0e-9)
            awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(3, name='blank4', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(4, name='blank5', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            
            cavity_I      = awg_sched.analog_channels['Cavity_I']
            #cavity_Q      = awg_sched.analog_channels['Cavity_Q']
            cavity_marker = awg_sched.digital_channels['Cavity_enable']


            # first measurement pulse
            cavity_I.add_pulse('pdh_I', position=cav_position,
                                amp_list = CAV_amp_list[0], phase_list = [0,0,0],
                                mod_freq = 0,
                                time = window_time)
            
            # SB(stimulation) pulse
            cavity_I.add_pulse('pdh_I', position=SB_position,
                                amp_list = SB_amp_list[pind], phase_list = [0,0,0],
                                mod_freq = mod_freq,
                                time = window_time)        
            
            # second measurement pulse
            cavity_I.add_pulse('pdh_I', position=cav2_position,
                        amp_list = CAV_amp_list[0], phase_list = [0,0,0],
                        mod_freq = 0,
                        time = window_time)

            # cavity marker on for whole time
            cavity_marker.add_window(start_time, start_time+cav2_delay+2*window_time)


            if first_it:
                awg_sched.plot_waveforms()
            
            [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['blank1', 'blank2'], ['blank3', 'Cavity_enable'])
            [ch3, ch4, marker2] = awg_sched.compile_schedule('HDAWG', ['Cavity_I', 'Cavity_Q'], ['blank4', 'blank5'])
            
            loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
            hdawg.AWGs[0].load_program(loadprog)
            hdawg.AWGs[0].stop()
            hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
            hdawg.AWGs[1].load_waveform('0', ch3, ch4, marker2)
            hdawg.AWGs[0].run_loop()
            time.sleep(0.1)
        
        else:
            awg_sched = scheduler_pdh(total_time=start_time+3*window_time+cav2_delay, sample_rate=2.4e9)

            awg_sched.add_analog_channel(1, name='Qubit_I')
            awg_sched.add_analog_channel(2, name='Qubit_Q')
            awg_sched.add_analog_channel(3, name='Cavity_I')
            awg_sched.add_analog_channel(4, name='Cavity_Q')
            
            awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=0e-9, HW_offset_off=0e-9)
            awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(3, name='blank1', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            awg_sched.add_digital_channel(4, name='blank2', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
            
            cavity_I      = awg_sched.analog_channels['Cavity_I']
            #cavity_Q      = awg_sched.analog_channels['Cavity_Q']
            cavity_marker = awg_sched.digital_channels['Cavity_enable']
            qubit_I       = awg_sched.analog_channels['Qubit_I']
            qubit_marker  = awg_sched.digital_channels['Qubit_enable']


            # Qubit pulse pi/2
            qubit_I.add_pulse('gaussian_square', position=qbit_position, 
                                amplitude=q_pulse['piAmp']/2, length = q_pulse['hold_time'], 
                                ramp_sigma=q_pulse['sigma'], num_sigma=q_pulse['num_sigma'])

            # first measurement pulse
            cavity_I.add_pulse('pdh_I', position=cav_position,
                                amp_list = CAV_amp_list[0], phase_list = [0,0,0],
                                mod_freq = 0,
                                time = window_time)
            
            # SB(stimulation) pulse
            cavity_I.add_pulse('pdh_I', position=SB_position,
                                amp_list = SB_amp_list[pind], phase_list = [0,0,0],
                                mod_freq = mod_freq,
                                time = window_time)
            
            # second measurement pulse
            cavity_I.add_pulse('pdh_I', position=cav2_position,
                        amp_list = CAV_amp_list[0], phase_list = [0,0,0],
                        mod_freq = 0,
                        time = window_time)

            # cavity marker on for whole time
            cavity_marker.add_window(start_time, start_time+cav2_delay+2*window_time)
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




        ## Acquire signal
        total_samples = card.samples

        Is_full  = np.zeros((segments, total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_full  = np.zeros((segments, total_samples))
        
        Is_back = np.zeros(Is_full.shape)
        Qs_back = np.zeros(Qs_full.shape)

        cavitygen.phase = 0
        qubitgen.phase = 0

    
        # read and process
        # I,Q : dim [segment, time trace]
        I_window, Q_window, I_window2, Q_window2, I_full, Q_full, xaxis = read_and_process_2meas(card, settings, 
                                                                        plot=first_it, 
                                                                        IQstorage = True)
        if exp_settings['subtract_background']:
            #Acquire background trace
            qubitgen.output='Off'
#                time.sleep(0.1)
            I_window_b, Q_window_b, I_window_b2, Q_window_b2, I_full_b, Q_full_b, xaxis_b = read_and_process_2meas(card, settings, 
                                                                plot=first_it, 
                                                                IQstorage = True)
        else:
            I_window_b, Q_window_b, I_window_b2, Q_window_b2, I_full_b, Q_full_b = 0,0,0,0,0,0
        
        ##Useful handles for variables
        #<I>, <Q> for signal time trace and get all segments
        I_sig, Q_sig   = [np.mean(I_window, axis=1), np.mean(Q_window, axis=1)] 
        I_back, Q_back = [np.mean(I_window_b), np.mean(Q_window_b)]
        I_sig2, Q_sig2   = [np.mean(I_window2, axis=1), np.mean(Q_window2, axis=1)] 
        I_back2, Q_back2 = [np.mean(I_window_b2), np.mean(Q_window_b2)]

        # theta_sig  = np.arctan2(Q_sig,I_sig)*180/np.pi #angle relative to x axis in IQ plane
        # theta_back = np.arctan2(Q_back, I_back)*180/np.pi #angle relative to x axis in IQ plane 
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back #compute <Q_net> in the data window
        I_final2 = I_sig2-I_back2 #compute <I_net> in the data window
        Q_final2 = Q_sig2-Q_back2 #compute <Q_net> in the data window
        I_full_net = I_full-I_full_b #full I data with background subtracted
        Q_full_net = Q_full-Q_full_b #full Q data with background subtracted

    
        # Store data
        Iblob_spec[pind] = I_final
        Qblob_spec[pind] = Q_final
        powerdat_spec[pind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_spec[pind] = np.arctan2(Q_final, I_final)*180/np.pi

        Iblob_spec2[pind] = I_final2
        Qblob_spec2[pind] = Q_final2
        powerdat_spec2[pind] = np.sqrt(I_final2**2 + Q_final2**2)
        phasedat_spec2[pind] = np.arctan2(Q_final2, I_final2)*180/np.pi

        
        # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
        Iblob_spec_avg = np.mean(Iblob_spec, axis=1, keepdims=True)
        Qblob_spec_avg = np.mean(Qblob_spec, axis=1, keepdims=True)
        Iblob_spec_avg2 = np.mean(Iblob_spec2, axis=1, keepdims=True)
        Qblob_spec_avg2 = np.mean(Qblob_spec2, axis=1, keepdims=True)

        # r_trans = np.zeros((power_points,1))        
        # for i in range(power_points):
        #     r_trans[i][0] = find_pointsincircle(Iblob_trans[i], Qblob_trans[i], Iblob_trans_avg[i], Qblob_trans_avg[i], segments)
            
        ##Packaging data for the plotting functions and saving 
        spec_data = {}
        spec_data['xaxis'] = SB_powers
        spec_data['mags'] = powerdat_spec
        spec_data['phases'] = phasedat_spec
        spec_data['mags2'] = powerdat_spec2
        spec_data['phases2'] = phasedat_spec2
        #trans_data['radius'] = r_trans


        if exp_settings['subtract_groundstate']:
            spec_data['Iblob'] = Iblob_spec - Iblob_trans_avg
            spec_data['Qblob'] = Qblob_spec - Qblob_trans_avg
            spec_data['center'] = [Iblob_spec_avg - Iblob_trans_avg, Qblob_spec_avg - Qblob_trans_avg]
            spec_data['Iblob2'] = Iblob_spec2 - Iblob_trans_avg2
            spec_data['Qblob2'] = Qblob_spec2 - Qblob_trans_avg2
            spec_data['center2'] = [Iblob_spec_avg2 - Iblob_trans_avg2, Qblob_spec_avg2 - Qblob_trans_avg2]
        else:
            spec_data['Iblob'] = Iblob_spec
            spec_data['Qblob'] = Qblob_spec
            spec_data['center'] = [Iblob_spec_avg, Qblob_spec_avg]
            spec_data['Iblob2'] = Iblob_spec2
            spec_data['Qblob2'] = Qblob_spec2
            spec_data['center2'] = [Iblob_spec_avg2, Qblob_spec_avg2]
        
        if first_it:
            first_it = False



    ## Save full data
    full_data = {}
    full_data['trans'] = trans_data
    full_data['spec'] = spec_data



    ############################################
    # Plot total
    ############################################
    # arb = 100
    # title = 'Mod freq(MHz): ' + str(np.abs(mod_freq)/1e6) 
    
    # for i in range(SB_power_points):
    #     fig = plt.figure(num=arb+i, figsize=(8,8))
    #     plt.clf()
    #     ax = plt.subplot(1,1,1)
    #     ax.scatter(Iblob_trans[i], Qblob_trans[i], s=4, c='k', label='1st meas')
    #     ax.scatter(Iblob_trans2[i], Qblob_trans2[i], s=4, c='r', label='2nd meas')
        
    #     # circle_trans = Circle((Iblob_trans_avg, Qblob_trans_avg), radius=r_trans, color ='black', fill =False)
    #     # plt.gca().add_patch(circle_trans)
    #     # circle_trans_path = circle_trans.get_path()

    #     # circle_spec = Circle((Iblob_spec_avg, Qblob_spec_avg), radius=r_spec, color ='red', fill =False)
    #     # plt.gca().add_patch(circle_spec)
    #     # circle_path = circle_spec.get_path()

        
    #     plt.title('Without sideband, {}, Single_shots: {}, SB power: {}'.format(title, segments, SB_powers[i]))
    #     plt.xlabel("I")
    #     plt.ylabel("Q")
    #     plt.axhline(y = 0, color ="black", linestyle ="-")
    #     plt.axvline(x = 0, color ="black", linestyle ="-")
    #     y_max = np.max(np.abs(ax.get_ylim()))
    #     x_max = np.max(np.abs(ax.get_xlim()))
    #     list = np.array([x_max, y_max])
    #     max = np.max(list)
    #     ax.set_ylim(ymin=-max, ymax=max)
    #     ax.set_xlim(xmin=-max, xmax=max)
    #     plt.legend()
    #     plt.show()

    #     save_name = filename + '_noSB_' + str(SB_powers[i])
    #     userfuncs.savefig(fig, save_name, saveDir, png=True)
        
    # for i in range(SB_power_points):
    #     fig = plt.figure(num=2*arb+i, figsize=(8,8))
    #     plt.clf()
    #     ax = plt.subplot(1,1,1)
    #     ax.scatter(Iblob_spec[i], Qblob_spec[i], s=4, c='k', label='1st meas')
    #     ax.scatter(Iblob_spec2[i], Qblob_spec2[i], s=4, c='r', label='2nd meas')
        
    #     # circle_trans = Circle((Iblob_trans_avg, Qblob_trans_avg), radius=r_trans, color ='black', fill =False)
    #     # plt.gca().add_patch(circle_trans)
    #     # circle_trans_path = circle_trans.get_path()

    #     # circle_spec = Circle((Iblob_spec_avg, Qblob_spec_avg), radius=r_spec, color ='red', fill =False)
    #     # plt.gca().add_patch(circle_spec)
    #     # circle_path = circle_spec.get_path()

        
    #     plt.title('With sideband, {}, Single_shots: {}, SB power: {}'.format(title, segments, SB_powers[i]))
    #     plt.xlabel("I")
    #     plt.ylabel("Q")
    #     plt.axhline(y = 0, color ="black", linestyle ="-")
    #     plt.axvline(x = 0, color ="black", linestyle ="-") 
    #     y_max = np.max(np.abs(ax.get_ylim()))
    #     x_max = np.max(np.abs(ax.get_xlim()))
    #     list = np.array([x_max, y_max])
    #     max = np.max(list)
    #     ax.set_ylim(ymin=-max, ymax=max)
    #     ax.set_xlim(xmin=-max, xmax=max)
    #     plt.legend()
    #     plt.show()

    #     save_name = filename + '_SB_' + str(SB_powers[i])
    #     userfuncs.savefig(fig, save_name, saveDir, png=True)
        
        
    
    lines = (SB_power_points-1)//4 + 1

    find_max_nosb = np.array([np.max(np.abs(Iblob_trans)), np.max(np.abs(Qblob_trans)), np.max(np.abs(Iblob_trans2)), np.max(np.abs(Qblob_trans2))])
    max_nosb = np.max(find_max_nosb)
    find_max_sb = np.array([np.max(np.abs(Iblob_spec)), np.max(np.abs(Qblob_spec)), np.max(np.abs(Iblob_spec2)), np.max(np.abs(Qblob_spec2))])
    max_sb = np.max(find_max_sb)
    
    
    # First, the no SB
    fig = plt.figure(num=300, figsize=(14,12))
    fig.set_figwidth(14)
    plt.clf()
    
    gs = gridspec.GridSpec(nrows=lines, ncols=4, left=0.05, right=0.95, hspace=0.15, wspace=0.2, top=0.9, bottom=0.1)
    plt.suptitle('Without sideband, Single_shots: {}'.format(segments))
    
    for i in range(SB_power_points):
        ax = fig.add_subplot(gs[i//4:i//4+1, i%4:i%4+1])
        h = ax.hist2d(Iblob_trans2[i], Qblob_trans2[i], bins=200, range = [[-max_nosb, max_nosb],[-max_nosb, max_nosb]], cmap = 'hot')
        
        ax.set_title('SB power (No SB): {}'.format(SB_powers[i]))
        ax.set_xlabel("I")
        ax.set_ylabel("Q")
        ax.axhline(y = 0, color ="black", linestyle ="-")
        ax.axvline(x = 0, color ="black", linestyle ="-")
        ax.set_ylim(ymin=-max_nosb, ymax=max_nosb)
        ax.set_xlim(xmin=-max_nosb, xmax=max_nosb)
        #fig.colorbar(h[3], ax=ax, fraction=0.1)
        plt.gca().set_aspect('equal', adjustable='box') 
    
    plt.tight_layout()
    plt.show()
    save_name = filename + '_NoSB'
    userfuncs.savefig(fig, save_name, saveDir, png=True)
        
    
    # Second, the SB
    fig = plt.figure(num=301, figsize=(14,12))
    fig.set_figwidth(14)
    plt.clf()
    
    gs = gridspec.GridSpec(nrows=lines, ncols=4, left=0.05, right=0.95, hspace=0.15, wspace=0.2, top=0.9, bottom=0.1)
    plt.suptitle('With sideband, Mod: {} MHz, Single_shots: {}'.format(mod_freq/1e6, segments))
    
    for i in range(SB_power_points):
        ax = fig.add_subplot(gs[i//4:i//4+1, i%4:i%4+1])
        h = ax.hist2d(Iblob_spec2[i], Qblob_spec2[i], bins=200, range = [[-max_sb, max_sb],[-max_sb, max_sb]], cmap = 'hot')
        
        ax.set_title('SB power: {}'.format(SB_powers[i]))
        ax.set_xlabel("I")
        ax.set_ylabel("Q")
        ax.axhline(y = 0, color ="black", linestyle ="-")
        ax.axvline(x = 0, color ="black", linestyle ="-")
        ax.set_ylim(ymin=-max_sb, ymax=max_sb)
        ax.set_xlim(xmin=-max_sb, xmax=max_sb)
        #fig.colorbar(h[3], ax=ax, fraction=0.1)
        plt.gca().set_aspect('equal', adjustable='box') 
    
    plt.tight_layout()
    plt.show()
    save_name = filename + '_SB'
    userfuncs.savefig(fig, save_name, saveDir, png=True)
        
        
        
    # Save the full data
    userfuncs.SaveFull(saveDir, filename, ['CAV_power', 'SB_powers', 'full_data'],
                        locals(), 
                        expsettings=settings, 
                        instruments=instruments, saveHWsettings=True)
    t2 = time.time()
    
    print('elapsed time = ' + str(t2-tstart0))
       
    cavitygen.Output = 'Off'
    qubitgen.Output = 'Off'
    LO.output = 'Off'
    
    return full_data