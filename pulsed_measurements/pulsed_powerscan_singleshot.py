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



def pulsed_powerscan_singleshot(instruments, settings):
    
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
    CAV_freq_list    = exp_settings['CAV_freq_list'] 
    CAV_freqs        = np.round(np.linspace(CAV_freq-10e6, CAV_freq+5e6, 101),-3)
    mod_freq         = exp_settings['mod_freq']


    ##Carrier and Sideband AWG settings    
    carrier_amp_list = np.zeros((CAV_power_points,3))
    SB_amp_list      = np.zeros((SB_power_points, 3))
    phase_list       = exp_settings['phase_list']
    
    for i in range(CAV_power_points):
        carrier_amp_list[i] = get_amp_comps([-100,CAV_power[i],-100], exp_settings['SGS_power'])
    
    for i in range(SB_power_points):
        SB_amp_list[i] = get_amp_comps([SB_power[i], CAV_power[0], SB_power[i]], exp_settings['SGS_power'])


    ##Qubit settings
    Qbit_freq  = exp_settings['Qbit_freq']
    freq_points = exp_settings['freq_points']
    freqs = np.array([Qbit_freq])

    Qbit_Attenuation = exp_globals['Qbit_Attenuation']
    Qbit_power = exp_settings['Qbit_power'] + Qbit_Attenuation
    
    ##Defining the number of single shots
    segments = int(exp_settings['segments'])
    
    ## Generator settings
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = SGS_power

    LO.power = 12
    LO.freq = CAV_freq - exp_globals['IF']    
    LO.output = 'On'
    
    cavitygen.enable_IQ()
    cavitygen.enable_pulse()
    cavitygen.output = 'On'




    #######################################################################
    # 0th, do cavity freq search for carrier power scan
    #######################################################################
    # if exp_settings['sidebandonoff'] == 'Off':
    #     cavitygen.freq   = CAV_freqs[0]
    #     cavitygen.power  = CAV_power[0]
        
    #     cavitygen.enable_pulse()
    #     cavitygen.output = 'On'
    
    #     LO.power  = 12
    #     LO.freq   = CAV_freqs[0] - exp_globals['IF']
    #     LO.output = 'On'
        
    #     ##Card settings
    #     configure_card(card, settings)
    
        
    #     progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\pulsed_measurements\HDAWG_sequencer_codes\hdawg_placeholder_4channels.cpp",'r')
    #     rawprog  = progFile.read()
    #     loadprog = rawprog
    #     progFile.close()
        
    #     m_pulse = exp_globals['measurement_pulse']
    #     start_time  = m_pulse['meas_pos']
    #     window_time = m_pulse['meas_window']
        
        
    #     awg_sched = scheduler(total_time=start_time+2*window_time, sample_rate=2.4e9)
        
    #     awg_sched.add_analog_channel(1, name='Qubit_I')
    #     awg_sched.add_analog_channel(2, name='Qubit_Q')
        
    #     awg_sched.add_digital_channel(1, name='Qubit_enable', polarity='Pos', HW_offset_on=55e-9, HW_offset_off=0e-9)
    #     awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        
    #     cavity_marker = awg_sched.digital_channels['Cavity_enable']
        
    #     cavity_marker.add_window(start_time, start_time+window_time)
        
    #     awg_sched.plot_waveforms()
        
    #     [ch1, ch2, marker] = awg_sched.compile_schedule('HDAWG', ['Qubit_I', 'Qubit_Q'], ['Qubit_enable', 'Cavity_enable'])
        
    #     loadprog = loadprog.replace('_samples_', str(awg_sched.samples))
    #     hdawg.AWGs[0].load_program(loadprog)
    #     hdawg.AWGs[0].load_waveform('0', ch1, ch2, marker)
    #     hdawg.AWGs[0].run_loop()
    #     time.sleep(0.1)
        
    #     ##create the digital down conversion filter if needed.
    #     if exp_globals['IF'] != 0:
    #         #create Chebychev type II digital filter
    #         filter_N = exp_globals['ddc_config']['order']
    #         filter_rs = exp_globals['ddc_config']['stop_atten']
    #         filter_cutoff = np.abs(exp_globals['ddc_config']['cutoff'])
    #         LPF = signal.cheby2(filter_N, filter_rs, filter_cutoff, btype='low', analog=False, output='sos', fs=card.sampleRate)
            
    #         xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
    #         digLO_sin = np.sin(2*np.pi*exp_globals['IF']*xaxis)
    #         digLO_cos = np.cos(2*np.pi*exp_globals['IF']*xaxis)
            
    #         #store in settings so that the processing functions can get to them
    #         settings['digLO_sin'] = digLO_sin 
    #         settings['digLO_cos'] = digLO_cos
    #         settings['LPF'] = LPF
    
    #     ## Start main measurement loop 
    #     powerdat = np.zeros((len(powers), len(freqs)))
    #     phasedat = np.zeros((len(powers), len(freqs)))
        
    #     tstart = time.time()
    #     first_it = True
    
    #     drive_powers_lin = 10**(powers/10)
    #     drive_amps_lin = np.sqrt(drive_powers_lin)
    
    #     print('asdf')
    #     for powerind in range(len(powers)):
    #         cavitygen.power = CAV_power[powerind]
            
    
    #         total_samples = card.samples 
            
    #         Is  = np.zeros((len(freqs), total_samples)) #the very rawest data is thrown away for heterodyne! Warning
    #         Qs  = np.zeros((len(freqs), total_samples))
            
    #         print('Current power:{}, max:{}'.format(CAV_power[powerind] - CAV_Attenuation, CAV_power[-1] - CAV_Attenuation))
        
    #         for find in range(0, len(freqs)):
    #             freq = CAV_freqs[find]
    
    #             cavitygen.freq = freq
    #             LO.freq = freq - exp_globals['IF']
    #             LO.output = 'On'
                
    #             cavitygen.phase = 0
    #             LO.phase = 0
    
    #             I_window, Q_window, I_full, Q_full, xaxis = read_and_process(card, settings, 
    #                                                                          plot=first_it, 
    #                                                                          IQstorage = True)
                
    #             I_final = np.mean(I_window) #compute <I> in the data window
    #             Q_final = np.mean(Q_window) #compute <Q> in the data window
                
    #             if first_it:
    #                 tstop = time.time()
    #                 estimate_time(tstart, tstop, len(powers)*len(freqs))
    #                 first_it = False
                    
                

    #             Is[find,:] = I_full
    #             Qs[find,:] = Q_full
    #             powerdat[powerind, find] = np.sqrt(I_final**2 + Q_final**2)
    #             phasedat[powerind, find] = np.arctan2(Q_final, I_final)*180/np.pi
                
    #         full_data_CAV = {}
    #         full_data_CAV['xaxis']  = freqs/1e9
    #         full_data_CAV['mags']   = powerdat[0:powerind+1]
    #         full_data_CAV['phases'] = phasedat[0:powerind+1]
            
    #         #rescale to fractional amplitude for the plots.
    #         plot_data = {}
    #         plot_data['xaxis']  = freqs/1e9
    #         plot_data['mags']   = np.transpose(np.transpose(powerdat[0:powerind+1])/drive_amps_lin[0:powerind+1])
    #         plot_data['phases'] = phasedat[0:powerind+1]
    
    #         single_data = {}
    #         single_data['xaxis'] = freqs/1e9
    #         single_data['mag']   = powerdat[powerind]
    #         single_data['phase'] = phasedat[powerind]
    
    #         yaxis  = powers[0:powerind+1] - CAV_Attenuation
    #         labels = ['Freq (GHz)', 'Power (dBm)']

    #         simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=1, IQdata = False) #normalized to drive level
    #         plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
    

    #         full_time = {}
    #         full_time['xaxis']  = xaxis*1e6
    #         full_time['Is']   = Is
    #         full_time['Qs']   = Qs
    

    #         single_time = {}
    #         single_time['xaxis'] = xaxis*1e6
    #         single_time['I']   = I_full
    #         single_time['Q'] = Q_full
    
    #         time_labels = ['Time (us)', 'Freq (GHz)']
    #         identifier = 'Power: {}dBm'.format(powers[powerind]-CAV_Attenuation)

    #         simplescan_plot(full_time, single_time, freqs/1e9, 
    #                         'Raw_time_traces\n'+filename, 
    #                         time_labels, 
    #                         identifier, 
    #                         fig_num=2,
    #                         IQdata = True)
    #         plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces.png'), dpi = 150)
    
    #         userfuncs.SaveFull(saveDir, filename, ['powers','freqs', 'powerdat', 'phasedat','xaxis','full_data_CAV', 'single_data', 'full_time', 'single_time'],
    #                              locals(), expsettings=settings, instruments=instruments, saveHWsettings=first_it)
    #     t2 = time.time()
        
    #     print('elapsed time = ' + str(t2-tstart))


    #     ## Calculate CAV_freq_list
    #     for i in range(CAV_power_points):
    #         CAV_freq_list[i] = 1e9 * full_data









    #######################################################################
    # First, do transmission scan
    #######################################################################
    ## Setting generators to some safe starting point before we turn it on
    print('Pulsed Trans')
    qubitgen.Output  = 'Off'
    
    cavitygen.Freq   = CAV_freq
    cavitygen.Power  = SGS_power
    cavitygen.Output = 'On'
    cavitygen.enable_pulse()
    cavitygen.enable_IQ()
    
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
    start_time  = m_pulse['meas_pos']
    window_time = m_pulse['meas_window']
    
    cav_position = start_time
    
    
    ## create the digital down conversion filter if needed.
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
    
    
    ## scanning power
    if exp_settings['sidebandonoff'] == 'On':
        powers = SB_power
        power_points = SB_power_points
        
    if exp_settings['sidebandonoff'] == 'Off':
        powers = CAV_power
        power_points = CAV_power_points
    
    
    first_it = True
    Iblob_trans = np.zeros((power_points, segments))
    Qblob_trans = np.zeros((power_points, segments))
    powerdat_trans = np.zeros((power_points, segments))
    phasedat_trans = np.zeros((power_points, segments))
    
    
    if first_it:
        tstart0 = time.time()

    for pind in range(power_points):
        print('Current power:{}'.format(powers[pind]))
        if exp_settings['sidebandonoff'] == 'Off':
            cavitygen.Freq = CAV_freq_list[pind]
        
        awg_sched = scheduler_pdh(total_time=start_time+2*window_time, sample_rate=2.4e9)

        awg_sched.add_analog_channel(1, name='blank1')
        awg_sched.add_analog_channel(2, name='blank2')
        awg_sched.add_analog_channel(3, name='Cavity_I')
        awg_sched.add_analog_channel(4, name='Cavity_Q')
        
        awg_sched.add_digital_channel(1, name='blank3', polarity='Pos', HW_offset_on=0e-9, HW_offset_off=0e-9)
        awg_sched.add_digital_channel(2, name='Cavity_enable', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        awg_sched.add_digital_channel(3, name='blank4', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        awg_sched.add_digital_channel(4, name='blank5', polarity='Pos', HW_offset_on=0, HW_offset_off=0)
        
        cavity_I      = awg_sched.analog_channels['Cavity_I']
        cavity_Q      = awg_sched.analog_channels['Cavity_Q']
        cavity_marker = awg_sched.digital_channels['Cavity_enable']

        
        
        if exp_settings['sidebandonoff'] == 'On':
            amp_list = SB_amp_list[pind]
        if exp_settings['sidebandonoff'] == 'Off':
            amp_list = carrier_amp_list[pind]
        
        cavity_I.add_pulse('pdh_I', position=cav_position,
                               amp_list = amp_list, phase_list = phase_list,
                               mod_freq = mod_freq,
                               time = window_time)

        cavity_Q.add_pulse('pdh_Q', position=cav_position,
                               amp_list = amp_list, phase_list = phase_list,
                               mod_freq = mod_freq,
                               time = window_time)
        
        cavity_marker.add_window(start_time, start_time+window_time)

        
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
        
        
        ##Acquire signal
        total_samples = card.samples

        Is_full  = np.zeros((segments, total_samples)) #the very rawest data is thrown away for heterodyne! Warning
        Qs_full  = np.zeros((segments, total_samples))
        
        Is_back = np.zeros(Is_full.shape)
        Qs_back = np.zeros(Qs_full.shape)

        cavitygen.phase = 0
        LO.phase = 0

    
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
        Iblob_trans[pind] = I_final
        Qblob_trans[pind] = Q_final
        powerdat_trans[pind] = np.sqrt(I_final**2 + Q_final**2)
        phasedat_trans[pind] = np.arctan2(Q_final, I_final)*180/np.pi

    
    # Park cavity to the normal freq
    cavitygen.Freq   = CAV_freq
    cavitygen.Output = 'On'
    
    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
    Iblob_trans_avg = np.mean(Iblob_trans, axis=1, keepdims=True)
    Qblob_trans_avg = np.mean(Qblob_trans, axis=1, keepdims=True)
    r_trans = np.zeros((power_points,1))
    for i in range(power_points):
        r_trans[i][0] = find_pointsincircle(Iblob_trans[i], Qblob_trans[i], Iblob_trans_avg[i], Qblob_trans_avg[i], segments)
        
    ##Packaging data for the plotting functions and saving 
    trans_data = {}
    trans_data['xaxis'] = powers
    trans_data['sidebandonoff'] = exp_settings['sidebandonoff']
    trans_data['mags'] = powerdat_trans
    trans_data['phases'] = phasedat_trans
    trans_data['radius'] = r_trans


    if exp_settings['subtract_groundstate']:
        trans_data['Iblob'] = Iblob_trans - Iblob_trans_avg
        trans_data['Qblob'] = Qblob_trans - Qblob_trans_avg
        trans_data['center'] = [Iblob_trans_avg - Iblob_trans_avg, Qblob_trans_avg - Qblob_trans_avg]
    else:
        trans_data['Iblob'] = Iblob_trans
        trans_data['Qblob'] = Qblob_trans
        trans_data['center'] = [Iblob_trans_avg, Qblob_trans_avg]








    #######################################################################
    # Second, do carrier+sideband spec scan with e state cav freq
    #######################################################################
    ## Setting qubit generator to some safe starting point before we turn it on
    print('Pulsed spec')
    qubitgen.Output  = 'On'
    qubitgen.Freq   = 4e9
    qubitgen.Power  = -20

    # Park cavity at e state
    cavitygen.Power = SGS_power
    cavitygen.Freq   = CAV_freq
    cavitygen.Output = 'On'
    cavitygen.enable_pulse()
    cavitygen.enable_IQ()
    
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

    delay = q_pulse['delay']
    sigma = q_pulse['sigma']
    num_sigma = q_pulse['num_sigma']
    cav_position = start_time
    qbit_position = start_time-delay-num_sigma*sigma-q_pulse['hold_time']
    
    
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
    
    
    
    if exp_settings['sidebandonoff'] == 'On':
        powers = SB_power
        power_points = SB_power_points
        
    if exp_settings['sidebandonoff'] == 'Off':
        powers = CAV_power
        power_points = CAV_power_points
    
    first_it = True
    Iblob_spec = np.zeros((power_points, segments))
    Qblob_spec = np.zeros((power_points, segments))
    powerdat_spec = np.zeros((power_points, segments))
    phasedat_spec = np.zeros((power_points, segments))

    if first_it:
        tstart1 = time.time()
        
    for pind in range(power_points):
        print('Current power:{}'.format(powers[pind]))
        if exp_settings['sidebandonoff'] == 'On':
            cavitygen.Freq = CAV_freq_list[pind]

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
        
        
        if exp_settings['sidebandonoff'] == 'On':
            amp_list = SB_amp_list[pind]
        if exp_settings['sidebandonoff'] == 'Off':
            amp_list = carrier_amp_list[pind]

    
        cavity_I.add_pulse('pdh_I', position=cav_position,
                            amp_list = amp_list, phase_list = phase_list,
                            mod_freq = mod_freq,
                            time = window_time)

        cavity_Q.add_pulse('pdh_Q', position=cav_position,
                            amp_list = amp_list, phase_list = phase_list,
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

    # Park cavity to the normal freq
    cavitygen.Freq   = CAV_freq
    cavitygen.Output = 'On'

    # Make FWHM gaussian circle, count by numbers of points in the gaussian circle. The normalized gaussian function contains erf(ln(2))~0.7609 points in the FWHM.
    Iblob_spec_avg = np.mean(Iblob_spec, axis=1, keepdims=True)
    Qblob_spec_avg = np.mean(Qblob_spec, axis=1, keepdims=True)
    r_spec = np.zeros((power_points,1))
    for i in range(power_points):
        r_spec[i][0] = find_pointsincircle(Iblob_spec[i], Qblob_spec[i], Iblob_spec_avg[i], Qblob_spec_avg[i], segments)
        
    ##Packaging data for the plotting functions and saving 
    spec_data = {}
    spec_data['xaxis'] = powers
    spec_data['sidebandonoff'] = exp_settings['sidebandonoff']
    spec_data['mags'] = powerdat_spec
    spec_data['phases'] = phasedat_spec
    spec_data['radius'] = r_spec

    if exp_settings['subtract_groundstate']:
        spec_data['Iblob'] = Iblob_spec - Iblob_trans_avg
        spec_data['Qblob'] = Qblob_spec - Qblob_trans_avg
        spec_data['center'] = [Iblob_spec_avg - Iblob_trans_avg, Qblob_spec_avg - Qblob_trans_avg]
    else:
        spec_data['Iblob'] = Iblob_spec
        spec_data['Qblob'] = Qblob_spec
        spec_data['center'] = [Iblob_spec_avg, Qblob_spec_avg]



    ## Save full data
    full_data = {}
    full_data['trans'] = trans_data
    full_data['spec'] = spec_data



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


    # Plot fstate plot
    arb = 100
    if exp_settings['sidebandonoff'] == 'On':
        letter = 'sideband'
    if exp_settings['sidebandonoff'] == 'Off':
        letter = 'carrier'
    
    for i in range(power_points):
        fig = plt.figure(num=arb+i, figsize=(8,8))
        plt.clf()
        ax = plt.subplot(1,1,1)
        ax.scatter(Iblob_trans[i], Qblob_trans[i], s=4, c='k', label='g')
        ax.scatter(Iblob_spec[i], Qblob_spec[i], s=4, c='r', label='e')
        
        # circle_trans = Circle((Iblob_trans_avg, Qblob_trans_avg), radius=r_trans, color ='black', fill =False)
        # plt.gca().add_patch(circle_trans)
        # circle_trans_path = circle_trans.get_path()

        # circle_spec = Circle((Iblob_spec_avg, Qblob_spec_avg), radius=r_spec, color ='red', fill =False)
        # plt.gca().add_patch(circle_spec)
        # circle_path = circle_spec.get_path()

        
        plt.title('Modulation Freq(MHz): {}, Single_shots: {}, Scanning {} power: {}'.format(exp_settings['mod_freq']/1e6, segments, letter, powers[i]))
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

        save_name = filename + letter + '_power scan_' + str(powers[i])
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