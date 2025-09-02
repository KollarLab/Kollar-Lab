# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:11:31 2024

@author: Kollarlab
"""
import os
import time
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, freqz, filtfilt, periodogram
import pickle
from datetime import datetime
from datetime import date

import userfuncs as uf


# Filter functions
def butter_lowpass(cutoff, fs, order):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order):
        b, a = butter_lowpass(cutoff, fs, order)
        y = lfilter(b, a, data)
        return y
    
def load_hdawg(I_data, Q_data, instruments):
    '''
    This should probably take in what hdawg is defined as.
    Here it assumes hdawg is loaded and defined as 'hdawg'
    '''
    hdawg = instruments["AWG"]
    I = I_data
    Q = Q_data
    progFile = open(r"K:\Data\Oxford\PDH\DataAcquisitionCode\hdawg_placeholder_hold.cpp", 'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    loadprog = loadprog.replace('_samples_', str(len(I)))
    hdawg.AWGs[0].load_program(loadprog) 
    hdawg.AWGs[0].load_waveform('0', I, Q, Q)  
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)
    
def modulating_IQ(carrier_amp, sidebandM_amp, sidebandP_amp, carrier_phase, sidebandM_phase, sidebandP_phase, modulation_freq, time):
    I = carrier_amp * np.cos(carrier_phase) + sidebandP_amp * np.cos(2*np.pi*modulation_freq*time + sidebandP_phase) + sidebandM_amp * np.cos(2*np.pi*modulation_freq*time - sidebandM_phase)
    Q = -carrier_amp * np.sin(carrier_phase) - sidebandP_amp * np.sin(2*np.pi*modulation_freq*time + sidebandP_phase) + sidebandM_amp * np.sin(2*np.pi*modulation_freq*time - sidebandM_phase)
    return I, Q

def card_data(rf_range, modulation_frequency, instruments, amp_list = [], phase_list = [], Normalize = True, detector = ''):
    gen = instruments["Rf_gen"]
    hdawg = instruments["AWG"]
    card = instruments["digitizer"]

    #Base settings for IQ waveform
    '''
    Move the following to exp globals
    '''
    sample_rate = 2.4e9
    total_time = 20e-6
    buffer = 1e-6
    hold = 5e-6

    #######################################################################
    modulation_freq = modulation_frequency 
    sidebandM_amp, carrier_amp, sidebandP_amp = [amp_list[0], amp_list[1], amp_list[2]]
    sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])
    if Normalize == True:
        sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])/(sidebandM_amp + carrier_amp + sidebandP_amp)
    sidebandM_phase, carrier_phase, sidebandP_phase = [phase_list[0]*np.pi, phase_list[1]*np.pi, phase_list[2]*np.pi]

    ##################################################################################

    # pulsed modulation setup
    time = np.linspace(0, total_time, int(total_time*sample_rate))
    time_hold = np.linspace(0, hold, int(hold*sample_rate))
    I_t, Q_t = modulating_IQ(carrier_amp, sidebandM_amp, sidebandP_amp, carrier_phase, sidebandM_phase, sidebandP_phase, modulation_freq, time_hold)
    I = np.zeros(len(time))
    I[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate))* I_t
    Q = np.zeros(len(time))
    Q[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate)) * Q_t 
    
    # Load IQ into HDAWG
    load_hdawg(I, Q, instruments)


    ########################## Get trimmed data
    trim_buffer = 0.15e-6
    trim_buffer_start = 3000
    trim_buffer_end = 6000
    #trim_buffer_start = int((buffer + trim_buffer) * card.sampleRate)
    #trim_buffer_end = int((hold - trim_buffer) * card.sampleRate)

    gen.output = 1 # THIS SHOULD BE IN EXP_GLOBALS or something
    RF = rf_range
    data = np.zeros((len(rf_range), int(card.samples)))
    data_trimmed = np.zeros((len(rf_range), trim_buffer_end - trim_buffer_start))
    
    ###########
    #data = []
    #data_trimmed = []
    ###########
    
    for ind in range(len(RF)):
        gen.freq = RF[ind] 
        card.ArmAndWait()
        I, Q = card.ReadAllData()
        if detector == 'envelope_detector':
            Q_trimmed = Q[0][trim_buffer_start:trim_buffer_end]
            I_trimmed = I[0][trim_buffer_start:trim_buffer_end]
            data_trimmed[ind] = Q_trimmed - I_trimmed
            data[ind] = Q[0]- I[0]
        elif detector == 'diode_detector':
            '''
            This assumes the diode detector is connected to the Q channel on card.
            Should find a way to make the right input channel something you give the code
            '''
            Q_trimmed = Q[0][trim_buffer_start:trim_buffer_end]
            data_trimmed[ind] = Q_trimmed 
            data[ind] = Q[0] 

    # digital downconversion on data from detector
    card_data = data
    card_data_trimmed = data_trimmed
    data_pdh_mean_im = np.zeros(len(card_data_trimmed))
    data_pdh_mean_re = np.zeros(len(card_data_trimmed))

    # analyze data for each RF-LO combination
    for k in range(len(card_data_trimmed)):
        reflected_power = card_data_trimmed[k]

        time_app = np.linspace(trim_buffer_start/card.sampleRate, trim_buffer_end/card.sampleRate, len(card_data_trimmed[k]))
        dc_re = np.cos(2*np.pi*modulation_freq*time_app) 
        dc_im = np.sin(2*np.pi*modulation_freq*time_app) 

        # downconvert signal 
        reflected_power_re = reflected_power * dc_re 
        reflected_power_im = reflected_power * dc_im 

        # second filter 
        order = 5
        cutoff = modulation_freq/10 # about 2MHz cutoff

        filter_re = butter_lowpass_filter(reflected_power_re, cutoff, card.sampleRate, order)
        filter_im = butter_lowpass_filter(reflected_power_im, cutoff, card.sampleRate, order)

    
        data_pdh_mean_im[k] = np.mean(filter_im)
        data_pdh_mean_re[k] = np.mean(filter_re)

    final_pdh_mean_im = data_pdh_mean_im
    final_pdh_mean_re = data_pdh_mean_re
    
    
    return card_data, card_data_trimmed, final_pdh_mean_im, final_pdh_mean_re, trim_buffer_start, trim_buffer_end

def modFreq_linewidth(settings, instruments):

    
    gen = instruments["Rf_gen"]
    hdawg = instruments["AWG"]
    card = instruments["digitizer"]
    
    
    mod_mult = settings['mod_mult']
    cavity_linewidth = settings['cavity_linewidth']
    cavity_freq = settings['cavity_freq']
    plot_data = settings['plot_data']
    power = settings['power']
    gen.power = power
    freq_points = settings['freq_points']
    # phase_factor = settings['phase_factor']
    
    

    amp_list = settings['amp_list']
    phase_list = settings['phase_list']
    Normalize = settings['Normalize']
    detector = settings['detector']
    span = settings['narrow_span']

    stamp = timestamp()
    path = saveDir(settings) 
    scanname = saveFilename(settings) 
    full_data = []

    col = 3
    row = len(mod_mult)
    fig = plt.figure(1)
    plt.clf()

    filename = scanname

    for k in range(len(mod_mult)):

        modulation_freq = mod_mult[k] * cavity_linewidth 
        if span:
            start_freq = cavity_freq - 15e6#1.5*modulation_freq #7.835e9
            stop_freq = cavity_freq + 15e6 # 1.5*modulation_freq #7.855e9
        else:
            start_freq = cavity_freq - 1.5*modulation_freq #7.835e9
            stop_freq = cavity_freq + 1.5*modulation_freq #7.855e9
        rf_range = np.linspace(start_freq, stop_freq, freq_points)

        init_data, init_trimmed, imag_comp, real_comp, trim_buffer_start, trim_buffer_end = card_data(rf_range, modulation_freq, instruments, amp_list, phase_list, Normalize, detector)
        amp_comp = np.sqrt(real_comp**2 + imag_comp**2)

        single_dat = {}
        single_dat['info'] = 'mod_freq {} data'.format(k+1)
        single_dat['raw_full'] = init_data
        single_dat['raw_trimmed'] = init_trimmed
        single_dat['imag_comp'] = imag_comp
        single_dat['real_comp'] = real_comp
        single_dat['amp_comp'] = amp_comp
        single_dat['mod_freq'] = modulation_freq
        single_dat['xaxis'] = rf_range
        single_dat['trim_buffer_start'] = trim_buffer_start
        single_dat['trim_buffer_end'] = trim_buffer_end
        single_dat['cavity_freq'] = cavity_freq

        full_data.insert(k, single_dat)

        imag_max = np.max(np.abs(imag_comp))
        real_max = np.max(np.abs(real_comp))

        # Plot
        if plot_data:

            ax = plt.subplot(row, col, (3*k)+1)
            plt.plot(rf_range/1e9,imag_comp, 'k', label='imag comp')
            plt.ylabel('V_pdh')
            plt.ylim(-1.2*imag_max, 1.2*imag_max)
            plt.xlabel('frequency (GHz)')
            plt.grid()
            plt.title("ModFreq={:.2f}MHz".format(modulation_freq/1e6))
            plt.legend()

            ax = plt.subplot(row, col, (3*k)+2)
            plt.plot(rf_range/1e9,real_comp, 'b', label = "real comp")
            plt.ylabel('V_pdh')
            plt.ylim(-1.2*real_max, 1.2*real_max)
            plt.xlabel('frequency (GHz)')
            plt.grid()
            plt.title("ModFreq={:.2f}MHz".format(modulation_freq/1e6))
            plt.legend()

            ax = plt.subplot(row, col, (3*k)+3)
            plt.plot(rf_range/1e9, amp_comp, 'b', label = "amp comp")
            plt.ylabel('V_pdh')
            plt.xlabel('frequency (GHz)')
            plt.grid()
            plt.title("ModFreq={:.2f}MHz".format(modulation_freq/1e6))
            plt.legend()
            
    if plot_data:
        plot_title = filename.split("_")
        plot_title = "_".join(plot_title) + stamp
        fig.suptitle(plot_title, fontsize=10)
        fig.set_size_inches([15,10])
        fig.tight_layout()
        plt.show()
     
        fig.savefig(os.path.join(path, filename + stamp + '.jpg'), dpi = 150)

    SaveFull(path, filename+stamp, full_data, expsettings=settings, instruments = instruments, saveHWsettings=True)
    gen.output = 0
    '''
    Should plot the last data points so we see if the data is trimmed right
    '''
    first_trimmed = np.zeros(len(init_data[0]))
    first_trimmed[trim_buffer_start:trim_buffer_end] = init_trimmed[0]
    
    last_trimmed = np.zeros(len(init_data[-1]))
    last_trimmed[trim_buffer_start:trim_buffer_end] = init_trimmed[-1]

    fig = plt.figure(2)
    plt.clf()
    ax = plt.subplot(1,2,1)
    plt.plot(init_data[0], color = 'deepskyblue', label='first-full')
    plt.plot(first_trimmed, color = 'orange', label='first-trimmed')
    ax.legend()
    ax = plt.subplot(1,2,2)
    plt.plot(init_data[-1], color = 'deepskyblue', label='last-full')
    plt.plot(last_trimmed, color = 'orange', label='last-trimmed')
    ax.legend()
    
    fig.suptitle('raw traces')
    fig.tight_layout()
    plt.show()
    
    return init_data, init_trimmed

def timestamp():
    date  = datetime.now()
    stamp = date.strftime('%Y%m%d_%H%M%S')
    return stamp

def saveDir(settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    
    today = date.today().strftime('%Y%m%d')
    
    root    = exp_globals['root_folder']
    project = exp_globals['project_name']
    device  = exp_globals['device_name']
    
    try:
        meas_type = exp_settings['meas_type']
    except:
        meas_type = exp_settings['spec']['meas_type']
    
    fullpath = os.path.join(root, project, device, meas_type, today)

    try:
        os.makedirs(fullpath)
    except:
        print('Dir already exists')
    return fullpath

def saveFilename(settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    power = str(settings['power']) + 'dB'
    avgs =  str(settings['averages']) + 'Avgs'
    
    today = date.today().strftime('%Y%m%d')
    
    project = exp_globals['project_name']
    device  = exp_globals['device_name']
    
    try:
        meas_type = exp_settings['meas_type']
        extra_comment = exp_settings['extra_comment']
    except:
        meas_type = exp_settings['spec']['meas_type']
    
    # fullpath = os.path.join(project, device, meas_type, extra_comment, power, avgs, today)
    fullpath = os.path.join(meas_type, extra_comment, power, avgs, today)

    try:
        os.makedirs(fullpath)
    except:
        print('Dir already exists')

    path_separator = os.sep

    filename = fullpath.replace(path_separator, '_')

    return filename

def SaveFull(path, name, data, expsettings={}, instruments={}, figures=[], saveHWsettings=True):
    
    if name[-4:] == '.pkl':
        saveName = name
    else:
        saveName = name + '.pkl'
        
    pathStr               = os.path.join(path, saveName)

    toSave                = {}
    toSave['Data']        = data
    toSave['ExpSettings'] = expsettings
    if saveHWsettings:
        toSave['HWSettings']  = SaveInst(instruments)
    toSave['Figures']     = figures

    pickle.dump(toSave, open(pathStr, 'wb'))
    
def SaveInst(instruments):
    HWsettings = {}
    for inst in instruments.keys():
        HWsettings[inst] = instruments[inst].settings
    return HWsettings

def NoDetector_data(modulation_frequency, power, instruments, amp_list = [], phase_list = [], Normalize = True):
    gen = instruments["Rf_gen"]
    hdawg = instruments["AWG"]
    card = instruments["digitizer"]
    lo_gen = instruments['lo']
    lo = lo_gen.ch1
    lo.power = 12
    lo.output = 1
    gen.output = 1
    gen.power = power

    #Base settings for IQ waveform
    '''
    Move the following to exp globals
    '''
    sample_rate = 2.4e9
    total_time = 20e-6
    buffer = 1e-6
    hold = 5e-6

    #######################################################################
    modulation_freq = modulation_frequency 
    sidebandM_amp, carrier_amp, sidebandP_amp = [amp_list[0], amp_list[1], amp_list[2]]
    sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])
    if Normalize == True:
        sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])/(sidebandM_amp + carrier_amp + sidebandP_amp)
    sidebandM_phase, carrier_phase, sidebandP_phase = [phase_list[0]*np.pi, phase_list[1]*np.pi, phase_list[2]*np.pi]

    ##################################################################################

    # pulsed modulation setup
    time = np.linspace(0, total_time, int(total_time*sample_rate))
    time_hold = np.linspace(0, hold, int(hold*sample_rate))
    I_t, Q_t = modulating_IQ(carrier_amp, sidebandM_amp, sidebandP_amp, carrier_phase, sidebandM_phase, sidebandP_phase, modulation_freq, time_hold)
    I = np.zeros(len(time))
    I[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate))* I_t
    Q = np.zeros(len(time))
    Q[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate)) * Q_t 
    
    # Load IQ into HDAWG
    load_hdawg(I, Q, instruments)


    ########################## Get trimmed data
    trim_buffer_start = 3000
    trim_buffer_end = 6000
    #trim_buffer_start = int((buffer + trim_buffer) * card.sampleRate)
    #trim_buffer_end = int((hold - trim_buffer) * card.sampleRate)

    gen.output = 1 # THIS SHOULD BE IN EXP_GLOBALS or something
    RF = np.linspace(6.07e9, 6.11e9, 501)
    rawdata_I = np.zeros((len(RF), int(card.samples)))
    rawdata_Q = np.zeros((len(RF), int(card.samples)))
    data_trimmed_I = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    data_trimmed_Q = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    
    
    for ind in range(len(RF)):
        gen.freq = RF[ind] 
        lo.freq = 6.04e9
        card.ArmAndWait()
        I, Q = card.ReadAllData()
        
        Q_trimmed = Q[0][trim_buffer_start:trim_buffer_end]
        I_trimmed = I[0][trim_buffer_start:trim_buffer_end]
        data_trimmed_Q[ind] = Q_trimmed 
        data_trimmed_I[ind] = I_trimmed
        rawdata_Q[ind] = Q[0]
        rawdata_I[ind] = I[0]
        
    lo.output = 0
    gen.output = 0
        
    # maybe plot the first and last data to see if things look normal
    first_trimmed_I = np.zeros(len(rawdata_I[0]))
    first_trimmed_I[trim_buffer_start:trim_buffer_end] = data_trimmed_I[0]
    first_trimmed_Q = np.zeros(len(rawdata_Q[0]))
    first_trimmed_Q[trim_buffer_start:trim_buffer_end] = data_trimmed_Q[0]
    
    last_trimmed_I = np.zeros(len(rawdata_I[-1]))
    last_trimmed_I[trim_buffer_start:trim_buffer_end] = data_trimmed_I[-1]
    last_trimmed_Q = np.zeros(len(rawdata_Q[-1]))
    last_trimmed_Q[trim_buffer_start:trim_buffer_end] = data_trimmed_Q[-1]
    
    fig = plt.figure(2)
    plt.clf()
    ax = plt.subplot(1,2,1)
    plt.plot(rawdata_I[0], color = 'deepskyblue', label='first-full_I')
    plt.plot(rawdata_Q[0], color = 'deepskyblue', linestyle='dotted', label='first-full_Q')
    plt.plot(first_trimmed_I, color = 'orange', label='first-trimmed_I')
    plt.plot(first_trimmed_Q, color = 'orange', linestyle = 'dotted', label='first-trimmed_Q')
    ax.legend()
    ax = plt.subplot(1,2,2)
    plt.plot(rawdata_I[-1], color = 'deepskyblue', label='last-full_I')
    plt.plot(rawdata_Q[-1], color = 'deepskyblue', linestyle='dotted', label='last-full_Q')
    plt.plot(last_trimmed_I, color = 'orange', label='last-trimmed_I')
    plt.plot(last_trimmed_Q, color = 'orange', linestyle = 'dotted', label='last-trimmed_Q')
    ax.legend()
    
    fig.suptitle('raw traces')
    fig.tight_layout()
    plt.show()
    
    return rawdata_I, rawdata_Q, data_trimmed_I, data_trimmed_Q, modulation_freq, trim_buffer_start, trim_buffer_end, RF

def NoDetector_dataNew(instruments, settings):
    
    # load needed instruments
    gen = instruments["Rf_gen"]
    hdawg = instruments["AWG"]
    card = instruments["digitizer"]
    lo_gen = instruments['lo']
    
    # load measurement settings
    modulation_frequency = settings['mod_freq']
    ddc_freq = settings['ddc_freq']
    amp_list = settings['amp_list']
    phase_list = settings['phase_list']
    Normalize = settings['normalize']
    lo_power = settings['lo_power']
    lo_freq = settings['lo_freq']
    gen_power = settings['gen_power']
    trim_buffer_start = settings['buffer_start']
    trim_buffer_end = settings['buffer_end']
    start_freq = settings['start_freq']
    stop_freq = settings['stop_freq']
    num_points = settings['num_points']
    
    ###########################################
    lo = lo_gen.ch1
    lo.power = lo_power # 12
    lo.output = 1
    gen.output = 1
    gen.power = gen_power

    #Base settings for IQ waveform
    '''
    Move the following to exp globals
    '''
    sample_rate = 2.4e9
    total_time = 20e-6
    buffer = 1e-6
    hold = 5e-6

    #######################################################################
    modulation_freq = modulation_frequency 
    sidebandM_amp, carrier_amp, sidebandP_amp = [amp_list[0], amp_list[1], amp_list[2]]
    sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])
    if Normalize == True:
        sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])/(sidebandM_amp + carrier_amp + sidebandP_amp)
    sidebandM_phase, carrier_phase, sidebandP_phase = [phase_list[0]*np.pi, phase_list[1]*np.pi, phase_list[2]*np.pi]

    ##################################################################################

    # pulsed modulation setup
    time = np.linspace(0, total_time, int(total_time*sample_rate))
    time_hold = np.linspace(0, hold, int(hold*sample_rate))
    I_t, Q_t = modulating_IQ(carrier_amp, sidebandM_amp, sidebandP_amp, carrier_phase, sidebandM_phase, sidebandP_phase, modulation_freq, time_hold)
    I = np.zeros(len(time))
    I[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate))* I_t
    Q = np.zeros(len(time))
    Q[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate)) * Q_t 
    
    # Load IQ into HDAWG
    load_hdawg(I, Q, instruments)


    ########################## Get trimmed data
    #trim_buffer_start = 3000
    #trim_buffer_end = 6000
    #trim_buffer_start = int((buffer + trim_buffer) * card.sampleRate)
    #trim_buffer_end = int((hold - trim_buffer) * card.sampleRate)

    # THIS SHOULD BE IN EXP_GLOBALS or something
    #RF = np.linspace(6.07e9, 6.11e9, 501)
    RF = np.linspace(start_freq, stop_freq, num_points)
    rawdata_I = np.zeros((len(RF), int(card.samples)))
    rawdata_Q = np.zeros((len(RF), int(card.samples)))
    data_trimmed_I = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    data_trimmed_Q = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    
    
    for ind in range(len(RF)):
        gen.freq = RF[ind] 
        lo.freq = lo_freq
        card.ArmAndWait()
        I, Q = card.ReadAllData()
        
        Q_trimmed = Q[0][trim_buffer_start:trim_buffer_end]
        I_trimmed = I[0][trim_buffer_start:trim_buffer_end]
        data_trimmed_Q[ind] = Q_trimmed 
        data_trimmed_I[ind] = I_trimmed
        rawdata_Q[ind] = Q[0]
        rawdata_I[ind] = I[0]
        
    lo.output = 0
    gen.output = 0
        
    # maybe plot the first and last data to see if things look normal
    first_trimmed_I = np.zeros(len(rawdata_I[0]))
    first_trimmed_I[trim_buffer_start:trim_buffer_end] = data_trimmed_I[0]
    first_trimmed_Q = np.zeros(len(rawdata_Q[0]))
    first_trimmed_Q[trim_buffer_start:trim_buffer_end] = data_trimmed_Q[0]
    
    last_trimmed_I = np.zeros(len(rawdata_I[-1]))
    last_trimmed_I[trim_buffer_start:trim_buffer_end] = data_trimmed_I[-1]
    last_trimmed_Q = np.zeros(len(rawdata_Q[-1]))
    last_trimmed_Q[trim_buffer_start:trim_buffer_end] = data_trimmed_Q[-1]
    
    fig = plt.figure(2)
    plt.clf()
    ax = plt.subplot(1,2,1)
    plt.plot(rawdata_I[0], color = 'deepskyblue', label='first-full_I')
    plt.plot(rawdata_Q[0], color = 'deepskyblue', linestyle='dotted', label='first-full_Q')
    plt.plot(first_trimmed_I, color = 'orange', label='first-trimmed_I')
    plt.plot(first_trimmed_Q, color = 'orange', linestyle = 'dotted', label='first-trimmed_Q')
    ax.legend()
    ax = plt.subplot(1,2,2)
    plt.plot(rawdata_I[-1], color = 'deepskyblue', label='last-full_I')
    plt.plot(rawdata_Q[-1], color = 'deepskyblue', linestyle='dotted', label='last-full_Q')
    plt.plot(last_trimmed_I, color = 'orange', label='last-trimmed_I')
    plt.plot(last_trimmed_Q, color = 'orange', linestyle = 'dotted', label='last-trimmed_Q')
    ax.legend()
    
    fig.suptitle('raw traces')
    fig.tight_layout()
    plt.show()
    
    
    ############## PDH LIKE ANALYSIS?
    sampleRate = card.sampleRate
    fig = plt.figure(12)
    plt.clf()
    data_trimmed_amp = np.sqrt(data_trimmed_I**2 + data_trimmed_Q**2)
    data_trimmed = [data_trimmed_I, data_trimmed_Q, data_trimmed_amp]
    label_list = ['I', 'Q', 'amp']
    
    for ind in range(len(data_trimmed)):
    
        # digital downconversion on data from detector
        card_data_trimmed = data_trimmed[ind]
        data_pdh_mean_im = np.zeros(len(card_data_trimmed))
        data_pdh_mean_re = np.zeros(len(card_data_trimmed))
        

    # analyze data for each RF-LO combination
        for k in range(len(card_data_trimmed)):
            reflected_power = card_data_trimmed[k] #trying to understand why this seems to work
    
            time_app = np.linspace(trim_buffer_start/sampleRate, trim_buffer_end/sampleRate, len(card_data_trimmed[k]))
            dc_re = np.cos(2*np.pi*ddc_freq*time_app) 
            dc_im = np.sin(2*np.pi*ddc_freq*time_app) 
            
    
            # downconvert signal 
            reflected_power_re = reflected_power * dc_re 
            reflected_power_im = reflected_power * dc_im 
    
            # second filter 
            order = 5
            cutoff = modulation_freq/10 # ###### check the diff is ddc_freq is used here
    
            filter_re = butter_lowpass_filter(reflected_power_re, cutoff, sampleRate, order)
            filter_im = butter_lowpass_filter(reflected_power_im, cutoff, sampleRate, order)
    
        
            data_pdh_mean_im[k] = np.mean(filter_im)
            data_pdh_mean_re[k] = np.mean(filter_re)
            

        final_pdh_mean_im = data_pdh_mean_im
        final_pdh_mean_re = data_pdh_mean_re
    
        imag_max = np.max(np.abs(final_pdh_mean_im))
        real_max = np.max(np.abs(final_pdh_mean_re))
    
        ax = plt.subplot(3,2,(2*ind)+1)
        plt.plot(RF, final_pdh_mean_im, color = 'orange', label='imag')
        plt.title("{} data pdh".format(label_list[ind]))
        plt.ylim(-1.2*imag_max, 1.2*imag_max)
        plt.grid()
        plt.legend()
        
        ax = plt.subplot(3,2,(2*ind)+2)
        plt.plot(RF, final_pdh_mean_re, color = 'orange', label='real')
        plt.ylim(-1.2*real_max, 1.2*real_max)
        plt.title("{} data pdh".format(label_list[ind]))
        plt.grid()
        plt.legend()
        

    fig.suptitle(f'NoDetector_{gen_power}dBm_{card.averages}Avgs_{modulation_freq/1e6}MHzMod_{ddc_freq/1e6}MHzDDC')
    fig.tight_layout()
    plt.show()
    
    ## mini data
    mini_data = {}
    mini_data['rawdata_I'] = rawdata_I
    mini_data['rawdata_Q'] = rawdata_Q
    mini_data['data_trimmed_I'] = data_trimmed_I
    mini_data['data_trimmed_Q'] = data_trimmed_Q
    
    
    return mini_data


def NoDetector_LOScan(instruments, settings):
    
    # load needed instruments
    gen = instruments["Rf_gen"]
    hdawg = instruments["AWG"]
    card = instruments["digitizer"]
    lo_gen = instruments['lo']
    
    # load measurement settings
    modulation_frequency = settings['mod_freq']
    ddc_freq = settings['ddc_freq']
    amp_list = settings['amp_list']
    phase_list = settings['phase_list']
    Normalize = settings['normalize']
    lo_power = settings['lo_power']
    IF_freq = settings['IF_freq']
    gen_power = settings['gen_power']
    trim_buffer_start = settings['buffer_start']
    trim_buffer_end = settings['buffer_end']
    start_freq = settings['start_freq']
    stop_freq = settings['stop_freq']
    num_points = settings['num_points']
    
    ###########################################
    lo = lo_gen
    lo.power = lo_power # 12
    lo.output = 1
    gen.output = 1
    gen.power = gen_power

    #Base settings for IQ waveform
    '''
    Move the following to exp globals
    '''
    sample_rate = 2.4e9
    total_time = 20e-6
    buffer = 1e-6
    hold = 5e-6

    #######################################################################
    modulation_freq = modulation_frequency 
    sidebandM_amp, carrier_amp, sidebandP_amp = [amp_list[0], amp_list[1], amp_list[2]]
    sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])
    if Normalize == True:
        sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])/(sidebandM_amp + carrier_amp + sidebandP_amp)
    sidebandM_phase, carrier_phase, sidebandP_phase = [phase_list[0]*np.pi, phase_list[1]*np.pi, phase_list[2]*np.pi]

    ##################################################################################

    # pulsed modulation setup
    time = np.linspace(0, total_time, int(total_time*sample_rate))
    time_hold = np.linspace(0, hold, int(hold*sample_rate))
    I_t, Q_t = modulating_IQ(carrier_amp, sidebandM_amp, sidebandP_amp, carrier_phase, sidebandM_phase, sidebandP_phase, modulation_freq, time_hold)
    I = np.zeros(len(time))
    I[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate))* I_t
    Q = np.zeros(len(time))
    Q[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate)) * Q_t 
    
    # Load IQ into HDAWG
    load_hdawg(I, Q, instruments)


    ########################## Get trimmed data
    #trim_buffer_start = 3000
    #trim_buffer_end = 6000
    #trim_buffer_start = int((buffer + trim_buffer) * card.sampleRate)
    #trim_buffer_end = int((hold - trim_buffer) * card.sampleRate)

    # THIS SHOULD BE IN EXP_GLOBALS or something
    #RF = np.linspace(6.07e9, 6.11e9, 501)
    RF = np.linspace(start_freq, stop_freq, num_points)
    rawdata_I = np.zeros((len(RF), int(card.samples)))
    rawdata_Q = np.zeros((len(RF), int(card.samples)))
    data_trimmed_I = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    data_trimmed_Q = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    
    
    for ind in range(len(RF)):
        gen.freq = RF[ind] 
        lo.freq = gen.freq - IF_freq
        card.ArmAndWait()
        I, Q = card.ReadAllData()
        
        Q_trimmed = Q[0][trim_buffer_start:trim_buffer_end]
        I_trimmed = I[0][trim_buffer_start:trim_buffer_end]
        data_trimmed_Q[ind] = Q_trimmed 
        data_trimmed_I[ind] = I_trimmed
        rawdata_Q[ind] = Q[0]
        rawdata_I[ind] = I[0]
        
    lo.output = 0
    gen.output = 0
        
    # maybe plot the first and last data to see if things look normal
    first_trimmed_I = np.zeros(len(rawdata_I[0]))
    first_trimmed_I[trim_buffer_start:trim_buffer_end] = data_trimmed_I[0]
    first_trimmed_Q = np.zeros(len(rawdata_Q[0]))
    first_trimmed_Q[trim_buffer_start:trim_buffer_end] = data_trimmed_Q[0]
    
    last_trimmed_I = np.zeros(len(rawdata_I[-1]))
    last_trimmed_I[trim_buffer_start:trim_buffer_end] = data_trimmed_I[-1]
    last_trimmed_Q = np.zeros(len(rawdata_Q[-1]))
    last_trimmed_Q[trim_buffer_start:trim_buffer_end] = data_trimmed_Q[-1]
    
    fig = plt.figure(109)
    plt.clf()
    ax = plt.subplot(1,2,1)
    plt.plot(rawdata_I[0], color = 'deepskyblue', label='first-full_I')
    plt.plot(rawdata_Q[0], color = 'deepskyblue', linestyle='dotted', label='first-full_Q')
    plt.plot(first_trimmed_I, color = 'orange', label='first-trimmed_I')
    plt.plot(first_trimmed_Q, color = 'orange', linestyle = 'dotted', label='first-trimmed_Q')
    ax.legend()
    ax = plt.subplot(1,2,2)
    plt.plot(rawdata_I[-1], color = 'deepskyblue', label='last-full_I')
    plt.plot(rawdata_Q[-1], color = 'deepskyblue', linestyle='dotted', label='last-full_Q')
    plt.plot(last_trimmed_I, color = 'orange', label='last-trimmed_I')
    plt.plot(last_trimmed_Q, color = 'orange', linestyle = 'dotted', label='last-trimmed_Q')
    ax.legend()
    
    fig.suptitle('raw traces')
    fig.tight_layout()
    plt.show()
    
    
    ############## PDH LIKE ANALYSIS?
    sampleRate = card.sampleRate
    fig = plt.figure(110)
    plt.clf()
    data_trimmed_amp = np.sqrt(data_trimmed_I**2 + data_trimmed_Q**2)
    data_trimmed = [data_trimmed_I, data_trimmed_Q, data_trimmed_amp]
    label_list = ['I', 'Q', 'amp']
    
    for ind in range(len(data_trimmed)):
    
        # digital downconversion on data from detector
        card_data_trimmed = data_trimmed[ind]
        data_pdh_mean_im = np.zeros(len(card_data_trimmed))
        data_pdh_mean_re = np.zeros(len(card_data_trimmed))
        

    # analyze data for each RF-LO combination
        for k in range(len(card_data_trimmed)):
            reflected_power = card_data_trimmed[k] #trying to understand why this seems to work
    
            time_app = np.linspace(trim_buffer_start/sampleRate, trim_buffer_end/sampleRate, len(card_data_trimmed[k]))
            dc_re = np.cos(2*np.pi*ddc_freq*time_app) 
            dc_im = np.sin(2*np.pi*ddc_freq*time_app) 
            
    
            # downconvert signal 
            reflected_power_re = reflected_power * dc_re 
            reflected_power_im = reflected_power * dc_im 
    
            # second filter 
            order = 5
            cutoff = modulation_freq/10 # ###### check the diff is ddc_freq is used here
    
            filter_re = butter_lowpass_filter(reflected_power_re, cutoff, sampleRate, order)
            filter_im = butter_lowpass_filter(reflected_power_im, cutoff, sampleRate, order)
    
        
            data_pdh_mean_im[k] = np.mean(filter_im)
            data_pdh_mean_re[k] = np.mean(filter_re)
            

        final_pdh_mean_im = data_pdh_mean_im
        final_pdh_mean_re = data_pdh_mean_re
    
        imag_max = np.max(np.abs(final_pdh_mean_im))
        real_max = np.max(np.abs(final_pdh_mean_re))
    
        ax = plt.subplot(3,2,(2*ind)+1)
        plt.plot(RF, final_pdh_mean_im, color = 'orange', label='imag')
        plt.title("{} data pdh".format(label_list[ind]))
        plt.ylim(-1.2*imag_max, 1.2*imag_max)
        plt.grid()
        plt.legend()
        
        ax = plt.subplot(3,2,(2*ind)+2)
        plt.plot(RF, final_pdh_mean_re, color = 'orange', label='real')
        plt.ylim(-1.2*real_max, 1.2*real_max)
        plt.title("{} data pdh".format(label_list[ind]))
        plt.grid()
        plt.legend()
        

    fig.suptitle(f'NoDetLOScan_{gen_power}dBm_{card.averages}Avgs_{modulation_freq/1e6}MHzMod_{ddc_freq/1e6}MHzDDC')
    fig.tight_layout()
    plt.show()
    
    ## mini data
    mini_data = {}
    mini_data['rawdata_I'] = rawdata_I
    mini_data['rawdata_Q'] = rawdata_Q
    mini_data['data_trimmed_I'] = data_trimmed_I
    mini_data['data_trimmed_Q'] = data_trimmed_Q
    
    # plot_fft
    fft_I = np.fft.fft(data_trimmed_I[0])
    fft_magnitude_I = np.abs(fft_I)
    fft_freq_I = np.fft.fftfreq(len(data_trimmed_I[0]), 1/card.sampleRate)
    
    fft_Q = np.fft.fft(data_trimmed_Q[0])
    fft_magnitude_Q = np.abs(fft_Q)
    fft_freq_Q = np.fft.fftfreq(len(data_trimmed_Q[0]), 1/card.sampleRate)
    
    fig = plt.figure(111)
    plt.subplot(1,2,1)
    plt.plot(fft_freq_I[:len(fft_freq_I)//2], fft_magnitude_I[:len(fft_magnitude_I)//2], label="I")
    plt.plot(fft_freq_Q[:len(fft_freq_Q)//2], fft_magnitude_Q[:len(fft_magnitude_Q)//2], label="Q")
    plt.legend()
    plt.show()
    
    ## plot periodogram
    xaxis_I, yaxis_I = periodogram(data_trimmed_I[0], card.sampleRate)
    xaxis_Q, yaxis_Q = periodogram(data_trimmed_Q[0], card.sampleRate)
    plt.subplot(1,2,2)
    plt.plot(xaxis_I, yaxis_I, label="I")
    plt.plot(xaxis_Q, yaxis_Q, label="Q")
    plt.legend()
    plt.show()
    
    
    
    
    
    
    return mini_data
  

def remove_IQ_ellipse(Is, Qs):
    
    #mix_corr = uf.fitEllipse(Is, Qs)
    
    center = [0,0] #mix_corr[1]#[0,0]
    axes   = [1,1] #mix_corr[0]#[1,1]
    phi    = 0 #mix_corr[2]#0

    Isprime = np.cos(phi)*(Is-center[0]) + np.sin(phi)*(Qs-center[1]) + center[0]
    Qsprime = -np.sin(phi)*(Is-center[0]) + np.cos(phi)*(Qs-center[1]) 
    Qsprime = Qsprime*axes[0]/axes[1] + center[1]
    return Isprime, Qsprime
    
def extract_DDC(rawdata, instruments, settings):
    
    card = instruments["digitizer"]
    ddc_freq = settings['ddc_freq']
    trim_buffer_start = settings['buffer_start']
    trim_buffer_end = settings['buffer_end']
    correction_angle = settings['correction_angle']
    cutoff_freq = settings['cutoff_freq']
    
    theta_lo = correction_angle * np.pi / 180
    
    time_ddc_full = np.linspace(0, card.samples, card.samples)/card.sampleRate
    dc_sin = np.sin(2*np.pi*ddc_freq*time_ddc_full + theta_lo)
    dc_cos = np.cos(2*np.pi*ddc_freq*time_ddc_full + theta_lo)
    
    ddc_data_sin = rawdata * dc_sin
    ddc_data_cos = rawdata * dc_cos
    
    # now filter
    order = 5
    cutoff = cutoff_freq#ddc_freq/10
    
    filtered_sin_full = butter_lowpass_filter(ddc_data_sin, cutoff, card.sampleRate, order)
    filtered_cos_full = butter_lowpass_filter(ddc_data_cos, cutoff, card.sampleRate, order)
    
    # extract trimmed ddc data
    time_ddc_trim = time_ddc_full[trim_buffer_start:trim_buffer_end]
    filtered_sin_trim = filtered_sin_full[trim_buffer_start:trim_buffer_end]
    filtered_cos_trim = filtered_cos_full[trim_buffer_start:trim_buffer_end]
    
    return filtered_sin_full, filtered_cos_full, time_ddc_full, filtered_sin_trim, filtered_cos_trim, time_ddc_trim
    
def NoDetector_LOScan_NEW(instruments, settings):
    
    # load needed instruments
    gen = instruments["Rf_gen"]
    hdawg = instruments["AWG"]
    card = instruments["digitizer"]
    LO = instruments['LO']
    LO_power = settings['LO_power']
    
    
    # load measurement settings
    modulation_frequency = settings['mod_freq']
    ddc_freq = settings['ddc_freq']
    amp_list = settings['amp_list']
    phase_list = settings['phase_list']
    Normalize = settings['normalize']
    IF_freq = settings['IF_freq']
    gen_power = settings['gen_power']
    start_freq = settings['start_freq']
    stop_freq = settings['stop_freq']
    num_points = settings['num_points']
    trim_buffer_start = settings['buffer_start']
    trim_buffer_end = settings['buffer_end']
    cutoff_freq = settings['cutoff_freq']
    
    # Base settingd for IQ waveform
    total_time = settings['total_IQtime'] #20e-6
    buffer = settings['data_buffer'] # 1e-6
    hold = settings['data_hold'] # 5e-6
    sample_rate = settings['HDAWG_sampleRate']  
    
    ###########################################
    LO.power = LO_power # 12
    LO.output = 1
    gen.output = 1
    gen.power = gen_power


    #######################################################################
    modulation_freq = modulation_frequency 
    sidebandM_amp, carrier_amp, sidebandP_amp = [amp_list[0], amp_list[1], amp_list[2]]
    sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])
    if Normalize == True:
        sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])/(sidebandM_amp + carrier_amp + sidebandP_amp)
    sidebandM_phase, carrier_phase, sidebandP_phase = [phase_list[0]*np.pi, phase_list[1]*np.pi, phase_list[2]*np.pi]

    ##################################################################################

    # pulsed modulation setup
    time = np.linspace(0, total_time, int(total_time*sample_rate))
    time_hold = np.linspace(0, hold, int(hold*sample_rate))
    I_t, Q_t = modulating_IQ(carrier_amp, sidebandM_amp, sidebandP_amp, carrier_phase, sidebandM_phase, sidebandP_phase, modulation_freq, time_hold)
    I = np.zeros(len(time))
    I[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate))* I_t
    Q = np.zeros(len(time))
    Q[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate)) * Q_t 
    
    # Load IQ into HDAWG
    load_hdawg(I, Q, instruments)


    # THIS SHOULD BE IN EXP_GLOBALS or something
    RF = np.linspace(start_freq, stop_freq, num_points)
    
    # STORAGE
    
    # 1) Store rawest data
    rawdata_I = np.zeros((len(RF), int(card.samples)))
    rawdata_Q = np.zeros((len(RF), int(card.samples)))
    
    rawdata_I_trimmed = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    rawdata_Q_trimmed = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    
    # 2) Store rawest data corrected for mixer imperfections (ellipticity))
    
    rawdata_I_mixer = np.zeros((len(RF), int(card.samples)))
    rawdata_Q_mixer = np.zeros((len(RF), int(card.samples)))
    
    rawdata_I_trimmed_mixer = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    rawdata_Q_trimmed_mixer = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    
    # 3) Store downconcerted IQ data 
    
    I_ddc_full = np.zeros((len(RF), int(card.samples)))
    Q_ddc_full = np.zeros((len(RF), int(card.samples)))
    
    I_ddc_trim = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    Q_ddc_trim = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    
    # 4) Store final IQ per freq
    
    I_final = np.zeros(len(RF))
    Q_final = np.zeros(len(RF))
    amp_final = np.zeros(len(RF))
    phase_final = np.zeros(len(RF))
    
    I_final2_sin = np.zeros(len(RF))
    Q_final2_sin = np.zeros(len(RF))
    
    I_final2_cos = np.zeros(len(RF))
    Q_final2_cos = np.zeros(len(RF))
    
    
    
    for ind in range(len(RF)):
        gen.freq = RF[ind] 
        LO.freq = gen.freq - IF_freq
        card.ArmAndWait()
        I, Q = card.ReadAllData()
        Iraw = np.mean(I, 0)
        Qraw = np.mean(Q, 0)
        
        # store raw IQ data
        rawdata_I[ind,:] = Iraw
        rawdata_Q[ind,:] = Qraw
        
        # store trimmed raw IQ data
        rawdata_I_trimmed[ind,:] = rawdata_I[ind][trim_buffer_start:trim_buffer_end]
        rawdata_Q_trimmed[ind,:] = rawdata_Q[ind][trim_buffer_start:trim_buffer_end]
        
        # correct for mixer imperfection
        Iraw_mixer, Qraw_mixer = remove_IQ_ellipse(Iraw, Qraw)
        
        # store mixer corrected IQ
        rawdata_I_mixer[ind,:] = Iraw_mixer
        rawdata_Q_mixer[ind,:] = Qraw_mixer
        
        rawdata_I_trimmed_mixer[ind,:] = rawdata_I_mixer[ind][trim_buffer_start:trim_buffer_end]
        rawdata_Q_trimmed_mixer[ind,:] = rawdata_Q_mixer[ind][trim_buffer_start:trim_buffer_end]
        
        # do digital downconversion
        
        Iddc_sin_full, Iddc_cos_full, Iddc_time_full, Iddc_sin_trim, Iddc_cos_trim, Iddc_time_trim, = extract_DDC(Iraw_mixer, instruments, settings)
        Qddc_sin_full, Qddc_cos_full, Qddc_time_full, Qddc_sin_trim, Qddc_cos_trim, Qddc_time_trim, = extract_DDC(Qraw_mixer, instruments, settings)
        
        ### ROTATE MIXER Q SIGNAL INTO I SO THEY CAN BE AVERAGED PROPERLY
        Qprime_sin_full = Qddc_cos_full
        Qprime_cos_full = -Qddc_sin_full
        
        Qprime_sin_trim = Qddc_cos_trim
        Qprime_cos_trim = -Qddc_sin_trim
        
        I_full = (Iddc_cos_full + Qprime_cos_full)/2
        Q_full = (Iddc_sin_full + Qprime_sin_full)/2
        
        I_trim = (Iddc_cos_trim + Qprime_cos_trim)/2
        Q_trim = (Iddc_sin_trim + Qprime_sin_trim)/2
        
        I_ddc_full[ind,:] = I_full
        Q_ddc_full[ind,:] = Q_full
        
        I_ddc_trim[ind,:] = I_trim
        Q_ddc_trim[ind,:] = Q_trim
        
        I_final[ind] = np.mean(I_trim)
        Q_final[ind] = np.mean(Q_trim)
        amp_final[ind] = np.sqrt(I_final[ind]**2 + Q_final[ind]**2)
        phase_final[ind] = np.arctan2(Q_final[ind], I_final[ind]) * 180 / np.pi
        
        '''
        at this point, I want to assume the ddc is perfect so I just have cosine and sine components of I and Q
        
        '''
        I_final2_cos[ind] = np.mean(Iddc_cos_trim)
        I_final2_sin[ind] = np.mean(Iddc_sin_trim)
        
        Q_final2_cos[ind] = np.mean(Qddc_cos_trim)
        Q_final2_sin[ind] = np.mean(Qddc_sin_trim)
        
        
        
         
    LO.output = 0
    gen.output = 0
        

    ## NOW PLOT
    # Plot a slice of raw IQ data
    fig = plt.figure(11)
    plt.clf()
    ax = plt.subplot(1,2,1)
    plt.plot(Iddc_time_full*1e6, rawdata_I[-1], color = 'deepskyblue', label='rawI_full')
    plt.plot(Iddc_time_trim*1e6, rawdata_I_trimmed[-1], color = 'orange', linestyle='dotted', label='rawI_trimmed')
    plt.xlabel('Time (us)')
    ax.legend()
    
    ax = plt.subplot(1,2,2)
    plt.plot(Qddc_time_full*1e6, rawdata_Q[-1], color = 'deepskyblue', label='rawQ_full')
    plt.plot(Qddc_time_trim*1e6, rawdata_Q_trimmed[-1], color = 'orange', linestyle='dotted', label='rawQ_trimmed')
    plt.xlabel('Time (us)')
    ax.legend()
    
    fig.suptitle('rawest IQ traces')
    fig.tight_layout()
    plt.show()
    
    # Plot a slice of the mixer corrected IQ data
    fig = plt.figure(12)
    plt.clf()
    ax = plt.subplot(1,2,1)
    plt.plot(Iddc_time_full*1e6, rawdata_I_mixer[-1], color = 'deepskyblue', label='rawI_full')
    plt.plot(Iddc_time_trim*1e6, rawdata_I_trimmed_mixer[-1], color = 'orange', linestyle='dotted', label='rawI_trimmed')
    plt.xlabel('Time (us)')
    ax.legend()
    
    ax = plt.subplot(1,2,2)
    plt.plot(Qddc_time_full*1e6, rawdata_Q_mixer[-1], color = 'deepskyblue', label='rawQ_full')
    plt.plot(Qddc_time_trim*1e6, rawdata_Q_trimmed_mixer[-1], color = 'orange', linestyle='dotted', label='rawQ_trimmed')
    plt.xlabel('Time (us)')
    ax.legend()
    
    fig.suptitle('Mixer corrected IQ traces')
    fig.tight_layout()
    plt.show()
    
    # Plot downconverted IQ signal at ddc freq
    fig = plt.figure(13)
    plt.clf()
    ax = plt.subplot(1,2,1)
    plt.plot(Iddc_time_full*1e6, I_ddc_full[-1], color = 'deepskyblue', label='Iddc_full')
    plt.plot(Iddc_time_trim*1e6, I_ddc_trim[-1], color = 'orange', linestyle='dotted', label='Iddc_trimmed')
    plt.xlabel('Time (us)')
    ax.legend()
    
    ax = plt.subplot(1,2,2)
    plt.plot(Qddc_time_full*1e6, Q_ddc_full[-1], color = 'deepskyblue', label='Qddc_full')
    plt.plot(Qddc_time_trim*1e6, Q_ddc_trim[-1], color = 'orange', linestyle='dotted', label='Qddc_trimmed')
    plt.xlabel('Time (us)')
    ax.legend()
    
    fig.suptitle('Down converted IQ traces')
    fig.tight_layout()
    plt.show()
    
    # Plot final IQ for each freq
    fig = plt.figure(14)
    plt.clf()
    ax = plt.subplot(2,2,1)
    plt.plot(RF*1e-9, I_final, color = 'deepskyblue', label='I')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    ax.legend()
    
    ax = plt.subplot(2,2,2)
    plt.plot(RF*1e-9, Q_final, color = 'deepskyblue', label='Q')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    ax.legend()
    
    ax = plt.subplot(2,2,3)
    plt.plot(RF*1e-9, amp_final, color = 'deepskyblue', label='amp')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    ax.legend()
    
    ax = plt.subplot(2,2,4)
    plt.plot(RF*1e-9, phase_final, color = 'deepskyblue', label='phase')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    ax.legend()

    fig.suptitle(f'NoDetLOScan_{gen_power}dBm_{card.averages}Avgs_{IF_freq/1e6}MHzIF_{modulation_freq/1e6}MHzMod_{ddc_freq/1e6}MHzDDC_{cutoff_freq/1e6}MHzCutoff')
    fig.tight_layout()
    plt.show()
    
    # plot periodograms
    
    # for raw IQs
    rawTrimmedI_x, rawTrimmedI_y = periodogram(rawdata_I_trimmed[-1], card.sampleRate)
    #rawTrimmedI_y = np.log10(rawTrimmedI_y)
    
    rawTrimmedQ_x, rawTrimmedQ_y = periodogram(rawdata_Q_trimmed[-1], card.sampleRate)
    #rawTrimmedQ_y = np.log10(rawTrimmedQ_y)
    
    ddcTrimmedI_x, ddcTrimmedI_y = periodogram(I_ddc_trim[-1], card.sampleRate)
    #ddcTrimmedI_y = np.log10(ddcTrimmedI_y)
    
    ddcTrimmedQ_x, ddcTrimmedQ_y = periodogram(Q_ddc_trim[-1], card.sampleRate)
    #ddcTrimmedQ_y = np.log10(ddcTrimmedQ_y)
    
    
    fig = plt.figure(15)
    plt.clf()
    ax = plt.subplot(2,2,1)
    plt.plot(rawTrimmedI_x/1e6, rawTrimmedI_y, label='raw_I')
    #plt.plot(rawTrimmedQ_x/1e6, rawTrimmedQ_y, label='raw_Q')
    plt.xlabel("Freq (MHz)")
    plt.xlim(0,50)
    ax.legend()
    
    ax = plt.subplot(2,2,2)
    plt.plot(rawTrimmedQ_x/1e6, rawTrimmedQ_y, label='raw_Q')
    #plt.plot(ddcTrimmedI_x/1e6, ddcTrimmedI_y, label='ddc_I')
    #plt.plot(ddcTrimmedQ_x/1e6, ddcTrimmedQ_y, label='ddc_Q')
    plt.xlabel("Freq (MHz)")
    plt.xlim(0,50)
    ax.legend()
    
    ax = plt.subplot(2,2,3)
    plt.plot(ddcTrimmedI_x/1e6, ddcTrimmedI_y, label='ddc_I')
    #plt.plot(rawTrimmedI_x/1e6, rawTrimmedI_y, label='raw_I')
    #plt.plot(rawTrimmedQ_x/1e6, rawTrimmedQ_y, label='raw_Q')
    plt.xlabel("Freq (MHz)")
    plt.xlim(0,50)
    ax.legend()
    
    ax = plt.subplot(2,2,4)
    #plt.plot(ddcTrimmedI_x/1e6, ddcTrimmedI_y, label='ddc_I')
    plt.plot(ddcTrimmedQ_x/1e6, ddcTrimmedQ_y, label='ddc_Q')
    plt.xlabel("Freq (MHz)")
    plt.xlim(0,50)
    ax.legend()
    
    cutoff_freq = settings['cutoff_freq']
    fig.suptitle(f'FFTs_{IF_freq/1e6}MHzIF_{modulation_freq/1e6}MHzMod_{ddc_freq/1e6}MHzDDC_{cutoff_freq/1e6}MHzCutoff')
    fig.tight_layout()
    plt.show()
    
    # plot IQ data with former processing idea
    fig = plt.figure(16)
    plt.clf()
    ax = plt.subplot(3,2,1)
    plt.plot(RF/1e9, I_final2_cos, label='Iddc_cos')
    #plt.plot(RF/1e9, I_final2_sin, label='Iddc_sin')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    plt.legend()
    
    ax = plt.subplot(3,2,2)
    plt.plot(RF/1e9, I_final2_sin, label='Iddc_sin')
    #plt.plot(RF/1e9, Q_final2_cos, label='Qddc_cos')
    #plt.plot(RF/1e9, Q_final2_sin, label='Qddc_sin')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    plt.legend()
    
    ax = plt.subplot(3,2,3)
    plt.plot(RF/1e9, Q_final2_cos, label='Qddc_cos')
    #plt.plot(RF/1e9, I_final2_cos - Q_final2_cos, label='ddc_cos')
    #plt.plot(RF/1e9, I_final2_sin - Q_final2_sin, label='ddc_sin')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    plt.legend()
    
    ax = plt.subplot(3,2,4)
    plt.plot(RF/1e9, Q_final2_sin, label='Qddc_sin')
    #plt.plot(RF/1e9, I_final2_cos, label='Iddc_cos')
    #plt.plot(RF/1e9, I_final2_sin, label='Iddc_sin')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    plt.legend()
    
    ax = plt.subplot(3,2,5)
    plt.plot(RF/1e9, I_final2_cos - Q_final2_cos, label='ddc_cos')
    #plt.plot(RF/1e9, Q_final2_cos, label='Qddc_cos')
    #plt.plot(RF/1e9, Q_final2_sin, label='Qddc_sin')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    plt.legend()
    
    ax = plt.subplot(3,2,6)
    #plt.plot(RF/1e9, I_final2_cos - Q_final2_cos, label='ddc_cos')
    plt.plot(RF/1e9, I_final2_sin - Q_final2_sin, label='ddc_sin')
    plt.xlabel('Freq (GHz)')
    plt.grid()
    plt.legend()
    
    fig.suptitle(f'PrevNoDetLOScan_{gen_power}dBm_{card.averages}Avgs_{IF_freq/1e6}MHzIF_{modulation_freq/1e6}MHzMod_{ddc_freq/1e6}MHzDDC_{cutoff_freq/1e6}MHzCutoff')
    plt.tight_layout()
    plt.show()
    
    # plot log form of the FFTS
    rawI_fft = rawdata_I_trimmed[-1]
    rawQ_fft = rawdata_Q_trimmed[-1]
    Iddc_fft = I_ddc_trim[-1]
    Qddc_fft = Q_ddc_trim[-1]

    ## process I & Q FFT 
    fft_mag_Iraw = np.abs(np.fft.fft(rawI_fft))
    fft_freq_Iraw = np.fft.fftfreq(len(rawI_fft),1/card.sampleRate)
    
    fft_mag_Iddc = np.abs(np.fft.fft(Iddc_fft))
    fft_freq_Iddc = np.fft.fftfreq(len(Iddc_fft),1/card.sampleRate)

    ## process Q FFT 
    fft_mag_Qraw = np.abs(np.fft.fft(rawQ_fft))
    fft_freq_Qraw = np.fft.fftfreq(len(rawQ_fft),1/card.sampleRate)
    
    fft_mag_Qddc = np.abs(np.fft.fft(Qddc_fft))
    fft_freq_Qddc = np.fft.fftfreq(len(Qddc_fft),1/card.sampleRate)
    
    
    fig = plt.figure(17)
    plt.clf()
    ax = plt.subplot(2,2,1)
    #plt.plot(rawTrimmedI_x/1e6, np.log10(rawTrimmedI_y), label='raw_I')
    plt.plot(fft_freq_Iraw[:len(fft_freq_Iraw)//2]/1e6, np.log10(fft_mag_Iraw[:len(fft_mag_Iraw)//2]), label='raw_I')
    plt.xlabel("Freq (MHz)")
    plt.xlim(0,50)
    ax.legend()
    
    ax = plt.subplot(2,2,2)
    #plt.plot(rawTrimmedQ_x/1e6, np.log10(rawTrimmedQ_y), label='raw_Q')
    plt.plot(fft_freq_Qraw[:len(fft_freq_Qraw)//2]/1e6, np.log10(fft_mag_Qraw[:len(fft_mag_Qraw)//2]), label='raw_Q')
    plt.xlabel("Freq (MHz)")
    plt.xlim(0,50)
    ax.legend()
    
    ax = plt.subplot(2,2,3)
    #plt.plot(ddcTrimmedI_x/1e6, np.log10(ddcTrimmedI_y), label='ddc_I')
    plt.plot(fft_freq_Iddc[:len(fft_freq_Iddc)//2]/1e6, np.log10(fft_mag_Iddc[:len(fft_mag_Iddc)//2]), label='ddc_I')
    plt.xlabel("Freq (MHz)")
    plt.xlim(0,50)
    ax.legend()
    
    ax = plt.subplot(2,2,4)
    #plt.plot(ddcTrimmedQ_x/1e6, np.log10(ddcTrimmedQ_y), label='ddc_Q')
    plt.plot(fft_freq_Qddc[:len(fft_freq_Qddc)//2]/1e6, np.log10(fft_mag_Qddc[:len(fft_mag_Qddc)//2]), label='ddc_Q')
    plt.xlabel("Freq (MHz)")
    plt.xlim(0,50)
    ax.legend()
    
    cutoff_freq = settings['cutoff_freq']
    fig.suptitle(f'FFTs_{IF_freq/1e6}MHzIF_{modulation_freq/1e6}MHzMod_{ddc_freq/1e6}MHzDDC_{cutoff_freq/1e6}MHzCutoff')
    fig.tight_layout()
    plt.show()
    
    
    ## mini data
    full_data = {}
    full_data['rawI'] = rawdata_I[-1]
    full_data['rawQ'] = rawdata_Q[-1]
    full_data['rawI_trimmed'] = rawdata_I_trimmed#[-1]
    full_data['rawQ_trimmed'] = rawdata_Q_trimmed#[-1]
    full_data['mixerI'] = rawdata_I_mixer[-1]
    full_data['mixerQ'] = rawdata_Q_mixer[-1]
    
    full_data['Iddc'] = I_ddc_full[-1]
    full_data['Qddc'] = Q_ddc_full[-1]
    full_data['I_final'] = I_final
    full_data['Q_final'] = Q_final

    
    
    return full_data



def NoDetector_TakeIQ(instruments, settings):
    # load needed instruments
    gen = instruments["Rf_gen"]
    hdawg = instruments["AWG"]
    card = instruments["digitizer"]
    lo_gen = instruments['lo']
    
    # load measurement settings
    modulation_frequency = settings['mod_freq']
    ddc_freq = settings['ddc_freq']
    amp_list = settings['amp_list']
    phase_list = settings['phase_list']
    Normalize = settings['normalize']
    lo_power = settings['lo_power']
    IF_freq = settings['IF_freq']
    gen_power = settings['gen_power']
    trim_buffer_start = settings['buffer_start']
    trim_buffer_end = settings['buffer_end']
    start_freq = settings['start_freq']
    stop_freq = settings['stop_freq']
    num_points = settings['num_points']
    plotIQraw = settings['plotIQraw']
    
    ###########################################
    lo = lo_gen
    lo.power = lo_power # 12
    lo.output = 1
    gen.output = 1
    gen.power = gen_power

    #Base settings for IQ waveform
    '''
    Move the following to exp globals
    '''
    sample_rate = 2.4e9
    total_time = 20e-6
    buffer = 1e-6
    hold = 5e-6

    #######################################################################
    modulation_freq = modulation_frequency 
    sidebandM_amp, carrier_amp, sidebandP_amp = [amp_list[0], amp_list[1], amp_list[2]]
    sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])
    if Normalize == True:
        sidebandM_amp, carrier_amp, sidebandP_amp = np.array([sidebandM_amp, carrier_amp, sidebandP_amp])/(sidebandM_amp + carrier_amp + sidebandP_amp)
    sidebandM_phase, carrier_phase, sidebandP_phase = [phase_list[0]*np.pi, phase_list[1]*np.pi, phase_list[2]*np.pi]

    ##################################################################################

    # pulsed modulation setup
    time = np.linspace(0, total_time, int(total_time*sample_rate))
    time_hold = np.linspace(0, hold, int(hold*sample_rate))
    I_t, Q_t = modulating_IQ(carrier_amp, sidebandM_amp, sidebandP_amp, carrier_phase, sidebandM_phase, sidebandP_phase, modulation_freq, time_hold)
    I = np.zeros(len(time))
    I[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate))* I_t
    Q = np.zeros(len(time))
    Q[int(buffer*sample_rate):int((buffer+hold)*sample_rate)] = np.ones(int(hold*sample_rate)) * Q_t 
    
    # Load IQ into HDAWG
    load_hdawg(I, Q, instruments)


    ########################## Get trimmed data
    RF = np.linspace(start_freq, stop_freq, num_points)
    rawdata_I = np.zeros((len(RF), int(card.samples)))
    rawdata_Q = np.zeros((len(RF), int(card.samples)))
    data_trimmed_I = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    data_trimmed_Q = np.zeros((len(RF), trim_buffer_end - trim_buffer_start))
    
    
    for ind in range(len(RF)):
        gen.freq = RF[ind] 
        lo.freq = gen.freq - IF_freq
        card.ArmAndWait()
        I, Q = card.ReadAllData()
        
        Q_trimmed = Q[0][trim_buffer_start:trim_buffer_end]
        I_trimmed = I[0][trim_buffer_start:trim_buffer_end]
        data_trimmed_Q[ind] = Q_trimmed 
        data_trimmed_I[ind] = I_trimmed
        rawdata_Q[ind] = Q[0]
        rawdata_I[ind] = I[0]
        
    lo.output = 0
    gen.output = 0
    
    if plotIQraw:
        # maybe plot the first and last data to see if things look normal
        first_trimmed_I = np.zeros(len(rawdata_I[0]))
        first_trimmed_I[trim_buffer_start:trim_buffer_end] = data_trimmed_I[0]
        first_trimmed_Q = np.zeros(len(rawdata_Q[0]))
        first_trimmed_Q[trim_buffer_start:trim_buffer_end] = data_trimmed_Q[0]
        
        last_trimmed_I = np.zeros(len(rawdata_I[-1]))
        last_trimmed_I[trim_buffer_start:trim_buffer_end] = data_trimmed_I[-1]
        last_trimmed_Q = np.zeros(len(rawdata_Q[-1]))
        last_trimmed_Q[trim_buffer_start:trim_buffer_end] = data_trimmed_Q[-1]
        
        fig = plt.figure(2)
        plt.clf()
        ax = plt.subplot(1,2,1)
        plt.plot(rawdata_I[0], color = 'deepskyblue', label='first-full_I')
        plt.plot(rawdata_Q[0], color = 'deepskyblue', linestyle='dotted', label='first-full_Q')
        plt.plot(first_trimmed_I, color = 'orange', label='first-trimmed_I')
        plt.plot(first_trimmed_Q, color = 'orange', linestyle = 'dotted', label='first-trimmed_Q')
        ax.legend()
        ax = plt.subplot(1,2,2)
        plt.plot(rawdata_I[-1], color = 'deepskyblue', label='last-full_I')
        plt.plot(rawdata_Q[-1], color = 'deepskyblue', linestyle='dotted', label='last-full_Q')
        plt.plot(last_trimmed_I, color = 'orange', label='last-trimmed_I')
        plt.plot(last_trimmed_Q, color = 'orange', linestyle = 'dotted', label='last-trimmed_Q')
        ax.legend()
        
        fig.suptitle('raw traces')
        fig.tight_layout()
        plt.show()
    
    data_dict = {}
    data_dict['rawdata_I'] = rawdata_I
    data_dict['rawdata_Q'] = rawdata_Q
    data_dict['data_trimmed_I'] = data_trimmed_I
    data_dict['data_trimmed_Q'] = data_trimmed_Q
    
    return data_dict

def compute_pdh(instruments, settings):
    

    card = instruments["digitizer"]
    
    # load measurement settings
    modulation_frequency = settings['mod_freq']
    ddc_freq = settings['ddc_freq']
    gen_power = settings['gen_power']
    trim_buffer_start = settings['buffer_start']
    trim_buffer_end = settings['buffer_end']
    start_freq = settings['start_freq']
    stop_freq = settings['stop_freq']
    num_points = settings['num_points']
    plot_ErrorSignal = settings['plot_ErrorSignal']
    correction_angle = settings['correction_angle']
    
    
    # get IQ data points
    data_dict = NoDetector_TakeIQ(instruments, settings)
    rawdata_I = data_dict['rawdata_I'] 
    rawdata_Q = data_dict['rawdata_Q']
    data_trimmed_I = data_dict['data_trimmed_I']
    data_trimmed_Q = data_dict['data_trimmed_Q']
    
    
    sampleRate = card.sampleRate
    RF = np.linspace(start_freq, stop_freq, num_points)
    fig = plt.figure(12)
    plt.clf()
    data_trimmed_amp = np.sqrt(data_trimmed_I**2 + data_trimmed_Q**2)
    data_trimmed = [data_trimmed_I, data_trimmed_Q, data_trimmed_amp]
    label_list = ['I', 'Q', 'amp']
    
    # save I, Q, amp PDH data
    final_pdh_mean_im = np.zeros((len(data_trimmed), len(data_trimmed[0])))
    final_pdh_mean_re = np.zeros((len(data_trimmed), len(data_trimmed[0])))
    
    for ind in range(len(data_trimmed)):
    
        # digital downconversion on data from detector
        card_data_trimmed = data_trimmed[ind]
        data_pdh_mean_im = np.zeros(len(card_data_trimmed))
        data_pdh_mean_re = np.zeros(len(card_data_trimmed))
        

    # analyze data for each RF-LO combination
        for k in range(len(card_data_trimmed)):
            reflected_power = card_data_trimmed[k] #trying to understand why this seems to work
    
            time_app = np.linspace(trim_buffer_start/sampleRate, trim_buffer_end/sampleRate, len(card_data_trimmed[k]))
            theta_lo = correction_angle * np.pi / 180
            dc_re = np.cos(2*np.pi*ddc_freq*time_app + theta_lo) 
            dc_im = np.sin(2*np.pi*ddc_freq*time_app + theta_lo) 
            
    
            # downconvert signal 
            reflected_power_re = reflected_power * dc_re 
            reflected_power_im = reflected_power * dc_im 
    
            # second filter 
            order = 5
            cutoff = modulation_frequency/10 # ###### check the diff is ddc_freq is used here
    
            filter_re = butter_lowpass_filter(reflected_power_re, cutoff, sampleRate, order)
            filter_im = butter_lowpass_filter(reflected_power_im, cutoff, sampleRate, order)
    
        
            data_pdh_mean_im[k] = np.mean(filter_im)
            data_pdh_mean_re[k] = np.mean(filter_re)
            

        final_pdh_mean_im[ind] = data_pdh_mean_im
        final_pdh_mean_re[ind] = data_pdh_mean_re
    
        imag_max = np.max(np.abs(data_pdh_mean_im))
        real_max = np.max(np.abs(data_pdh_mean_re))
    
        if plot_ErrorSignal:
            ax = plt.subplot(3,2,(2*ind)+1)
            plt.plot(RF, data_pdh_mean_im, color = 'orange', label='imag')
            plt.title("{} data pdh".format(label_list[ind]))
            plt.ylim(-1.2*imag_max, 1.2*imag_max)
            plt.grid()
            plt.legend()
            
            ax = plt.subplot(3,2,(2*ind)+2)
            plt.plot(RF, data_pdh_mean_re, color = 'orange', label='real')
            plt.ylim(-1.2*real_max, 1.2*real_max)
            plt.title("{} data pdh".format(label_list[ind]))
            plt.grid()
            plt.legend()
        
    if plot_ErrorSignal:
        fig.suptitle(f'NoDetLOScan_{gen_power}dBm_{card.averages}Avgs_{modulation_frequency/1e6}MHzMod_{ddc_freq/1e6}MHzDDC')
        fig.tight_layout()
        plt.show()
        
    
    data_dict2 = {}
    data_dict2['IQ_raw'] = data_dict
    data_dict2['real'] = final_pdh_mean_re
    data_dict2['imag'] = final_pdh_mean_im
    
    return data_dict2

def find_minimum_in_range(frequencies, values, start_freq, end_freq):
    # Create a mask for the desired frequency range
    mask = (frequencies >= start_freq) & (frequencies <= end_freq)
    # Apply the mask to both frequencies and values
    filtered_frequencies = frequencies[mask]
    filtered_values = values[mask]
    # Find the index of the minimum value within the filtered range
    min_index_in_range = np.argmin(filtered_values)
    max_index_in_range = np.argmax(filtered_values)
    
    # Get the actual index in the original array
    original_index_min = np.where(frequencies == filtered_frequencies[min_index_in_range])[0][0]
    original_index_max = np.where(frequencies == filtered_frequencies[max_index_in_range])[0][0]
    
    return original_index_min, original_index_max

def stability_pdh(instruments, settings, num_measure):
    
    card = instruments["digitizer"]
    
    # load measurement settings
    modulation_frequency = settings['mod_freq']
    ddc_freq = settings['ddc_freq']
    gen_power = settings['gen_power']
    start_freq = settings['start_freq']
    stop_freq = settings['stop_freq']
    num_points = settings['num_points']
    center = settings['cavity_freq']
    
    #get first PDH error signal data for stability points calibration
    data_dict2 = compute_pdh(instruments, settings)
    final_pdh_mean_re = data_dict2['real']
    final_pdh_mean_im = data_dict2['imag']

    
    label_list = ['I', 'Q', 'combo']
    lower_bound = center - 0.25*modulation_frequency
    upper_bound = center + 0.25*modulation_frequency
    RF = np.linspace(start_freq, stop_freq, num_points)
    
    minInd = np.zeros(2)
    maxInd = np.zeros(2)
    
    fig = plt.figure(14)
    plt.clf()
    for ind in range(3):
        if ind!=2:
            real_comp = final_pdh_mean_re[ind] 
            #get min and max index for the real comp
            min_index_real, max_index_real = find_minimum_in_range(RF, real_comp, lower_bound, upper_bound)
            
            imag_comp = final_pdh_mean_im[ind]
            #get min and max index for the imag comp
            min_index_imag, max_index_imag = find_minimum_in_range(RF, imag_comp, lower_bound, upper_bound)
            
            imag_max = np.max(np.abs(imag_comp))
            real_max = np.max(np.abs(real_comp))
            
            # get the frequencies
            freq_at_min_imag = RF[min_index_imag]
            freq_at_max_imag = RF[max_index_imag]
    
            freq_at_min_real = RF[min_index_real]
            freq_at_max_real = RF[max_index_real]
            
            ax = plt.subplot(3,2, (2*ind)+1)
            plt.plot(RF/1e9,imag_comp, 'k', label=f'imag comp {label_list[ind]}')
            plt.scatter(freq_at_min_imag/1e9, imag_comp[min_index_imag], marker='o', label='Min_point')
            plt.scatter(freq_at_max_imag/1e9, imag_comp[max_index_imag], marker='x', label='Max_point')
            plt.ylim(-1.2*imag_max, 1.2*imag_max)
            plt.xlabel('frequency (GHz)')
            plt.grid()
            plt.legend()

            ax = plt.subplot(3,2, (2*ind)+2)
            plt.plot(RF/1e9,real_comp, 'b', label = f'real comp {label_list[ind]}')
            plt.scatter(freq_at_min_real/1e9, real_comp[min_index_real], marker='o', label='Min_point')
            plt.scatter(freq_at_max_real/1e9, real_comp[max_index_real], marker='x', label='Max_point')
            plt.ylim(-1.2*real_max, 1.2*real_max)
            plt.xlabel('frequency (GHz)')
            plt.grid()
            plt.legend()
        else:
            I_pdh_re = final_pdh_mean_re[0]
            I_pdh_im = final_pdh_mean_im[0]
            
            Q_pdh_re = final_pdh_mean_re[1]
            Q_pdh_im = final_pdh_mean_im[1]
            
            combo_pdh_re = I_pdh_re - Q_pdh_re
            combo_pdh_im = I_pdh_im - Q_pdh_im
            
            real_comp = combo_pdh_re 
            #get min and max index for the real comp
            min_index_real, max_index_real = find_minimum_in_range(RF, real_comp, lower_bound, upper_bound)
            
            imag_comp = combo_pdh_im
            #get min and max index for the imag comp
            min_index_imag, max_index_imag = find_minimum_in_range(RF, imag_comp, lower_bound, upper_bound)
            
            imag_max = np.max(np.abs(imag_comp))
            real_max = np.max(np.abs(real_comp))
            
            #save the min and max indices 
            minInd[0] = min_index_imag
            minInd[1] = min_index_real
            
            maxInd[0] = max_index_imag
            maxInd[1] = max_index_real
            
            # get the frequencies
            freq_at_min_imag = RF[min_index_imag]
            freq_at_max_imag = RF[max_index_imag]
    
            freq_at_min_real = RF[min_index_real]
            freq_at_max_real = RF[max_index_real]
            
            ax = plt.subplot(3,2, (2*ind)+1)
            plt.plot(RF/1e9,imag_comp, 'k', label=f'imag comp {label_list[ind]}')
            plt.scatter(freq_at_min_imag/1e9, imag_comp[min_index_imag], marker='o', label='Min_point')
            plt.scatter(freq_at_max_imag/1e9, imag_comp[max_index_imag], marker='x', label='Max_point')
            plt.ylim(-1.2*imag_max, 1.2*imag_max)
            plt.xlabel('frequency (GHz)')
            plt.grid()
            plt.legend()
            
            ax = plt.subplot(3,2, (2*ind)+2)
            plt.plot(RF/1e9,real_comp, 'b', label = f'real comp {label_list[ind]}')
            plt.scatter(freq_at_min_real/1e9, real_comp[min_index_real], marker='o', label='Min_point')
            plt.scatter(freq_at_max_real/1e9, real_comp[max_index_real], marker='x', label='Max_point')
            plt.ylim(-1.2*real_max, 1.2*real_max)
            plt.xlabel('frequency (GHz)')
            plt.grid()
            plt.legend()
    fig.suptitle(f'NoDetLOScan_{gen_power}dBm_{card.averages}Avgs_{modulation_frequency/1e6}MHzMod_{ddc_freq/1e6}MHzDDC')
    fig.tight_layout()
    plt.show()
        
    # now do repeated measurements
    wait_time = 1
    ts = np.zeros(num_measure)
    combo_reals_atMax = np.zeros(num_measure)
    combo_reals_atMin = np.zeros(num_measure) 

    combo_imags_atMax = np.zeros(num_measure)
    combo_imags_atMin = np.zeros(num_measure) 
    
    start_time = time.time()
    for ind in range(num_measure):
        data_dict2 = compute_pdh(instruments, settings)
        final_pdh_mean_re = data_dict2['real']
        final_pdh_mean_im = data_dict2['imag']
        
        # the final_pdh_mean_re/im should have pdh error signal data from I, Q, and AMP
        I_pdh_re = final_pdh_mean_re[0]
        I_pdh_im = final_pdh_mean_im[0]
        
        Q_pdh_re = final_pdh_mean_re[1]
        Q_pdh_im = final_pdh_mean_im[1]
        
        combo_pdh_re = I_pdh_re - Q_pdh_re
        combo_pdh_im = I_pdh_im - Q_pdh_im
        
        imag_maxInd = int(maxInd[0])
        combo_imags_atMax[ind] = combo_pdh_im[imag_maxInd]
        
        imag_minInd = int(minInd[0])
        combo_imags_atMin[ind] = combo_pdh_im[imag_minInd]
        
        real_maxInd = int(maxInd[1])
        combo_reals_atMax[ind] = combo_pdh_re[real_maxInd]
        
        real_minInd = int(minInd[1])
        combo_reals_atMin[ind] = combo_pdh_re[real_minInd]
        
        
        currT = time.time() - start_time
        ts[ind] = currT
        print(f'{ind+1} of {num_measure}, time elapsed: {currT} secs')
    
    # plot data 
    fig = plt.figure(19)
    plt.clf()
    ax = plt.subplot(1, 2, 1)
    plt.plot(ts, combo_imags_atMax, color = 'mediumblue', label = 'imag_Max')
    plt.plot(ts, combo_imags_atMin, color = 'darkorange', label = 'imag_Min')
    plt.xlabel('time (s)')
    plt.grid()
    plt.legend()

    ax = plt.subplot(1, 2, 2)
    plt.plot(ts, combo_reals_atMax, color = 'mediumblue', label = 'real_Max')
    plt.plot(ts, combo_reals_atMin, color = 'darkorange', label = 'real_Min')
    plt.xlabel('time (s)')
    plt.grid()
    plt.legend()
    
    fig.suptitle(f'Stability_NoDetLOScan_{gen_power}dBm_{card.averages}Avgs_{modulation_frequency/1e6}MHzMod_{ddc_freq/1e6}MHzDDC')
    fig.tight_layout()
    plt.show()
    
