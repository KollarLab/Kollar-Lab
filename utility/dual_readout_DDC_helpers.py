# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:42:51 2024

@author: kollarlab
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
def configure_card_dual(card, card_config, settings):
    '''
    Helper function to configure the card from a set of settings. This will
    force a standard definition of what the digitizer timing should be. Computes
    the total acquisition time and converts it to samples (also configures the 
    rest of the digitizer but this is the main part that requires logic)
    Inputs:
        settings dictionary with the following keys (should be initialized from
        the get_default_settings method)
        meas_window: width of measurment tone
        meas_pos: position of measurement tone relative to the rise edge of the 
            trigger
        emp_delay: line delay and other delays accumulated between the 
            AWG and the digitizer
        init_buffer: buffer to collect data before the pulse
        post_buffer: buffer after measurment tone to wait out the ringing before
            background subtraction
        averages: number of averages for the card to perform in HW
        segments: number of segments (will be averaged together), useful when
            number of averages gets very high ~>50e3
        sampleRate: sampling rate of digitizer (make this lower for longer acquisitions)
        activeChannels: active channels on digitizer (both by default)
        timeout: timeout (in s) for VISA communication (make this longer for 
                         longer acquisitions)
        channelRange: full scale (peak to peak) of each channel, 0.5 or 2.5 V
    '''
    default_meas = settings['base']
    boost = settings['boost']
    readout = settings['readout']
    
    meas_positions = [boost['position'], readout['position']]
    
    start_time = min(meas_positions)
    total_time = max([boost['length'], abs(readout['position']-boost['position'])+readout['length']])+readout['length']
    
    delay = np.round(boost['position']-readout['position'], 9)
    if delay>0:
        boost['shift'] = delay
        readout['shift'] = 0
    else:
        boost['shift'] = 0
        readout['shift'] = -delay
    #Compute total acquisition time

    emp_delay   = default_meas['emp_delay']
    init_buffer = default_meas['init_buffer']
    post_buffer = default_meas['post_buffer']
    card_time = total_time + emp_delay + init_buffer + post_buffer

    meas_samples = card_config['sampleRate']*card_time

    #Compute trigger delay, has to be multiple of 32ns for Acquiris
    trigger_delay = np.floor((start_time - init_buffer)/32e-9)*32e-9
    
    card.channelRange   = card_config['channelRange']
    card.timeout        = card_config['timeout']
    card.sampleRate     = card_config['sampleRate']
    card.activeChannels = card_config['activeChannels']
    card.triggerSlope   = card_config['triggerSlope']
    card.averages       = settings['averages']
    card.segments       = settings['segments']
    card.triggerDelay   = trigger_delay
    card.samples        = int(meas_samples)
    card.SetParams()    

def generate_filter(card, settings):
    DDC_settings = settings['DDC']
    #create Chebychev type II digital filter
    filter_N = DDC_settings['order']
    filter_rs = DDC_settings['stop_atten']
    filter_cutoff = np.abs(DDC_settings['cutoff'])
    LPF = signal.cheby2(filter_N, filter_rs, filter_cutoff, btype='low', analog=False, output='sos', fs=card.sampleRate)
    
    xaxis = np.arange(0, card.samples, 1) * 1/card.sampleRate
    digLO_sin = np.sin(2*np.pi*DDC_settings['IF']*xaxis)
    digLO_cos = np.cos(2*np.pi*DDC_settings['IF']*xaxis)
    
    #store in settings so that the processing functions can get to them
    settings['digLO_sin'] = digLO_sin 
    settings['digLO_cos'] = digLO_cos
    settings['LPF'] = LPF
    
def extract_data_heterodyne(raw_data, xaxis, settings, channel):
    '''
    Extracts the data from a trace. Finds the indices of start/ stop for both
    data and background (assumes that they are the same length, which should
    be enforced by the card samples in the actual exp script) and splices the 
    correct ranges from the raw data. Returns a time axis for the data window 
    and subtracts mean of background if specified
    Key input params:
    ALL VALUES SHOULD BE IN TIME UNITS (base seconds)
        init_buffer: buffer specified to card before the measurement tone
        emp_delay: combination of line delay and other delays that correctly
                   shifts the digitizer to match the HDAWG time axis
        meas_window: width of the measurement pulse
        post_buffer: time to wait after the pulse before measuring background
        
    modified from extract_data so that it will work with heterodyne and digital 
    down conversion
        
    '''
    
    init_buffer = settings['base']['init_buffer']+settings[channel]['shift']
    emp_delay   = settings['base']['emp_delay']
    meas_window = settings[channel]['length']
    post_buffer = settings['base']['post_buffer']
    
    timestep   = xaxis[1] - xaxis[0]
    data_start = int((init_buffer + emp_delay)/timestep)
    back_start = int((init_buffer + emp_delay + meas_window + post_buffer)/timestep)
    window_width = int(meas_window/timestep)
    
    data_x = xaxis[data_start:data_start+window_width]
    background = raw_data[back_start:back_start+window_width]
    back_val = np.mean(background)

    raw_data = raw_data-back_val

    #now I am ready to filter the data.
    LPF = settings['LPF']
    digLO_sin = settings['digLO_sin']
    digLO_cos = settings['digLO_cos']
    
    filtered_sin_full = signal.sosfilt(LPF, raw_data*digLO_sin)
    filtered_cos_full = signal.sosfilt(LPF, raw_data*digLO_cos)
    
    filtered_sin = filtered_sin_full[data_start:data_start+window_width]
    filtered_cos = filtered_cos_full[data_start:data_start+window_width]
    
    return filtered_cos, filtered_sin, data_x, filtered_cos_full, filtered_sin_full, xaxis

def read_and_process_dual_readout(card, settings, plot):
    card.ArmAndWait()
    I, Q = card.ReadAllData()
    readout = np.mean(I, 0)
    boost = np.mean(Q, 0)
    
    boost = boost-np.mean(boost)
    readout = readout-np.mean(readout)
    
    xaxis = np.linspace(0, card.samples, card.samples, endpoint=False)/card.sampleRate
    
    I, Q, time, I_full, Q_full, time_full = extract_data_heterodyne(boost, xaxis,settings, 'boost')
    boost_I, boost_Q = I, Q
    boost_I_full, boost_Q_full = I_full, Q_full
    boost_time = time
    boost_time_full = time_full
    I, Q, time, I_full, Q_full, time_full = extract_data_heterodyne(readout, xaxis, settings, 'readout')
    readout_I, readout_Q = I, Q
    readout_I_full, readout_Q_full = I_full, Q_full
    readout_time = time
    readout_time_full = time_full
    if plot:
        fig0 = plt.figure(97)
        plt.clf()
        ax = plt.subplot(1,3,1)
        plt.plot(xaxis*1e6, boost, label = 'raw boost')
        plt.plot(xaxis*1e6, readout, label = 'raw readout')
        plt.xlabel('Time (us)')
        plt.ylabel('Voltage')
        plt.title('Card Voltages')
        ax.legend(loc = 'upper right')
        
        ax = plt.subplot(132)
        plt.plot(boost_time_full, boost_I_full)
        plt.plot(boost_time_full, boost_Q_full)
        plt.plot(boost_time, boost_I)
        plt.plot(boost_time, boost_Q)
        plt.legend(['Full I', 'Full Q', 'Filt I', 'Filt Q'])
        plt.title('Boost')
        
        ax = plt.subplot(133)
        plt.plot(readout_time_full, readout_I_full)
        plt.plot(readout_time_full, readout_Q_full)
        plt.plot(readout_time, readout_I)
        plt.plot(readout_time, readout_Q)
        plt.legend(['Full I', 'Full Q', 'Filt I', 'Filt Q'])
        plt.title('Readout')
        
        fig0.canvas.draw()
        fig0.canvas.flush_events()
    
    return boost_I, boost_Q, readout_I, readout_Q, xaxis           