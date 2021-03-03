# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 09:20:52 2021

@author: Kollarlab
"""
import numpy as np

def check_inputs(inputs, defaults):
    '''
    Checks that the given input dictionary has the correct settings. It's easy
    to accidentally mistype something or us an old naming convention so this 
    function will throw an error if the keys in the dictionaries don't match
    Input params:
        inputs: input dictionary that is modified by user
        defaults: default dictionary specified by script
    '''
    diff1 = set(inputs) - set(defaults)
    diff2 = set(defaults) - set(inputs)
    if len(diff1) !=0 or len(diff2) !=0:
        print('Differing keys:')
        if len(diff1)!=0:
            print(diff1)
        else:
            print(diff2)
        raise ValueError('Incorrect number of inputs, please check the default settings for variable names')

def remove_IQ_ellipse(Is, Qs, axes, center, phi):
    '''
    Removes IQ imbalances from the mixer. This is critical for signals with 
    small SNR since the differences between high and low can be wiped out by
    the ellipse eccentricity. Expects axes center phi in the convention defined
    by the 'fitEllipse' function in the ellipse_fitting file
    Input params:
        Is, Qs: raw I, Q signals 
        axes: major and minor axes of the ellipse
        center: center of ellipse  
        phi: phase rotation of ellipse
    '''
    Isprime = np.cos(phi)*(Is-center[0]) + np.sin(phi)*(Qs-center[1]) + center[0]
    Qsprime = -np.sin(phi)*(Is-center[0]) + np.cos(phi)*(Qs-center[1]) 
    Qsprime = Qsprime*axes[0]/axes[1] + center[1]
    return Isprime, Qsprime
   
def extract_data(raw_data, xaxis, settings):
    '''
    Extracts the data from a trace. Finds the indices of start/ stop for both
    data and background (assumes that they are the same length, which should
    be enforced by the card samples in the actual exp script) and splices the 
    corrct ranges from the raw data. Returns a time axis for the data window 
    and subtracts mean of background if specified
    Key input params:
    ALL VALUES SHOULD BE IN TIME UNITS (base seconds)
        init_buffer: buffer specified to card before the measurement tone
        emp_delay: combination of line delay and other delays that correctly
                   shifts the digitizer to match the HDAWG time axis
        meas_window: width of the measurement pulse
        pulse_buffer: time to wait after the pulse before measuring background
    '''
    init_buffer = settings['init_buffer']
    emp_delay   = settings['emp_delay']
    meas_window = settings['meas_window']
    pulse_buffer= settings['pulse_buffer']
    
    timestep   = xaxis[1] - xaxis[0]
    data_start = int((init_buffer + emp_delay)/timestep)
    back_start = int((init_buffer + emp_delay + meas_window + pulse_buffer)/timestep)
    window_width = int(meas_window/timestep)
    
    print(timestep)
    print(data_start)
    print(back_start)
    print(window_width)
    
    data_x = xaxis[data_start:data_start+window_width]
    back_x = xaxis[back_start:back_start+window_width]
    
    data    = raw_data[data_start:data_start+window_width]
    background = raw_data[back_start:back_start+window_width]
    back_val = np.mean(background)
    
    if settings['subtract_background']:
        return data - back_val, data_x, background, back_x
    else:
        return data, data_x, background, back_x
    
    