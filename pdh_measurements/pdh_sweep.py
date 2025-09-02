# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:05:51 2024

@author: Kollarlab
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, freqz, filtfilt
import pickle
from datetime import datetime
from datetime import date
from pdh_measurements.pdh_funcs import modFreq_linewidth, NoDetector_data, NoDetector_pdh2, NoDetector_pdh1

# instruments used
settings = {}

exp_globals = {}
exp_globals['root_folder'] = r"K:\Data\Oxford"
exp_globals['project_name'] = 'PDH'
exp_globals['device_name'] = 'PDH_v5'
settings['exp_globals'] = exp_globals

exp_settings = {}
exp_settings['meas_type'] = "PDH" 
exp_settings['extra_comment'] = "NoDetectporPDH_test2"
settings['exp_settings'] = exp_settings

instruments = {}
instruments["Rf_gen"] = gen
instruments["AWG"] = hdawg
instruments["digitizer"] = card
instruments['lo'] = holz

###################### ONLY THINGS TO CHANGE ###########################
averages = 1000
card.averages = averages
card.SetParams()
norm_amp = [1,2,1]
set1 = [0.38,0.24,0.38]
set2 = [0.12,0.76,0.12]
set3 = [2,1,2]
set4 = [0.5, 4, 0.5]
amp_list = norm_amp #set4:[0.5, 4, 0.5]#set3:[2,1,2]#[1, 2, 1] #[1, 2, 1] set1=[0.38,0.24,0.38], set2=[0.12,0.76,0.12]
phase_list = [1, 0, 0]
mod_list = [5, 15, 40]#[5,15,40]#[5, 10, 15] # [40, 60, 80]
Normalize = True 
settings['power'] = -5
settings['cavity_linewidth'] = 1e6 #0.7e6 
settings['cavity_freq'] = 6.09e9 #7.28e9 #6.808e9
settings['freq_points']  = 501
settings['detector'] = 'envelope_detector'  # change to 'envelope_detector' if envelope detector is used
settings['averages'] = card.averages
settings['mod_mult'] = mod_list 
settings['plot_data'] = True
settings['narrow_span'] = False # if true, 60MHz span
settings['amp_list'] = amp_list 
settings['phase_list'] = phase_list
settings['Normalize'] = Normalize

########################################################################

start_time = time.time()
# measurement
#card_data, card_data_trimmed = modFreq_linewidth(settings, instruments)
rawdata_I, rawdata_Q, data_trimmed_I, data_trimmed_Q, modulation_freq, trim_buffer_start, trim_buffer_end, RF = NoDetector_data(10e6, settings['power'],
                                                                                                                                instruments, 
                                                                                                                                amp_list = amp_list,
                                                                                                                                phase_list = phase_list, 
                                                                                                                                Normalize = True)
im, re = NoDetector_pdh2(data_trimmed_I, data_trimmed_Q, 10e6, trim_buffer_start, trim_buffer_end, RF)
end_time = time.time() - start_time
print(f'time_elasped:{end_time}')


################################################################################
from scipy.signal import butter, lfilter, freqz, filtfilt
def butter_lowpass(cutoff, fs, order):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order):
        b, a = butter_lowpass(cutoff, fs, order)
        y = lfilter(b, a, data)
        return y

def pdh_data(data, mod_freq, sampleRate):
    
    data_pdh_mean_im = np.zeros(len(data))
    data_pdh_mean_re = np.zeros(len(data))
    for k in range(len(data)):
        reflected_power = data[k]  
        time_app = np.linspace(trim_buffer_start/sampleRate, trim_buffer_end/sampleRate, len(data[k]))
        dc_re = np.cos(2*np.pi*mod_freq*time_app) 
        dc_im = np.sin(2*np.pi*mod_freq*time_app) 
        

        # downconvert signal 
        reflected_power_re = reflected_power * dc_re 
        reflected_power_im = reflected_power * dc_im 

        # second filter 
        order = 5
        cutoff = mod_freq/10 # about 2MHz cutoff

        filter_re = butter_lowpass_filter(reflected_power_re, cutoff, sampleRate, order)
        filter_im = butter_lowpass_filter(reflected_power_im, cutoff, sampleRate, order)

    
        data_pdh_mean_im[k] = np.mean(filter_im)
        data_pdh_mean_re[k] = np.mean(filter_re)
        

    final_pdh_mean_im = data_pdh_mean_im
    final_pdh_mean_re = data_pdh_mean_re


    return final_pdh_mean_im, final_pdh_mean_re

def run_stability(num_points):
    
    
    #I_data = np.zeros((num_points, 501, 3000))
    #Q_data = np.zeros((num_points, 501, 3000))
    
    pdh_I_im = np.zeros(num_points)
    pdh_I_re = np.zeros(num_points)
    
    pdh_Q_im = np.zeros(num_points)
    pdh_Q_re = np.zeros(num_points)
    #I_data = []
    #Q_data = []
    ts = np.zeros(num_points)
    start_time = time.time()
    
    for ind in range(num_points):
        rawdata_I, rawdata_Q, data_trimmed_I, data_trimmed_Q, modulation_freq, trim_buffer_start, trim_buffer_end, RF = NoDetector_data(10e6, settings['power'],
                                                                                                                                        instruments, 
                                                                                                                                        amp_list = amp_list,
                                                                                                                                        phase_list = phase_list, 
                                                                                                                                        Normalize = True)
        #I_data[ind] = data_trimmed_I
       # Q_data[ind] = data_trimmed_Q
    
        
        I_im, I_re = pdh_data(data_trimmed_I, 10e6, 1e9)
        pdh_I_im[ind], pdh_I_re[ind] = np.max(I_im), np.max(I_re)
      
        Q_im, Q_re = pdh_data(data_trimmed_Q, 10e6, 1e9)
        pdh_Q_im[ind], pdh_Q_re[ind] = np.max(Q_im), np.max(Q_re)
        
        currT = time.time() - start_time
        ts[ind] = currT
        
        time.sleep(1)
    
    fig = plt.figure(13)
    plt.clf()
    
    ax = plt.subplot(1,2,1)
    plt.plot(ts, pdh_I_im, label = "I_im")
    plt.plot(ts, pdh_I_re, label = "I_re")
    plt.grid()
    plt.legend()
    
    ax = plt.subplot(1,2,2)
    plt.plot(ts, pdh_Q_im, label = "Q_im")
    plt.plot(ts, pdh_Q_re, label = "Q_re")
    plt.grid()
    plt.legend()
    
    fig.tight_layout()
    plt.show()
    
    print(f'time elasped: {currT}')
