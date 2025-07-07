# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab

Modified by: jhyang
"""

from CW_spec import cw_spec, get_cw_spec_settings

# gain_vals = [1000, 2500, 5000] #cav only, qubit fixed at 500
# windows = [10,100]
# for g_vals in gain_vals:
#     for w_vals in windows:
    
settings = get_cw_spec_settings()
settings['scanname'] = 'Two_tone_w_amp_no_sidebands'

settings['cav_freq']  = 7.230796e9  #4.389545e9 + 1e6 #7.08979e9 #7.1995e9 #7.8392e9 original
settings['cav_gain'] = 300
settings['meas_window'] = 1000000
settings['cav_pulse_len'] = 10

qub_center = 3.6117e9 #5.33866e9
span = 50e6 
freq_points = 51

settings['qub_gain']     = 0 #500
settings['qub_len'] = 10 #us

#Sweep Parameters
settings['freq_start']      = qub_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 1


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = cw_spec(soc,soccfg,instruments,fullsettings)