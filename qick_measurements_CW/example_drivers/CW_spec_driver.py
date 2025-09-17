# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab

Modified by: jhyang
"""

from qick_measurements_CW.CW_spec import cw_spec, get_cw_spec_settings

# gain_vals = [1000, 2500, 5000] #cav only, qubit fixed at 500
# windows = [10,100]
# for g_vals in gain_vals:
#     for w_vals in windows:
    
settings = get_cw_spec_settings()
settings['scanname'] = 'Cavity_1_Test'

settings['cav_freq']  = 6.52426e9#8.0684e9#7.2309e9  #4.389545e9 + 1e6 #7.08979e9 #7.1995e9 #7.8392e9 original
settings['cav_gain'] = 800#6000
settings['meas_window'] = 10000 #1e6
settings['cav_pulse_len'] = 10 # LEGACY DOES NOT AFFECT OPERATION

qub_center = 5e9#4.275e9#9.46e9#3.6117e9 #5.33866e9
span = 5000e6 
freq_points = 5001

settings['qub_gain_start']     = 5000#8000 #7000
settings['qub_gain_step']      = 5000
settings['qub_gain_points']    = 5

settings['qub_len'] = 10 #us #LEGACY DOES NOT AFFECT OPERATION

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

instruments = {}
full_data = cw_spec(soc,soccfg,instruments,fullsettings)