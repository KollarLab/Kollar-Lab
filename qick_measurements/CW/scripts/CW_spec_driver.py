# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab

Modified by: jhyang
"""

from CW_spec import cw_spec, get_cw_spec_settings

settings = get_cw_spec_settings()
settings['scanname'] = 'Jon_qubit_hi_pow'

settings['cav_freq']  = 100e6 #4.389545e9 + 1e6 #7.08979e9 #7.1995e9 #7.8392e9 original
settings['cav_gain'] = 600
settings['meas_window'] = 900

qub_center = 5.33866e9
span = 50e6 
freq_points = 26

settings['qub_gain']     = 1
settings['qub_len'] = 10 #us

#Sweep Parameters
settings['freq_start']      = qub_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 500


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = cw_spec(soc,soccfg,instruments,fullsettings)