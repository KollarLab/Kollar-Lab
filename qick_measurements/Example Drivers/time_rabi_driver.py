# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from time_rabi import time_Rabi_sweep, get_time_Rabi_settings

instruments = {}
instruments['LO'] = logen

settings = get_time_Rabi_settings()
settings['scanname'] = 'Hold_time_test_11'

settings['cav_freq']  = 7.31266e9
settings['cav_gain'] = 1000

# qub_center = 5.33e9
# span = 30e6 
# freq_points = 16

settings['qub_freq']       = 5.87781e9
settings['qub_gain']       = 6000

settings['hold_start']     = 0.01e-6
settings['hold_step']      = 0.002e-6
settings['hold_points']    = 30

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 1500

# settings['freq_start']      = qub_center-span/2
# settings['freq_step']       = span/(freq_points-1)
# settings['freq_points']     = freq_points 

#settings['subtract_background'] = True 
#Currently no background subtraction

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = time_Rabi_sweep(soc,soccfg,instruments,fullsettings)