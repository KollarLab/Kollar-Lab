# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from qick_measurements.quasi_cw import quasi_cw, get_quasi_cw_settings

instruments = {}
instruments['LO'] = logen

settings = get_quasi_cw_settings()
settings['scanname'] = 'Zach_qubit_hi_pow'

settings['cav_freq']  = 7.8392e9 
settings['cav_gain'] = 6000

qub_center = 5.33866e9
span = 50e6 
freq_points = 26

settings['qub_gain']     = 1000
settings['quasi_CW_len'] = 10 #us

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

full_data = quasi_cw(soc,soccfg,instruments,fullsettings)