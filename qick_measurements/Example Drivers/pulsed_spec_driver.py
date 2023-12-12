# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from pulsed_spec import pulsed_spec, get_spec_settings

instruments = {}
instruments['LO'] = logen

settings = get_spec_settings()
settings['scanname'] = 'Zach_qubit_hi_pow'

settings['cav_freq']  = 7.8392e9 
settings['cav_gain'] = 6000

qub_center = 5.33866e9
span = 50e6 
freq_points = 26

settings['gain_start']     = 1000
settings['gain_step']      = 3000
settings['gain_points']    = 11

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 500

settings['freq_start']      = qub_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 

#settings['subtract_background'] = True 
#Currently no background subtraction

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = pulsed_spec(soc,soccfg,instruments,fullsettings)