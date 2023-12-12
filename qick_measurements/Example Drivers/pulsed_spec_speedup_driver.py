# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from pulsed_spec_speedup import pulsed_spec_sweep, get_spec_settings

instruments = {}
instruments['LO'] = logen

settings = get_spec_settings()
settings['scanname'] = 'test_scan_fine_4'

settings['cav_freq']  = 7.31266e9
settings['cav_gain'] = 1000

qub_center = 5.87115e9 + 1.537e6 + 0.658e6
span = 40e6
freq_points = 81

settings['freq_start']      = qub_center-span/2
settings['freq_stop']       = qub_center+span/2
settings['freq_points']     = freq_points 

settings['gain_start']     = 7000
settings['gain_stop']      = 8000
settings['gain_points']    = 3

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 1000




settings['subtract_background'] = False 
#Currently no background subtraction

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = pulsed_spec_sweep(soc,soccfg,instruments,fullsettings)