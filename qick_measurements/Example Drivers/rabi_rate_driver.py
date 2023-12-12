# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from rabi_rate import Rabi_sweep, get_Rabi_settings

instruments = {}
instruments['LO'] = logen

settings = get_Rabi_settings()
settings['scanname'] = 'Longer_sigma_Rabi'

settings['cav_freq']  = 7.683e9 
settings['cav_gain'] = 5000

# qub_center = 5.33e9
# span = 30e6 
# freq_points = 16

settings['qub_freq']       = 5.332e9

settings['gain_start']     = 1000
settings['gain_step']      = 1000
settings['gain_points']    = 30

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 500

# settings['freq_start']      = qub_center-span/2
# settings['freq_step']       = span/(freq_points-1)
# settings['freq_points']     = freq_points 

#settings['subtract_background'] = True 
#Currently no background subtraction

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = Rabi_sweep(soc,soccfg,instruments,fullsettings)