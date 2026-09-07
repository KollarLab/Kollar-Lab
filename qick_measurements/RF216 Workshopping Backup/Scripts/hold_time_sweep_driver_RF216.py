# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from hold_time_sweep_RF216 import hold_time_sweep, get_hold_time_settings

instruments = {}


settings = get_hold_time_settings()
settings['scanname'] = 'Zach_qubit_hi_pow'

settings['cav_freq']  = 5e9 #7.8392e9
settings['cav_mixer_detuning'] = -250
settings['cav_gain'] = 0.2

settings['qub_freq']  = 5e9
settings['qub_mixer_detuning'] = -250
settings['qub_gain']     = 0.5


settings['filter'] = True

#Sweep Parameters
settings['hold_start']      = 0.1
settings['hold_stop']       = 4
settings['hold_points']     = 20 

#ADC settings
settings['reps']      = 500
settings['rounds']  = 1


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = hold_time_sweep(soc,soccfg,instruments,fullsettings)