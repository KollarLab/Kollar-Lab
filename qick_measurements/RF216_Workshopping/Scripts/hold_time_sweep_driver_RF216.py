# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from hold_time_sweep_RF216 import hold_time_sweep, get_hold_time_settings

instruments = {}


settings = get_hold_time_settings()
settings['scanname'] = 'flux_qubit_hold_sweep'

offset=-0.0035
SRS2.voltage_ramp(0.093+offset)
SRS3.voltage_ramp(0.21)

settings['cav_freq']  = 6.1013e9 #7.8392e9
settings['cav_mixer_detuning'] = 200
settings['cav_gain'] = 0.5

settings['qub_freq']  = 3.5e9
settings['qub_mixer_detuning'] = -250
settings['qub_gain']     = 0.5

settings['filter'] = 'not_qubit'

#Sweep Parameters
settings['hold_start']      = 0.01#0.01
settings['hold_stop']       = 0.2#0.1
settings['hold_points']     = 30 

#ADC settings
settings['reps']      = 80000
settings['rounds']  = 1


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = hold_time_sweep(soc,soccfg,instruments,fullsettings)