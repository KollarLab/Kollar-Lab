# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from quasi_cw_RF216 import quasi_cw, get_quasi_cw_settings

instruments = {}


settings = get_quasi_cw_settings()
settings['scanname'] = 'Zach_qubit_hi_pow'

settings['cav_freq']  = 5e9 #7.8392e9
settings['cav_mixer_detuning'] = -250
settings['qub_mixer_detuning'] = -250
settings['cav_gain'] = 0.1

settings['fit'] = False

qub_center = 5e9
span = 50e6 
freq_points = 51

settings['qub_gain']     = 0.5
settings['quasi_CW_len'] = 10 #us


settings['filter'] = True

#Sweep Parameters
settings['freq_start']      = qub_center-span/2
settings['freq_stop']       = qub_center+span/2
settings['freq_points']     = freq_points 

#ADC settings
settings['reps']      = 1
settings['rounds']  = 500


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = quasi_cw(soc,soccfg,instruments,fullsettings)