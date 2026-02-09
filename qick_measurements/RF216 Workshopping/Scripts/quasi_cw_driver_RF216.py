# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from quasi_cw_RF216 import quasi_cw, get_quasi_cw_settings

instruments = {}


settings = get_quasi_cw_settings()
settings['scanname'] = 'finding_flux_qubit_no_amp'
offset=-0.0035
SRS2.voltage_ramp(0.093+offset)
SRS3.voltage_ramp(0.21)

settings['cav_freq'] = 6.1013e9#6.10208e9 
settings['cav_mixer_detuning'] = 200
settings['qub_mixer_detuning'] = -250
settings['cav_gain'] = 0.5

settings['fit'] = False

qub_center = 3.418e9#3.411e9#
span = 50e6 
freq_points = 51

settings['qub_gain']     = 0.005
settings['quasi_CW_len'] = 10 #us


settings['filter'] = 'not_qubit'

#Sweep Parameters
settings['freq_start']      = qub_center-span/2
settings['freq_stop']       = qub_center+span/2
settings['freq_points']     = freq_points 

#ADC settings
settings['reps']      = 500000
settings['rounds']  =1


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = quasi_cw(soc,soccfg,instruments,fullsettings)