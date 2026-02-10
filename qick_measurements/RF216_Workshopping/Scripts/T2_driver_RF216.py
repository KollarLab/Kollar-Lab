# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 11:32:50 2023

@author: kollarlab
"""

from T2_RF216 import meas_T2, get_T2_settings

instruments = {}



settings = get_T2_settings()
settings['scanname'] = 'flux_qubit_T2_test'

offset=-0.0035
SRS2.voltage_ramp(0.093+offset)
SRS3.voltage_ramp(0.21)

settings['cav_freq']  = 6.10174e9#6.1014e9
settings['cav_mixer_detuning'] = 200
settings['cav_gain'] = 0.5

settings['qub_freq']  = 3.411e9
settings['qub_mixer_detuning'] = -250
#settings['qub_gain']  = 1.0 #Set in experiment_globals as piGain instead

settings['Tau_min']    = 100e-9
settings['Tau_max']    = 3e-6
settings['Tau_points'] = 26
settings['T2_guess']   = 1e-6
settings['detuning']   = 2e6

#ADC settings
settings['reps']       = 400000
settings['rounds']  = 1

settings['phase_reset'] = False #Not currently implemented
settings['filter']      = True

settings['subtract_background'] = False #Don't even try it man

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = meas_T2(soc,soccfg,instruments,fullsettings)