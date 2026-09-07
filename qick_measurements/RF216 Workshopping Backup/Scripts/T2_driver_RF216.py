# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 11:32:50 2023

@author: kollarlab
"""

from T2_RF216 import meas_T2, get_T2_settings

instruments = {}



settings = get_T2_settings()
settings['scanname'] = 'T2_test_6'

settings['cav_freq']  = 7.31266e9 
settings['cav_mixer_detuning'] = -250
settings['cav_gain'] = 1.0

settings['qub_freq']  = 5.87781e9
settings['qub_mixer_detuning'] = -250
#settings['qub_gain']  = 1.0 #Set in experiment_globals as piGain instead

settings['Tau_min']    = 100e-9
settings['Tau_max']    = 1e-6
settings['Tau_points'] = 26
settings['T2_guess']   = 1e-6
settings['detuning']   = 4e6

#ADC settings
settings['reps']       = 1000
settings['rounds']  = 1

settings['phase_reset'] = False #Not currently implemented
settings['filter']      = True

settings['subtract_background'] = False #Don't even try it man

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = meas_T2(soc,soccfg,instruments,fullsettings)