# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 11:32:50 2023

@author: kollarlab
"""

from T1_RF216 import meas_T1, get_T1_settings

#scans = 80

T1s = []


instruments = {}

settings = get_T1_settings()
 
settings['scanname'] = 'T1_test_corrected_2'

settings['cav_freq']  = 5e9
settings['cav_mixer_detuning'] = -250
settings['cav_gain'] = 0.1

settings['qub_freq']  = 5e9
settings['qub_mixer_detuning'] = -250
settings['qub_gain']  = 0.5

settings['Tau_min']    = 200e-9
settings['Tau_max']    = 40e-6
settings['Tau_points'] = 21
settings['T1_guess']   = 5e-6

#ADC settings
settings['reps']      = 100
settings['rounds']  = 1

settings['phase_reset'] = False #Not currently implemented
settings['filter']      = True


settings['subtract_background'] = False # Don't even try it man

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings
    
full_data = meas_T1(soc,soccfg,instruments,fullsettings) 
T1s.append(full_data[0])
    

#Currently no background subtraction

