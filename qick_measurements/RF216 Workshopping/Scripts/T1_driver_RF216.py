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
 
settings['scanname'] = 'T1_test_tmon_qubit'

offset=-0.0035
SRS2.voltage_ramp(1.449)#0.093+offset)
SRS3.voltage_ramp(0.21)

settings['cav_freq']  = 7.1765e9 #7.8392e9
settings['cav_mixer_detuning'] = 200
settings['cav_gain'] = 0.5

settings['qub_freq']  = 9.3165e9
settings['qub_mixer_detuning'] = -250
#settings['qub_gain']  = 0.5 #Set in experiment_globals as piGain instead

settings['Tau_min']    = 10e-9
settings['Tau_max']    = 1e-6
settings['Tau_points'] = 11
settings['T1_guess']   = 1e-6

#ADC settings
settings['reps']      = 400000
settings['rounds']  = 1

settings['phase_reset'] = False #Not currently implemented
settings['filter']      = 'not_qubit' 


settings['subtract_background'] = False # Don't even try it man

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings
    
full_data = meas_T1(soc,soccfg,instruments,fullsettings) 
T1s.append(full_data[0])
    

#Currently no background subtraction

