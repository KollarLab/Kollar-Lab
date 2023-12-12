# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 11:32:50 2023

@author: kollarlab
"""

from T2 import meas_T2, get_T2_settings

instruments = {}
instruments['LO'] = logen


settings = get_T2_settings()
settings['scanname'] = 'T2_test_6'

settings['cav_freq']  = 7.31266e9 
settings['meas_gain'] = 1000

settings['qub_freq']  = 5.87781e9
settings['qub_gain']  = 6000

settings['Tau_min']    = 100e-9
settings['Tau_max']    = 1e-6
settings['Tau_points'] = 26
settings['T2_guess']   = 1e-6
settings['detuning']   = 4e6

#ADC settings
settings['reps']       = 1
settings['soft_avgs']  = 1000


settings['subtract_background'] = True

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = meas_T2(soc,soccfg,instruments,fullsettings)