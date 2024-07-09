# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from qick_measurements.quasi_cw import quasi_cw, get_quasi_cw_settings

instruments = {}
instruments['LO'] = logen

settings = get_quasi_cw_settings()
settings['scanname'] = 'Zach_qubit_hi_pow'

settings['cav_freq']  = 7.8392e9 
settings['cav_gain'] = 6000

settings['qub_freq']     = 4e9

#Sweep Parameters
settings['gain_start']      = 1000
settings['gain_stop']       = 30000
settings['gain_points']     = 11

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 500


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = quasi_cw(soc,soccfg,instruments,fullsettings)