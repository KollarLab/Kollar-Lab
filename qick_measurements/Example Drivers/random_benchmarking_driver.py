# -*- coding: utf-8 -*-
"""
Created on Mon July 1 10:28:28 2024

@author: Ruthie Vogel
"""

from qick_measurements.random_benchmarking import random_bench, get_RB_settings

instruments = {}
instruments['LO'] = logen

settings = get_RB_settings()
settings['scanname'] = 'Zach_qubit_hi_pow'

settings['cav_freq']  = 6.8054e9 
settings['cav_gain'] = 5000

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 500
settings['qub_freq'] = 8.4808e9
settings['qub_gain'] = 5300

# RB settings
settings['pulse_schedule'] = ['X']
settings['pulse_data']     = [] 
settings['buffer']         = .05 #[us]
settings['num_gates']      = 10
settings['used_gates']     = 1
settings['gen_new']        = False
settings['pulse_plot']     = []



#settings['subtract_background'] = True 
#Currently no background subtraction

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = random_bench(soc,soccfg,instruments,fullsettings)