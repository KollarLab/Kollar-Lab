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

###test sequences
test_t1 = 2*(['X'] + 15*['I'])
test_toggle = 5*['I', 'X']
test_toggle2 = 5*['X', 'I']
test_pipi = ['I', 'X', 'X/2'] + 10*['X']

# RB settings
settings['gen_new']        = False # True when you want to generate new sequences
settings['run_full']       = False # True takes all IQ output data in a run
settings['full_seq_dict']  = {}
settings['pulse_plot']     = []
settings['pulse_data']     = []
settings['buffer']         = .01 #[us]

settings['num_gates']      = 10
settings['used_gates']     = 1 
# for testing, use np.arange(1,len(settings['full_seq_dict'][0])+1, 1)

#settings['subtract_background'] = True 
#Currently no background subtraction

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = random_bench(soc,soccfg,instruments,fullsettings)