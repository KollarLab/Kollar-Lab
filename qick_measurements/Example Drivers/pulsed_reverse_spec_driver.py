# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 12:27:18 2025

@author: KollarLab
"""

from qick_measurements.Other_Scripts.pulsed_reverse_spec import pulsed_reverse_spec, get_reverse_spec_settings

#scans = 80

#from T1
instruments = {}
#instruments['LO'] = logen

settings = get_reverse_spec_settings()

 
settings['scanname'] = 'Reverse_Spec_Test_qubit_gain_sweep'
cav_center = 6.534627e9#6.53235e9 #7.2309e9 #7.08979e9 #7.1995e9 #7.19962e9 #q2 7.11736e9 #q3 7.1409e9
span = 5e6 
freq_points = 101

settings['cav_freq']  =   6.534622e9
settings['freq_start']      = cav_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 
settings['meas_gain'] =  2200

settings['qub_freq']  =  1.4416e9
settings['qub_gain']  =  1#14400

settings['Tau_min']    = 0e-9
settings['Tau_max']    = 1e-9
settings['Tau_points'] = 1
settings['T1_guess']   = 1e-6

settings['gain_start']     = 0#350#20000#8000 #7000
settings['gain_step']      = 14400
settings['gain_points']    = 2

#ADC settings
settings['reps']      = 501
settings['soft_avgs']  = 1
settings['dphi_df']   = 0 #5.96965198e7

settings['debug'] = False
settings['debug_time'] = 4 # 0 is defined as the start of the measurement pulse

settings['subtract_background'] = False

settings['phase_reset'] = False
settings['sweep_device'] = 'qubit' #can gain sweep 'cavity' or 'qubit'

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = pulsed_reverse_spec(soc,soccfg,instruments,fullsettings) 
