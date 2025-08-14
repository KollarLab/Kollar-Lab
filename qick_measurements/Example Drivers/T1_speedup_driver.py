# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 11:32:50 2023

@author: kollarlab
"""

from qick_measurements.T1 import meas_T1, get_T1_settings
from kollar_instruments.SGS import SGS100A

#scans = 80

T1s = []

logen = SGS100A('TCPIP0::rssgs100a110738::inst0::INSTR')

logen.Ref.Source = 'INT'
logen.Ref.Frequency = 10e6

instruments = {}
#instruments['LO'] = logen

settings = get_T1_settings()

 
settings['scanname'] = 'T1_with_new_amp'
settings['debug'] = False
settings['debug_time'] = 0

settings['cav_freq']  = 7.230796e9 #7.31266e9
settings['meas_gain'] = 5000 #1000

settings['qub_freq']  = 3.6117e9#5.87781e9
settings['qub_gain']  = 500#6000

settings['Tau_min']    = 200e-9
settings['Tau_max']    = 40e-6
settings['Tau_points'] = 21
settings['T1_guess']   = 5e-6

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 1000


settings['subtract_background'] = True

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings
    
full_data = meas_T1(soc,soccfg,instruments,fullsettings) 
T1s.append(full_data[0])
    

#Currently no background subtraction

