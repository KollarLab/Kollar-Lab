# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 19:47:15 2020

@author: Kollarlab
"""
import numpy as np
from pdh_measurements.pulsed_2Readouts_Stability import get_default_settings, pulsed_2Readouts_Stability
instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = holz.ch2
instruments['card'] = card
instruments['AWG'] = hdawg

exp_settings = get_default_settings()
exp_settings['scanname'] = 'NewCodeTest_2Readouts_Stability_Ground'

#Cavity parameters

exp_settings['CAV_power']       =  -50
exp_settings['CAV_freq']        = 6.770532e9 #6.77067e9

#Qubit paramters
exp_settings['Qbit_freq'] = 4.85045e9
exp_settings['Qbit_power'] = 5

# readout IQ settings
exp_settings['amp_list'] = [0,1,1]
exp_settings['phase_list'] = [0,0,0]
exp_settings['mod_freq'] = 18e6 #22e6 #5e6

# scan settings
exp_settings['start_freq']      = exp_settings['CAV_freq'] #- 1.5*exp_settings['mod_freq'] #cav_center-span/2
exp_settings['stop_freq']       = exp_settings['CAV_freq'] #+ 1.5*exp_settings['mod_freq'] #cav_center+span/2
exp_settings['freq_points']     = 1 #01
exp_settings['stability_points'] = 200
exp_settings['Qstate'] = 'Ground'  # should be Ground or Excited

# exp_settings['start_power']     = -2
# exp_settings['stop_power']      = 5
# exp_settings['power_points']    = 6

##Pulse 
exp_settings['Quasi_CW'] = False#True#False
exp_settings['reverse'] = False

exp_settings['segments']         = 1
exp_settings['reads']            = 1
exp_settings['averages']         = 5e3

exp_settings['subtract_background'] = False

settings = {}
settings['exp_globals'] = exp_globals
settings['exp_settings'] = exp_settings

spec_data = pulsed_2Readouts_Stability(instruments, settings)
