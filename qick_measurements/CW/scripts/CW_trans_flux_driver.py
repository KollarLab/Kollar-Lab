# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: jhyang
"""

from CW_trans_flux import cw_trans_flux, get_cw_trans_flux_settings
    
settings = get_cw_trans_flux_settings()
settings['scanname'] = 'Trans_flux_testing'

settings['start_voltage']  = 0
settings['stop_voltage']   = 0.1
settings['voltage_points'] = 5
settings['stability'] = False

cav_center = 3.6117e9 #7.230796e9 #4.389545e9 + 1e6 #7.08979e9 #7.1995e9 #7.8392e9 original
span = 50e6 
freq_points = 51

settings['cav_gain'] = 2000
settings['meas_window'] = 10000
settings['cav_pulse_len'] = 10
settings['initial_phase'] = -0.001041 -6.860519651869944e-07 # 0.13 # rad

#Sweep Parameters
settings['freq_start']      = cav_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 1


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = cw_trans_flux(soc,soccfg,instruments,fullsettings)