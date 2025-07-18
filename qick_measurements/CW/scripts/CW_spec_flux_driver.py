# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: jhyang
"""

from CW_spec_flux import cw_spec_flux, get_cw_spec_flux_settings
    
settings = get_cw_spec_flux_settings()
settings['scanname'] = 'Spec_flux_testing'

settings['start_voltage']  = 0
settings['stop_voltage']   = 0.1
settings['voltage_points'] = 2
settings['stability'] = False

autoscan = settings['autoscan']
autoscan['freq_start']     = 7.2e9
autoscan['freq_stop']      = 7.25e9
autoscan['freq_points']    = 101
autoscan['reps']           = 10
autoscan['soft_avgs']      = 1

#settings['cav_freq']  = 7.230796e9  #4.389545e9 + 1e6 #7.08979e9 #7.1995e9 #7.8392e9 original
settings['cav_gain'] = 1000
settings['cav_meas_window'] = 10000
settings['cav_pulse_len'] = 10

qub_center = 3.6117e9 #5.33866e9
span = 50e6 
freq_points = 51

settings['qub_gain'] = 500
settings['qub_len'] = 10 #us
settings['qub_readout_length'] = 5000

#Sweep Parameters
settings['freq_start']      = qub_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 

#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 1


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

instruments = {}
#instruments['DCSupply'] = 'SRS'
full_data = cw_spec_flux(soc,soccfg,instruments,fullsettings)