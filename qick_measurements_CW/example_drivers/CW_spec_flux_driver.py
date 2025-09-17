# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: jhyang
"""

from qick_measurements_CW.CW_spec_flux import cw_spec_flux, get_cw_spec_flux_settings

#SRS2.output = 'Off'
SRS2.voltage_ramp(-0.15)
SRS3.voltage_ramp(0)

instruments = {}
instruments['DCsupply'] = SRS3

full = get_cw_spec_flux_settings()

settings = full['spec']
settings['scanname'] = 'CFQsweep'
settings['start_voltage']  = 0.18
settings['stop_voltage']   = 0.2
settings['voltage_points'] = 10
# settings['stability'] = False #Add later

#settings['cav_freq']  = 7.230796e9  #4.389545e9 + 1e6 #7.08979e9 #7.1995e9 #7.8392e9 original
settings['cav_gain'] = 400
settings['meas_window'] = 10000
settings['cav_pulse_len'] = 10

qub_center = 7.8e9#5.25e9 #5.33866e9
span = 2e9 
freq_points = 1001#201

settings['qub_gain'] = 3000

#Sweep Parameters
settings['freq_start']      = 2.0e9#4.3e9#qub_center-span/2
settings['freq_stop']       = 3.9e9#4.8e9#qub_center+span/2 
settings['freq_points']     = freq_points 

settings['reps']      = 1
settings['soft_avgs']  = 1


autoscan = full['autoscan']
autoscan['freq_start']     = 6.53e9
autoscan['freq_stop']      = 6.54e9#6.55
autoscan['freq_points']    = 251
autoscan['reps']           = 1
autoscan['soft_avgs']      = 1

fullsettings = {}

fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = full

#instruments['DCSupply'] = 'SRS'
full_data = cw_spec_flux(soc,soccfg,instruments,fullsettings)

SRS2.voltage_ramp(0)
SRS3.voltage_ramp(0)