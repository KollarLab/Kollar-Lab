# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 20:06:44 2021

@author: Kollarlab
"""
import numpy as np

from cw_measurements.vna_trans_flux_scan_Keysight_source import get_default_settings, vna_trans_flux_scan

instruments = {}
instruments['VNA'] = vna
instruments['DCsupply'] = Dual_gen

settings = get_default_settings()

#Save location

settings['scanname']    = 'Cav_2_FBL_2_Keysight_source_take_1'
settings['meas_type']   = 'trans_flux_scan'

#Sweep parameter

settings['start_voltage']  = -3
settings['stop_voltage']   = 3
settings['voltage_points'] = 31
settings['start_power'] = -64
settings['stop_power'] = -64
settings['power_points'] = 1

#settings['avg_times'] = np.array([30, 30, 22, 20, 20, 13, 13, 11, 9, 8, 8, 6, 6, 6, 4, 4, 4, 4, 10, 10, 10, 10, 10])

settings['avg_times'] = np.array([30])

#VNA settings
center = 6.54531e9
span = 40e6
settings['channel']  = 1
#settings['avg_time'] = 12 #looks like this is ignored
settings['measurement'] = 'S21'
settings['start_freq']  = center-span/2
settings['stop_freq']   = center+span/2
settings['freq_points'] = 4001
settings['ifBW'] = 1e3
settings['mode'] = 'MOV'

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings


Dual_gen.Ch2_dc_voltage_ramp(0)
Dual_gen.Ch2_output = 1



vna_trans_flux_scan(instruments, fullsettings, channel=2)




