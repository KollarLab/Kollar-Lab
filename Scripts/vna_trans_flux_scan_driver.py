# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 20:06:44 2021

@author: Kollarlab
"""
import numpy as np

from vna_trans_flux_scan import get_default_settings, vna_trans_flux_scan

instruments = {}
instruments['VNA'] = vna
instruments['DCsupply'] = SRS

settings = get_default_settings()

#Save location
settings['scanname']    = 'cavity_wide_scan'
settings['meas_type']   = 'trans_flux_scan'
settings['project_dir'] = r'Z:\Data\HouckDualHangerFluxonium'


#Sweep parameter
settings['CAV_Attenuation'] = 30

settings['start_voltage']  = -2
settings['stop_voltage']   = 2
settings['voltage_points'] = 401

settings['start_power'] = -55
settings['stop_power'] = -55
settings['power_points'] = 1

settings['avg_times'] = np.array([10])

center = 7.58e9
span = 60e6
#VNA settings
settings['channel']  = 1
settings['avg_time'] = 20
settings['measurement'] = 'S21'
settings['start_freq']  = center-span/2 
settings['stop_freq']   = center+span/2 
settings['freq_points'] = 201
settings['ifBW'] = 1e3

vna_trans_flux_scan(instruments, settings)