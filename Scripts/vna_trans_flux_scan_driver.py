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
settings['scanname']    = 'flux_Scan'
settings['meas_type']   = 'trans_flux_scan'
settings['project_dir'] = r'Z:\Data'

#Sweep parameter
settings['CAV_Attenuation'] = 30

settings['start_volt']  = 0.4
settings['stop_volt']   = 0.9
settings['volt_points'] = 15

settings['start_power'] = -50
settings['stop_power'] = -48
settings['power_points'] = 1

settings['avg_times'] = np.array([0.8])

#VNA settings
settings['channel']  = 1
settings['avg_time'] = 1
settings['measurement'] = 'S21'
settings['start_freq']  = 7.576e9-40e6 
settings['stop_freq']   = 7.576e9+40e6 
settings['freq_points'] = 501
settings['ifBW'] = 4e3

vna_trans_flux_scan(instruments, settings)