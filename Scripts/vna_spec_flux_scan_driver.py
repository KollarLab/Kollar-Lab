# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 20:22:55 2021

@author: Kollarlab
"""
import numpy as np

from vna_spec_flux_scan import get_default_settings, vna_spec_flux_scan

instruments = {}
instruments['VNA'] = vna
instruments['DCsupply'] = SRS

fullsettings = get_default_settings()
settings = fullsettings['spec']
settings['scanname']    = 'flux_Scan'
settings['meas_type']   = 'spec_flux_scan'
settings['project_dir'] = r'Z:\Data'

settings['start_voltage'] = -2
settings['stop_voltage'] =  2
settings['voltage_points'] = 251

settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1

settings['channel'] = 1
settings['avg_time'] = 10
settings['measurement'] = 'S21'
settings['start_freq'] = 6.5e9
settings['stop_freq'] = 9.5e9
settings['freq_points'] = 3001
settings['CAVpower'] = -55
settings['RFpower'] = -15
settings['ifBW'] = .2e3

settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10

autoscan_settings = fullsettings['autoscan']
autoscan_settings['freq_points'] = 501
autoscan_settings['ifBW'] = settings['ifBW']
autoscan_settings['avg_time'] = 15
autoscan_settings['start_freq'] = 7.6e9
autoscan_settings['stop_freq'] = 7.7e9
autoscan_settings['RFpower'] = settings['CAVpower']

fullsettings['spec'] = settings
fullsettings['autoscan'] = autoscan_settings

vna_spec_flux_scan(instruments, settings)