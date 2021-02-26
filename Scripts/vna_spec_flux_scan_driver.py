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
settings['scanname']    = 'plasmon_check'
settings['meas_type']   = 'spec_flux_scan'
settings['project_dir'] = r'Z:\Data\HouckDualHangerFluxonium'

settings['start_voltage'] = 0.112
settings['stop_voltage'] = 0.116
settings['voltage_points'] = 5

settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1

settings['channel'] = 1
settings['avg_time'] = 30
settings['measurement'] = 'S21'
settings['start_freq'] = 8.75e9
settings['stop_freq'] = 8.95e9
settings['freq_points'] = 201
settings['CAVpower'] = -55
settings['RFpower'] = -30
settings['ifBW'] = 1e3

settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10

autoscan_settings = fullsettings['autoscan']
autoscan_settings['freq_points'] = 201
autoscan_settings['ifBW'] = 1e3
autoscan_settings['avg_time'] = 20
autoscan_settings['start_freq'] = 7.55e9
autoscan_settings['stop_freq'] = 7.6e9
autoscan_settings['RFpower'] = settings['CAVpower']
autoscan_settings['background_subtract'] = False

fullsettings['spec'] = settings
fullsettings['autoscan'] = autoscan_settings

vna_spec_flux_scan(instruments, fullsettings)
