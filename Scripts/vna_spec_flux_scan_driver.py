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

spec_powers = [5]
cavity_powers = [-65]

for c_power in cavity_powers:
    for s_power in spec_powers:
        fullsettings = get_default_settings()
        settings = fullsettings['spec']
        settings['scanname']    = 'High_freq_sweep_{}_spec_{}_cavity'.format(s_power, c_power)
        settings['meas_type']   = 'spec_flux_scan'
        settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'
        
        settings['start_voltage'] = 0.1
        settings['stop_voltage'] = 0.35
        settings['voltage_points'] = 126
        
        settings['RFport'] = 3
        settings['Mport'] = 2
        settings['CAVport'] = 1
        
        settings['channel'] = 1
        settings['avg_time'] = 30
        settings['measurement'] = 'S21'
        settings['start_freq'] = 12.5e9
        settings['stop_freq'] = 14e9
        settings['freq_points'] = 751
        settings['CAVpower'] = c_power
        settings['RFpower'] = s_power
        settings['ifBW'] = 1e3
        
        settings['CAV_Attenuation'] = 30
        settings['Qbit_Attenuation'] = 10
        
        autoscan_settings = fullsettings['autoscan']
        autoscan_settings['freq_points'] = 201
        autoscan_settings['ifBW'] = 1e3
        autoscan_settings['avg_time'] = 10
        autoscan_settings['start_freq'] = 6.52e9
        autoscan_settings['stop_freq'] = 6.57e9
        autoscan_settings['RFpower'] = settings['CAVpower']
        autoscan_settings['background_subtract'] = False
        
        fullsettings['spec'] = settings
        fullsettings['autoscan'] = autoscan_settings
        
        vna_spec_flux_scan(instruments, fullsettings)
        
spec_powers = [5]
cavity_powers = [-65]

for c_power in cavity_powers:
    for s_power in spec_powers:
        fullsettings = get_default_settings()
        settings = fullsettings['spec']
        settings['scanname']    = 'Integer_flux_below_cavity_{}_spec_{}_cavity'.format(s_power, c_power)
        settings['meas_type']   = 'spec_flux_scan'
        settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'
        
        settings['start_voltage'] = 0.1
        settings['stop_voltage'] = 0.35
        settings['voltage_points'] = 126
        
        settings['RFport'] = 3
        settings['Mport'] = 2
        settings['CAVport'] = 1
        
        settings['channel'] = 1
        settings['avg_time'] = 30
        settings['measurement'] = 'S21'
        settings['start_freq'] = 4e9
        settings['stop_freq'] = 5.5e9
        settings['freq_points'] = 751
        settings['CAVpower'] = c_power
        settings['RFpower'] = s_power
        settings['ifBW'] = 1e3
        
        settings['CAV_Attenuation'] = 30
        settings['Qbit_Attenuation'] = 10
        
        autoscan_settings = fullsettings['autoscan']
        autoscan_settings['freq_points'] = 201
        autoscan_settings['ifBW'] = 1e3
        autoscan_settings['avg_time'] = 10
        autoscan_settings['start_freq'] = 6.52e9
        autoscan_settings['stop_freq'] = 6.57e9
        autoscan_settings['RFpower'] = settings['CAVpower']
        autoscan_settings['background_subtract'] = False
        
        fullsettings['spec'] = settings
        fullsettings['autoscan'] = autoscan_settings
        
        vna_spec_flux_scan(instruments, fullsettings)