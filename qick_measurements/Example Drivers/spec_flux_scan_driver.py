# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:28:28 2022

@author: Kollarlab
"""

from qick_measurements.spec_flux_scan import spec_flux_scan, get_spec_flux_scan_settings

instruments = {} #I don't remember which supply
instruments['DCsupply'] = SRS2

settings = get_spec_flux_scan_settings()
settings['scanname'] = 'Power_Mod_Settings'
settings['meas_type'] = 'SpecFluxScan'
settings['stability'] = False

settings['start_voltage']  = -0.078#-0.053
settings['stop_voltage']   = -0.082#-0.056
settings['voltage_points'] = 3

settings['meas_gain'] = 2500

settings['qub_gain']     =  4000 #1000 #30000
settings['quasi_CW_len'] = 40 #us

#Sweep Parameters
qub_center = 5.87785e9 #5.87785e9 #2.9814e9 #3.02e9 #3.53754e9 + 0.01e6
span = 50e6 
freq_points = 201

settings['freq_start']      = 5.15e9#qub_center-span/2
settings['freq_stop']       = 5.35e9#qub_center+span/2
settings['freq_points']     = freq_points 

#ADC settings
settings['reps']      = 1500
settings['soft_avgs']  = 1

autoscan = {}
autoscan['freq_start']     = 4.3893e9
autoscan['freq_stop']      = 4.3903e9
autoscan['freq_points']    = 51
autoscan['reps']           = 301
autoscan['soft_avgs']      = 1

settings['autoscan'] = autoscan

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

data = spec_flux_scan(soc,soccfg,instruments,fullsettings)
#stuff = spec_flux_scan(soc,soccfg,instruments,fullsettings)


# Is = data[1]['mags'][0]*np.cos(data[1]['phases'][0])
# Qs = data[1]['mags'][0]*np.sin(data[1]['phases'][0])

# plt.figure(7)
# plt.clf()
# plt.plot(Is,Qs,'x')
