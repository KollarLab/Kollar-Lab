# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 19:47:15 2020

@author: Kollarlab
"""
import numpy as np
from pulsed_spec import get_default_settings, pulsed_spec

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

vna.SweepType = 'CW'



settings = get_default_settings()
settings['scanname'] = 'q2_fine_hold_script'
settings['project_dir']  = r'Z:\Data\HouckQuadTransmon'

settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 0

#Cavity parameters
settings['CAVpower']        = -45
settings['CAV_freq']        = 7.89725e9

#Qubit parameters
q_freq = 5.363e9
span = 100e6
settings['start_freq']      = q_freq - span/2
settings['stop_freq']       = q_freq + span/2
settings['freq_points']     = 101

settings['start_power']     = -35
settings['stop_power']      = 5
settings['power_points']    = 81

#triggergen.Freq       = '100 kHz'

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 5e3
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['meas_window']      = 1e-6
settings['empirical_delay']  = np.floor(375e-9/32e-9) * 32e-9  #1e-6
settings['timeout']          = 30

##Pulse 
settings['Quasi_CW'] = False
settings['meas_pos']    = np.floor(7e-6/32e-9) * 32e-9 #15e-6
settings['pulse_delay'] = 100e-9
settings['pulse_width'] = 80e-9

pulsed_spec(instruments, settings)




#settings = get_default_settings()
#settings['scanname'] = 'plasmon_recheck'
#settings['project_dir']  = r'Z:\Data\HouckDualHangerFluxonium'
#
#settings['CAV_Attenuation'] = 30
#settings['Qbit_Attenuation'] = 0
#
##Cavity parameters
#settings['CAVpower']        = -45
#settings['CAV_freq']        = 7.57630e9
#
##Qubit parameters
#q_freq = 4.615e9
#span = 100e6
#settings['start_freq']      = 8.81e9
#settings['stop_freq']       = 8.84e9
#settings['freq_points']     = 9
#
#settings['start_power']     = -30
#settings['stop_power']      = -30
#settings['power_points']    = 1
#
##Card settings
#settings['segments']         = 1
#settings['reads']            = 1
#settings['averages']         = 60e3
#settings['activeChannels']   = [1,2]
#settings['channelRange']     = 0.5
#settings['sampleRate']       = 2e9/8
#settings['meas_window']      = 5e-6
#settings['empirical_delay']  = np.floor(375e-9/32e-9) * 32e-9  #1e-6
#settings['timeout']          = 30
#
###Pulse 
#settings['Quasi_CW'] = True
#settings['meas_pos']    = np.floor(7e-6/32e-9) * 32e-9 #15e-6
#settings['pulse_delay'] = 200e-9
#settings['pulse_width'] = 1.4e-7
#
#pulsed_spec(instruments, settings)

