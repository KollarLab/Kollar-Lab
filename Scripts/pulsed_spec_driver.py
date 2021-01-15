# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 19:47:15 2020

@author: Kollarlab
"""

from pulsed_spec import get_default_settings, pulsed_spec

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

vna.SweepType = 'CW'

settings = get_default_settings()
settings['scanname'] = 'testing'
settings['saveDir']  = r'Z:\Data\deleteme'

settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10

#Cavity parameters
settings['CAVpower']        = -18
settings['CAV_freq']        = 8.126e9

#Qubit parameters
q_freq = 4.204e9
span = 100e6
settings['start_freq']      = q_freq-span/2 
settings['stop_freq']       = q_freq+span/2
settings['freq_points']     = 31

settings['start_power']     = -20
settings['stop_power']      = 0
settings['power_points']    = 11

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 5e3
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['meas_window']      = 10e-6
settings['empirical_delay']  = 1e-6
settings['timeout']          = 30

##Pulse settings
settings['meas_pos'] = 79e-6
settings['pulse_delay'] = 200e-9
settings['pulse_width'] = 80e-9

pulsed_spec(instruments, settings)