# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 19:47:15 2020

@author: Kollarlab
"""

from pulsed_spec import get_default_settings, pulsed_spec

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = SMB
instruments['card'] = card
instruments['AWG'] = hdawg

settings = get_default_settings()
settings['scanname'] = 'q4_PulseLengthSweep_5ns_widersweep'
settings['saveDir']  = r'Z:\Data\HouckQuadTransmon\PulsedSpec\20201208'

#Cavity parameters
settings['CAVpower']        = -18
settings['CAV_freq']        = 8.126e9

#Qubit parameters
q_freq = 4.20554796e9
span = 200e6
settings['start_freq']      = q_freq-span/2 
settings['stop_freq']       = q_freq+span/2
settings['freq_points']     = 30

settings['start_power']     = 5
settings['stop_power']      = 20
settings['power_points']    = 31

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 5e3
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['meas_window']      = 10e-6
settings['timeout']          = 30

##Pulse settings
settings['meas_pos'] = 80e-6
settings['pulse_delay'] = 200e-9
settings['pulse_width'] = 10e-9

pulsed_spec(instruments, settings)