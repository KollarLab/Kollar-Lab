# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:42:04 2021

@author: Kollarlab
"""
from rabi_chevron import get_default_settings, rabi_chevron
import numpy as np

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

vna.SweepType = 'CW'

settings = get_default_settings()
settings['scanname'] = 'rabi_chevron_fine'
settings['meas_type'] = 'rabi_chevron'
settings['project_dir'] = r'Z:\Data\HouckQuadTransmon'

settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 0

#Cavity parameters
settings['CAVpower']        = -45
settings['CAV_freq']        = 7.89725e9
settings['Q_power']         = -19

#Qubit parameters
settings['start_freq']      = 5.345*1e9  
settings['stop_freq']       = 5.375*1e9 
settings['freq_points']     = 61

settings['start_time']     = 10e-9
settings['stop_time']      = 1e-6
settings['time_points']    = 161

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 5e3
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['empirical_delay']  = np.floor(375e-9/32e-9) * 32e-9
settings['meas_window']      = 10e-6
settings['timeout']          = 30

##Pulse settings
settings['meas_pos']    = np.floor(5e-6/32e-9) * 32e-9
settings['pulse_delay'] = 200e-9
settings['pulse_width'] = 80e-9

rabi_chevron(instruments, settings)