# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:29:44 2020

@author: Kollarlab
"""
import numpy as np
from pulsed_trans import get_default_settings, pulsed_trans

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

vna.SweepType = 'CW'

#qubitgen.Power = -30
#qubitgen.Freq = 8.824e9
#qubitgen.Output = 'On'

settings = get_default_settings()
settings['scanname'] = 'transmissionRecheck'
settings['project_dir']  = r'Z:\Data\HouckQuadTransmon'

settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10

cav_center = 7.897e9
span = 10e6
settings['start_freq']      = cav_center-span/2
settings['stop_freq']       = cav_center+span/2
settings['freq_points']     = 31

settings['start_power']     = -45
settings['stop_power']      = -45
settings['power_points']    = 1

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 5e3
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['empirical_delay']  = np.floor(375e-9/32e-9) * 32e-9  #1e-6
settings['timeout']          = 30

##Pulse settings
settings['meas_pos']    = np.floor(1e-6/32e-9) * 32e-9 #15e-6
settings['meas_window'] = 2.5e-6

pulsed_trans(instruments, settings)