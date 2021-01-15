# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:29:44 2020

@author: Kollarlab
"""

from pulsed_trans import get_default_settings, pulsed_trans

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

vna.SweepType = 'CW'

settings = get_default_settings()
settings['scanname'] = 'testingjunk'
settings['project_dir']  = r'Z:\Data\deleteme'

settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10

settings['start_freq']      = 7.7e9
settings['stop_freq']       = 7.8e9
settings['freq_points']     = 30

settings['start_power']     = -20
settings['stop_power']      = -18
settings['power_points']    = 2

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 1e3
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['empirical_delay']  = 1e-6
settings['timeout']          = 30

##Pulse settings
settings['meas_pos']    = 15e-6
settings['meas_window'] = 10e-6

pulsed_trans(instruments, settings)