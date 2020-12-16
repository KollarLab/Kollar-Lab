# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:29:44 2020

@author: Kollarlab
"""

from CleanPulsedtrans import GetDefaultSettings, PulsedTrans

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

settings = GetDefaultSettings()
settings['scanname'] = 'F4_138p4mT_cavity_30mK'
settings['project_dir']  = r'C:\Users\Kollarlab\Desktop\Data\PhosporusDopedSilicon'

settings['start_freq']      = 3.878e9
settings['stop_freq']       = 3.884e9
settings['freq_points']     = 151

settings['start_power']     = 0
settings['stop_power']      = 0
settings['power_points']    = 1

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 20e3
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['timeout']          = 30

##Pulse settings
settings['meas_pos']    = 16e-6
settings['meas_window'] = 140e-6

PulsedTrans(instruments, settings)