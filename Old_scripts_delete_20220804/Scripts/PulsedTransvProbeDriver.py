# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:29:44 2020

@author: Kollarlab
"""

from PulseTransVProbe import GetDefaultSettings,  PulsedTransVProbe

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

settings = GetDefaultSettings()
settings['scanname'] = 'F4_137mT_100mK'
settings['project_dir']  = r'C:\Users\Kollarlab\Desktop\Data\PhosporusDopedSilicon'

settings['start_freq']      = 3.881e9
settings['stop_freq']       = 3.883e9
settings['freq_points']     = 201

settings['cavity_power']    = -30
    
settings['probe_freq']      = 3.85e9
settings['start_power']     = -30
settings['stop_power']      = 0
settings['power_points']    = 2

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 1e3
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['timeout']          = 30

##Pulse settings
settings['meas_pos']    = 16.0e-6
settings['meas_window'] = 140e-6

PulsedTransVProbe(instruments, settings)

# Turn outputs off at the end of measurement
qubitgen.Output = 'Off'
cavitygen.Output = 'Off'
vna.output = 'Off'