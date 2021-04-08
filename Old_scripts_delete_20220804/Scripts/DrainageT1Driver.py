# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 16:38:37 2020

@author: Kollarlab
"""
import numpy as np
from DrainageT1 import T1_Drainage

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

settings = {}
    
settings['scanname'] = 'F4_B137mT_T1_100mK_overnight'
settings['project_dir']  = r'C:\Users\Kollarlab\Desktop\Data\PhosporusDopedSilicon'

#Cavity parameters
center = 3.8817e9
span = 2e6
num_points = 31

settings['CAVpower']        = -30
settings['CAV_freq']        = center

freqs = np.round(np.linspace(center-span/2, center+span/2, num_points),-3)
#Qubit parameters
settings['probe_freq']      = 3.835e9
settings['probe_power']     = 0

#data acquisition settings
settings['segments']         = 500
settings['reads']            = 250
settings['averages']         = 30

#Card settings
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['meas_window']      = 140e-6
settings['empirical_delay']  = 1e-6
settings['timeout']          = 30

settings['drive_time']       = 1
settings['meas_type'] = 'DrainageT1'

settings['meas_pos'] = 15.0e-6
settings['trigger_freq'] = 3e3 

for freq in freqs:
    print('Freq: {}'.format(freq))
    settings['CAV_freq'] = freq
    T1_Drainage(instruments, settings)

qubitgen.Output = 'Off'
cavitygen.Output = 'Off'
vna.output = 'Off'
