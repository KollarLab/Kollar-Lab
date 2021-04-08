# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 16:42:56 2020

@author: Kollarlab
"""
import time
import numpy as np

from T1 import GetDefaultSettings, meas_T1

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

settings = GetDefaultSettings()
settings['scanname'] = 'q2_T1'
settings['project_dir'] = r'Z:\Data\HouckQuadTransmon'
settings['meas_type'] = 'Tmeas'

settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 0

settings['Q_Freq'] = 5.360e9
settings['Q_Power'] = -19
settings['CAV_Freq'] = 7.89725e9
settings['CAV_Power'] = -45

Meas_pos = 12e-6
#Meas_pos = np.floor(7e-6/32e-9) * 32e-9 #15e-6
#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 5e3
settings['activeChannels']   = [1,2]
settings['sampleRate']       = 2e9/8
settings['trigger_buffer']   = np.floor((Meas_pos+350e-9)/32e-9) * 32e-9
settings['meas_window']      = 5e-6

settings['Tau_min'] = 1e-8
settings['Tau_max'] = 10e-6
settings['Tau_points'] = 25
settings['pulse_width'] = 80e-9

settings['Measurement_pos'] = Meas_pos


T1 = meas_T1(instruments, settings)