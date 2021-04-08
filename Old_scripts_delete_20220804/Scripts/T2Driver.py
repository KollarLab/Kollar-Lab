# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 17:38:30 2021

@author: Kollarlab
"""

import time
import numpy as np

from T2 import GetDefaultSettings, meas_T2

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

settings = GetDefaultSettings()
settings['scanname'] = 'fluxon_T2_112mV_short_T2'
settings['project_dir'] = r'Z:\Data\HouckDualHangerFluxonium'
settings['meas_type'] = 'Tmeas'

settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10

detuning = 1.5e6

settings['Q_Freq'] = 4.62582e9
settings['Q_Power'] = 9
settings['CAV_Freq'] = 7.57630e9
settings['CAV_Power'] = -45
settings['detuning'] = detuning

Meas_pos = 219e-6
#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 5e3
settings['activeChannels']   = [1,2]
settings['sampleRate']       = 2e9/8
settings['trigger_buffer']   = Meas_pos
settings['meas_window']      = 30e-6

settings['Tau_min'] = 0e-9
settings['Tau_max'] = 50e-9
settings['Tau_points'] = 21

settings['Measurement_pos'] = Meas_pos
settings['wait_time'] = 1e-7
settings['pulse_width'] = 5e-7

T2, delta = meas_T2(instruments, settings)