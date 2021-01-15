# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 16:42:56 2020

@author: Kollarlab
"""
import time

from T1WIP import GetDefaultSettings, meas_T1

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['qubitgen'] = qubitgen
instruments['LO'] = vna
instruments['card'] = card
instruments['AWG'] = hdawg

settings = GetDefaultSettings()
settings['scanname'] = 'q4_T1_test'
settings['project_dir'] = r'Z:\Data\HouckQuadTransmon'
settings['meas_type'] = 'Tmeas'


settings['Q_Freq'] = 4.20431e9
settings['Q_Power'] = -11
settings['CAV_Freq'] = 8.126e9
settings['CAV_Power'] = -18

Meas_pos = 80e-6
#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 20e3
settings['activeChannels']   = [1,2]
settings['sampleRate']       = 2e9/8
settings['trigger_buffer']   = Meas_pos
settings['meas_window']      = 20e-6

settings['Tau_min'] = 200e-9
settings['Tau_max'] = 50e-6
settings['Tau_points'] = 101

settings['Measurement_pos'] = Meas_pos

meas_repeats = 80

T1vec = numpy.zeros(meas_repeats)

for ind in range(meas_repeats):
    T1 = meas_T1(instruments, settings)
    T1vec[ind] = T1
    
plt.plot(T1vec)