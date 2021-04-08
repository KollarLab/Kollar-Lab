# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:20:43 2021

@author: Kollarlab
"""
card.triggerDelay = 20e-6
qubitgen.Freq = 8.8246e9
qubitgen.Power = -30
qubitgen.Output = 0

cavitygen.Power = -15
cavitygen.Freq = 5.57630e9
vna.Freq = 7.57630e9

meas_window = 20e-6
card.samples = meas_window*card.sampleRate
card.averages = 10e3
card.SetParams()
xaxis = 1e6*np.linspace(0, card.samples, card.samples)/card.sampleRate
card.ArmAndWait()
I, Q = card.ReadAllData()

plt.plot(xaxis, np.sqrt(I[0]**2+Q[0]**2))
