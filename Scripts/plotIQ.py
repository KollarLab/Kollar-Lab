# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 12:13:02 2020

@author: Kollarlab
"""

#7.168413793103
freq = 7.765e9
cavitygen.IQ.Mod = 'on'
#vna.inst.write('SOUR:FREQ:CW {}'.format(freq))
cavitygen.Freq = freq
qubitgen.Freq = freq
SMB.Freq = freq
cavitygen.Power = -18
cavitygen.output = 'On'

card.averages = 1e4
card.samples  = 4e3
card.triggerSlope = 'Rising'
card.SetParams()

xaxis = (numpy.array(range(card.samples))/card.sampleRate)
xaxis_us = xaxis*1e6

card.ArmAndWait()
I,Q = card.ReadAllData()

plt.figure()
plt.clf()
plt.plot(xaxis*1e6, I[0])
plt.plot(xaxis*1e6, Q[0])



#cavitygen.output = 'Off'