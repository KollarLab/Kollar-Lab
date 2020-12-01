# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 12:13:02 2020

@author: Kollarlab
"""

#7.168413793103
freq = 7.165e9
cavitygen.IQ.Mod = 'on'
vna.inst.write('SOUR:FREQ:CW {}'.format(freq))
cavitygen.Freq = freq
cavitygen.Power = 0
cavitygen.output = 'On'

card.averages = 1e3
card.triggerSlope = 'Rising'
card.SetParams()

card.ArmAndWait()
I,Q = card.ReadAllData()

plt.figure(1)
plt.clf()
plt.plot(xaxis*1e6, I[0])
plt.plot(xaxis*1e6, Q[0])



cavitygen.output = 'Off'