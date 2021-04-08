# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 15:33:46 2020

@author: Kollarlab
"""

import numpy as np
import time
import matplotlib.pyplot as plt
import userfuncs as uf

pi = np.pi

def remove_IQ_ellipse(Is, Qs, axes, center, phi):
   
    Isprime = np.cos(phi)*(Is-center[0]) + np.sin(phi)*(Qs-center[1]) + center[0]
    Qsprime = -np.sin(phi)*(Is-center[0]) + np.cos(phi)*(Qs-center[1]) 
    Qsprime = Qsprime*axes[0]/axes[1] + center[1]
    return Isprime, Qsprime

center = [0.027, -0.034]
phi = 0.044*2*pi/180
axes = [0.024, 0.018]

vna.freq = 3.8816e9
cavitygen.Freq = 3.8816e9
    
time.sleep(0.1)
      
card.ArmAndWait()
rawI,rawQ = card.ReadAllData()

Ip, Qp = remove_IQ_ellipse(rawI, rawQ, axes, center, phi)
DC_I = np.mean(Ip[-50:])
DC_Q = np.mean(Qp[-50:])
Idat = Ip-DC_I
Qdat = Qp-DC_Q
            
        
            
amp = np.sqrt(Idat**2+Qdat**2)
phase = np.arctan2(Qdat, Idat)*180/np.pi

   
plt.figure(1)
plt.clf()
plt.plot(Idat)
plt.plot(Qdat)

plt.show()