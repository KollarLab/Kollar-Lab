# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 12:13:02 2020

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
#plt.figure()
#plt.clf()

numpoints = 50
#freqs = numpy.linspace(3.88, 3.886, numpoints)*1e9
freqs = [3.89*1e9]*numpoints

Is = numpy.zeros(numpoints)
Qs = numpy.zeros(numpoints)

Ips = numpy.zeros(numpoints)
Qps = numpy.zeros(numpoints)

for ind, freq in enumerate(freqs):

    vna.freq = 3.8e9
    cavitygen.Freq = 3.8e9
    
    vna.freq = freq
    cavitygen.Freq = freq
    
    time.sleep(0.1)
      
    card.ArmAndWait()
    I,Q = card.ReadAllData()
    
#    plt.plot(xaxis*1e6, np.sqrt(I[0]**2+Q[0]**2))
    Idat = np.mean(I[0:data_window])
    Qdat = np.mean(Q[0:data_window])
    Is[ind] = Idat
    Qs[ind] = Qdat
    
    Ip, Qp = remove_IQ_ellipse(Idat, Qdat, axes, center, phi)
    
    Ips[ind] = Ip
    Qps[ind] = Qp


plt.figure(1)
plt.clf()
ax = plt.subplot(1,1,1)
plt.plot(Ips, Qps, 'b.')
plt.plot(Is, Qs, 'r.')
#ax.set_xlim([-0.15,0.15])
#ax.set_ylim([-0.15,0.15])
ax.set_aspect('equal')

plt.show()


#plt.figure(7)
#plt.clf()
#
#axes, center, phi = uf.fitEllipse(Is,Qs, verbose = True)
#xx, yy = uf.make_elipse(axes,  center, phi, 150)
#        
#ax = plt.subplot(1,1,1)
#plt.plot(Is, Qs, linestyle = '', marker = 'o', markersize = 5, color = 'mediumblue')
#plt.plot(xx, yy, color = 'firebrick')
#
#
#    
#
#
#Qsprime = Qsprime*axes[0]/axes[1]
#axesprime, centerprime, phiprime = uf.fitEllipse(Isprime,Qsprime, verbose = True)
#xx, yy = uf.make_elipse(axesprime,  centerprime, phiprime, 150)
#plt.plot(xx +  center[0] , yy+ center[1], color = 'dodgerblue')
#plt.plot(Isprime + center[0], Qsprime+ center[1], linestyle = '', marker = 'o', markersize = 5, color = 'darkgoldenrod')
#
#ax.set_aspect('equal')
#
#plt.show()





