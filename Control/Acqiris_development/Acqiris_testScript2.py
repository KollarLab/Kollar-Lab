# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:48:02 2020

@author: Kollarlab
"""

# import comtypes
import os
import time
#import subprocess

#import re
import scipy
import pylab
#import tarfile
#import struct
#import glob
import numpy
import time

#import pickle
#import datetime
#import itertools
import sys



from Acqiris import Acqiris




hardwareAddress = "PXI23::0::0::INSTR"

IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
sys.path.append(IVIbinPath)


############
#minimum acquisition
###############

card = Acqiris(hardwareAddress)

#card.activeChannels = [1]
card.activeChannels = [1,2]

card.timeout = 10

segs = 2

card.samples = 1024*600 #too long for averaging mode, but fine for regular
card.segments = segs
card.averages = 1
card.triggerDelay = 25*10**-6
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
if len(card.activeChannels) == 1:
    data1 = card.ReadAllData() #read data for the active channels.
else:
    data1, data2 = card.ReadAllData() #read data for the active channels.
ts = 10**6* scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate
#
print('Took regular data')
#
#card.ReInitialize()

card.samples = 1024*512
card.segments = segs
card.averages = 100
card.triggerDelay = 25*10**-6
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
if len(card.activeChannels) == 1:
    avData1 = card.ReadAllData() #read data for the active channels.
else:
    avData1, avData2 = card.ReadAllData() #read data for the active channels.
avTs = 10**6* scipy.arange(0, avData1.shape[1],1.)*1/card.sampleRate

print('Took regular averaged data')


card.samples = 1024*605 #too long for averaging mode, but fine for regular
card.segments = segs
card.averages = 1
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
if len(card.activeChannels) == 1:
    data1 = card.ReadAllData() #read data for the active channels.
else:
    data1, data2 = card.ReadAllData() #read data for the active channels.
ts = 10**6* scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate
#
print('Took regular data again')

pylab.figure(4)
pylab.clf()
ax = pylab.subplot(2,2,1)
pylab.plot(ts, data1[0,:], label = 'single')
pylab.plot(avTs, avData1[0,:], label = 'averaged')
pylab.title('Channel 1, Segment 1')
pylab.ylabel('Voltage')
pylab.xlabel('t ($\mu$s)')
ax.legend(loc = 'upper left')


if len(card.activeChannels) == 2:
    ax = pylab.subplot(2,2,2)
    pylab.plot(ts,data2[0,:], label = 'single')
    pylab.plot(avTs, avData2[0,:], label = 'averaged')
    ax.set_ylim([-0.04,0.04])
    pylab.title('Channel 2, Segment 1')
    pylab.ylabel('Voltage')
    pylab.xlabel('t ($\mu$s)')
    ax.legend(loc = 'upper left')


ax = pylab.subplot(2,2,3)
pylab.plot(ts, data1[1,:], label = 'single')
pylab.plot(avTs, avData1[1,:], label = 'averaged')
pylab.title('Channel 1, Segment 2')
pylab.ylabel('Voltage')
pylab.xlabel('t ($\mu$s)')
ax.legend(loc = 'upper left')


if len(card.activeChannels) == 2:
    ax = pylab.subplot(2,2,4)
    pylab.plot(ts, data2[1,:], label = 'single')
    pylab.plot(avTs, avData2[1,:], label = 'averaged')
    ax.set_ylim([-0.04,0.04])
    pylab.title('Channel 2, Segment 2')
    pylab.ylabel('Voltage')
    pylab.xlabel('t ($\mu$s)')
    ax.legend(loc = 'upper left')

#pylab.suptitle('Tah-Dah!')
pylab.tight_layout()
pylab.show()











print ('single seg data')

card.samples = 1024*600 #too long for averaging mode, but fine for regular
card.segments = 1
card.averages = 1
card.triggerDelay = 25*10**-6
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
if len(card.activeChannels) == 1:
    d1data1 = card.ReadAllData() #read data for the active channels.
else:
    d1data1, d1data2 = card.ReadAllData() #read data for the active channels.
d1ts = 10**6* scipy.arange(0, d1data1.shape[0],1.)*1/card.sampleRate
#

#


card.samples = 1024*512
card.segments = 1
card.averages = 100
card.triggerDelay = 25*10**-6
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
if len(card.activeChannels) == 1:
    d1avData1 = card.ReadAllData() #read data for the active channels.
else:
    d1avData1, d1avData2 = card.ReadAllData() #read data for the active channels.
d1avTs = 10**6* scipy.arange(0, d1avData1.shape[0],1.)*1/card.sampleRate



pylab.figure(5)
pylab.clf()
ax = pylab.subplot(1,2,1)
pylab.plot(d1ts, d1data1[:], label = 'single')
pylab.plot(d1avTs, d1avData1[:], label = 'averaged')
pylab.title('Channel 1, Only Segment')
pylab.ylabel('Voltage')
pylab.xlabel('t ($\mu$s)')
ax.legend(loc = 'upper left')


if len(card.activeChannels) == 2:
    ax = pylab.subplot(1,2,2)
    pylab.plot(d1ts,d1data2[:], label = 'single')
    pylab.plot(d1avTs, d1avData2[:], label = 'averaged')
#    ax.set_ylim([-0.04,0.04])
    pylab.title('Channel 2, Only Segment')
    pylab.ylabel('Voltage')
    pylab.xlabel('t ($\mu$s)')
    ax.legend(loc = 'upper left')




















#card.samples = 1024*500 #going down to something short enogh for regular
#card.segments = 100
#card.averages = 1
#card.triggerDelay = 25*10**-6
#card.SetParams() #pushes default to to card if the fields haven't been edited
#card.ArmAndWait() #initiates aquisition and calibrates if need be
#data1, data2 = card.ReadAllData() #read data for the active channels.
#ts2 = 10**6* scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate
#
#card.samples = 1024*500 #going down to something short enogh for regular
#card.segments = 1
#card.averages = 100
#card.triggerDelay = 25*10**-6
#card.SetParams() #pushes default to to card if the fields haven't been edited
#card.ArmAndWait() #initiates aquisition and calibrates if need be
#avData1, avData2 = card.ReadAllData() #read data for the active channels.
#
#
#
#pylab.figure(5)
#pylab.clf()
#ax = pylab.subplot(2,2,1)
#pylab.title('Ch1: 100 single segments')
#for seg in range(0, card.segments):
#    pylab.plot(ts2, data1[seg,:])
#
#pylab.subplot(2,2,3)
#pylab.plot(ts2, avData1)
#pylab.title('Ch2: average of 100 segments')    
#
#
#
#ax = pylab.subplot(2,2,2)
#pylab.title('Ch2: 100 single segments')
#for seg in range(0, card.segments):
#    pylab.plot(ts2, data2[seg,:])
#    
#pylab.subplot(2,2,4)
#pylab.plot(ts2, avData2)
#pylab.title('Ch2: average of 100 segments')    
#
#pylab.show()





#card.close() #terminate connection with the card. which we probably don't want to do.





