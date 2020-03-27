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



card.samples = 1024*515 + 2 #too long for averaging mode, but fine for regular
card.segments = 2
card.averages = 1
card.triggerDelay = 25*10**-6
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
data1, data2 = card.ReadAllData() #read data for the active channels.
ts = 10**6* scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate
#
print('Took regular data')

card.ReInitialize()

card.samples = 1024*513 + 2 #going down to something short enogh for regular
card.segments = 2
card.averages = 100
card.triggerDelay = 25*10**-6
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
avData1, avData2 = card.ReadAllData() #read data for the active channels.
avTs = 10**6* scipy.arange(0, avData1.shape[1],1.)*1/card.sampleRate

print('Took regular averaged data')

#card.samples = 1024*513 + 2 #too long for averaging mode, but fine for regular
#card.segments = 2
#card.averages = 1
#card.triggerDelay = 25*10**-6
#card.SetParams() #pushes default to to card if the fields haven't been edited
#card.ArmAndWait() #initiates aquisition and calibrates if need be
#data1, data2 = card.ReadAllData() #read data for the active channels.
#ts = 10**6* scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate
##
#print('Took regular data again')


pylab.figure(4)
pylab.clf()
ax = pylab.subplot(2,2,1)
pylab.plot(ts, data1[0,:], label = 'single')
pylab.plot(avTs, avData1[0,:], label = 'averaged')
pylab.title('Channel 1, Segment 1')
pylab.ylabel('Voltage')
pylab.xlabel('t ($\mu$s)')
ax.legend(loc = 'upper left')


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





#card.close()



