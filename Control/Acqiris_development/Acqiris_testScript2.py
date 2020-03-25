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

card.samples = 1024*513 #too long for averaging mode, but fine for regular
card.segments = 2
card.averages = 1
card.triggerDelay = 25*10**-6
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
data1, data2 = card.ReadAllData() #read data for the active channels.
#
print('Took regular data')
#
#card.ReInitialize()

card.samples = 1024*500 #going down to something short enogh for regular
card.segments = 2
card.averages = 100
card.triggerDelay = 25*10**-6
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
avData1, avData2 = card.ReadAllData() #read data for the active channels.

print('Took regular averaged data')


#card.samples = 1024*513 #too long for averaging mode, but fine for regular
#card.segments = 2
#card.averages = 1
#card.SetParams() #pushes default to to card if the fields haven't been edited
#card.ArmAndWait() #initiates aquisition and calibrates if need be
#data1, data2 = card.ReadAllData() #read data for the active channels.
##
#print('Took regular data again')

pylab.figure(4)
pylab.clf()
ax = pylab.subplot(2,2,1)
pylab.plot(data1[0,:], label = 'single')
pylab.plot(avData1[0,:], label = 'averaged')
pylab.title('Channel 1, Segment 1')
ax.legend(loc = 'upper left')


ax = pylab.subplot(2,2,2)
pylab.plot(data2[0,:], label = 'single')
pylab.plot(avData2[0,:], label = 'averaged')
pylab.title('Channel 2, Segment 1')
ax.legend(loc = 'upper left')


ax = pylab.subplot(2,2,3)
pylab.plot(data1[1,:], label = 'single')
pylab.plot(avData1[1,:], label = 'averaged')
pylab.title('Channel 1, Segment 2')
ax.legend(loc = 'upper left')


ax = pylab.subplot(2,2,4)
pylab.plot(data2[1,:], label = 'single')
pylab.plot(avData2[1,:], label = 'averaged')
pylab.title('Channel 2, Segment 2')
ax.legend(loc = 'upper left')

pylab.suptitle('Tah-Dah!')
pylab.show()

#card.close() #terminate connection with the card. which we probably don't want to do.





