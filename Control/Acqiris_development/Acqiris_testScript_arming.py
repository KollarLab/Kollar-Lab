# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 10:49:30 2020

@author: Kollarlab
"""


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



card.verbose = True
card.timeout = 2
card.samples = 1024*600 #too long for averaging mode, but fine for regular
card.segments = 2
card.averages = 1
card.triggerDelay = 25*10**-6
card.SetParams() #pushes default to to card if the fields haven't been edited




print('armed = ' + str(card.armed))
print('finished = ' + str(card.acquisitionFinished))
card.Arm() #initiates aquisition and calibrates if need be
print('armed = ' + str(card.armed))
print('finished = ' + str(card.acquisitionFinished))
if len(card.activeChannels) == 1:
    data1 = card.ReadAllData() #read data for the active channels.
else:
    data1, data2 = card.ReadAllData() #read data for the active channels.
ts = 10**6* scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate
print('armed = ' + str(card.armed))
print('finished = ' + str(card.acquisitionFinished))





#print('armed = ' + str(card.armed))
#print('finished = ' + str(card.acquisitionFinished))
#card.ArmAndWait() #initiates aquisition and calibrates if need be
#print('armed = ' + str(card.armed))
#print('finished = ' + str(card.acquisitionFinished))
#if len(card.activeChannels) == 1:
#    data1 = card.ReadAllData() #read data for the active channels.
#else:
#    data1, data2 = card.ReadAllData() #read data for the active channels.
#ts = 10**6* scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate
#print('armed = ' + str(card.armed))
#print('finished = ' + str(card.acquisitionFinished))












