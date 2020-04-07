# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 15:50:37 2020

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

card = Acqiris(hardwareAddress)

t0 = time.time()

dt = 0
rind = 0
T = 0
maxT = 15*60
runStep = 100
while T < maxT:
    
    
    card.activeChannels = [1,2]
    card.timeout = 2
    segs = 8
    card.verbose = False
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


    if numpy.mod(rind, runStep) == 99:
        t1 = time.time()
        dt = (t1-t0) -T #elapsed time for last set of runs
        T = t1 - t0 #new absolute time
        print('Round : ' + str(rind))
        print('time per ' +str(runStep) + ' runs = ' + str(numpy.round(dt, 2)) + ' seconds')
    rind = rind+1

print()
print('samples : ' + str(card.samples))
print('segments : ' + str(card.segments))
print('averages : ' + str(card.averages))
print('reads : ' + str(rind))
print('total time : ' + str(numpy.round(T/60,2)) + ' minutes')

card.close()    







#rinds = scipy.arange(0,10,1.)
#for rind in rinds:
#    print()
#    print('Round : ' + str(rind))
#    
#    card.activeChannels = [1,2]
#    card.timeout = 2
#    segs = 2
#    card.verbose = True
#    card.samples = 1024*600 #too long for averaging mode, but fine for regular
#    card.segments = segs
#    card.averages = 1
#    card.triggerDelay = 25*10**-6
#    card.SetParams() #pushes default to to card if the fields haven't been edited
#    card.ArmAndWait() #initiates aquisition and calibrates if need be
#    if len(card.activeChannels) == 1:
#        data1 = card.ReadAllData() #read data for the active channels.
#    else:
#        data1, data2 = card.ReadAllData() #read data for the active channels.
#    ts = 10**6* scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate
#    #
#    print('Took regular data')
#    
#
#card.close()    
    
    
    