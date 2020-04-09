# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 10:53:31 2020

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

card = Acqiris(hardwareAddress)




params = {}
params['timeout'] = 2
params['verbose'] = False

params['samples'] = 1024*600
params['sampleRate'] = 10**9

params['averages'] = 1
params['averageMode'] = 0

params['segments'] = 1

params['activeChannels'] = [1,2]
params['channelRange'] = 0.5
params['channelOffset'] = 0



params['triggerSource'] = 'External1'
params['triggerLevel'] = 0.1
params['triggerSlope'] = 'Falling'
params['triggerDelay'] = 25*10**-6

params['clockSource'] = 'Internal'
params['clockFrequency'] = 10**7



card.SetParams(params)



#card.activeChannels = [1,2]
#card.timeout = 2
#segs = 1
#card.verbose = False
#card.samples = 1024*600 #too long for averaging mode, but fine for regular
#card.segments = segs
#card.averages = 1
#card.triggerDelay = 25*10**-6
#card.SetParams() #pushes default to to card if the fields haven't been edited
#card.ArmAndWait() #initiates aquisition and calibrates if need be
#if len(card.activeChannels) == 1:
#    data1 = card.ReadAllData() #read data for the active channels.
#else:
#    data1, data2 = card.ReadAllData() #read data for the active channels.
#ts = 10**6* scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate
#    
    
    
    
    
    
    
    
    
    
    
    
    
    