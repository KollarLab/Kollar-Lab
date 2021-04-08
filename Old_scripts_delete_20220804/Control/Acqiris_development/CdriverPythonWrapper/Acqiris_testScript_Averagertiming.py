# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 18:27:58 2020

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

def get_color(ind):
    colorlist = ['firebrick',  'mediumblue', 'deepskyblue', 'darkgoldenrod', 'forestgreen', 'indigo', 'dodgerblue']
    
    nind = numpy.mod(ind, len(colorlist))
    return colorlist[nind]

hardwareAddress = "PXI23::0::0::INSTR"

IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
sys.path.append(IVIbinPath)


card.activeChannels = [1,2]
card.timeout = 120


avs = 1
segs = 1



card.samples = 1024*125
card.segments = segs
card.averages = 1
#
#card.SetParams() #here is the danger of not using properties to set everything.
##without this here, card.samples isn't 

dataMat1 = numpy.zeros((len(delays), card.samples)) 
dataMat1_av = numpy.zeros((len(delays), card.samples)) 
dataMat2 = numpy.zeros((len(delays), card.samples)) 
dataMat2_av = numpy.zeros((len(delays), card.samples)) 
tMat = numpy.zeros((len(delays), card.samples)) 



#pretake data to set everything up for test
card.averages = 1
card.triggerDelay = 0
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be
if len(card.activeChannels) == 1:
    data1 = card.ReadAllData() #read data for the active channels.
else:
    data1, data2 = card.ReadAllData() #read data for the active channels.



t0 = time.time()

for ind in range(0, avs):
    card.averages = 1
    card.triggerDelay = 0
    card.SetParams() #pushes default to to card if the fields haven't been edited
    card.ArmAndWait() #initiates aquisition and calibrates if need be
    if len(card.activeChannels) == 1:
        data1 = card.ReadAllData() #read data for the active channels.
    else:
        data1, data2 = card.ReadAllData() #read data for the active channels.
    
    ts =  ( delay + scipy.arange(0, len(data1),1.)*1/card.sampleRate) 

t1 = time.time()


card.averages = 50
card.triggerDelay = 0
card.SetParams() #pushes default to to card if the fields haven't been edited
card.ArmAndWait() #initiates aquisition and calibrates if need be


if len(card.activeChannels) == 1:
    avData1 = card.ReadAllData() #read data for the active channels.
else:
    avData1, avData2 = card.ReadAllData() #read data for the active channels.
t2 = time.time()

d1 = numpy.round(t1-t0, 3)
d2 = numpy.round(t2-t1, 3)


print('segments = ' + str(segs))
print('averages = ' + str(avs))
print('time for ' + str(avs) + ' single (possibly multiseg) runs = ' + str(d1) )
print('time for ' + str(avs) + ' averages on card = ' + str(d2) )

















