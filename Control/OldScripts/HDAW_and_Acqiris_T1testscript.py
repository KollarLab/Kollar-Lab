# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 15:14:36 2020

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



AcqirisPath = "C:\\Users\\Kollarlab\\Desktop\\Kollar-Lab\\Control\\Acqiris_development\\"
sys.path.append(AcqirisPath)


IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
sys.path.append(IVIbinPath)



from Acqiris import Acqiris
from Instruments.HDAWG import HDAWG

hardwareAddress = "PXI23::0::0::INSTR"



#fullInit = True
fullInit = False
if fullInit:
    card = Acqiris(hardwareAddress)
    
    AWG=HDAWG('dev8163')
    AWG.enable_channels([0])
    AWG.set_AWGamp([1.],[0])
    AWG.enable_markers([0])


##########
#digitizer card
###########
card.activeChannels = [1,2]
card.timeout = 5

#digitizerTime = 200*10**-6
digitizerTime = 20*10**-6
minSamples = digitizerTime*card.sampleRate
card.samples = int(minSamples)
card.segments = 5
card.averages = 10
card.triggerSource = 'External1'

#trigger of AWG marker
card.triggerSlope = 'Rising' 
card.triggerLevel = 1

##trigger of trigger generator
#card.triggerSlope = 'Falling' 
#card.triggerLevel = 0.25

#card.triggerDelay = 0.100*10**-6
card.triggerDelay = 0.*10**-6


card.SetParams() #pushes default to to card if the fields haven't been edited



##########
#AWG
#############

#progFile = open("HDAWG_sequencer_codes/T1Measurement.cpp",'r')
#rawprog  = progFile.read()
#AWG.load_program(rawprog)
#AWG.AWG_run()

#card.ArmAndWait() #initiates aquisition and calibrates if need be


card.Arm() #tell the card to look for data
AWG.AWG_run()   #have the AWG send data
card.WaitForAcquisition() #get the card to realize that it is done.


data1, data2 = card.ReadAllData() #read data for the active channels.
if card.segments == 1:
    ts = scipy.arange(0, len(data1),1.)*1/card.sampleRate
else:
    ts = scipy.arange(0, data1.shape[1],1.)*1/card.sampleRate



awgSampleRate = 2.4*10**9
awgSamples = 800
awgTicks = scipy.arange(0,800,1.)
awgTs = awgTicks*1/awgSampleRate
awgWaveForm = numpy.exp( -(awgTicks- 400)**2/(2*(100**2))    )  #factors of two etc


pylab.figure(1)
pylab.clf()

waterfall = 0.25
ax = pylab.subplot(1,2,1)
#pylab.plot(ts*10**6, data1, '.-',label = 'Channel 1')
if card.segments == 1:
    pylab.plot(ts*10**6, data2, '-', label = '0')
else:
    for seg in range(0, card.segments):
        labelstr = str(seg)
        pylab.plot(ts*10**6, data2[seg,:] + waterfall*seg, '-', label = labelstr)
pylab.title('Digitizer')
pylab.xlabel('time ($\mu$s)')
pylab.ylabel('Voltage')
ax.legend()
ax.set_ylim([-1,2])



ax = pylab.subplot(1,2,2)


pylab.plot(awgTs*10**6, awgWaveForm, label = 'guess')
pylab.title('AWG')
pylab.xlabel('time ($\mu$s)')
pylab.ylabel('Voltage')
ax.legend()


pylab.show()


