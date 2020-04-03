# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 17:57:12 2020

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


#card = Acqiris(hardwareAddress)


##########
#set settings
#########

card.activeChannels = [1,2]
card.timeout = 3


avs = 50
delays = numpy.linspace(0.,125,6)*10**-6



card.samples = 1024*125
card.segments = 1
card.averages = 1
#
#card.SetParams() #here is the danger of not using properties to set everything.
##without this here, card.samples isn't 

dataMat1 = numpy.zeros((len(delays), card.samples)) 
dataMat1_av = numpy.zeros((len(delays), card.samples)) 
dataMat2 = numpy.zeros((len(delays), card.samples)) 
dataMat2_av = numpy.zeros((len(delays), card.samples)) 
tMat = numpy.zeros((len(delays), card.samples)) 


for dind in range(0, len(delays)):
    delay = delays[dind]
    
    
    card.averages = 1
    card.triggerDelay = delay
    card.SetParams() #pushes default to to card if the fields haven't been edited
    card.ArmAndWait() #initiates aquisition and calibrates if need be
    if len(card.activeChannels) == 1:
        data1 = card.ReadAllData() #read data for the active channels.
    else:
        data1, data2 = card.ReadAllData() #read data for the active channels.
    
    ts =  ( delay + scipy.arange(0, len(data1),1.)*1/card.sampleRate) 
    
    #
    print('Took regular data')
    
    
    card.averages = avs
    card.SetParams() #pushes default to to card if the fields haven't been edited
    card.ArmAndWait() #initiates aquisition and calibrates if need be
    if len(card.activeChannels) == 1:
        avData1 = card.ReadAllData() #read data for the active channels.
    else:
        avData1, avData2 = card.ReadAllData() #read data for the active channels.
    
    print('Took regular averaged data')

    dataMat1[dind, :] = data1
    dataMat1_av[dind, :] = avData1
    dataMat2[dind, :] = data2
    dataMat2_av[dind, :] = avData2
    tMat[dind,:] = ts
    
    

waterfall = 0.1
pylab.figure(7)
pylab.clf()
ax = pylab.subplot(1,2,1)
for dind in range(0, len(delays)):
    delay = delays[dind]
    labelstr = str(numpy.round(delay*10**6, 1)) + ' us'
    pylab.plot(tMat[dind,:]*10**6, dataMat1[dind,:] + waterfall*dind, zorder = 1, label = labelstr, color=get_color(dind))
    
    labelstr = str(numpy.round(delay*10**6, 1)) + ' us : Av'
    pylab.plot(tMat[dind,:]*10**6, dataMat1_av[dind,:] + waterfall*dind+ 0.01, zorder = 2, label = labelstr, color=get_color(dind+2))
pylab.xlabel('t ($\mu$s)')
pylab.ylabel('Voltage (waterfall)')
pylab.title('Channel 1')
ax.legend(loc = 'lower right')
ax.set_ylim([-0.04, 0.08 + waterfall*dind])
ax.set_xlim([0,numpy.max(ts*10**6) + 75])


if len(card.activeChannels) == 2:
    ax = pylab.subplot(1,2,2)
    for dind in range(0, len(delays)):
        delay = delays[dind]
        labelstr = str(numpy.round(delay*10**6, 1)) + ' us'
        pylab.plot(tMat[dind,:]*10**6, dataMat2[dind,:] + waterfall*dind, zorder = 1, label = labelstr, color=get_color(dind))
        
        labelstr = str(numpy.round(delay*10**6, 1)) + ' us : Av'
        pylab.plot(tMat[dind,:]*10**6, dataMat2_av[dind,:] + waterfall*dind + 0.01, zorder = 2, label = labelstr, color=get_color(dind+2))

    pylab.xlabel('t ($\mu$s)')
    pylab.ylabel('Voltage (waterfall)')
    pylab.title('Channel 2')
    ax.legend(loc = 'lower right')
    ax.set_ylim([-0.04, 0.08 + waterfall*dind])
    ax.set_xlim([0,numpy.max(ts*10**6) + 75])


pylab.suptitle('Trigger Delay Test')






