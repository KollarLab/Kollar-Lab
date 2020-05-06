# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:39:58 2020

@author: Kollarlab
"""

import pylab
import sys
import numpy
import scipy
import time

from Instruments.HDAWG import HDAWG
from Acqiris_development.Acqiris import Acqiris
from Instruments.SGS import RFgen

#Card Settings
hardwareAddress = "PXI23::0::0::INSTR"
IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
if not IVIbinPath in sys.path:
    sys.path.append(IVIbinPath)

card = Acqiris(hardwareAddress)
card.triggerSlope = 'Rising'
card.triggerLevel = 0.1

###HDAWG Configuration
hdawg = HDAWG('dev8163')
hdawg.AWGs[0].samplerate = '2.4GHz'
hdawg.channelgrouping = '1x4'
hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
hdawg.Channels[1].configureChannel(marker_out='Trigger', hold='True')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
hdawg.OSCs[1].freq = 10e6
hdawg.Channels[2].analog_outs = [0.5,0]
hdawg.Channels[3].analog_outs = [0,1.0]
hdawg.Channels[2].configureChannel(amp=1.0)
hdawg.Channels[3].configureChannel(amp=2.0)

#Generators
logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')