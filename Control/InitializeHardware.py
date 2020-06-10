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
from mpldatacursor import datacursor

from Instruments.HDAWG import HDAWG
from Instruments.Acqiris import Acqiris
from Instruments.SGS import RFgen
from Instruments.Generator import Generator
from SGShelper import SGS_coupling

#Digitizer
hardwareAddress = "PXI23::0::0::INSTR"
card = Acqiris(hardwareAddress)

#HDAWG
hdawg = HDAWG('dev8163')

#Generators
logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')
SGS_coupling(logen,rfgen)

triggergen  = Generator('')
syncgen     = Generator('')

triggergen.reference = 'Ext'
triggergen.waveform  = 'SQU'
triggergen.frequency = '500 Hz'
triggergen.volts     = '0.5 V'

syncgen.reference = 'Ext'
syncgen.waveform  = 'SIN'
syncgen.frequency = '10 MHz'
syncgen.volts     = '3 dBm'