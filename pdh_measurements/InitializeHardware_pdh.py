# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:09:24 2024

@author: Kollarlab
"""

import pyvisa

from kollar_instruments.HDAWG import HDAWG
from kollar_instruments.Acqiris import Acqiris
from kollar_instruments.SGS import SGS100A
from kollar_instruments.Generator import Keysight33500B
from kollar_instruments.HS9000B import HS9000B

#Digitizer
hardwareAddress = "PXI23::0::0::INSTR" 
card = Acqiris(hardwareAddress)

#HDAWG
hdawg = HDAWG('dev8687')
hdawg.channelgrouping='1x4'
hdawg.Channels[0].configureChannel(amp=0.5, marker_out='Marker', hold='False')
hdawg.Channels[1].configureChannel(amp=0.5, marker_out='Marker', hold='False')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising', channel='Trigger in 1')


#Generators
# Using two holzworths 
holz = HS9000B('ASRL6::INSTR') 
holz.ref.source = 'Ext'
holz.ref.frequency = 10e6

holz2 = HS9000B('ASRL4::INSTR') 
holz2.ref.source = 'Ext'
holz2.ref.frequency = 10e6

rm = pyvisa.ResourceManager()
rm.close()

# cavity gen
gen = SGS100A('TCPIP0::rssgs100a110738::inst0::INSTR')
gen.enable_IQ()
gen.power = -60
gen.Ref.Source = 'Ext'
gen.Ref.Frequency = 10e6
gen.RefOut.Frequency = 1e9
gen.output = 'Off'

# TriggerGen
triggergen  = Keysight33500B('USB0::0x0957::0x2507::MY58000660::0::INSTR')

triggergen.Ref.Source = 'Ext'
triggergen.Freq       = 40e3
triggergen.Waveform   = 'SQU'
triggergen.Volt       = '2 V'
triggergen.Output     = 'ON'
triggergen.Duty_cycle = '0.1'

# # Acqiris Card
card.channelRange = 2.5
card.activeChannels = [1,2]
card.timeout = 30
card.triggerSlope = 'Falling'
card.averages = 1000
card.segments = 1
card.triggerLevel = 0.1
card.verbose = False
card.sampleRate = 1e9
card.samples = 20e3
card.SetParams()


