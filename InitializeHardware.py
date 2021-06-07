# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:39:58 2020

@author: Kollarlab
"""
import pyvisa

from kollar_instruments.HDAWG import HDAWG
from kollar_instruments.Acqiris import Acqiris
from kollar_instruments.SGS import SGS100A
from kollar_instruments.Generator import Keysight33500B
from kollar_instruments.DC205 import DC205
from kollar_instruments.VNA import VNA

from kollar_instruments.SGShelper import SGS_coupling

#Digitizer
hardwareAddress = "PXI23::0::0::INSTR" 
card = Acqiris(hardwareAddress)

#HDAWG
hdawg = HDAWG('dev8163')
  
#Generators

rm = pyvisa.ResourceManager()
rm.close()

qubitgen = SGS100A('TCPIP0::rssgs100a110425::inst0::INSTR')
cavitygen = SGS100A('TCPIP0::rssgs100a110739::inst0::INSTR')
#SGS_coupling(qubitgen,cavitygen)
cavitygen.Ref.Source = 'Ext'
cavitygen.Ref.Frequency = 10e6

cavitygen.RefOut.Source = 'Ref'
cavitygen.RefOut.Frequency = 1e9

qubitgen.Ref.Source = 'Ext'
qubitgen.Ref.Frequency = 1e9

qubitgen.IQ.Imp = 'On'
qubitgen.IQ.Ileak = exp_globals['qubitgen_config']['Ileak']
qubitgen.IQ.Qleak = exp_globals['qubitgen_config']['Qleak']

cavitygen.IQ.Imp = 'On'
cavitygen.IQ.Ileak = exp_globals['cavitygen_config']['Ileak']
cavitygen.IQ.Qleak = exp_globals['cavitygen_config']['Qleak']

triggergen  = Keysight33500B('USB0::0x0957::0x2507::MY58000681::0::INSTR')

triggergen.Ref.Source = 'Ext'
triggergen.Waveform   = 'SQU'
triggergen.Freq       = '50 kHz'
triggergen.Volt       = '2 V'
triggergen.Output     = 'ON'

SRS = DC205('ASRL3::INSTR', reset = False)
#
vna = VNA('TCPIP0::192.168.1.7::inst0::INSTR', reset = False)
#
#vars_to_save = dir()
#vars_to_save += ['vars_to_save']