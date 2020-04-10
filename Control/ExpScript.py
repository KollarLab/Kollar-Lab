from Instruments.HDAWG import HDAWG
from Instruments.SGS import RFgen
import numpy as np
import time

#RFsource=RFgen('TCPIP0::rssgs100a110739::inst0::INSTR') #SGS visa resource name

## Connect to HDAWG and initialize it 
hdawg = HDAWG('dev8163') #HDAWG device name
hdawg.AWGs[0].samplerate = '2.4GHz'
hdawg.channelgrouping = '1x4'
hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')

## Read in the sequencer program we want to use
progFile = open("HDAWG_sequencer_codes/T1Measurement.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close

## Configure program to our parameters
piAmp           = 1.0 #pi pulse amplitude, relative to output amplitude
piTime          = 0.5e-6 #pi pulse length (s)
piWidth         = piTime/8 #pi pulse width (gaussian) (s)
qubit_lifetime  = 10e-6 #qubit lifetime (s)
wait_increment  = 1e-6 #timestep between different iterations (s)
marker_position = 0 #timing of marker signal (s)
averages        = 10 #number of averages per setting
number_configs  = 5 #number of different configurations

loadprog = loadprog.replace('_piAmp_', str(piAmp))
loadprog = loadprog.replace('_piTime_', str(piTime))
loadprog = loadprog.replace('_piWidth_', str(piWidth))
loadprog = loadprog.replace('_qlifetime_', str(qubit_lifetime))
loadprog = loadprog.replace('_waitInc_', str(wait_increment))
loadprog = loadprog.replace('_markerPos_', str(marker_position))
loadprog = loadprog.replace('_averages_', str(averages))
loadprog = loadprog.replace('_numConfig_', str(number_configs))
hdawg.AWGs[0].load_program(loadprog)

hdawg.AWGs[0].run()