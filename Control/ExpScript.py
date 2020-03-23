from Instruments.HDAWG import HDAWG
from Instruments.SGS import RFgen
import numpy as np
import time

#RFsource=RFgen('TCPIP0::rssgs100a110739::inst0::INSTR') #SGS visa resource name

## Connect to HDAWG and initialize it 
AWG = HDAWG('dev8163') #HDAWG device name
AWG.enable_channels([0,1]) #turn on output for channels 0 and 1
AWG.set_AWGamp([1.,1.],[0,1]) #set output amplitude
AWG.enable_markers([0,1]) #enable markers for our channels 

## Read in the sequencer program we want to use
progFile = open("HDAWG_sequencer_codes/T1Measurement",'r')
rawprog  = progFile.read()
loadprog = rawprog

## Configure program to our parameters
NumSamples = 800 #sets number of samples to use for the waveform (try to make it a multiple of 16 otherwise AWG pads with zeros)
waitInc    = 1000 #sets increment of wait time between two pulses (in units of clock cycles, approx 3.3 ns)
clearTime  = 1000 #time between two measurements (clock cycles)
markerPos  = 300 #sets position of trigger relative to signal (0 is coicident)

loadprog = loadprog.replace('_NumSamples_', str(NumSamples))
loadprog = loadprog.replace('_waitInc_', str(waitInc))
loadprog = loadprog.replace('_clearTime_', str(clearTime))
loadprog = loadprog.replace('_markerPos_', str(markerPos))
AWG.load_program(loadprog)

## Create waveforms we want to use for channels 1 and 2
time    = np.linspace(-np.pi,np.pi, NumSamples)
Isignal = np.sin(time)
Qsignal = np.cos(time)

#AWG.load_waveform(0,Isignal,Qsignal)

## Run AWG program
AWG.AWG_run()