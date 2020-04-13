from Instruments.HDAWG import HDAWG
from Instruments.SGS import RFgen
import numpy
import time

#RFsource=RFgen('TCPIP0::rssgs100a110739::inst0::INSTR') #SGS visa resource name

## Connect to HDAWG and initialize it 
hdawg = HDAWG('dev8163') #HDAWG device name
hdawg.AWGs[0].samplerate = '2.4GHz'
hdawg.channelgrouping = '1x4'
hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker')
hdawg.Channels[1].configureChannel(marker_out='Trigger')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')

## Read in the sequencer program we want to use
progFile = open("HDAWG_sequencer_codes/T1Measurement.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close

## Configure program to our parameters
piAmp     = 1.0 #pi pulse amplitude, relative to output amplitude
piTime    = 0.5e-6 #pi pulse length (s)
piWidth   = piTime/8 #pi pulse width (gaussian) (s)
meas_wait = 0.5e-6 #wait time between last pulse and measurement (s)
meas_time = 10e-6 #measurement window size (s)
max_time  = 200e-6 #maximum separation between pulses (s)
min_time  = piTime*2

loadprog = loadprog.replace('_piAmp_', str(piAmp))
loadprog = loadprog.replace('_piTime_', str(piTime))
loadprog = loadprog.replace('_piWidth_', str(piWidth))
loadprog = loadprog.replace('_meas_wait_', str(meas_wait))
loadprog = loadprog.replace('_meas_time_', str(meas_time))
loadprog = loadprog.replace('_max_time_', str(max_time))

# Create a logarithmically spaced array for exponential measurements 
measure_points = 50
taus = numpy.logspace(min_time,max_time,num=measure_points)

for time in taus:
    finalprog = loadprog
    finalprog = finalprog.replace('_tau_',str(time))
    hdawg.AWGs[0].load_program(finalprog)
    hdawg.AWGs[0].run()
    #if Digitizer Happy:
    #    hdawg.AWGs[0].stop()