import pylab
import sys
import numpy
import scipy
import time
from datetime import datetime

from mplcursors import cursor as datacursor
from SGShelper import SGS_coupling, HDAWG_clock
import userfuncs as uf

def GetDefaultSettings():
    settings = {}
    settings['lopower']       = 12
    settings['rfpower']       = 0
    settings['one_shot_time'] = 1e-6
    settings['carrier_freq']  = 8e9
    settings['IF_freq']       = 40e6
    settings['trigger_rate']  = 500
    #Card settings
    settings['segments']         = 1
    settings['averages']         = 1
    settings['activeChannels']   = [1,2]
    settings['sampleRate']       = 2e9/8
    settings['trigger_buffer']   = 0e-6
    settings['digitizer_buffer'] = 10e-6

    return settings

#def removeIQleak(instruments, settings):

settings = GetDefaultSettings()
instruments = {}
instruments['AWG'] = hdawg
instruments['Digitizer'] = card
instruments['LO'] = logen
instruments['RFgen'] = rfgen

## Instruments used
hdawg = instruments['AWG']
logen = instruments['LO']
rfgen = instruments['RFgen']
card  = instruments['Digitizer']

## General reference configuration
carrier_freq       = settings['carrier_freq']
IF_freq            = settings['IF_freq']
rfpower            = settings['rfpower']
lopower            = settings['lopower']    
measDur            = settings['one_shot_time']         

## HDAWG settings
hdawg.AWGs[0].samplerate = '2.4GHz'
hdawg.channelgrouping = '1x4'
hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')

## Generator settings
freq_GHz    = carrier_freq/1e9
IF_freq_GHz = IF_freq/1e9
logen.set_Freq(freq_GHz)
logen.set_Amp(lopower)
logen.mod_Off()

rfgen.set_Freq(freq_GHz+IF_freq_GHz)
rfgen.set_Amp(rfpower)
rfgen.mod_On()
rfgen.imp_On()

logen.power_On() 
rfgen.power_On()

## Digitizer settings
card.clockSource  = 'External'
card.triggerSlope = 'Rising'
card.triggerLevel = 0.1
card.verbose      = False

card.averages       = settings['averages']
card.segments       = settings['segments']
card.sampleRate     = settings['sampleRate']
card.activeChannels = settings['activeChannels']
card.triggerDelay   = settings['trigger_buffer']

## Read in the sequencer program we want to use and configure it
#progFile = open("HDAWG_sequencer_codes/Blank.cpp",'r')
#rawprog  = progFile.read()
#loadprog = rawprog
#progFile.close()
#loadprog = loadprog.replace('_Time_',str(measDur))
#hdawg.AWGs[0].load_program(loadprog)
#hdawg.AWGs[0].run_loop()

card.samples = numpy.ceil(10e-6*card.sampleRate)
card.SetParams() #warning. this may round the number of smaples to multiple of 1024

cardTicks = scipy.arange(0, card.samples, 1.) # time in sample clocks from the digitizers point of view
cardXaxis = cardTicks /card.sampleRate        # time in seconds from the digitizers point of view

    ################################################
    #Array Initialization
    ################################################
def collect_data(card):
    ch1 = numpy.zeros(card.samples) #store all the raw data from a single readS
    ch2 = numpy.zeros(card.samples) #store all the raw data from a single read
    
    card.ArmAndWait()
    ch1, ch2 = card.ReadAllData()
    
    ch1 = ch1-numpy.mean(ch1)
    ch2 = ch2-numpy.mean(ch2)
    
    amp1 = max(ch1[0])-min(ch1[0])
    amp2 = max(ch2[0])-min(ch2[0])
    return max(amp1,amp2)

def run_me(card, rfgen, ipoints, qpoints, irange = 1, qrange = 1):
    amps   = numpy.zeros((ipoints, qpoints))
    ileaks = numpy.linspace(-irange,irange,ipoints)
    qleaks = numpy.linspace(-qrange,qrange,qpoints)
    for i in range(ipoints):
        rfgen.leak_I(ileaks[i])
        for j in range(qpoints):
            rfgen.leak_Q(qleaks[j])
            amps[i,j] = collect_data(card)
    pylab.imshow(amps, extent=[-qrange,qrange,-irange,irange])
    pylab.colorbar()
    
    bests = numpy.where(amps == amps.min())
    ibest = ileaks[bests[0][0]]
    qbest = qleaks[bests[1][0]]
    return ibest, qbest
