# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 16:21:14 2020

@author: Kollarlab
"""

settings = {}

settings['lopower']       = 12
settings['rfpower']       = 0
settings['one_shot_time'] = 1e-6
settings['carrier_freq']  = 8e9
settings['trigger_rate']  = 500

#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 1
settings['activeChannels']   = [1,2]
settings['sampleRate']       = 2e9/8
settings['trigger_buffer']   = 0e-6
settings['digitizer_buffer'] = 10e-6
#Misc settings
settings['savepath']      = r'C:\Users\Kollarlab\Desktop'

## General reference configuration
carrier_freq       = settings['carrier_freq']
rfpower            = settings['rfpower']
lopower            = settings['lopower']             

## HDAWG settings
hdawg.AWGs[0].samplerate = '2.4GHz'
hdawg.channelgrouping = '1x4'
hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')

## Generator settings
logen.Freq   = carrier_freq
logen.Power  = lopower
logen.IQ.Mod = 'Off'

rfgen.Freq = carrier_freq
rfgen.Power = rfpower
rfgen.IQ.Mod = 'On'

logen.Output = 'On'
rfgen.Output = 'On'

## Digitizer settings
card.clockSource  = 'External'
card.triggerSlope = 'Rising'
card.triggerLevel = 0.1
card.verbose      = False

reads               = settings['reads']
digitizer_buffer    = settings['digitizer_buffer']
card.averages       = settings['averages']
card.segments       = settings['segments']
card.sampleRate     = settings['sampleRate']
card.activeChannels = settings['activeChannels']
card.triggerDelay   = settings['trigger_buffer']

## Read in the sequencer program we want to use and configure it
progFile = open("spectest.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close()
hdawg.AWGs[0].load_program(loadprog)
hdawg.AWGs[0].run_loop()
time.sleep(0.1)

#set up some houskeeping for the timing
digitizer_max_time = P2_time + pulse_width*(1+ramp_frac) + digitizer_buffer
digitizer_min_time = 0
digitizer_time     = digitizer_max_time - digitizer_min_time

card.samples = numpy.ceil(digitizer_time*card.sampleRate)
card.SetParams() #warning. this may round the number of smaples to multiple of 1024

empiricalFudgeFactor = 0.137e-6   #this is a very exact number to be off by!!!!!!
digitizerTimeOffset  = P2_time + empiricalFudgeFactor
cardTicks = scipy.arange(0, card.samples, 1.) # time in sample clocks from the digitizers point of view
cardXaxis = cardTicks /card.sampleRate - digitizerTimeOffset #time in seconds from the digitizers point of view

################################################
#Array Initialization
################################################

raw_read1 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single readS
raw_read2 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single read


################################################
#Data Collection
################################################
    for rind in range(0, reads):
        card.ArmAndWait()
        data1, data2 = card.ReadAllData() #return a matrix which segments x samples

        raw_read1 = data1
        raw_read2 = data2

        dataVec1 = data1[0,:] #first segments is representative data vector
        dataVec2 = data2[0,:]

        for sind in range(0, card.segments):

            Is1 = data1[sind, pulse1_range[0]:pulse1_range[1]]
            Qs1 = data2[sind, pulse1_range[0]:pulse1_range[1]]

            Is2 = data1[sind, pulse2_range[0]:pulse2_range[1]]
            Qs2 = data2[sind, pulse2_range[0]:pulse2_range[1]]

            I1 = numpy.mean(Is1)
            Q1 = numpy.mean(Qs1)

            I2 = numpy.mean(Is2)
            Q2 = numpy.mean(Qs2)

            phase1 = numpy.arctan2(Q1, I1)*180/numpy.pi
            phase2 = numpy.arctan2(Q2, I2)*180/numpy.pi