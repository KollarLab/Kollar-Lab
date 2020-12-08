# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:52:18 2020

@author: Kollarlab
"""
import time
from mplcursors import cursor as datatip
import numpy
import os
from VNAplottingTools import base_power_plot_imshow
import matplotlib.pyplot as plt
import userfuncs

settings = {}

settings['scanname'] = 'T1_meas'

Q_Freq = 4.20431e9
Q_Power = -11
CAV_Freq = 8.126e9
CAV_Power = -18

Meas_pos = 80e-6
#Card settings
settings['segments']         = 1
settings['reads']            = 1
settings['averages']         = 25e3
settings['activeChannels']   = [1,2]
settings['sampleRate']       = 2e9/8
settings['trigger_buffer']   = Meas_pos
settings['meas_window']      = 20e-6

settings['Tau_min'] = 200e-9
settings['Tau_max'] = 30e-6
settings['Tau_points'] = 30

settings['Measurement_pos'] = Meas_pos

#load up the settings
data_type_name = 'PulseSpec'
stamp = userfuncs.timestamp()
filename = data_type_name + '_' + settings['scanname'] + '_' + stamp

## Generator settings
SMB.Ref.Source = 'EXT'
SMB.Power = 12
SMB.Freq = CAV_freq
SMB.Output = 'On'

cavitygen.Freq   = CAV_Freq
cavitygen.Power  = CAV_Power
cavitygen.IQ.Mod = 'On'

qubitgen.Freq   = Q_Freq
qubitgen.Power  = Q_Power
qubitgen.IQ.Mod = 'On'

cavitygen.Output = 'On'
qubitgen.Output = 'On'

## Card Settings
meas_samples = settings['sampleRate']*settings['meas_window']

card.averages       = settings['averages']
card.segments       = settings['segments']
card.sampleRate     = settings['sampleRate']
card.activeChannels = settings['activeChannels']
card.triggerDelay   = settings['trigger_buffer']
card.timeout        = 30
card.samples        = int(meas_samples*2.5)
card.channelRange   = 0.5
card.SetParams()

data_window = int(meas_samples) #for 100us measurement pulses 

xaxis = (numpy.array(range(card.samples))/card.sampleRate)
xaxis_us = xaxis*1e6 + settings['trigger_buffer']

## HDAWG
hdawg.AWGs[0].samplerate = '2.4GHz'
hdawg.channelgrouping = '1x4'
hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='True')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')

progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\Control\HDAWG_sequencer_codes\T1.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close()

loadprog = loadprog.replace('_max_time_', str(settings['Measurement_pos']))
loadprog = loadprog.replace('_meas_window_', str(settings['meas_window']))
taus = numpy.linspace(settings['Tau_min'],settings['Tau_max'] , settings['Tau_points'] )

plt.figure()

amp_int = numpy.zeros(len(taus))
amps    = numpy.zeros((len(taus),card.samples))

tstart = time.time()

for tind in range(len(taus)):
        
    tau = taus[tind]
    print('Tau: {}'.format(tau))
    finalprog = loadprog
    finalprog = finalprog.replace('_tau_',str(tau))
    hdawg.AWGs[0].load_program(finalprog)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)
    
    card.ArmAndWait()
    I,Q = card.ReadAllData()
    
    DC_I = numpy.mean(I[0][-data_window:])
    DC_Q = numpy.mean(Q[0][-data_window:])
    
    Idat = I[0]-DC_I
    Qdat = Q[0]-DC_Q
    
    amp = numpy.sqrt(Idat**2+Qdat**2)
    
    amps[tind] = amp
    amp_int[tind] = numpy.mean(amp[0:int(data_window)])
    if tind == 0:
        tstop = time.time()
        singlePointTime = tstop-tstart
        
        estimatedTime = singlePointTime*len(taus)
        print('    ')
        print('estimated time for this scan : ' + str(numpy.round(estimatedTime/60, 1)) + ' minutes')
        print('estimated time for this scan : ' + str(numpy.round(estimatedTime/60/60, 2)) + ' hours')
        print('    ')

t2 = time.time()


print('Elapsed time: {}'.format(t2-tstart))
plt.plot(taus, amp_int)

fig = plt.figure(figsize=(13,8))
plt.clf()
    
ax = plt.subplot(1,1,1)
base_power_plot_imshow(fig, ax, xaxis_us, taus*1e6, amps, ['Time (us)', 'Tau (us)', 'Amp'], attenuation=0)