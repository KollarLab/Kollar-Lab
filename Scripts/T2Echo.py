# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 15:06:42 2020

@author: Kollarlab
"""

import time
from mplcursors import cursor as datatip
import numpy
import os
from VNAplottingTools import base_power_plot, base_raw_time_plot_spec
import matplotlib.pyplot as plt
import userfuncs

saveDir = r'Z:\Data\HouckQuadTransmon\Tmeas\20201208'

settings = {}

settings['scanname'] = 'T2_echo_meas_2Pulses'

extra_shift = 500e3
detuning = 1.23e6
final_shift = 492.04e3
Q_Freq = 4.20431e9 + detuning + extra_shift - final_shift

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
settings['wait_time']        = 1e-6

settings['Tau_min'] = 200e-9
settings['Tau_max'] = 40e-6
settings['Tau_points'] = 52

settings['Measurement_pos'] = Meas_pos

#load up the settings
data_type_name = 'Fine_Data'
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

progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\Control\HDAWG_sequencer_codes\T2.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close()

loadprog = loadprog.replace('_max_time_', str(settings['Measurement_pos']))
loadprog = loadprog.replace('_meas_window_', str(settings['meas_window']))
loadprog = loadprog.replace('_wait_time_', str(settings['wait_time']))

taus = numpy.linspace(settings['Tau_min'],settings['Tau_max'] , settings['Tau_points'] )

plt.figure(1)
plt.clf()
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

tau0 = 5e-6
offset0 = numpy.mean(amp_int[-10])
amp0 = numpy.mean(amp_int-offset0)

fit_guess = [tau0, amp0, offset0]
fit_T_exp(taus, amp_int, fit_guess)

plt.suptitle(filename)
plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)

fig = plt.figure(figsize=(13,8))
plt.clf()
    
ax = plt.subplot(1,1,1)
base_power_plot_imshow(fig, ax, xaxis_us, taus*1e6, amps, ['Time (us)', 'Tau (us)', 'Amp'], attenuation=0)

plt.suptitle(filename)
plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
userfuncs.SaveFull(saveDir, filename, ['taus','xaxis_us', 'amps', 'amp_int'], locals(), expsettings=settings)