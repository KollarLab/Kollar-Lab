# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:20:43 2021

@author: Kollarlab
"""
import matplotlib.pyplot as plt
import time
import numpy as np
import userfuncs
import os

settings = {}
#Misc
settings['scanname'] = 'testing'
settings['project_dir']  = r'Z:\Data\HouckDualHangerFluxonium'
settings['meas_type'] = 'drainage_T1'
settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 0

#Cavity parameters
settings['CAV_power'] = -45
settings['CAV_freq']  = 7.57630e9
settings['Q_power']   = -45
#settings['Q_freq']    = 8.806e9
#settings['Q_freq']    = 8.8246e9
settings['Q_freq']    = 7.57630e9

#Card settings
settings['segments']         = 5
settings['reads']            = 1
settings['averages']         = 60e3
settings['activeChannels']   = [1,2]
settings['channelRange']     = 0.5
settings['sampleRate']       = 2e9/8
settings['empirical_delay']  = np.floor(375e-9/32e-9) * 32e-9  #1e-6
settings['timeout']          = 30

##Pulse settings
settings['meas_pos']    = np.floor(2e-6/32e-9) * 32e-9 #15e-6
settings['meas_window'] = 3e-6
settings['overlap'] = 1e-6

saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
stamp = userfuncs.timestamp()
filename = settings['scanname'] + '_' + stamp
    
qubitgen.Freq   = settings['Q_freq']
qubitgen.Power  = settings['Q_power'] + settings['Qbit_Attenuation']
qubitgen.IQ.Mod = 'On'
qubitgen.Output = 'On'

cavitygen.Freq   = settings['CAV_freq']
cavitygen.Power  = settings['CAV_power'] + settings['CAV_Attenuation']
cavitygen.IQ.Mod = 'On'
cavitygen.Output = 'On'

vna.Ref.Source = 'EXT'
vna.Power = 12
vna.Freq = settings['CAV_freq']
vna.Output = 'On'

##Card setings

meas_samples = settings['sampleRate']*(settings['meas_window']+settings['meas_pos'])

card.averages       = settings['averages']
card.segments       = settings['segments']
card.sampleRate     = settings['sampleRate']
card.activeChannels = settings['activeChannels']
#card.triggerDelay   = settings['meas_pos'] + settings['empirical_delay']
card.triggerDelay = 0
card.timeout        = settings['timeout']
card.samples        = int(meas_samples*1.5)
card.channelRange   = settings['channelRange']
card.SetParams()

xaxis = (np.array(range(card.samples))/card.sampleRate)
xaxis_us = xaxis*1e6

##HDAWG settings
hdawg.AWGs[0].samplerate = '2.4GHz'
hdawg.channelgrouping = '1x4'
hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')

#progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\Control\HDAWG_sequencer_codes\pulsedtrans.cpp",'r')
progFile = open(r"C:\Users\Kollarlab\Desktop\Kollar-Lab\Control\HDAWG_sequencer_codes\overlapped.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close()

loadprog = loadprog.replace('_max_time_', str(settings['meas_pos']))
loadprog = loadprog.replace('_meas_window_', str(settings['meas_window']))
loadprog = loadprog.replace('_overlap_', str(settings['overlap']))

hdawg.AWGs[0].load_program(loadprog)
hdawg.AWGs[0].run_loop()
time.sleep(0.1)
    
    
card.ArmAndWait()
I, Q = card.ReadAllData()

Ip = np.mean(I, 0)
Qp = np.mean(Q, 0)

DC_I = np.mean(Ip[-50:])
DC_Q = np.mean(Qp[-50:])
Idat = Ip-DC_I
Qdat = Qp-DC_Q 
amp = np.sqrt(Idat**2+Qdat**2)
phase = np.arctan2(Qdat, Idat)*180/np.pi
            
plt.figure(1)
plt.plot(xaxis_us, amp)
plt.xlabel('Time (us)')
plt.ylabel('Voltage (V)')
plt.title('Drainage T1, overlap: {}'.format(settings['overlap']))
#plt.title('Normal transmission')

plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)

userfuncs.SaveFull(saveDir, filename, ['amp', 'phase','xaxis','xaxis_us'], locals(), expsettings=settings)