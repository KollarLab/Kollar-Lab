# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:14:33 2020

@author: Kollarlab
"""

import pylab
import sys
import numpy
import scipy
import time
import threading

from mpldatacursor import datacursor

from Instruments.HDAWG import HDAWG
from Acqiris_development.Acqiris import Acqiris
from Instruments.SGS import RFgen

import userfuncs as uf


#############
#initialize stuff
############

#Card Settings
hardwareAddress = "PXI23::0::0::INSTR"
IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
if not IVIbinPath in sys.path:
    sys.path.append(IVIbinPath)

#card = Acqiris(hardwareAddress)
card.triggerSlope = 'Rising'
card.triggerLevel = 0.1

####HDAWG Configuration
#hdawg = HDAWG('dev8163')
#hdawg.AWGs[0].samplerate = '2.4GHz'
#hdawg.channelgrouping = '1x4'
#hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
#hdawg.Channels[1].configureChannel(marker_out='Trigger', hold='True')
#hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
#hdawg.OSCs[1].freq = 10e6
#hdawg.Channels[2].analog_outs = [0.5,0]
#hdawg.Channels[3].analog_outs = [0,1.0]
#hdawg.Channels[2].configureChannel(amp=1.0)
#hdawg.Channels[3].configureChannel(amp=2.0)

#Generators
#logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
#rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')
freq_GHz = 7;
logen.set_Freq(freq_GHz)
logen.set_Amp(12)
logen.mod_Off()
logen.set_External_Reference()
logen.power_On() 

rfgen.set_Freq(freq_GHz)
rfgen.set_Amp(2)
rfgen.mod_On()
rfgen.set_External_Reference()
rfgen.power_On()

time.sleep(0.05)


##############
## Experimental parameters, that we change
###########
timeSpacings = numpy.ones(7200)
#timeSpacings = numpy.ones(60)

pulseamp = 0.5
pulselength = 200e-9
overflowBuffer = 2e-6




##############
## Experimental parameters, that wedon't change
## or are determined by the ones we've already chozen
###########
#T2 measurement protocol parameters
measure_points = 1
tau_max = 2e-6
tau_min = 1e-6 
trig_buffer = 0e-6
taus = numpy.linspace(tau_min, tau_max, measure_points)
maxtime = tau_max + overflowBuffer


#################
#Push the final harware settings
##################
#Data acquisition parameters
card.averages = 1 #on-board averages
card.segments = 1
card.triggerDelay = trig_buffer
card.activeChannels = [1,2]
reads = 1  #reads of the card
card.verbose = False
card.sampleRate = 2e9
card.clockSource = 'External'

## Read in the sequencer program we want to use and configure it
progFile = open("HDAWG_sequencer_codes/HomodynePulse.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close
loadprog = loadprog.replace('_Amp_', str(pulseamp))
loadprog = loadprog.replace('_Time_', str(pulselength))
loadprog = loadprog.replace('_max_time_', str(maxtime))
tau = taus[0]
finalprog = loadprog
finalprog = finalprog.replace('_tau_',str(tau))
hdawg.AWGs[0].load_program(finalprog)
hdawg.AWGs[0].run_loop()
time.sleep(0.1)


#set up some houskeeping for the timing
digitizer_max_time = maxtime + 1e-6
digitizer_min_time = 0
digitizer_time = digitizer_max_time - digitizer_min_time

card.samples = numpy.ceil(digitizer_time*card.sampleRate)
card.SetParams() #warning. this may round the number of smaples to multiple of 1024

empiricalFudgeFactor = -0.096e-6
digitizerTimeOffset = tau_max + overflowBuffer + pulselength + empiricalFudgeFactor
cardTicks = scipy.arange(0, card.samples, 1.) # time in sample clocks from the digitizers point of view
cardXaxis = cardTicks /card.sampleRate - digitizerTimeOffset #time in seconds from the digitizers point of view



##############################
#data taking

################
#cailbrate the mixer elipse
###############
numPoints = 35

Idata = numpy.zeros(card.samples)
Qdata = numpy.zeros(card.samples)

Amps = numpy.zeros(numPoints)
Angles = numpy.zeros(numPoints)
Is = numpy.zeros(numPoints)
Qs = numpy.zeros(numPoints)

phases = numpy.linspace(0,360,numPoints)

#find the second pulse
cut1 = numpy.where(cardXaxis > -pulselength/2)[0][0]
cut2 = numpy.where(cardXaxis < pulselength/2)[0][-1] -1 
pulse2_range = [cut1,cut2]


for tind in range(0, numPoints):
    rfgen.set_Phase(phases[tind])
    time.sleep(0.05)
    
    card.ArmAndWait()
    Idata, Qdata = card.ReadAllData()
    
    Is2 = Idata[0,pulse2_range[0]:pulse2_range[1]]
    Qs2 = Qdata[0,pulse2_range[0]:pulse2_range[1]]
    
    Iav = numpy.mean(Is2)
    Qav = numpy.mean(Qs2)
    
    Amp = numpy.sqrt(Iav**2 + Qav**2)
    Angle = numpy.arctan2(Iav, Qav)*180/numpy.pi
    
    Amps[tind] = Amp
    Angles[tind] = Angle
    
    Is[tind] = Iav
    Qs[tind] = Qav
    
mixerAxes, mixerCenter, mixerPhi = uf.fitEllipse(Is,Qs, verbose = True)

xx, yy = uf.make_elipse(mixerAxes,  mixerCenter, mixerPhi, 150)
##########################



##################
#measure long time stability
###############
actualTimes = numpy.zeros(len(timeSpacings))

Is = numpy.zeros(len(timeSpacings))
Qs = numpy.zeros(len(timeSpacings))
Amps = numpy.zeros(len(timeSpacings))
Angles = numpy.zeros(len(timeSpacings))


plotSpacing = 30
fig1 = pylab.figure(1)
pylab.clf()



def plot_main_fig(fig):
    fig.clf()
    pylab.clf()
    ax = pylab.subplot(1,2,1)
    
    cVec = actualTimes[0:tind]
    cmap = 'cividis'
    pylab.scatter(Is[0:tind], Qs[0:tind], c = cVec, cmap = cmap , marker = 'o', s = 75,
                  vmin = 0, vmax = sum(timeSpacings), edgecolors = 'midnightblue', zorder = 2)
    
    
    pylab.plot(xx, yy, color = 'firebrick', zorder = 0)
    cbar = pylab.colorbar()
    cbar.set_label('elapsed time (s)', rotation=270)

    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    # Eliminate upper and right axes
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    # Show ticks in the left and lower axes only
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_aspect('equal')
    
    ax.set_xlim([-0.16, 0.16])
    ax.set_ylim([-0.16, 0.16])

    ax = pylab.subplot(1,2,2)
    pylab.plot(actualTimes[0:tind], Angles[0:tind], color = 'mediumblue', linestyle = '', marker = 'd', markersize = 3)
    pylab.xlabel('Time (s)')
    pylab.ylabel('Homodyne Phase (degrees)')
    
    pylab.suptitle('Pulse Homodyne Stability v. Time')
    
#        pylab.tight_layout()
    fig.canvas.draw()
    fig.canvas.flush_events()



t0 = time.time()
fig1 = pylab.figure(1)
pylab.show()
for tind in range(0, len(timeSpacings)):
    spacing = timeSpacings[tind]
#    print("current spacing = ", numpy.round(spacing,3))
    time.sleep(spacing)
    
    card.ArmAndWait()
    currT = time.time()
    print("current time = ", numpy.round(currT-t0,3))
    Idata, Qdata = card.ReadAllData()
    
    Is2 = Idata[0,pulse2_range[0]:pulse2_range[1]]
    Qs2 = Qdata[0,pulse2_range[0]:pulse2_range[1]]
    
    Iav = numpy.mean(Is2)
    Qav = numpy.mean(Qs2)
    
    Amp = numpy.sqrt(Iav**2 + Qav**2)
    Angle = numpy.arctan2(Qav, Iav)*180/numpy.pi
    
    Amps[tind] = Amp
    Angles[tind] = Angle
    
    Is[tind] = Iav
    Qs[tind] = Qav
    
    actualTimes[tind] = currT - t0
    
    
#    if tind == 0:
#        ############
#        #diagnostic figures
#        
#        fig4 = pylab.figure(4)
#        pylab.clf()
#        ax = pylab.subplot(1,1,1)
#        pylab.plot(cardXaxis*1e6,Idata[0,:] )
#        pylab.plot(cardXaxis*1e6,Qdata[0,:] )
#        
#        ta = cardXaxis[pulse2_range[0]]
#        tb = cardXaxis[pulse2_range[1]]
#        pylab.plot([ta*1e6,ta*1e6], [-0.02,0.02])
#        pylab.plot([tb*1e6,tb*1e6], [-0.03,0.03])
#        
#        #find the first pulse
#        cut1 = numpy.where(cardXaxis > -pulselength/2 - tau)[0][0] 
#        cut2 = numpy.where(cardXaxis < pulselength/2-tau)[0][-1] -1  
#        pulse1_range = [cut1,cut2]
#        
#        ta = cardXaxis[pulse1_range[0]]
#        tb = cardXaxis[pulse1_range[1]]
#        pylab.plot([ta*1e6,ta*1e6], [-0.02,0.02])
#        pylab.plot([tb*1e6,tb*1e6], [-0.03,0.03])
#        
#        pylab.xlabel('time (us)')
#        pylab.ylabel('voltage')
#        pylab.title('Single raw trace')
#        pylab.show()
#        fig4.canvas.draw()
#        fig4.canvas.flush_events()


    if numpy.mod(tind, plotSpacing) == 0:
        plot_main_fig(fig1)
        
#        thr = threading.Thread(target=plot_main_fig, kwargs = {'fig': fig1})
#        thr.start() 
        
#        fig1 = pylab.figure(1)
#        pylab.clf()
#        ax = pylab.subplot(1,2,1)
#        
#        cVec = actualTimes[0:tind]
#        cmap = 'cividis'
#        pylab.scatter(Is[0:tind], Qs[0:tind], c = cVec, cmap = cmap , marker = 'o', s = 75,
#                      vmin = 0, vmax = sum(timeSpacings), edgecolors = 'midnightblue', zorder = 2)
#        
#        
#        pylab.plot(xx, yy, color = 'firebrick', zorder = 0)
#        cbar = pylab.colorbar()
#        cbar.set_label('elapsed time (s)', rotation=270)
#    
#        # Move left y-axis and bottim x-axis to centre, passing through (0,0)
#        ax.spines['left'].set_position('center')
#        ax.spines['bottom'].set_position('center')
#        # Eliminate upper and right axes
#        ax.spines['right'].set_color('none')
#        ax.spines['top'].set_color('none')
#        # Show ticks in the left and lower axes only
#        ax.xaxis.set_ticks_position('bottom')
#        ax.yaxis.set_ticks_position('left')
#        ax.set_aspect('equal')
#        
#        ax.set_xlim([-0.16, 0.16])
#        ax.set_ylim([-0.16, 0.16])
#
#        ax = pylab.subplot(1,2,2)
#        pylab.plot(actualTimes[0:tind], Angles[0:tind], color = 'mediumblue', linestyle = '', marker = 'd', markersize = 3)
#        pylab.xlabel('Time (s)')
#        pylab.ylabel('Homodyne Phase (degrees)')
#        
#        pylab.suptitle('Pulse Homodyne Stability v. Time')
#        
##        pylab.tight_layout()
#        fig1.canvas.draw()
#        fig1.canvas.flush_events()


plot_main_fig(fig1)

#fig1 = pylab.figure(1)
#pylab.clf()
#ax = pylab.subplot(1,2,1)
#cVec = actualTimes
#cmap = 'cividis'
#pylab.scatter(Is, Qs, c = cVec, cmap = cmap , marker = 'o', s = 75,
#              vmin = 0, vmax = sum(timeSpacings), edgecolors = 'midnightblue', zorder = 2)
#
#
#pylab.plot(xx, yy, color = 'firebrick', zorder = 0)
#cbar = pylab.colorbar()
#cbar.set_label('elapsed time (s)', rotation=270)
#
## Move left y-axis and bottim x-axis to centre, passing through (0,0)
#ax.spines['left'].set_position('center')
#ax.spines['bottom'].set_position('center')
## Eliminate upper and right axes
#ax.spines['right'].set_color('none')
#ax.spines['top'].set_color('none')
## Show ticks in the left and lower axes only
#ax.xaxis.set_ticks_position('bottom')
#ax.yaxis.set_ticks_position('left')
#ax.set_aspect('equal')
#
#ax.set_xlim([-0.16, 0.16])
#ax.set_ylim([-0.16, 0.16])
#
#
#ax = pylab.subplot(1,2,2)
#pylab.plot(actualTimes, Angles, color = 'mediumblue', linestyle = '', marker = 'd', markersize = 3)
#pylab.xlabel('Time (s)')
#pylab.ylabel('Homodyne Phase (degrees)')
#
#pylab.suptitle('Pulse Homodyne Stability v. Time')
#
##        pylab.tight_layout()
#fig1.canvas.draw()
#fig1.canvas.flush_events()


print("std(angles) = " , numpy.std(Angles))



rfgen.power_Off()
logen.power_Off()  


























