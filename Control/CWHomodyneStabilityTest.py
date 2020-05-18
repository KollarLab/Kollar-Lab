# -*- coding: utf-8 -*-
"""
Created on Fri May  1 11:21:50 2020

@author: Kollarlab
"""

from Instruments.HDAWG import HDAWG
from Instruments.SGS import RFgen
import numpy
import time
import sys
import scipy
import pylab
import scipy.optimize
from mplcursors import cursor as datacursor


from userfuncs import freeze
import userfuncs as uf


from Acqiris_development.Acqiris import Acqiris


IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
if not IVIbinPath in sys.path:
    sys.path.append(IVIbinPath)
    
    
def calibrate_mixer_IQ(freq, power, numPoints, measDur = 5e-6, verbose = False, showFig = True):
    freq_GHz = freq/1e9
    phases = numpy.linspace(0,360,numPoints)
    
    card.triggerSlope = 'Rising'
    card.triggerLevel = 0.1
    card.averages = 1 #on-board averages
    card.segments = 1
    card.triggerDelay = 0
    card.activeChannels = [1,2]
    card.verbose = False
    card.sampleRate = 2e9
    card.clockSource = 'External'
    card.channelRange = 0.5
    card.samples = numpy.ceil(measDur*card.sampleRate)
    card.SetParams() #warning. this may round the number of smaples to multiple of 1024
    
    logen.set_Freq(freq_GHz)
    logen.set_Amp(12)
    logen.mod_Off()
    logen.power_On() 
    
    rfgen.set_Freq(freq_GHz)
    rfgen.set_Amp(power)
    rfgen.mod_Off()
    rfgen.power_On()

    time.sleep(0.05)
    
    Idata = numpy.zeros(card.samples)
    Qdata = numpy.zeros(card.samples)
    
    Amps = numpy.zeros(numPoints)
    Angles = numpy.zeros(numPoints)
    Is = numpy.zeros(numPoints)
    Qs = numpy.zeros(numPoints)
    
    for tind in range(0, numPoints):
#        logen.set_Phase(phases[tind])
        time.sleep(0.05)
        
        card.ArmAndWait()
        Idata, Qdata = card.ReadAllData()
        
        Iav = numpy.mean(Idata)
        Qav = numpy.mean(Qdata)
        
        Amp = numpy.sqrt(Iav**2 + Qav**2)
        Angle = numpy.arctan2(Iav, Qav)*180/numpy.pi
        
        Amps[tind] = Amp
        Angles[tind] = Angle
        
        Is[tind] = Iav
        Qs[tind] = Qav
        
    axes, center, phi = uf.fitEllipse(Is,Qs, verbose = True)
    
    if showFig:
        xx, yy = uf.make_elipse(axes,  center, phi, 150)
        
        fig = pylab.figure(10)
        pylab.clf()
        ax = pylab.subplot(1,1,1)
        pylab.plot(Is, Qs, linestyle = '', marker = 'o', markersize = 5, color = 'mediumblue')
        pylab.plot(xx, yy, color = 'firebrick')
        
        
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
        titleStr = 'Mixer performance '
        pylab.title(titleStr)
    #    pylab.show(block = False)
        
        fig.canvas.draw()
        fig.canvas.flush_events()
    
    
    return axes, center, phi
    
    
global card
global rfgen
global logen

#################################
#Configuration settings
#################################
reference_signal   = 'None'
reference_freq_MHz = 100
coupling_type      = 'Ref'
#Timing
short = numpy.ones(1800)

#initializeHardware = True
#initializeHardware = False
try:
    card.samples == 1024
except:
#if initializeHardware:
    
    hardwareAddress = "PXI23::0::0::INSTR"

    
    card = Acqiris(hardwareAddress)
#       
    ##set up the HDAWG. 
    ##in this case, we just need channels 3,4 for our fake clock
    #### Connect to HDAWG and initialize it 
    hdawg = HDAWG('dev8163') #HDAWG device name
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[2].analog_outs = [0.5,0]
    hdawg.Channels[3].analog_outs = [0,1]
    hdawg.Channels[2].configureChannel(amp=1.0)
    hdawg.Channels[3].configureChannel(amp=2.0)
    
    
    logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
    rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')
    
    logen.power_Off()
    rfgen.power_Off()
    

hdawg.OSCs[1].configure_sine(0,10e6)
hdawg.OSCs[1].configure_sine(1,reference_freq_MHz*1e6)

############
#measurement params
###########
measDur = 1e-6
freq = 8e9
power = 0


rfgen.set_Internal_Reference() #RF gen by itself, to it's own drum
logen.set_Internal_Reference() #LO gen following RF gen
time.sleep(0.5)



#cailbrate the mixer for reference
mixerAxes, mixerCenter, mixerPhi = calibrate_mixer_IQ(freq, power, 35, measDur = 5e-6, verbose = False)
xx, yy = uf.make_elipse(mixerAxes,  mixerCenter, mixerPhi, 150)

#########
#now look at homodyne versus time.
freq_GHz = freq/1e9

card.triggerSlope = 'Rising'
card.triggerLevel = 0.1
card.averages = 1 #on-board averages
card.segments = 1
card.triggerDelay = 0
card.activeChannels = [1,2]
card.verbose = False
card.sampleRate = 2e9
card.clockSource = 'External'
card.channelRange = 0.5
card.samples = numpy.ceil(measDur*card.sampleRate)
card.SetParams() #warning. this may round the number of smaples to multiple of 1024

logen.set_Freq(freq_GHz)
logen.set_Amp(12)
logen.mod_Off()

rfgen.set_Freq(freq_GHz)
rfgen.set_Amp(power)
rfgen.mod_Off()
if reference_signal == 'HDAWG':
    rfgen.set_External_Reference(freq = reference_freq_MHz)
else:
    rfgen.set_Internal_Reference()
if coupling_type == 'LO':
    rfgen.set_RefLO_output(output='LO')
else:
    rfgen.set_RefLO_output(output='Ref', freq = reference_freq_MHz)
rfgen.power_On()

if coupling_type == 'Ref':
    logen.set_External_Reference(freq = reference_freq_MHz)
    logen.set_Internal_LO()
if coupling_type == 'LO':
    logen.set_External_LO()
    
logen.power_On() 

time.sleep(30)

timeSpacings = short


actualTimes = numpy.zeros(len(timeSpacings))

Is = numpy.zeros(len(timeSpacings))
Qs = numpy.zeros(len(timeSpacings))
Amps = numpy.zeros(len(timeSpacings))
Angles = numpy.zeros(len(timeSpacings))


plotSpacing = 30
fig1 = pylab.figure(1)
pylab.clf()

t0 = time.time()
for tind in range(0, len(timeSpacings)):
    spacing = timeSpacings[tind]
#    print("current spacing = ", numpy.round(spacing,3))
    time.sleep(spacing)
    
    
    
    card.ArmAndWait()
    currT = time.time()
    print("current time = ", numpy.round(currT-t0,3))
    Idata, Qdata = card.ReadAllData()
    
    Iav = numpy.mean(Idata)
    Qav = numpy.mean(Qdata)
    
    Amp = numpy.sqrt(Iav**2 + Qav**2)
    Angle = numpy.arctan2(Qav, Iav)*180/numpy.pi
    
    Amps[tind] = Amp
    Angles[tind] = Angle
    
    Is[tind] = Iav
    Qs[tind] = Qav
    
    actualTimes[tind] = currT - t0


    if numpy.mod(tind, plotSpacing) == 0:
        fig1 = pylab.figure(1)
        pylab.clf()
        ax = pylab.subplot(1,2,1)
#        pylab.plot(Is[0:tind], Qs[0:tind], linestyle = '', marker = 'o', markersize = 5, color = 'mediumblue')
        
#        cVec = actualTimes[0:tind]/sum(timeSpacings)
#        cmap = 'cividis'
#        pylab.scatter(Is[0:tind], Qs[0:tind], c = cVec, cmap = cmap , marker = 'o', s = 75, edgecolors = 'midnightblue', zorder = 2)
        
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
        
        pylab.suptitle('Homodyne Stability v. Time')
        
#        pylab.tight_layout()
        fig1.canvas.draw()
        fig1.canvas.flush_events()



fig1 = pylab.figure(1)
pylab.clf()
ax = pylab.subplot(1,2,1)
cVec = actualTimes
cmap = 'cividis'
pylab.scatter(Is, Qs, c = cVec, cmap = cmap , marker = 'o', s = 75,
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
pylab.plot(actualTimes, Angles, color = 'mediumblue', linestyle = '', marker = 'd', markersize = 3)
pylab.xlabel('Time (s)')
pylab.ylabel('Homodyne Phase (degrees)')

pylab.suptitle('Homodyne Stability v. Time\n Reference Signal: {}, Ref frequency: {}, Coupling:{}'.format(reference_signal,reference_freq_MHz,coupling_type))

#        pylab.tight_layout()
fig1.canvas.draw()
fig1.canvas.flush_events()


print("std(angles) = " , numpy.std(Angles))

rfgen.power_Off()
logen.power_Off()    
