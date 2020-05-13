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



################
##shamelessly stolen elipse fitting functions
##############
#import numpy as np
#from numpy.linalg import eig, inv
#    
#def fitEllipse(x,y):
#    x = x[:,np.newaxis]
#    y = y[:,np.newaxis]
#    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
#    S = np.dot(D.T,D)
#    C = np.zeros([6,6])
#    C[0,2] = C[2,0] = 2; C[1,1] = -1
#    E, V =  eig(numpy.dot(inv(S), C))
#    n = np.argmax(np.abs(E))
#    a = V[:,n]
#    return a    
#    
#def ellipse_center(a):
#    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
#    num = b*b-a*c
#    x0=(c*d-b*f)/num
#    y0=(a*f-b*d)/num
#    return np.array([x0,y0])
#
#
##this one is maybe wrong
##def ellipse_angle_of_rotation( a ):
##    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
##    return 0.5*np.arctan(2*b/(a-c))
#
#
#def ellipse_axis_length( a ):
#    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
#    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
#    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
#    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
#    res1=np.sqrt(up/down1)
#    res2=np.sqrt(up/down2)
#    return np.array([res1, res2])
#    
#def ellipse_angle_of_rotation2( a ):
#    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
#    if b == 0:
#        if a > c:
#            return 0
#        else:
#            return np.pi/2
#    else:
#        if a > c:
#            return np.arctan(2*b/(a-c))/2
#        else:
#            return np.pi/2 + np.arctan(2*b/(a-c))/2    
#        
###########################
            
        
#def make_elipse(axes, phi, center,  numPoints):
#    a, b = axes
#    thetas = np.linspace(0,2*np.pi, numPoints)
#    
#    xx = center[0] + a*np.cos(thetas)*np.cos(phi) - b*np.sin(thetas)*np.sin(phi)
#    yy = center[1] + a*np.cos(thetas)*np.sin(phi) + b*np.sin(thetas)*np.cos(phi)
#    return xx, yy
    
    
def calibrate_mixer_IQ(freq, power, numPoints, measDur = 5e-6, verbose = False):
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
    #logen.set_Internal_Reference()
#    logen.set_External_Reference()
    logen.power_On() 
    
    rfgen.set_Freq(freq_GHz)
    rfgen.set_Amp(power)
    rfgen.mod_Off()
#    rfgen.set_External_Reference()
    rfgen.power_On()

    time.sleep(0.05)
    
    Idata = numpy.zeros(card.samples)
    Qdata = numpy.zeros(card.samples)
    
    Amps = numpy.zeros(numPoints)
    Angles = numpy.zeros(numPoints)
    Is = numpy.zeros(numPoints)
    Qs = numpy.zeros(numPoints)
    
    for tind in range(0, numPoints):
        rfgen.set_Phase(phases[tind])
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
        
    
#    #stolen online, lest squares elipse fit
#    a = fitEllipse(Is,Qs)
#    center = ellipse_center(a)
#    phi = ellipse_angle_of_rotation2(a)
#    axes = ellipse_axis_length(a)
#    
#    
#    #fix the angle of the stig because this fit function 
#    # is having domain issues
#    #This mixer seems to have major axis along pi/4, roughly
#    if phi > numpy.pi/4:
#        phi = phi - numpy.pi/2
#        
#    if axes[1] > axes[0]:
#        #major axis is second
#        axes = [axes[1], axes[0]]
#        phi = phi + numpy.pi/2
#        
#    if phi < 0:
#        phi = phi +numpy.pi
#    
##    xOffset = center[0]
##    yOffset = center[1]
#    stigAngle =  180*phi/numpy.pi #degrees
#    stig = (axes[1]-axes[0])/numpy.mean(axes)
#    
#    if verbose:
#        print("    ")
#        print("frequency = ", freq_GHz)
#        print("center = ",  numpy.round(center,3))
#        print("angle of rotation = " + str(numpy.round( stigAngle, 3)) + ' degrees')
#        print("axes = ", numpy.round(axes,3))
#        print("stig = ", numpy.round(stig,3))
        
    axes, center, phi = uf.fitEllipse(Is,Qs, verbose = True)
    
    return axes, center, phi
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
global card
global rfgen
global logen

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
    #hdawg = HDAWG('dev8163') #HDAWG device name
    ##hdawg.AWGs[0].samplerate = '2.4GHz'
    ##hdawg.channelgrouping = '1x4'
    ##hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
    ##hdawg.Channels[1].configureChannel(marker_out='Trigger', hold='True')
    ##hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    ###hdawg.daq.setInt('/dev8163/awgs/0/outputs/0/hold',1)
    ###hdawg.daq.setInt('/dev8163/awgs/0/outputs/1/hold',1)
    #hdawg.OSCs[1].freq = 10e6
    #hdawg.Channels[2].analog_outs = [0.5,0]
    #hdawg.Channels[3].analog_outs = [0,0.5]
    #hdawg.Channels[2].configureChannel(amp=1.0)
    #hdawg.Channels[3].configureChannel(amp=1.0)

    
    
    logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
    rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')
    
    logen.power_Off()
    rfgen.power_Off()
    


############
#measurement params
###########
measDur = 1e-6
freq = 8e9
power = 2


#rfgen.set_Internal_Reference() #RF gen by itself, to it's own drum
#logen.set_External_Reference() #LO gen following RF gen
#time.sleep(0.5)



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
#logen.set_External_Reference()
logen.power_On() 

rfgen.set_Freq(freq_GHz)
rfgen.set_Amp(power)
rfgen.mod_Off()
#rfgen.set_External_Reference()
rfgen.power_On()

time.sleep(0.05)


#temp1 = numpy.linspace(0, 10, 10)
#temp2 = numpy.linspace(10, 60*1, 4)
#timeSpacings = numpy.concatenate((temp1, temp2))

#temp3 = numpy.linspace(0, 1, 25)
#timeSpacings = temp3

temp4 = numpy.ones(25)
temp5 = numpy.ones(25)*5
temp6 = numpy.ones(10)*30
#timeSpacings = numpy.concatenate((temp4, temp5))
#timeSpacings = numpy.concatenate((temp4, temp5, temp6)) #ultralong

#timeSpacings = temp4

timeSpacings = numpy.ones(10)


actualTimes = numpy.zeros(len(timeSpacings))

Is = numpy.zeros(len(timeSpacings))
Qs = numpy.zeros(len(timeSpacings))
Amps = numpy.zeros(len(timeSpacings))
Angles = numpy.zeros(len(timeSpacings))


plotSpacing = 30
fig1 = pylab.figure(1)
pylab.clf()


#rfgen.set_Internal_Reference() #RF gen by itself, to it's own drum
#logen.set_External_Reference() #LO gen following RF gen
#time.sleep(0.5)

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
#    Angle = numpy.arctan2(Iav, Qav)*180/numpy.pi
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

pylab.suptitle('Homodyne Stability v. Time')

#        pylab.tight_layout()
fig1.canvas.draw()
fig1.canvas.flush_events()


print("std(angles) = " , numpy.std(Angles))



rfgen.power_Off()
logen.power_Off()    














     
        
        
        