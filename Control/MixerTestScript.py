# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:33:16 2020

@author: Kollarlab

First version of test script to look at our mixer. Will need to be fixed later with programtic control of the SGSs

"""

from Instruments.HDAWG import HDAWG
from Instruments.SGS import RFgen
import numpy
import time
import sys
import scipy
import pylab
import scipy.optimize


from Acqiris_development.Acqiris import Acqiris


hardwareAddress = "PXI23::0::0::INSTR"

IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
if not IVIbinPath in sys.path:
    sys.path.append(IVIbinPath)


###############
#shamelessly stolen elimpse fitting functions
#############
import numpy as np
from numpy.linalg import eig, inv
    
def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(numpy.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a    
    
def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


#this one is maybe wrong
#def ellipse_angle_of_rotation( a ):
#    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
#    return 0.5*np.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])
    
def ellipse_angle_of_rotation2( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return np.pi/2
    else:
        if a > c:
            return np.arctan(2*b/(a-c))/2
        else:
            return np.pi/2 + np.arctan(2*b/(a-c))/2    
    
    
    
    
    
######################    







#measurement parameters
measDur = 5e-6



#setup the digitizer
#card = Acqiris(hardwareAddress)
card.triggerSlope = 'Rising'
card.triggerLevel = 0.1
card.averages = 1 #on-board averages
card.segments = 1
card.timeDelay = 0
card.activeChannels = [1,2]
card.verbose = False
card.sampleRate = 2e9
card.clockSource = 'External'
card.channelRange = 0.5
card.samples = numpy.ceil(measDur*card.sampleRate)
card.SetParams() #warning. this may round the number of smaples to multiple of 1024



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



#lo generator
#(upper, 110738)
#freq = 8 GHz
#level = 12 dBm
#rf on
#mod off
#ext ref on (for good phase), or ext ref off for random phase


#rf generator
#(lower, 110739)
#freq = 8 GHz
#level = 0 dBm
#rf on
#mod off
#ext ref on



#mixer setup
#mixer I  = digitizer channel 1
#mixer Q = digitizer channel 2



numTries = 50

Idata = numpy.zeros(card.samples)
Qdata = numpy.zeros(card.samples)

Amps = numpy.zeros(numTries)
Angles = numpy.zeros(numTries)
Is = numpy.zeros(numTries)
Qs = numpy.zeros(numTries)
for tind in range(0, numTries):
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
    
    
    
#stolen online, lest squares version

a = fitEllipse(Is,Qs)
center = ellipse_center(a)
phi = ellipse_angle_of_rotation2(a)
axes = ellipse_axis_length(a)
thetas = np.arange(0,2*np.pi, 0.01)

#print("center = ",  center)
#print("angle of rotation = ",  phi)
#print("axes = ", axes)

#fit the angle of the stig because this fit function 
# is having domain issues
#This mixer seems to have major axis along pi/4, roughly
if phi > numpy.pi/4:
    phi = phi - numpy.pi/2


xOffset = center[0]
yOffset = center[1]
stigAngle =  180*phi/numpy.pi #degrees
stig = (axes[1]-axes[0])/numpy.mean(axes)
print("center = ",  numpy.round(center,3))
print("angle of rotation = " + str(numpy.round( stigAngle, 3)) + ' degrees')
print("axes = ", axes)
print("stig = ", numpy.round(stig,3))


a, b = axes
xx = center[0] + a*np.cos(thetas)*np.cos(phi) - b*np.sin(thetas)*np.sin(phi)
yy = center[1] + a*np.cos(thetas)*np.sin(phi) + b*np.sin(thetas)*np.cos(phi)



#a = 0.041
#b = 0.037
#phi = numpy.pi/4
#xx = center[0] + a*np.cos(thetas)*np.cos(phi) - b*np.sin(thetas)*np.sin(phi)
#yy = center[1] + a*np.cos(thetas)*np.sin(phi) + b*np.sin(thetas)*np.cos(phi)



#compute residuals, something is sneakily wrong here, and stuff doesn;t land right, I don't know why.
newIs = Is - center[0]
newQs = Qs - center[1]
dataThetas = numpy.arctan2(newQs, newIs)
dataMags = numpy.sqrt(newIs**2 + newQs**2)
#fitXs = center[0] + a*np.cos(dataThetas-phi)*np.cos(phi) - b*np.sin(dataThetas-phi)*np.sin(phi)
#fitYs = center[1] + a*np.cos(dataThetas-phi)*np.sin(phi) + b*np.sin(dataThetas-phi)*np.cos(phi)
fitXs =  a*np.cos(dataThetas-phi)*np.cos(phi) - b*np.sin(dataThetas-phi)*np.sin(phi)
fitYs =  a*np.cos(dataThetas-phi)*np.sin(phi) + b*np.sin(dataThetas-phi)*np.cos(phi)
#fitXs = 0.33*np.cos(dataThetas-phi)*np.cos(phi) - 0.33*np.sin(dataThetas-phi)*np.sin(phi)
#fitYs = 0.33*np.cos(dataThetas-phi)*np.sin(phi) + 0.33*np.sin(dataThetas-phi)*np.cos(phi)
#fitXs = fitXs/4
#fitYs = fitYs/4

fitThetas = numpy.arctan2(fitYs, fitXs)
fitMags  = numpy.sqrt(fitXs**2 + fitYs**2)

residuals = numpy.sqrt( (Is - fitXs)**2 + (Qs - fitYs)**2)
meanResidual = numpy.mean(residuals)/numpy.mean(axes)
print("average residual = ", numpy.round(meanResidual, 3))


pylab.figure(2)
pylab.clf()
ax = pylab.subplot(1,1,1)
pylab.polar( dataThetas, dataMags, 'b.')
pylab.polar( fitThetas, fitMags, 'r.')
pylab.show()



pylab.figure(1)
pylab.clf()
ax = pylab.subplot(1,1,1)
pylab.plot(Is, Qs, linestyle = '', marker = 'o', markersize = 5, color = 'mediumblue')
pylab.plot(xx, yy, color = 'firebrick')
pylab.plot(fitXs, fitYs, color = 'goldenrod', marker = '+', linestyle = '')


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
#ax.set_xlim([-0.3,0.3])
#ax.set_ylim([-.3,0.3])
#pylab.xlabel('I (Voltage)')
#pylab.ylabel('Q (voltage)')
titleStr = 'Mixer performance at 8 GHz' ######!!!!!!!!!!!!! fix this
pylab.title(titleStr)
pylab.show()

    
    
    
    

######make test elipse
#elipticity = 0.9
#Eangle = -45*numpy.pi/180.  
#
##circle  
#thetas = numpy.linspace(0, 2*numpy.pi, 100)
#ampFactor = 1.05
#xs = ampFactor*numpy.mean(Amps) *numpy.cos(thetas)
#ys = ampFactor*numpy.mean(Amps) *numpy.sin(thetas)
#
##distort
#ys = elipticity*ys
#
##rotate and shift
##dx = +0.003
##dy = -0.001
#dx = numpy.mean(Is) ###!!!! works if you have enough points.
#dy = numpy.mean(Qs)
#newXs = numpy.cos(Eangle)*xs + numpy.sin(Eangle)*ys  + dx
#newYs = -numpy.sin(Eangle)*xs + numpy.cos(Eangle)*ys +  dy
    
  
    
##fit with curve_fit
####################    
##fit parameters
##stig, angle, ampFactor, xOffset, yOffset 
#
#def elipse_ff(theta, a,e, theta_0):
##    a, e = p
#    return a * (1. - e**2)/(1. - e*numpy.cos(theta - theta_0))
#
#Ioffset = numpy.mean(Is)
#Qoffset = numpy.mean(Qs)
#
#thetas = numpy.arctan2(Is - Ioffset, Qs- Qoffset )
#rs = numpy.sqrt((Is - Ioffset)**2 + (Qs- Qoffset)**2)
#
#fit_guess = [numpy.mean(Amps), 0.9, numpy.pi/4]
#fit_out, pcov = scipy.optimize.curve_fit(elipse_ff, thetas, rs, p0= fit_guess)
#
#
#pylab.figure(1)
#pylab.clf()
#ax = pylab.subplot(1,1,1)
##pylab.plot(xs, ys, color = 'dodgerblue')
##pylab.plot(newXs, newYs, color = 'firebrick')
##pylab.plot(Is, Qs, linestyle = '', marker = 'o', markersize = 5)
#
#pylab.polar(thetas, rs, 'k.')
#pylab.polar(thetas, elipse_ff(thetas, fit_guess[0], fit_guess[1], fit_guess[2]), 'r.')
##pylab.polar(thetas, elipse_ff(thetas, fit_out[0], fit_out[1], fit_out[2]), 'b.')
#
#
#ax.set_aspect('equal')
##ax.set_xlim([-0.3,0.3])
##ax.set_ylim([-.3,0.3])
#pylab.xlabel('I (Voltage)')
#pylab.ylabel('Q (voltage)')
#titleStr = 'Mixer performance at 8 GHz' ######!!!!!!!!!!!!! fix this
#pylab.title(titleStr)
#pylab.show()













