# -*- coding: utf-8 -*-
"""
Created on Fri May  8 14:11:28 2020

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

import threading

from userfuncs import freeze
import userfuncs as uf


from Acqiris_development.Acqiris import Acqiris


hardwareAddress = "PXI23::0::0::INSTR"

IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
if not IVIbinPath in sys.path:
    sys.path.append(IVIbinPath)
    
    
    
######################    
#measurement parameters
measDur = 5e-6

numFreqs = 20
freqs = numpy.linspace(4e9,10e9, numFreqs)
freqs = numpy.flipud(freqs)

numPoints = 15
#phases = numpy.linspace(0, numpy.pi,numPoints)
phases = numpy.linspace(0, 360,numPoints)



#setup the digitizer
#card = Acqiris(hardwareAddress)
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
logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
logen.set_Freq(8)
logen.set_Amp(12)
logen.mod_Off()
#logen.set_Internal_Reference()
logen.set_External_Reference()
logen.power_On()


#rf generator
#(lower, 110739)
#freq = 8 GHz
#level = 0 dBm
#rf on
#mod off
#ext ref on
rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')
rfgen.set_Freq(8)
rfgen.set_Amp(-4)
rfgen.mod_Off()
rfgen.set_External_Reference()
rfgen.power_On()



def plot_fig1():
    fig = pylab.figure(1)
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
    titleStr = 'Mixer performance at ' + str(numpy.round(freq_GHz, 3)) + ' GHz'
    pylab.title(titleStr)
#    pylab.show(block = False)
    
    datacursor()
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return


def plot_main_fig(fig):
    fig.clf()
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
    titleStr = 'Mixer performance at ' + str(numpy.round(freq_GHz, 3)) + ' GHz'
    pylab.title(titleStr)
#    pylab.show(block = False)
    
    datacursor()
    
    pylab.title('thread test figure')
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return

def thread_test():
    print(1)
    time.sleep(1)
    print(2)
    time.sleep(1)
    print(3)
    time.sleep(1)
    return

def thread_fig(fig):
    ax= pylab.subplot(1,1,1)
    xs = numpy.linspace(-5,5,50)
    ys = xs**2
    pylab.plot(xs, ys)
    datacursor()
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return



stigVec = numpy.zeros(len(freqs))
phiVec = numpy.zeros(len(freqs))

for find in range(0, len(freqs)):
    freq = freqs[find]
    freq_GHz = freq/1e9
    rfgen.set_Freq(freq_GHz)
    logen.set_Freq(freq_GHz)
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
        
        
        
    
    mixerAxes, mixerCenter,  mixerPhi = uf.fitEllipse(Is,Qs, verbose = True)
    xx, yy = uf.make_elipse(mixerAxes,  mixerCenter, mixerPhi, 150)
    
    
    stig = (mixerAxes[1]-mixerAxes[0])/numpy.mean(mixerAxes)
    stigVec[find] = stig
    phiVec[find] = mixerPhi
    
    
    
#    fig = pylab.figure(1)
#    pylab.clf()
#    ax = pylab.subplot(1,1,1)
#    pylab.plot(Is, Qs, linestyle = '', marker = 'o', markersize = 5, color = 'mediumblue')
#    pylab.plot(xx, yy, color = 'firebrick')
#    
#    
#    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
#    ax.spines['left'].set_position('center')
#    ax.spines['bottom'].set_position('center')
#    
#    # Eliminate upper and right axes
#    ax.spines['right'].set_color('none')
#    ax.spines['top'].set_color('none')
#    
#    # Show ticks in the left and lower axes only
#    ax.xaxis.set_ticks_position('bottom')
#    ax.yaxis.set_ticks_position('left')
#    
#    
#    ax.set_aspect('equal')
#    titleStr = 'Mixer performance at ' + str(numpy.round(freq_GHz, 3)) + ' GHz'
#    pylab.title(titleStr)
##    pylab.show(block = False)
#    
#    datacursor()
#    
#    fig.canvas.draw()
#    fig.canvas.flush_events()
    
    
    plot_fig1()
#    thr = threading.Thread(target=thread_test)
    
#    if numpy.mod(find,4) == 0:
    if find == 0:
        fig8 = pylab.figure(8)
        ax = pylab.subplot(1,1,1)
        pylab.plot([1,2], [3,4])
        pylab.show()
#        thr = threading.Thread(target=thread_fig, kwargs = {'fig': fig8})
        thr = threading.Thread(target=plot_main_fig, kwargs = {'fig': fig8})
        thr.start() 

 
    
stigVec_dB = numpy.log10(stigVec+1)*10
    
fig2 = pylab.figure(2)
pylab.clf()

ax = pylab.subplot(2,2,1)
pylab.plot(freqs/1e9, stigVec, 'b.')
pylab.xlabel('Frequency (GHz)')
pylab.ylabel('Astigmatism (linear)')
pylab.title('Linear Astigmatism')

ax = pylab.subplot(2,2,2)
pylab.plot(freqs/1e9, stigVec_dB, 'r.')
pylab.xlabel('Frequency (GHz)')
pylab.ylabel('Astigmatism (dB)')
pylab.title('Log Astigmatism')



ax = pylab.subplot(2,2,3)
pylab.plot(freqs/1e9, 180*phiVec/numpy.pi, 'b.')
pylab.xlabel('Frequency (GHz)')
pylab.ylabel('Astigmatism Angle (degrees)')
pylab.title('Absolute Astigmatism Angle')

ax = pylab.subplot(2,2,4)
pylab.plot(freqs/1e9, 180*phiVec/numpy.pi - 45, 'r.')
pylab.xlabel('Frequency (GHz)')
pylab.ylabel('Astigmatism Angle (degrees) - 45')
pylab.title('IQ Angle Imbalance')


pylab.suptitle('Mixer Calibration')
pylab.tight_layout()
pylab.show()
    
    
    

rfgen.power_Off()
logen.power_Off()   