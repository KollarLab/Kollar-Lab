# -*- coding: utf-8 -*-
"""
Created on Fri May  1 11:21:50 2020

@author: Kollarlab
"""

import time
import numpy
import pylab
import scipy
from datetime import datetime
from mplcursors import cursor as datacursor

import userfuncs as uf
import Calibration.mixerIQcal as mixer
from SGShelper import HDAWG_clock, SGS_coupling


def GetDefaultSettings():
    settings = {}
    settings['measure_time']  = 900
    settings['savepath']      = r'C:\Users\Kollarlab\Desktop\CWHeterodyneData'
    settings['lopower']       = 12
    settings['rfpower']       = 0
    settings['one_shot_time'] = 1e-6
    settings['frequency_lo']  = 8e9
    settings['frequency_rf']  = 8.04e9
    settings['save']          = True

    #Mixer settings
    settings['mixerIQcal']    = mixer.GetDefaultSettings()

    #Card settings
    settings['triggerDelay']   = 0
    settings['activeChannels'] = [1,2]
    settings['channelRange ']  = 0.5
    settings['sampleRate']     = 2e9
    
    settings['averages'] = 1 
    settings['segments'] = 1
    
    return settings

def CWHomodyneStabilityTest(instruments, settings):

    ## Instruments used
    hdawg = instruments['AWG']
    logen = instruments['LO']
    rfgen = instruments['RFgen']
    card  = instruments['Digitizer']

    ## Measurement parameters
    measure_time = settings['measure_time']
    measDur      = settings['one_shot_time']
    freq_lo      = settings['frequency_lo']
    freq_rf      = settings['frequency_rf']
    rfpower      = settings['rfpower']
    lopower      = settings['lopower']             
    freq_GHz_lo  = freq_lo/1e9
    freq_GHz_rf  = freq_rf/1e9

    ## Misc settings
    savepath = settings['savepath']
    
    logen.power_Off()
    rfgen.power_Off()
    time.sleep(0.5) #Making sure that the settings have been applied
    
    ## Calibrate the mixer for reference
    mixerAxes, mixerCenter, mixerPhi = mixer.calibrate_mixer_IQ(instruments, settings['mixerIQcal'])
    xx, yy = uf.make_elipse(mixerAxes,  mixerCenter, mixerPhi, 150)
    xmax = 1.1 * max(abs(xx))
    ymax = 1.1 * max(abs(yy))
    
    ## Digitizer card settings 
    card.triggerSlope = 'Rising'
    card.triggerLevel = 0.1
    card.clockSource = 'External'
    card.verbose = False
    
    card.triggerDelay   = settings['triggerDelay']
    card.activeChannels = settings['activeChannels']
    card.channelRange   = settings['channelRange']
    card.sampleRate     = settings['sampleRate']
    
    card.averages = settings['averages']
    card.segments = settings['segments']

    card.samples = numpy.ceil(measDur*card.sampleRate)
    card.SetParams() #warning. this may round the number of smaples to multiple of 1024
    
    ## SGS settings
    logen.set_Freq(freq_GHz_lo)
    logen.set_Amp(lopower)
    logen.mod_Off()
    
    rfgen.set_Freq(freq_GHz_rf)
    rfgen.set_Amp(rfpower)
    rfgen.mod_Off()

    rfgen.power_On()
    logen.power_On()

    #################################################
    # Measurement and analysis
    #################################################
    
    short = numpy.ones(measure_time) #multiply by an int to get larger time spacing (1 second at the moment)
    timeSpacings = short #if desired, can concatenate multiple different time ranges/ spacings for especially large sets
    actualTimes = numpy.zeros(len(timeSpacings)) #array to hold the 'real' wall time ellapsed in measurement
    
    Is = numpy.zeros(len(timeSpacings)) #array to hold I data
    Qs = numpy.zeros(len(timeSpacings)) #array to hold Q data
    Amps = numpy.zeros(len(timeSpacings)) #array to hold amplitude data
    Angles = numpy.zeros(len(timeSpacings)) #array to hold phase data
    
    plotSpacing = 30
    figure_size = (14,8)
    fig1 = pylab.figure(1, figsize=figure_size)
    pylab.clf()

    date  = datetime.now()
    stamp = date.strftime('%Y%m%d_%H%M%S')
    filename = 'CWHeterodyne{}_{}'.format(measure_time, stamp)
    fullpath = savepath+'\\'+filename
    figtitle = 'Heterodyne Stability v. Time\n Total Measure Time: {}s\n Path:{}'.format(measure_time, fullpath)
    
    t0 = time.time()

    #Make a synthetic cos/sin signal and digitally mix down to DC
    digitizertime = scipy.arange(0, card.samples,1.)/card.sampleRate
    freq = numpy.abs(freq_lo-freq_rf)
    CosPC = numpy.cos(2*numpy.pi()*freq*digitizertime)
    SinPC = numpy.sin(2*numpy.pi()*freq*digitizertime)
    Idata = numpy.zeros(card.segments)
    Qdata = numpy.zeros(card.segments)

    for tind in range(0, len(timeSpacings)):
        spacing = timeSpacings[tind]
        time.sleep(spacing)
        
        #Read data from digitizer
        card.ArmAndWait()
        currT = time.time()
        print("current time = ", numpy.round(currT-t0,3))
        Ch1, Ch2 = card.ReadAllData()
        
        for i in range(card.segments):
            Idata[i] = CosPC*Ch1[i,:]
            Qdata[i] = SinPC*Ch1[i,:]

        Iav = numpy.mean(Idata)
        Qav = numpy.mean(Qdata)
        
        Amp = numpy.sqrt(Iav**2 + Qav**2)
        Angle = numpy.arctan2(Qav, Iav)*180/numpy.pi
        
        Amps[tind] = Amp
        Angles[tind] = Angle
        
        Is[tind] = Iav
        Qs[tind] = Qav
        
        actualTimes[tind] = currT - t0
    
        # Every plotspacing iterations, we want to update the figure and on the last iteration
        if numpy.mod(tind, plotSpacing) == 0 or tind==len(timeSpacings-1):
            #fig1 = pylab.figure(1)
            pylab.clf()
            ax = pylab.subplot(1,2,1)
            
            #Plot I/Q data on parameteric plot to show IQ vector
            cVec = actualTimes[0:tind]
            cmap = 'cividis'
            pylab.scatter(Is[0:tind], Qs[0:tind], c = cVec, cmap = cmap , marker = 'o', s = 75,
                          vmin = 0, vmax = sum(timeSpacings), edgecolors = 'midnightblue', zorder = 2)
            
            #Plot mixer ellipse on behind data for reference 
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
            
            ax.set_xlim([-xmax, xmax])
            ax.set_ylim([-ymax, ymax])
    
            #Plot phase data on standard phase vs t plot
            ax = pylab.subplot(1,2,2)
            pylab.plot(actualTimes[0:tind], Angles[0:tind], color = 'mediumblue', linestyle = '', marker = 'd', markersize = 3)
            pylab.xlabel('Time (s)')
            pylab.ylabel('Homodyne Phase (degrees)')
            pylab.suptitle(figtitle)
            
            fig1.canvas.draw()
            fig1.canvas.flush_events()
    
    #Save png of figure for easy viewing
    uf.savefig(fig1,filename,savepath, png = True)
    
    print("std(angles) = " , numpy.std(Angles))

    rfgen.power_Off()
    logen.power_Off()    

    dataTosave = ['Is','Qs','Amps','Angles', 'actualTimes','xx','yy','Idata', 'Qdata', 'Iav', 'Qav']
    figsTosave = [fig1]

    if settings['save']:    
        uf.SaveFull(savepath, filename, dataTosave, locals(), settings, instruments, figsTosave) 