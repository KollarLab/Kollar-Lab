# -*- coding: utf-8 -*-
"""
Created on Fri May  1 11:21:50 2020

@author: Kollarlab
"""

import time
from datetime import datetime

import numpy
#import sys
import pylab
from mplcursors import cursor as datacursor

#from Instruments.Acqiris import Acqiris
import userfuncs as uf
#from Instruments.HDAWG import HDAWG
#from Instruments.SGS import RFgen
from SGShelper import HDAWG_clock, SGS_coupling


def CWHomodyneStabilityTest(instruments, settings):
    '''
    General useage function that performs a homodyne phase stability measurement. Produces a figure showing mixer
    ellipticity and phase over time
    Arguments:
        ref (str) : Sets the reference used by the rf generator ('HDAWG' or 'None')
        ref_freq (int) : Sets the external reference frequency (10,100,1000 MHz)
        SGS_ref_freq (int) : Sets the reference frequency between the two SGS units (1000 MHz gives cleanest results)
        coupling (str) : Sets coupling type between SGS units ('LO': direct coupling or 'Ref' uses reference signal to lock PLL)
        measure_time (int) : Number of measurements (1 measurement~1s)
        savepath (str) : Path where figure will be saved
    '''
    ## Instruments used
    hdawg = instruments['AWG']
    logen = instruments['LO']
    rfgen = instruments['RFgen']
    card  = instruments['Digitizer']

    ## General reference configuration
    reference_signal   = settings['ref'] 
    reference_freq_MHz = settings['ref_freq'] 
    SGS_ref_freq       = settings['SGS_ref_freq'] 
    coupling_type      = settings['coupling'] 
    
    ## Measurement parameters
    measure_time = settings['measure_time']
    measDur      = settings['one_shot_time']
    freq         = settings['frequency']
    rfpower      = settings['rfpower']
    lopower      = settings['lopower']             
    freq_GHz     = freq/1e9

    ## Misc settings
    savepath = settings['savepath']

    ## HDAWG Settings
    freqs = [10e6, reference_freq_MHz*1e6]
    channels = [2,3]
    amps = [1,2]
    HDAWG_clock(hdawg, freqs, channels, amps)
    
    logen.power_Off()
    rfgen.power_Off()
    time.sleep(0.5) #Making sure that the settings have been applied
    
    ## Calibrate the mixer for reference
    mixerAxes, mixerCenter, mixerPhi = calibrate_mixer_IQ(freq, rfpower, 35, measDur = 5e-6, verbose = False, showFig=False)
    xx, yy = uf.make_elipse(mixerAxes,  mixerCenter, mixerPhi, 150)
    xmax = 1.1 * max(abs(xx))
    ymax = 1.1 * max(abs(yy))
    
    ## Digitizer card settings 
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
    
    ## SGS settings
    logen.set_Freq(freq_GHz)
    logen.set_Amp(lopower)
    logen.mod_Off()
    
    rfgen.set_Freq(freq_GHz)
    rfgen.set_Amp(rfpower)
    rfgen.mod_Off()

    SGS_coupling(logen, rfgen, reference_signal, reference_freq_MHz, coupling_type, SGS_ref_freq)

    rfgen.power_On()
    logen.power_On()
    
    # Wait and make sure that the settings have been applied
    #TO DO: write SGS function to check PLL status to avoid these long waits
    time.sleep(30)

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
    filename = 'CWHomodyne{}{}{}{}_{}'.format(reference_signal, reference_freq_MHz, coupling_type,SGS_ref_freq, stamp)
    fullpath = savepath+'\\'+filename
    figtitle = 'Homodyne Stability v. Time\n Reference Signal: {}, Ref frequency: {}, Coupling:{}, SGS ref freq:{}\n Path:{}'.format(reference_signal,reference_freq_MHz,coupling_type, SGS_ref_freq, fullpath)
    
    t0 = time.time()
    for tind in range(0, len(timeSpacings)):
        spacing = timeSpacings[tind]
        time.sleep(spacing)
        
        #Read data from digitizer
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
    
    ## Don't need this any more in theory:
    '''
    fig1 = pylab.figure(1, figsize=figure_size)
    pylab.clf()
    ax = pylab.subplot(1,2,1)

    #Plot I/Q data on parametric plot with mixer ellipse trace underneath
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
    
    ax.set_xlim([-xmax, xmax])
    ax.set_ylim([-ymax, ymax])
    
    #Phase vs time plot 
    ax = pylab.subplot(1,2,2)
    pylab.plot(actualTimes, Angles, color = 'mediumblue', linestyle = '', marker = 'd', markersize = 3)
    pylab.xlabel('Time (s)')
    pylab.ylabel('Homodyne Phase (degrees)')
    
    pylab.suptitle(figtitle)
    
    fig1.canvas.draw()
    fig1.canvas.flush_events()
    '''

    #Save figure in both .pkl and .png format for easy viewing 
    uf.savefig(fig1,filename,savepath)
    uf.savefig(fig1,filename,savepath, png = True)
    
    print("std(angles) = " , numpy.std(Angles))
    
    rfgen.power_Off()
    logen.power_Off()    

def GetDefaultSettings():
    settings = {}
    settings['ref'] = 'HDAWG'
    settings['ref_freq'] = 10
    settings['SGS_ref_freq'] = 1000
    settings['coupling'] = 'Ref'
    settings['measure_time'] = 900
    settings['savepath'] = r'C:\Users\Kollarlab\Desktop'
    settings['lopower'] = 12
    settings['rfpower'] = 0
    settings['one_shot_time'] = 1e-6
    settings['frequency'] = 8e9
    
    return settings

###########################################################################################################
# Example usage of CWHomodyneStabilityTest
###########################################################################################################

if __name__ == '__main__':
    global card
    global hdawg
    global rfgen
    global logen

    # Check and see whether instruments have already been defined, if not, create new instances
    try:
        card.samples = 1024
    except:
        card  = Acqiris('PXI23::0::0::INSTR')
        hdawg = HDAWG('dev8163')
        logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
        rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')

    #################################
    #Configuration settings
    #################################
    reference_signal   = 'HDAWG'
    coupling_type      = 'Ref'
    savepath           = r'C:\Users\Kollarlab\Desktop\CWHomodyneData'
    measure_time       = 1800

    reference_freq_MHz = 10
    SGS_ref_freq       = 1000
    CWHomodyneStabilityTest(ref = reference_signal, ref_freq = reference_freq_MHz, SGS_ref_freq = SGS_ref_freq, coupling='Ref',measure_time = measure_time, savepath=savepath)
    SGS_ref_freq       = 10
    CWHomodyneStabilityTest(ref = reference_signal, ref_freq = reference_freq_MHz, SGS_ref_freq = SGS_ref_freq, coupling='Ref',measure_time = measure_time, savepath=savepath)
