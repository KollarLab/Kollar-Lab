import numpy, time, pylab
import userfuncs as uf

def GetDefaultSettings():
    settings = {}
    settings['frequency']      = 8e9
    settings['lopower']        = 12
    settings['rfpower']        = 0
    settings['NumPoints']      = 30
    settings['SingleShotTime'] = 5e-6
    settings['verbose']        = False
    settings['showFig']        = False
    return settings

def calibrate_SGS_IQ(instruments, settings):

    logen = instruments['LO']
    rfgen = instruments['RFgen']
    card  = instruments['Digitizer'] 
    hdawg = instruments['AWG']

    freq      = settings['frequency']
    lopower   = settings['lopower']
    rfpower   = settings['rfpower']
    numPoints = settings['NumPoints']
    measDur   = settings['SingleShotTime']
    verbose   = settings['verbose']
    showFig   = settings['showFig']

    freq_GHz = freq/1e9
    phases = numpy.linspace(0,2*numpy.pi,numPoints)

    ## HDAWG settings
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    progFile = open("HDAWG_sequencer_codes/IQcontrol.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    loadprog = loadprog.replace('_time_', str(measDur))

    ## Digitizer card settings 
    card.triggerSlope = 'Rising'
    card.triggerLevel = 0.1
    card.clockSource = 'External'
    card.verbose = False
    
    card.triggerDelay = 0
    card.activeChannels = [1,2]
    card.channelRange = 0.5
    card.sampleRate = 2e9
    
    card.averages = 1 #on-board averages
    card.segments = 1
    card.samples = numpy.ceil(measDur*card.sampleRate)
    card.SetParams() #warning. this may round the number of smaples to multiple of 1024
    
    ## SGS unit settings
    logen.set_Freq(freq_GHz)
    logen.set_Amp(lopower)
    logen.mod_Off()
    logen.power_On() 
    
    rfgen.set_Freq(freq_GHz)
    rfgen.set_Amp(rfpower)
    rfgen.mod_Off()
    rfgen.power_On()

    ## Wait for settings to percolate
    time.sleep(0.5)
    
    Idata = numpy.zeros(card.samples)
    Qdata = numpy.zeros(card.samples)
    
    Amps = numpy.zeros(numPoints)
    Angles = numpy.zeros(numPoints)
    Is = numpy.zeros(numPoints)
    Qs = numpy.zeros(numPoints)
    
    for tind in range(0, numPoints):
        finalprog = loadprog.replace('_phase_', str(phases[tind]))
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)
        
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
        xx, yy = uf.make_elipse(axes, center, phi, 150)
        
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