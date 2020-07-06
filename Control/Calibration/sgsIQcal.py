import numpy, time, pylab
import userfuncs as uf

def GetDefaultSettings():
    settings = {}
    settings['frequency']      = 8e9
    settings['lopower']        = 12
    settings['rfpower']        = 0
    settings['NumPoints']      = 30
    settings['SingleShotTime'] = 10e-6
    settings['AngleImp']       = 0
    settings['AmpImp']         = 1.
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
    angle_impairment = settings['AngleImp']
    amp_impairment   = settings['AmpImp']
    measDur   = settings['SingleShotTime']
    showFig   = settings['showFig']
    verbose   = settings['verbose']

    freq_GHz = freq/1e9
    #phases = numpy.linspace(0,2*numpy.pi,numPoints)
    amp_impairment_Q = 1.
    amp_impairment_I = 1.
    
    if amp_impairment>1:
        amp_impairment_Q = 1./amp_impairment
        if verbose:
            print('Q impairment: {}'.format(amp_impairment_Q))
    else:
        amp_impairment_I = amp_impairment
        if verbose:
            print('I impairment: {}'.format(amp_impairment_I))
    if verbose:
        print('Angle impairment: {}'.format(angle_impairment))
    ## HDAWG settings
        print('HDAWG config')
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    if verbose:
        print('HDAWG config ended')
    progFile = open("HDAWG_sequencer_codes/IQcontrol.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    loadprog = loadprog.replace('_time_', str(1.1*measDur))
    loadprog = loadprog.replace('_npoints_', str(numPoints))
    loadprog = loadprog.replace('_angimp_', str(angle_impairment))
    loadprog = loadprog.replace('_ampimpI_', str(amp_impairment_I))
    loadprog = loadprog.replace('_ampimpQ_', str(amp_impairment_Q))
    if verbose:
        print('Program loading')
    hdawg.AWGs[0].load_program(loadprog)
    hdawg.AWGs[0].run_loop()
    if verbose:
        print('Program loaded')
    
    ## Digitizer card settings 
    card.triggerSlope = 'Rising'
    card.triggerLevel = 0.1
    card.clockSource = 'External'
    card.verbose = False
    
    card.triggerDelay = 0.2e-6
    card.activeChannels = [1,2]
    card.channelRange = 2.5
    card.sampleRate = 2e9
    
    card.averages = 1 #on-board averages
    card.segments = numPoints
    card.samples = numpy.ceil(measDur*card.sampleRate)
    card.SetParams() #warning. this may round the number of smaples to multiple of 1024
    
    ## SGS unit settings
    logen.set_Freq(freq_GHz)
    logen.set_Amp(lopower)
    logen.mod_Off()
    logen.power_On() 
    
    rfgen.set_Freq(freq_GHz)
    rfgen.set_Amp(rfpower)
    rfgen.mod_On()
    rfgen.power_On()

    ## Wait for settings to percolate
    time.sleep(0.5)
    
    Idata = numpy.zeros((card.segments, card.samples))
    Qdata = numpy.zeros((card.segments, card.samples))
    
    Amps = numpy.zeros(numPoints)
    Angles = numpy.zeros(numPoints)
    Is = numpy.zeros(numPoints)
    Qs = numpy.zeros(numPoints)
    if verbose:
        print('Acquisition started')
    card.ArmAndWait()
    Idata, Qdata = card.ReadAllData()
    if verbose:
        print('Acquisition ended')
    for tind in range(0, numPoints):

        Iav = numpy.mean(Idata[tind])
        Qav = numpy.mean(Qdata[tind])
        
        Amp = numpy.sqrt(Iav**2 + Qav**2)
        Angle = numpy.arctan2(Iav, Qav)*180/numpy.pi
        
        Amps[tind] = Amp
        Angles[tind] = Angle
        
        Is[tind] = Iav
        Qs[tind] = Qav

    axes, center, phi, ecc = uf.fit_ell_martin(Is,Qs, verbose)
    
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
        titleStr = 'SGS IQ performance'
        pylab.title(titleStr)
    #    pylab.show(block = False)
        
        fig.canvas.draw()
        fig.canvas.flush_events()
        
    return axes, center, phi, ecc