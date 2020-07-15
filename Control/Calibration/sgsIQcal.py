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

def calibrate_SGS_IQ_fast(instruments, settings):

    logen   = instruments['LO']
    rfgen   = instruments['RFgen']
    card    = instruments['Digitizer'] 
    hdawg   = instruments['AWG']
    trigger = instruments['Trigger']

    freq      = settings['frequency']
    lopower   = settings['lopower']
    rfpower   = settings['rfpower']
    numPoints = settings['NumPoints']
    angmin = settings['MinAngle']
    angmax = settings['MaxAngle']
    angnum = settings['AngleNum']
    measDur   = settings['SingleShotTime']
    amp_ramp_points = settings['Amp_points']
    verbose   = settings['verbose']

    freq_GHz = freq/1e9
    
    ## HDAWG config
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    # Program config
    progFile = open("HDAWG_sequencer_codes/IQcontrol.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    loadprog = loadprog.replace('_time_', str(1.1*measDur))
    loadprog = loadprog.replace('_npoints_', str(numPoints))
    loadprog = loadprog.replace('_amppoints_', str(amp_ramp_points))
    
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
    card.segments = numPoints*amp_ramp_points
    card.samples = numpy.ceil(measDur*card.sampleRate)
    card.SetParams() #warning. this may round the number of smaples to multiple of 1024

    ## Wait for settings to percolate
    time.sleep(0.5)
    
    ## SGS unit settings
    logen.set_Freq(freq_GHz)
    logen.set_Amp(lopower)
    logen.mod_Off()
    logen.power_On() 
    
    rfgen.set_Freq(freq_GHz)
    rfgen.set_Amp(rfpower)
    rfgen.mod_On()
    rfgen.power_On()
    
    trigger.output    = 'OFF'
    
    angles = numpy.linspace(angmin, angmax, angnum)
    ecc    = numpy.zeros((angnum, amp_ramp_points))
    
    for angle in range(len(angles)):
        print('Angle {} out of {}'.format(angle, len(angles)))
        if verbose:
            print('Program loading')
        hdawg.AWGs[0].stop()
        finalprog = loadprog
        finalprog = finalprog.replace('_angimp_', str(angles[angle]))
        hdawg.AWGs[0].load_program(finalprog)
        
        if verbose:
            print('Program loaded')
        
        Idata = numpy.zeros((numPoints*amp_ramp_points, card.samples))
        Qdata = numpy.zeros((numPoints*amp_ramp_points, card.samples))

        Is = numpy.zeros(numPoints)
        Qs = numpy.zeros(numPoints)
        if verbose:
            print('Acquisition started')
        trigger.output    = 'ON'
        hdawg.AWGs[0].run_loop()
        card.ArmAndWait()
        Idata, Qdata = card.ReadAllData()
        trigger.output    = 'OFF'
        if verbose:
            print('Acquisition ended')
        for amp_point in range(0,amp_ramp_points):
            for tind in range(0, numPoints):
        
                Iav = numpy.mean(Idata[tind+amp_point*numPoints])
                Qav = numpy.mean(Qdata[tind+amp_point*numPoints])
                
                Is[tind] = Iav
                Qs[tind] = Qav
    
            axes, center, phi, ecc[angle, amp_point] = uf.fit_ell_martin(Is,Qs)
            pylab.figure()
            pylab.plot(Is,Qs)
        
    return ecc

def calibrate_SGS_IQ(instruments, settings, corrections):

    logen = instruments['LO']
    rfgen = instruments['RFgen']
    card  = instruments['Digitizer'] 
    hdawg = instruments['AWG']
    trigger = instruments['Trigger']

    freq      = settings['frequency']
    lopower   = settings['lopower']
    rfpower   = settings['rfpower']
    numPoints = settings['NumPoints']
#    angle_impairment = settings['AngleImp']
#    amp_impairment   = settings['AmpImp']
    measDur   = settings['SingleShotTime']
    showFig   = settings['showFig']
    verbose   = settings['verbose']
    
    freq_GHz = freq/1e9
    #phases = numpy.linspace(0,2*numpy.pi,numPoints)
#    amp_impairment_Q = 1.
#    amp_impairment_I = 1.
#    
#    if amp_impairment>1:
#        amp_impairment_Q = 1./amp_impairment
#        if verbose:
#            print('Q impairment: {}'.format(amp_impairment_Q))
#    else:
#        amp_impairment_I = amp_impairment
#        if verbose:
#            print('I impairment: {}'.format(amp_impairment_I))
#    if verbose:
#        print('Angle impairment: {}'.format(angle_impairment))
    ## HDAWG settings
    if verbose:
        print('HDAWG config')
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    if verbose:
        print('HDAWG config ended')
    progFile = open("HDAWG_sequencer_codes/IQcorrected.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    loadprog = loadprog.replace('_npoints_', str(numPoints))
    loadprog = loadprog.replace('_time_', str(measDur))
    
    loadprog = loadprog.replace('_a_', str(corrections[0][0]))
    loadprog = loadprog.replace('_b_', str(corrections[0][1]))
    loadprog = loadprog.replace('_c_', str(corrections[1][0]))
    loadprog = loadprog.replace('_d_', str(corrections[1][1]))

    trigger.output = 'OFF'
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
    
    card.triggerDelay = 0
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
        
    
    trigger.output = 'ON'
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
     
    return axes, center, phi, ecc, Is, Qs

def calibrate_SGS_IQ_basic(instruments, settings, corrections):

    logen = instruments['LO']
    rfgen = instruments['RFgen']
    card  = instruments['Digitizer'] 
    hdawg = instruments['AWG']

    freq      = settings['frequency']
    lopower   = settings['lopower']
    rfpower   = settings['rfpower']
    numPoints = settings['NumPoints']
    measDur   = settings['SingleShotTime']
    showFig   = settings['showFig']
    verbose   = settings['verbose']
    
    freq_GHz = freq/1e9
    phases = numpy.linspace(0,2*numpy.pi,numPoints)

    ## HDAWG settings
    if verbose:
        print('HDAWG config')
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
    if verbose:
        print('HDAWG config ended')
    progFile = open("HDAWG_sequencer_codes/BasicIQ.cpp",'r')
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
    card.channelRange = 2.5
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
    rfgen.mod_On()
    rfgen.power_On()

    ## Wait for settings to percolate
    time.sleep(0.5)
    
    Idata = numpy.zeros(card.samples)
    Qdata = numpy.zeros(card.samples)
    
    Amps = numpy.zeros(numPoints)
    Angles = numpy.zeros(numPoints)
    Is = numpy.zeros(numPoints)
    Qs = numpy.zeros(numPoints)
    
    traces = numpy.zeros((2, numPoints, card.samples))
    
    init = [corrections@numpy.array([numpy.cos(th), numpy.sin(th)]) for th in phases]
    maxval = max([max(pt) for pt in init])
    cal = init/maxval
    
    for tind in range(0, numPoints):
        print('{} out of {}'.format(tind, numPoints))
        
        finalprog = loadprog
        finalprog = finalprog.replace('_Iamp_', str(cal[tind][0]))
        finalprog = finalprog.replace('_Qamp_', str(cal[tind][1]))
    
        hdawg.AWGs[0].stop()
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        
        time.sleep(0.5)
        
        card.ArmAndWait()
        Idata, Qdata = card.ReadAllData()
    
        traces[0][tind] = Idata
        traces[1][tind] = Qdata
        
        Iav = numpy.mean(Idata[0])
        Qav = numpy.mean(Qdata[0])
        
        Amp = numpy.sqrt(Iav**2 + Qav**2)
        Angle = numpy.arctan2(Iav, Qav)*180/numpy.pi
        
        Amps[tind] = Amp
        Angles[tind] = Angle
        
        Is[tind] = Iav
        Qs[tind] = Qav
        
    
    return Is, Qs