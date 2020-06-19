import pylab
import sys
import numpy
import scipy
import time
from datetime import datetime

from mplcursors import cursor as datacursor
from SGShelper import SGS_coupling, HDAWG_clock
import userfuncs as uf

def GetDefaultSettings():
    settings = {}
    settings['lopower']       = 12
    settings['rfpower']       = 0
    settings['one_shot_time'] = 1e-6
    settings['carrier_freq']  = 8e9
    settings['trigger_rate']  = 500
    #Pulse settings
    settings['pulse_width']   = 200e-9
    settings['pulse_amp']     = 0.5
    settings['ramp_frac']     = 0.1
    settings['IQangle1']      = 0
    settings['IQangle2']      = 0
    #Tau range settings
    settings['tau_min']       = 1e-6
    settings['tau_max']       = 1000e-6
    settings['num_points']    = 25
    settings['spacing']       = 'Linear'
    #Card settings
    settings['segments']         = 100
    settings['reads']            = 2
    settings['averages']         = 1
    settings['activeChannels']   = [1,2]
    settings['sampleRate']       = 2e9/8
    settings['trigger_buffer']   = 0e-6
    settings['digitizer_buffer'] = 10e-6
    #Misc settings
    settings['savepath']      = r'C:\Users\Kollarlab\Desktop'
    settings['save']          = True

    return settings

def PulsedHomodyneStability(instruments, settings):

    ## Instruments used
    hdawg = instruments['AWG']
    logen = instruments['LO']
    rfgen = instruments['RFgen']
    card  = instruments['Digitizer']

    ## General reference configuration
    carrier_freq       = settings['carrier_freq']
    rfpower            = settings['rfpower']
    lopower            = settings['lopower']             

    ## Misc settings
    savepath = settings['savepath']
    periodSteps = 30  #how many periods to compare

    ##Experimental parameters
    #Pulse settings
    pulse_width  = settings['pulse_width']
    pulse_amp    = settings['pulse_amp']
    ramp_frac    = settings['ramp_frac']
    IQangle1     = settings['IQangle1']
    IQangle2     = settings['IQangle2']
    #Tau range settings
    tau_min      = settings['tau_min']
    tau_max      = settings['tau_max']
    num_points   = settings['num_points']
    spacing      = settings['spacing']
    P2_time      = tau_max + 3e-6

    ## Create array of taus and check validity of user provided ranges
    T_min = pulse_width*(1+ramp_frac)
    T_max = P2_time-T_min/2

    if tau_min < T_min:
        print('tau_min is too small, using {}'.format(T_min))
        tau_min = tau_min
    if tau_max > T_max:
        print('tau_max is too large, using {}'.format(T_max))

    if spacing == 'Linear':
        taus = numpy.linspace(tau_min, tau_max, num_points)
    if spacing == 'Log':
        taus = numpy.logspace(tau_min, tau_max, num_points)

    ## HDAWG settings
    hdawg.AWGs[0].samplerate = '2.4GHz'
    hdawg.channelgrouping = '1x4'
    hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.Channels[1].configureChannel(amp=1.0,marker_out='Marker', hold='False')
    hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')

    ## Generator settings
    freq_GHz = carrier_freq/1e9
    logen.set_Freq(freq_GHz)
    logen.set_Amp(lopower)
    logen.mod_Off()

    rfgen.set_Freq(freq_GHz)
    rfgen.set_Amp(rfpower)
    rfgen.mod_On()

    logen.power_On() 
    rfgen.power_On()

    ## Digitizer settings
    card.clockSource  = 'External'
    card.triggerSlope = 'Rising'
    card.triggerLevel = 0.1
    card.verbose      = False

    reads               = settings['reads']
    digitizer_buffer    = settings['digitizer_buffer']
    card.averages       = settings['averages']
    card.segments       = settings['segments']
    card.sampleRate     = settings['sampleRate']
    card.activeChannels = settings['activeChannels']
    card.triggerDelay   = settings['trigger_buffer']

    ## Read in the sequencer program we want to use and configure it
    progFile = open("HDAWG_sequencer_codes/PulsePair.cpp",'r')
    rawprog  = progFile.read()
    loadprog = rawprog
    progFile.close()
    loadprog = loadprog.replace('_Amp_', str(pulse_amp))
    loadprog = loadprog.replace('_Time_', str(pulse_width))
    loadprog = loadprog.replace('_frac_', str(ramp_frac))
    loadprog = loadprog.replace('_max_time_', str(P2_time))
    loadprog = loadprog.replace('_IQangle1_', str(IQangle1))
    loadprog = loadprog.replace('_IQangle2_', str(IQangle2))

    #set up some houskeeping for the timing
    digitizer_max_time = P2_time + pulse_width*(1+ramp_frac) + digitizer_buffer
    digitizer_min_time = 0
    digitizer_time     = digitizer_max_time - digitizer_min_time

    card.samples = numpy.ceil(digitizer_time*card.sampleRate)
    card.SetParams() #warning. this may round the number of smaples to multiple of 1024

    empiricalFudgeFactor = 0.137e-6   #this is a very exact number to be off by!!!!!!
    digitizerTimeOffset  = P2_time + empiricalFudgeFactor
    cardTicks = scipy.arange(0, card.samples, 1.) # time in sample clocks from the digitizers point of view
    cardXaxis = cardTicks /card.sampleRate - digitizerTimeOffset #time in seconds from the digitizers point of view

    ################################################
    #Array Initialization
    ################################################

    raw_read1 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single readS
    raw_read2 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single read

    phaseMat1 = numpy.zeros((len(taus), reads, card.segments))
    phaseMat2 = numpy.zeros((len(taus), reads, card.segments))

    IMat1 = numpy.zeros((len(taus), reads,  card.segments))
    QMat1 = numpy.zeros((len(taus), reads, card.segments))
    IMat2 = numpy.zeros((len(taus), reads, card.segments))
    QMat2 = numpy.zeros((len(taus), reads, card.segments))


    readTimes = numpy.zeros(len(taus))

    ################################################
    #Data Collection
    ################################################


    t0 = time.time()
    for tind in range(0, len(taus)):
        tau = taus[tind]
        print('Tau: {}'.format(tau))
        finalprog = loadprog
        finalprog = finalprog.replace('_tau_',str(tau))
        hdawg.AWGs[0].load_program(finalprog)
        hdawg.AWGs[0].run_loop()
        time.sleep(0.1)

        ## Find the pulses on the x-axis
        #find the second pulse
        cut1 = numpy.where(cardXaxis > -pulse_width/2)[0][0]
        cut2 = numpy.where(cardXaxis < pulse_width/2)[0][-1] -1 
        pulse2_range = [cut1,cut2]

        #find the first pulse
        cut1 = numpy.where(cardXaxis > -pulse_width/2 - tau)[0][0] 
        cut2 = numpy.where(cardXaxis < pulse_width/2-tau)[0][-1] -1  
        pulse1_range = [cut1,cut2]

        for rind in range(0, reads):
            card.ArmAndWait()
            data1, data2 = card.ReadAllData() #return a matrix which segments x samples

            currT = time.time()
            readTimes[tind] = currT - t0

            raw_read1 = data1
            raw_read2 = data2

            dataVec1 = data1[0,:] #first segments is representative data vector
            dataVec2 = data2[0,:]

            for sind in range(0, card.segments):

                Is1 = data1[sind, pulse1_range[0]:pulse1_range[1]]
                Qs1 = data2[sind, pulse1_range[0]:pulse1_range[1]]

                Is2 = data1[sind, pulse2_range[0]:pulse2_range[1]]
                Qs2 = data2[sind, pulse2_range[0]:pulse2_range[1]]

                I1 = numpy.mean(Is1)
                Q1 = numpy.mean(Qs1)

                I2 = numpy.mean(Is2)
                Q2 = numpy.mean(Qs2)

                phase1 = numpy.arctan2(Q1, I1)*180/numpy.pi
                phase2 = numpy.arctan2(Q2, I2)*180/numpy.pi

                #store data
                IMat1[tind,rind, sind] = I1
                QMat1[tind,rind, sind] = Q1

                IMat2[tind,rind, sind] = I2
                QMat2[tind,rind, sind] = Q2

                phaseMat1[tind,rind, sind] = phase1
                phaseMat2[tind,rind, sind] = phase2

    ################################################
    #Data Processing
    ################################################

    #Look at phase difference vs tau
    mean_v_tau = numpy.zeros(len(taus))
    std_v_tau = numpy.zeros(len(taus))
    for tind in range(0, len(taus)):
        phaseDiff = phaseMat2[tind,:,:] - phaseMat1[tind,:,:]
        std_v_tau[tind] = numpy.std(phaseDiff)   
        mean_v_tau[tind] = numpy.mean(phaseDiff)

    #Look at phase difference across repetitions (triggers)
    phaseDiff_v_periods = numpy.zeros(( len(taus), reads, periodSteps, card.segments-periodSteps))
    for tind in range(0, len(taus)):
        #step thtrough the tau values
        for rind in range(0,reads):
            #step through the different card reads
            for diff in range(0,periodSteps):
                phaseDiff = phaseMat2[tind,rind, 0:-periodSteps] - phaseMat2[tind,rind, diff:(-periodSteps+diff)]
                phaseDiff_v_periods[tind, rind, diff,:] = phaseDiff

    trigPeriod = 1./settings['trigger_rate']
    period_ints = scipy.arange(0, periodSteps,1.)*trigPeriod
    std_v_period = numpy.zeros(periodSteps)
    for diff in range(0,periodSteps):
        std_v_period[diff] = numpy.std(phaseDiff_v_periods[:,:,diff,:])

    angle_std_vs_time_fig = pylab.figure(2)
    pylab.clf()

    #Phase scatter vs tau
    ax = pylab.subplot(2,3,1)
    for tind in range(0, len(taus)):
        ts = taus[tind]*numpy.ones(card.segments*reads)
        phaseDiff = phaseMat2[tind,:,:] - phaseMat1[tind,:,:]
        pylab.scatter(ts*1e6, phaseDiff, c = 'dodgerblue', s = 15)
    pylab.xlabel('time (us)')
    pylab.ylabel('phase difference between pulses (degrees)')
    pylab.title('individual phase separations \n single read')

    #Phase scatter vs reps 
    ax = pylab.subplot(2,3,2)
    for tind in range(0, len(taus)):
        ts = taus[tind]*numpy.ones(card.segments*reads)
        phaseDiff = phaseMat2[tind,:,:] - phaseMat1[tind,:,:]
        pylab.scatter(ts*1e3, phaseDiff, c = 'deepskyblue', s = 7)

    for diff in range(0,periodSteps):
        vals =  phaseDiff_v_periods[:,:,diff,:]
        ts = period_ints[diff]*numpy.ones(vals.size)
        pylab.scatter(ts*1e3, vals, c = 'mediumblue', s = 15)
    pylab.xlabel('time (ms)')
    pylab.ylabel('phase difference between pulses (degrees)')
    pylab.title('individual separations')

    #Pulse phases vs absolute time
    ax = pylab.subplot(2,3,3)
    pylab.plot(readTimes, phaseMat1[:,1,1], 'r', linestyle = '', marker = '.', label = 'pulse 1')
    pylab.plot(readTimes, phaseMat2[:,1,1], 'b', linestyle = '', marker = '.', label = 'pulse 2')
    pylab.xlabel('~time from measurement start (S)')
    pylab.ylabel('phase (degrees)')
    pylab.title('absolute phases of first segments')
    ax.legend(loc = 'upper right')

    #Phase stds vs tau
    ax = pylab.subplot(2,3,4)
    pylab.scatter(taus*1e6,std_v_tau, c = 'deepskyblue',s = 15)
    pylab.xlabel('time (us)')
    pylab.ylabel('std of pulse diffs (degrees) \n single read')
    pylab.title('phase separations uncertainties')

    #Phase stds vs rep rate
    ax = pylab.subplot(2,3,5)
    pylab.scatter(taus*1e3,std_v_tau, c = 'deepskyblue',s = 7)
    pylab.scatter(period_ints[1:]*1e3,std_v_period[1:], c = 'mediumblue',s = 15)
    pylab.xlabel('time (ms)')
    pylab.ylabel('std of pulse diffs (degrees)')
    pylab.title('phase separations uncertainties')

    pylab.tight_layout()
    pylab.show()

    ##Histograms of phase differences (vs tau and longer time scales)
    histograms_phase_separation = pylab.figure(3)
    pylab.clf()

    #long time phase separations 
    temp = phaseDiff_v_periods[:, :, 6,:]
    temp2= temp.reshape(temp.shape[0]*temp.shape[1]*temp.shape[2])
    dev= numpy.std(temp)
    dev2 = numpy.std(temp2)

    #short time phase separations
    temp3 = phaseMat2[7,:,:] - phaseMat1[7,:,:]
    temp4 = temp3.reshape(temp3.shape[0]*temp3.shape[1])
    dev3 = numpy.std(temp4)

    longHist, longBins = numpy.histogram(temp2, bins = 50, range = (-10,10))
    shortHist, shortBins = numpy.histogram(temp4, bins = 50, range = (-10,10))

    plotBins_long = (longBins[0:-1] + longBins[1:])/2
    binWidth_long = longBins[1]- longBins[0]

    plotBins_short = (shortBins[0:-1] + shortBins[1:])/2
    binWidth_short = shortBins[1]- shortBins[0]

    longHist = longHist/(len(taus)*reads*card.segments)
    shortHist = shortHist/(reads*card.segments)

    ax = pylab.subplot(1,1,1)
    pylab.bar(plotBins_long, longHist, binWidth_long, color = 'k', facecolor = 'mediumblue', alpha = 0.8, label = 'rep rate')
    pylab.bar(plotBins_short, shortHist, binWidth_short, color = 'k', facecolor = 'deepskyblue', alpha = 0.5, label = 'tau')
    pylab.xlabel('phase diff (degrees)')
    pylab.ylabel('density of counts')
    pylab.title('raw histograms')
    ax.legend(loc = 'upper left')

    pylab.show()

    #Raw trace for debugging
    raw_trace_figure = pylab.figure(4)
    pylab.clf()
    ax = pylab.subplot(1,1,1)
    pylab.plot(cardXaxis*1e6,dataVec1 )
    pylab.plot(cardXaxis*1e6,dataVec2 )

    t0 = cardXaxis[pulse2_range[0]]
    t1 = cardXaxis[pulse2_range[1]]
    pylab.plot([t0*1e6,t0*1e6], [-0.05,0.05])
    pylab.plot([t1*1e6,t1*1e6], [-0.06,0.06])

    t0 = cardXaxis[pulse1_range[0]]
    t1 = cardXaxis[pulse1_range[1]]
    pylab.plot([t0*1e6,t0*1e6], [-0.05,0.05])
    pylab.plot([t1*1e6,t1*1e6], [-0.06,0.06])

    pylab.xlabel('time (us)')
    pylab.ylabel('voltage')
    pylab.title('Single raw trace')
    pylab.show()

    rfgen.power_Off()
    logen.power_Off()

    dataTosave = ['IMat1', 'IMat2', 'QMat1', 'QMat2', 'phaseMat1', 'phaseMat2', 'dataVec1', 'dataVec2',
        'cardXaxis', 'pulse1_range', 'pulse2_range', 'readTimes', 'mean_v_tau', 'std_v_tau', 'std_v_period',
        'phaseDiff_v_periods', 'taus']

    date  = datetime.now()
    stamp = date.strftime('%Y%m%d_%H%M%S')
    filename = 'PulsedHomodyne{}'.format(stamp)

    figsTosave = [angle_std_vs_time_fig, histograms_phase_separation, raw_trace_figure]
    if settings['save']:
        uf.SaveFull(savepath, filename, dataTosave, locals(), settings, instruments, figsTosave)
