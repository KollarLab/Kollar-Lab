import pylab
import sys
import numpy
import scipy

from Instruments.HDAWG import HDAWG
from Acqiris_development.Acqiris import Acqiris
from Instruments.SGS import RFgen

#Card Settings
hardwareAddress = "PXI23::0::0::INSTR"
IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
if not IVIbinPath in sys.path:
    sys.path.append(IVIbinPath)

#card = Acqiris(hardwareAddress)
card.triggerSlope = 'Rising'
card.triggerLevel = 0.1

####HDAWG Configuration
#hdawg = HDAWG('dev8163')
#hdawg.AWGs[0].samplerate = '2.4GHz'
#hdawg.channelgrouping = '1x4'
#hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
#hdawg.Channels[1].configureChannel(marker_out='Trigger', hold='True')
#hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
#hdawg.OSCs[1].freq = 10e6
#hdawg.Channels[2].analog_outs = [0.5,0]
#hdawg.Channels[3].analog_outs = [0,1.0]
#hdawg.Channels[2].configureChannel(amp=1.0)
#hdawg.Channels[3].configureChannel(amp=2.0)

#Generators
#logen = RFgen('TCPIP0::rssgs100a110738::inst0::INSTR')
#rfgen = RFgen('TCPIP0::rssgs100a110739::inst0::INSTR')
freq_GHz = 7;
logen.set_Freq(freq_GHz)
logen.set_Amp(12)
logen.mod_Off()
logen.set_External_Reference()
logen.power_On() 

rfgen.set_Freq(freq_GHz)
rfgen.set_Amp(2)
rfgen.mod_On()
rfgen.set_External_Reference()
rfgen.power_On()

time.sleep(0.05)

## Experimental parameters
pulseamp = 0.5
pulselength = 200e-9
overflowBuffer = 10e-6
#maxtime = 110e-6
#tau = 1e-6
#T2 measurement protocol parameters
measure_points = 25
tau_max = 500e-6
tau_min = 1e-6
trig_buffer = 0e-6
taus = numpy.linspace(tau_min, tau_max, measure_points)
maxtime = tau_max + overflowBuffer

#Data acquisition parameters
card.averages = 1 #on-board averages
card.segments = 25
card.triggerDelay = trig_buffer
card.activeChannels = [1,2]
reads = 1  #reads of the card
card.verbose = False
card.sampleRate = 2e9
card.clockSource = 'External'

## Read in the sequencer program we want to use and configure it
progFile = open("HDAWG_sequencer_codes/HomodynePulse.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close
loadprog = loadprog.replace('_Amp_', str(pulseamp))
loadprog = loadprog.replace('_Time_', str(pulselength))
loadprog = loadprog.replace('_max_time_', str(maxtime))


#set up some houskeeping for the timing
digitizer_max_time = maxtime + 1e-6
digitizer_min_time = 0
digitizer_time = digitizer_max_time - digitizer_min_time

card.samples = numpy.ceil(digitizer_time*card.sampleRate)
card.SetParams() #warning. this may round the number of smaples to multiple of 1024

empiricalFudgeFactor = 0.0e-6
digitizerTimeOffset = tau_max + overflowBuffer + pulselength + empiricalFudgeFactor
cardTicks = scipy.arange(0, card.samples, 1.) # time in sample clocks from the digitizers point of view
cardXaxis = cardTicks /card.sampleRate - digitizerTimeOffset #time in seconds from the digitizers point of view


total_averages = card.averages*card.segments*reads


raw_read1 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single read
#raw_data1 = numpy.zeros((measure_points, card.samples)) #an averaged time trace for every value of tau

raw_read2 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single read
#raw_data2 = numpy.zeros((measure_points, card.samples)) #an averaged time trace for every value of tau

phaseMat1 = numpy.zeros((len(taus), card.segments))
phaseMat2 = numpy.zeros((len(taus), card.segments))

IMat1 = numpy.zeros((len(taus), card.segments))
QMat1 = numpy.zeros((len(taus), card.segments))
IMat2 = numpy.zeros((len(taus), card.segments))
QMat2 = numpy.zeros((len(taus), card.segments))

readTimes = numpy.zeros(len(taus))



pylab.figure(1)
pylab.clf()
#ax1 = pylab.subplot(1,1,1)
##ax2 = pylab.subplot(1,2,2)

ax1 = pylab.subplot(2,1,1)
ax2 = pylab.subplot(2,1,2)

waterfall = 0.1
xshift = 25
t0 = time.time()
for tind in range(0, len(taus)):
    tau = taus[tind]
    finalprog = loadprog
    finalprog = finalprog.replace('_tau_',str(tau))
    hdawg.AWGs[0].load_program(finalprog)
    hdawg.AWGs[0].run_loop()
    time.sleep(0.1)
    for rind in range(0, reads):
        card.ArmAndWait()
        data1, data2 = card.ReadAllData() #return a matrix which segments x samples
        
        currT = time.time()
        readTimes[tind] = currT - t0
        
        raw_read1 = data1
        raw_read2 = data2
        
#        #if card.segments == 1:
#        #    dataVec1 = data1
#        #    dataVec2 = data2
#        #else:
#        dataVec1 = numpy.mean(data1,0)
#        dataVec2 = numpy.mean(data2,0)
#        
#        raw_data1[tind,:] = raw_data1[tind,:] + dataVec1
#        raw_data2[tind,:] = raw_data2[tind,:] + dataVec2 
        
        
        dataVec1 = data1[0,:] #first segements is representative data vector
        dataVec2 = data2[0,:]
  
        
        #find the second pulse
        cut1 = numpy.where(cardXaxis > -pulselength/2)[0][0]
        cut2 = numpy.where(cardXaxis < pulselength/2)[0][-1] -1 
        pulse2_range = [cut1,cut2]
        
        #find the first pulse
        cut1 = numpy.where(cardXaxis > -pulselength/2 - tau)[0][0] 
        cut2 = numpy.where(cardXaxis < pulselength/2-tau)[0][-1] -1  
        pulse1_range = [cut1,cut2]
        
        for sind in range(0, card.segments):
            Is1 = data1[sind, pulse1_range[0]:pulse1_range[1]]
            Qs1 = data2[sind, pulse1_range[0]:pulse1_range[1]]
            
            Is2 = data1[sind, pulse2_range[0]:pulse2_range[1]]
            Qs2 = data2[sind, pulse2_range[0]:pulse2_range[1]]
            
            I1 = numpy.mean(Is1)
            Q1 = numpy.mean(Qs1)
            
            I2 = numpy.mean(Is2)
            Q2 = numpy.mean(Qs2)
            
            phase1 = numpy.arctan2(Q1, I1)
            phase2 = numpy.arctan2(Q2, I2)
        
            #store data
            IMat1[tind,sind] = I1
            QMat1[tind,sind] = Q1
            
            IMat2[tind,sind] = I2
            QMat2[tind,sind] = Q2
            
            phaseMat1[tind,sind] = phase1
            phaseMat2[tind,sind] = phase2
        
    
    
    #done reading. Time to plot
    pylab.sca(ax1)
    pylab.plot(cardXaxis*1e6, dataVec1 + waterfall*tind)
#    pylab.plot(cardXaxis*1e6, dataVec2)
    pylab.xlabel('Time (us)')
    pylab.ylabel('Voltage')
    
    
    pylab.sca(ax2)
    pylab.plot(cardXaxis*1e6, dataVec1)
    pylab.plot(cardXaxis*1e6, dataVec2)
    pylab.xlabel('Time (us)')
    pylab.ylabel('Voltage')
    
    pylab.show()
    
    
    
pylab.figure(2)
pylab.clf()

ax = pylab.subplot(1,3,1)
pylab.plot(readTimes, phaseMat1[:,1]*180/numpy.pi, 'r', linestyle = '', marker = '.', label = 'pulse 1')
pylab.plot(readTimes, phaseMat2[:,1]*180/numpy.pi, 'b', linestyle = '', marker = '.', label = 'pulse 2')
pylab.xlabel('~time from measurement start (S)')
pylab.ylabel('phase (degrees)')
pylab.title('absolute phases')
ax.legend(loc = 'upper right')


ax = pylab.subplot(1,3,2)
for tind in range(0, len(taus)):
    ts = taus[tind]*numpy.ones(card.segments)
    phaseDiff = 180*(phaseMat2[tind,:] - phaseMat1[tind,:])/numpy.pi
    pylab.scatter(ts*1e6, phaseDiff, c = 'b', s = 15)
pylab.xlabel('time (us)')
pylab.ylabel('phase difference between pulses (degrees)')
pylab.title('individual phase separations')
    
ax = pylab.subplot(1,3,3)
stdVec = numpy.zeros(len(taus))
for tind in range(0, len(taus)):
    phaseDiff = 180*(phaseMat2[tind,:] - phaseMat1[tind,:])/numpy.pi
    stdVec[tind] = numpy.std(phaseDiff)
pylab.scatter(taus*1e6, stdVec, c = 'b',s = 15)
pylab.xlabel('time (us)')
pylab.ylabel('std of pulse diffs (degrees)')
pylab.title('phase separations uncertainties')
pylab.show()
    
    
    
    
    