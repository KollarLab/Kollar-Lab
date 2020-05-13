import pylab
import sys
import numpy
import scipy
import time

from mplcursors import cursor as datacursor

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
measure_points = 10
tau_max = 1000e-6
tau_min = 1e-6
trig_buffer = 0e-6
taus = numpy.linspace(tau_min, tau_max, measure_points)
maxtime = tau_max + overflowBuffer

#Data acquisition parameters
card.averages = 1 #on-board averages
card.segments = 25
card.triggerDelay = trig_buffer
card.activeChannels = [1,2]
reads = 2  #reads of the card
card.verbose = False
card.sampleRate = 2e9/8
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

empiricalFudgeFactor = -0.096e-6   #this is a very exact number to be off by!!!!!!
digitizerTimeOffset = tau_max + overflowBuffer + pulselength + empiricalFudgeFactor
cardTicks = scipy.arange(0, card.samples, 1.) # time in sample clocks from the digitizers point of view
cardXaxis = cardTicks /card.sampleRate - digitizerTimeOffset #time in seconds from the digitizers point of view


total_averages = card.averages*card.segments*reads


raw_read1 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single read
#raw_data1 = numpy.zeros((measure_points, card.samples)) #an averaged time trace for every value of tau

raw_read2 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single read
#raw_data2 = numpy.zeros((measure_points, card.samples)) #an averaged time trace for every value of tau


phaseMat1 = numpy.zeros((len(taus), reads, card.segments))
phaseMat2 = numpy.zeros((len(taus), reads, card.segments))

IMat1 = numpy.zeros((len(taus), reads,  card.segments))
QMat1 = numpy.zeros((len(taus), reads, card.segments))
IMat2 = numpy.zeros((len(taus), reads, card.segments))
QMat2 = numpy.zeros((len(taus), reads, card.segments))


readTimes = numpy.zeros(len(taus))



#pylab.figure(1)
#pylab.clf()
##ax1 = pylab.subplot(1,1,1)
###ax2 = pylab.subplot(1,2,2)
#
#ax1 = pylab.subplot(2,1,1)
#ax2 = pylab.subplot(2,1,2)

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
#            dind = card.segments*rind + sind
            
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
        
    
    
#    #done reading. Time to plot
#    pylab.sca(ax1)
#    pylab.plot(cardXaxis*1e6, dataVec1 + waterfall*tind)
##    pylab.plot(cardXaxis*1e6, dataVec2)
#    pylab.xlabel('Time (us)')
#    pylab.ylabel('Voltage')
#    
#    
#    pylab.sca(ax2)
#    pylab.plot(cardXaxis*1e6, dataVec1)
#    pylab.plot(cardXaxis*1e6, dataVec2)
#    pylab.xlabel('Time (us)')
#    pylab.ylabel('Voltage')
#    
#    pylab.show()
            
#process the data
##############
            
#first versus tau
mean_v_tau = numpy.zeros(len(taus))
std_v_tau = numpy.zeros(len(taus))
for tind in range(0, len(taus)):
    phaseDiff = phaseMat2[tind,:,:] - phaseMat1[tind,:,:]
    std_v_tau[tind] = numpy.std(phaseDiff)   
    mean_v_tau[tind] = numpy.mean(phaseDiff)
            
#next versus rep rate multiples
periodSteps = 20           

phaseDiff_v_periods = numpy.zeros(( len(taus), reads, periodSteps, card.segments-periodSteps))
for tind in range(0, len(taus)):
    #step thtrough the tau values
    for rind in range(0,reads):
        #step through the different card reads
        for diff in range(0,periodSteps):
            phaseDiff = phaseMat2[tind,rind, 0:-periodSteps] - phaseMat2[tind,rind, diff:(-periodSteps+diff)]
            
#            noise = numpy.std(phaseDiff)
            phaseDiff_v_periods[tind, rind, diff,:] = phaseDiff

trigPeriod = 1/500. ####!!!!!! eventually this should be determined adaptively, not hard coded
period_ints = scipy.arange(0, periodSteps,1.)*trigPeriod
std_v_period = numpy.zeros(periodSteps)
for diff in range(0,periodSteps):
    std_v_period[diff] = numpy.std(phaseDiff_v_periods[:,:,diff,:])
    
    
    
    
pylab.figure(2)
pylab.clf()

ax = pylab.subplot(2,3,1)
for tind in range(0, len(taus)):
    ts = taus[tind]*numpy.ones(card.segments*reads)
    phaseDiff = phaseMat2[tind,:,:] - phaseMat1[tind,:,:]
    pylab.scatter(ts*1e6, phaseDiff, c = 'dodgerblue', s = 15)
pylab.xlabel('time (us)')
pylab.ylabel('phase difference between pulses (degrees)')
pylab.title('individual phase separations \n single read')
    

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


ax = pylab.subplot(2,3,3)
pylab.plot(readTimes, phaseMat1[:,1,1], 'r', linestyle = '', marker = '.', label = 'pulse 1')
pylab.plot(readTimes, phaseMat2[:,1,1], 'b', linestyle = '', marker = '.', label = 'pulse 2')
pylab.xlabel('~time from measurement start (S)')
pylab.ylabel('phase (degrees)')
pylab.title('absolute phases of first segments')
ax.legend(loc = 'upper right')


ax = pylab.subplot(2,3,4)
pylab.scatter(taus*1e6,std_v_tau, c = 'deepskyblue',s = 15)
pylab.xlabel('time (us)')
pylab.ylabel('std of pulse diffs (degrees) \n single read')
pylab.title('phase separations uncertainties')


ax = pylab.subplot(2,3,5)
pylab.scatter(taus*1e3,std_v_tau, c = 'deepskyblue',s = 7)
pylab.scatter(period_ints[1:]*1e3,std_v_period[1:], c = 'mediumblue',s = 15)
pylab.xlabel('time (ms)')
pylab.ylabel('std of pulse diffs (degrees)')
pylab.title('phase separations uncertainties')

 
#ax = pylab.subplot(2,3,6)
##stdVec = numpy.zeros(len(taus))
##for tind in range(0, len(taus)):
##    phaseDiff = 180*(phaseMat2[tind,:,:] - phaseMat1[tind,:,:])/numpy.pi
##    stdVec[tind] = numpy.std(phaseDiff)
#pylab.scatter(period_ints[1:]*1e3,std_v_period[1:], c = 'b',s = 15)
#pylab.xlabel('time (ms)')
#pylab.ylabel('std of pulse diffs (degrees)')
#pylab.title('phase separations uncertainties')


pylab.tight_layout()
pylab.show()
        
    
    
    
    
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

pylab.figure(3)
pylab.clf()

ax = pylab.subplot(1,1,1)
pylab.bar(plotBins_long, longHist, binWidth_long, color = 'k', facecolor = 'mediumblue', alpha = 0.8, label = 'rep rate')
pylab.bar(plotBins_short, shortHist, binWidth_short, color = 'k', facecolor = 'deepskyblue', alpha = 0.5, label = 'tau')
pylab.xlabel('phase diff (degrees)')
pylab.ylabel('density of counts')
pylab.title('raw histograms')
ax.legend(loc = 'upper left')


#ax = pylab.subplot(1,2,2)
#pylab.bar(plotBins_long, longHist*25, binWidth_long, color = 'k', facecolor = 'mediumblue', alpha = 0.8, label = 'rep rate')
#pylab.bar(plotBins_short, shortHist, binWidth_short, color = 'k', facecolor = 'deepskyblue', alpha = 0.5, label = 'tau')
#pylab.xlabel('phase diff (degrees)')
#pylab.ylabel('scaled counts')
#pylab.title('amplitudes scaled')
#ax.set_ylim([0, numpy.max(shortHist)*1.2])
#ax.legend(loc = 'upper left')

pylab.show()

    




pylab.figure(4)
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









##plot testing
#from mpldatacursor import datacursor
#
#pylab.figure(5)
#pylab.clf()
#ax = pylab.subplot(1,1,1)
#pylab.plot(cardXaxis*1e6,dataVec1 )
#pylab.plot(cardXaxis*1e6,dataVec2 )
#
#t0 = cardXaxis[pulse2_range[0]]
#t1 = cardXaxis[pulse2_range[1]]
#pylab.plot([t0*1e6,t0*1e6], [-0.05,0.05])
#pylab.plot([t1*1e6,t1*1e6], [-0.06,0.06])
#
#t0 = cardXaxis[pulse1_range[0]]
#t1 = cardXaxis[pulse1_range[1]]
#pylab.plot([t0*1e6,t0*1e6], [-0.05,0.05])
#pylab.plot([t1*1e6,t1*1e6], [-0.06,0.06])
#
#ax.set_xlim([-0.5,0.5])
#
#pylab.xlabel('time (us)')
#pylab.ylabel('voltage')
#pylab.title('Single raw trace')
#
#datacursor()
#pylab.show()



































    
#pylab.hist(temp2, bins = 50, range = (-10,10), color = 'b', alpha = 0.5)
#pylab.hist(temp4, bins = 50, range = (-10,10), color = 'r', alpha = 0.5)



#pylab.figure(2)
#pylab.clf()
#
#ax = pylab.subplot(1,3,1)
#pylab.plot(readTimes, phaseMat1[:,1,1]*180/numpy.pi, 'r', linestyle = '', marker = '.', label = 'pulse 1')
#pylab.plot(readTimes, phaseMat2[:,1,1]*180/numpy.pi, 'b', linestyle = '', marker = '.', label = 'pulse 2')
#pylab.xlabel('~time from measurement start (S)')
#pylab.ylabel('phase (degrees)')
#pylab.title('absolute phases')
#ax.legend(loc = 'upper right')
#
#
#ax = pylab.subplot(1,3,2)
#for tind in range(0, len(taus)):
#    ts = taus[tind]*numpy.ones(card.segments*reads)
#    phaseDiff = 180*(phaseMat2[tind,:,:] - phaseMat1[tind,:,:])/numpy.pi
#    pylab.scatter(ts*1e6, phaseDiff, c = 'b', s = 15)
#pylab.xlabel('time (us)')
#pylab.ylabel('phase difference between pulses (degrees)')
#pylab.title('individual phase separations')
#    
#ax = pylab.subplot(1,3,3)
##stdVec = numpy.zeros(len(taus))
##for tind in range(0, len(taus)):
##    phaseDiff = 180*(phaseMat2[tind,:,:] - phaseMat1[tind,:,:])/numpy.pi
##    stdVec[tind] = numpy.std(phaseDiff)
#pylab.scatter(taus*1e6,std_v_tau, c = 'b',s = 15)
#pylab.xlabel('time (us)')
#pylab.ylabel('std of pulse diffs (degrees)')
#pylab.title('phase separations uncertainties')
#pylab.tight_layout()
#pylab.show()
#    
    
    
    
    