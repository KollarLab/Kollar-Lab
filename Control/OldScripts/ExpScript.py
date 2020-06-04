from Instruments.HDAWG import HDAWG
from Instruments.SGS import RFgen
import numpy
import time
import sys
import scipy
import pylab


from Acqiris_development.Acqiris import Acqiris


hardwareAddress = "PXI23::0::0::INSTR"

IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
if not IVIbinPath in sys.path:
    sys.path.append(IVIbinPath)

#card = Acqiris(hardwareAddress)
card.triggerSlope = 'Rising'
card.triggerLevel = 0.1



#RFsource=RFgen('TCPIP0::rssgs100a110739::inst0::INSTR') #SGS visa resource name

### Connect to HDAWG and initialize it 
hdawg = HDAWG('dev8163') #HDAWG device name
hdawg.AWGs[0].samplerate = '2.4GHz'
hdawg.channelgrouping = '1x4'
hdawg.Channels[0].configureChannel(amp=1.0,marker_out='Marker', hold='True')
hdawg.Channels[1].configureChannel(marker_out='Trigger', hold='True')
hdawg.AWGs[0].Triggers[0].configureTrigger(slope='rising',channel='Trigger in 1')
#hdawg.daq.setInt('/dev8163/awgs/0/outputs/0/hold',1)
#hdawg.daq.setInt('/dev8163/awgs/0/outputs/1/hold',1)
hdawg.Channels[2].configureAnalogOut([0.5,0])
hdawg.Channels[3].configureAnalogOut([0,0.5])
hdawg.Channels[2].configureChannel(amp=1.0)
hdawg.Channels[3].configureChannel(amp=1.0)

## Read in the sequencer program we want to use
progFile = open("HDAWG_sequencer_codes/TestingIQ.cpp",'r')
rawprog  = progFile.read()
loadprog = rawprog
progFile.close

## Configure program to our parameters
#pi pulse parameters
piAmp     = 1.0 #pi pulse amplitude, relative to output amplitude
piTime    = 10.0e-6 #pi pulse length (s)
piWidth   = piTime/8 #pi pulse width (gaussian) (s)


#T2 measurement protocol parameters
measure_points = 10
tau_max = 200e-6
tau_min = 50e-6
trig_buffer = 0e-6
taus = numpy.linspace(tau_min, tau_max, measure_points)
#taus = numpy.asarray([tau_min, tau_min])
#taus = numpy.logspace(min_time,max_time,num=measure_points)
meas_wait = 0.5e-6 #wait time between last pulse and measurement (s)
meas_time = 10e-6 #measurement window size (s)

#implict AWG things
max_time  = tau_max + piTime + trig_buffer #maximum separation between pulses (s)
min_time = tau_min #shortest time, checking the pulses don't overlap



##T2 measurement protocol parameters
#meas_wait = 0.5e-6 #wait time between last pulse and measurement (s)
#meas_time = 10e-6 #measurement window size (s)
#max_time  = 200e-6 #maximum separation between pulses (s)
#min_time = 10e-6 #shortest time
#measure_points = 50
#taus = numpy.linspace(min_time, max_time, measure_points)
##taus = numpy.logspace(min_time,max_time,num=measure_points)


#Data acquisition parameters
card.averages = 1 #on-board averages
card.segments = 1
card.timeDelay = trig_buffer
card.activeChannels = [1,2]
reads = 1  #reads of the card
card.verbose = False
card.sampleRate = 2e9

#card.clockSource = 'Internal'
card.clockSource = 'External'




#set up some houskeeping for the timing
digitizer_max_time = max_time + 1e-6 + meas_wait + meas_time
digitizer_min_time = 0
digitizer_time = digitizer_max_time - digitizer_min_time

card.samples = numpy.ceil(digitizer_time*card.sampleRate)
card.SetParams() #warning. this may round the number of smaples to multiple of 1024

empiricalFudgeFactor = 0.2e-6
digitizerTimeOffset = tau_max + piTime + empiricalFudgeFactor
cardTicks = scipy.arange(0, card.samples, 1.) # time in sample clocks from the digitizers point of view
cardXaxis = cardTicks /card.sampleRate - digitizerTimeOffset #time in seconds from the digitizers point of view

#ch1_waveform, ch2_waveform, markers =  hdawg.AWGs[0].read_waveform(0, channels=1, markers_present=True)
#awgTicks = scipy.arange(0, len(ch1_waveform), 1.)
#awgXaxis = awgTicks/(2.4*10-9) + trig_buffer


#implict settings
if min_time < piTime*2:
    print('BAD!!!!')
#min_time  = piTime*2

total_averages = card.averages*card.segments*reads


raw_read1 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single read
raw_data1 = numpy.zeros((measure_points, card.samples)) #an averaged time trace for every value of tau

raw_read2 = numpy.zeros((card.segments, card.samples)) #store all the raw data from a single read
raw_data2 = numpy.zeros((measure_points, card.samples)) #an averaged time trace for every value of tau




loadprog = loadprog.replace('_piAmp_', str(piAmp))
loadprog = loadprog.replace('_piTime_', str(piTime))
loadprog = loadprog.replace('_piWidth_', str(piWidth))
loadprog = loadprog.replace('_meas_wait_', str(meas_wait))
loadprog = loadprog.replace('_meas_time_', str(meas_time))
loadprog = loadprog.replace('_max_time_', str(max_time))
    
# Create a logarithmically spaced array for exponential measurements 
#measure_points = 50
##taus = numpy.logspace(min_time,max_time,num=measure_points)
#

pylab.figure(1)
pylab.clf()
#ax1 = pylab.subplot(1,1,1)
##ax2 = pylab.subplot(1,2,2)

ax1 = pylab.subplot(2,1,1)
ax2 = pylab.subplot(2,1,2)

waterfall = 0
xshift = 25
for tind in range(0, len(taus)):
    tau = taus[tind]
    finalprog = loadprog
    finalprog = finalprog.replace('_tau_',str(tau))
    hdawg.AWGs[0].load_program(finalprog)
#    hdawg.AWGs[0].run()
    hdawg.AWGs[0].run_loop()
    #if Digitizer Happy:
    #    hdawg.AWGs[0].stop()
    
    for rind in range(0, reads):
        card.ArmAndWait()
        data1, data2 = card.ReadAllData() #return a matrix which segments x samples
        
        raw_read1 = data1
        raw_read2 = data2
        
        if card.segments == 1:
            dataVec1 = data1
            dataVec2 = data2
        else:
            dataVec1 = numpy.mean(data1,0)
            dataVec2 = numpy.mean(data2,0)
        
        raw_data1[tind,:] = raw_data1[tind,:] + dataVec1
        raw_data2[tind,:] = raw_data2[tind,:] + dataVec2 
    
    
    #done reading. Time to plot
    pylab.sca(ax1)
    pylab.plot(cardXaxis*1e6 +xshift*tind, dataVec1 + waterfall*tind)
#    pylab.plot(cardXaxis*1e6, dataVec2)
    pylab.xlabel('Time (us)')
    pylab.ylabel('Voltage')
    
    
    pylab.sca(ax2)
    pylab.plot(cardXaxis*1e6, dataVec1)
    pylab.plot(cardXaxis*1e6, dataVec2/3.)
    pylab.xlabel('Time (us)')
    pylab.ylabel('Voltage')
    
    pylab.show()
    

 
#pulseInd = numpy.where(raw_data1[0,:] == numpy.max(raw_data1[0,:]))[0][0] 
#width = 20
#pulseRange = scipy.arange(pulseInd - width, pulseInd+ width, 1)
#        
    

#cut3 = numpy.where( (cardXaxis < 1e-6) and (cardXaxis > -1e-6))[0] 
#cut3 = numpy.where( -1e-6 < cardXaxis < 1e-6)[0]
    
#finalprog=loadprog
#hdawg.AWGs[0].load_program(finalprog.replace('_tau_','1e-6'))
#hdawg.AWGs[0].run()

cut1 = numpy.where(cardXaxis > -1e-6)[0]  
cut2 = numpy.where(cardXaxis < 1e-6)[0]  
limit = numpy.where(cut1 == cut2[-1])[0][0]
cut3 = cut1[0:limit]

thetas_sgs = numpy.zeros(len(taus))
thetas_hdawg = numpy.zeros(len(taus))
for tind in range(0, len(taus)):
    fft1 = numpy.fft.fft(raw_data1[tind,cut3])
    peakInd = numpy.where(numpy.abs(fft1) == numpy.max(numpy.abs(fft1)))[0][0]
    theta1 = numpy.angle(fft1[peakInd])
    thetas_sgs[tind] = theta1
    
    fft1 = numpy.fft.fft(raw_data2[tind,cut3])
    peakInd = numpy.where(numpy.abs(fft1) == numpy.max(numpy.abs(fft1)))[0][0]
    theta1 = numpy.angle(fft1[peakInd])
    thetas_hdawg[tind] = theta1


pylab.figure(2)
pylab.clf()

ax = pylab.subplot(1,1,1)
#pylab.plot(cut3, raw_data1[0,cut3])
pylab.plot(thetas_sgs*180/numpy.pi)
offset =( -thetas_hdawg[0] + thetas_sgs[0])*180/numpy.pi
pylab.plot(thetas_hdawg*180/numpy.pi + offset)

pylab.xlabel('run number - 1')
pylab.ylabel('Ch1 v Ch2 phase difference (degrees)')
pylab.show()












