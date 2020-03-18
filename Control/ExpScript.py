from Instruments.HDAWG import HDAWG
from Instruments.SGS import RFgen
import numpy as np
import time

AWG=HDAWG('dev8163') #HDAWG device name
#RFsource=RFgen('TCPIP0::rssgs100a110739::inst0::INSTR') #SGS visa resource name

AWG.enable_channels([0,1])
AWG.set_AWGamp([1.,1.],[0,1])

#Most basic program for AWG
initprog="""
const NumSamples=_samplesparam_;
wave w1=ones(NumSamples);
wave w2=zeros(NumSamples);

for(var i=1; i<5;i=i+1){
playWave(w1,w2);
waitWave();
wait(1000*i); //Sets the amount of time between the two pulses (is multiples of sequencer clock, around 3.3 ns)
playWave(w2,w1);
waitWave();
wait(1000);
}
"""
#Program using markers
markerprog="""
const marker_pos = 1000;
wave w_gauss = gauss(8000, 4000, 1000);
wave w_left = marker(marker_pos, 0);
wave w_right = marker(8000-marker_pos, 1);
wave w_marker = join(w_left, w_right);
wave w_gauss_marker = w_gauss + w_marker;
playWave(w_gauss_marker,w_gauss_marker);
"""

NumSamples=800 #sets number of samples to use for the waveform (try to make it a multiple of 16 otherwise AWG pads with zeros)
x=np.linspace(-np.pi,np.pi, NumSamples)
Isignal=np.sin(x)
Qsignal=np.cos(x)

loadprog=initprog.replace('_samplesparam_', str(NumSamples))
AWG.load_program(loadprog)
#AWG.load_waveform(0,Isignal,Qsignal)
t1=time.time()

AWG.AWG_run()