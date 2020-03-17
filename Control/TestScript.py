from HDAWG import HDAWG
from SGS import RFgen
import numpy as np

AWG=HDAWG('dev8163') #HDAWG device name
RFsource=RFgen('TCPIP0::rssgs100a110739::inst0::INSTR') #SGS visa resource name

AWG.enable_channels([0,1])
AWG.set_AWGamp([1.,1.],[0,1])

initprog="""
const NumSamples=8000;
wave w1=zeros(NumSamples);
wave w2=zeros(NumSamples);

setTrigger(0b0001);
playWave(w1,w2);
"""

AWG.load_program(initprog)