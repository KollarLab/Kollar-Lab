# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 12:24:11 2020

@author: Kollarlab
"""


import numpy as np

import userfuncs
from VNAplottingTools import power_plot

saveDir = r'Z:\Data\HouckTaTransmon\20201110'

stamp = userfuncs.timestamp()

name = 'cavitytrans_qbitpowersweep'

scanname = name + '_' + stamp

qfreqs = [4.94247e9, 5.126e9]

settings = vna.trans_default_settings()

settings['averages'] = 30
settings['start'] = 7.165e9
settings['stop'] = 7.177e9
settings['sweep_points'] = 2001
settings['ifBW'] = 1e3
settings['RFpower'] = -20

logen.Output = 'Off'  
logen.Freq = qfreqs[0]

powers = numpy.linspace(-15, 0, 30)

mags = []
phases = []

for power in powers:
    logen.Output = 'Off'
    logen.Power = power
    logen.Output = 'On'
    m, p, f = vna.trans_meas(settings)
    mags.append(m)
    phases.append(p)
    
power_plot(f, mags, phases, powers, scanname+'_q0', HWattenuation = -10)

userfuncs.SaveFull(saveDir, scanname+'_q0', ['mags', 'phases', 'freqs', 'powers', 'HWattenuation'], locals(), expsettings=settings)
plt.savefig(os.path.join(saveDir, scanname+'_q0.png'), dpi = 150)

logen.Output = 'Off'
logen.Freq = qfreqs[1]

powers = numpy.linspace(-15, 0, 30)

mags1 = []
phases1 = []

for power in powers:
    logen.Output = 'Off'
    logen.Power = power
    logen.Output = 'On'
    m, p, f = vna.trans_meas(settings)
    mags1.append(m)
    phases1.append(p)
    
power_plot(f, mags1, phases1, powers, scanname+'_q1', HWattenuation = -10)

userfuncs.SaveFull(saveDir, scanname+'_q1', ['mags', 'phases', 'freqs', 'powers', 'HWattenuation'], locals(), expsettings=settings)
plt.savefig(os.path.join(saveDir, scanname+'_q1.png'), dpi = 150)