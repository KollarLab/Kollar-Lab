# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 14:07:57 2020

@author: Kollarlab
"""
import userfuncs
import numpy as np
import time
import os
import VNAplottingTools as VNAplots
import matplotlib.pyplot as plt

saveDir = r'Z:\Data\HouckTaTransmon\Trans\20201120'

stamp = userfuncs.timestamp()

name = 'cavityscan'

scanname = name + '_' + stamp

settings = vna.trans_default_settings()

settings['channel'] = 1
settings['avg_time'] = 30
settings['measurement'] = 'S21'
settings['start'] = 7.165e9
settings['stop'] = 7.174e9
settings['sweep_points'] = 1001
settings['RFpower'] = -20
settings['ifBW'] = 1e3

HWattenuation = -30
numPowers = 10
startpower = -18
stoppower = -18.1

powers = np.linspace(startpower, stoppower, numPowers)

mags = np.zeros((len(powers), settings['sweep_points']))
phases = np.zeros((len(powers), settings['sweep_points']))

t0 = time.time()
for powerind in range(len(powers)):
    power = powers[powerind]
    print('Power: {}, final power: {}'.format(power, powers[-1]))
    settings['RFpower'] = power
    
    data = vna.trans_meas(settings)
    
    vna.autoscale()
    
    mags[powerind] = data['mag']
    phases[powerind] = data['phase']
    
    if powerind==0:
        t1=time.time()
        tdiff = t1-t0
        ttotal = tdiff*len(powers)
        print('Single run time: {}, estimated total time: {}'.format(tdiff, ttotal))

t2 = time.time()
print('Elapsed time: {}'.format(t2-t0))

freqs = data['xaxis']   

VNAplots.power_plot(freqs, mags, phases, powers, scanname, HWattenuation)

userfuncs.SaveFull(saveDir, scanname, ['mags', 'phases', 'freqs', 'powers', 'HWattenuation'], locals(), expsettings=settings)
plt.savefig(os.path.join(saveDir, scanname+'.png'), dpi = 150)