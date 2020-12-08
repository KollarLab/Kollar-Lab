# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 18:31:45 2020

@author: Kollarlab
"""
import userfuncs
import time
import os
import numpy as np
import VNAplottingTools as VNAplots
import matplotlib.pyplot as plt

saveDir = r'Z:\Data\HouckQuadTransmon\Spec\20201203'

stamp = userfuncs.timestamp()

name = 'q4specscan'

scanname = name + '_' + stamp

settings = vna.spec_default_settings()

CAV_Attenuation = 30

settings['channel'] = 1
settings['avg_time'] = 30
settings['measurement'] = 'S21'
settings['start'] = 3.5e9
settings['stop'] = 4.5e9
settings['sweep_points'] = 501
settings['RFpower'] = -25
settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1
settings['CAVpower'] = -55 + CAV_Attenuation
settings['CAVfreq'] = 8.12555e9
settings['ifBW'] = 2e2

HWattenuation = -10
numPowers = 21
startpower = -40
stoppower = 0

powers = np.linspace(startpower, stoppower, numPowers)

mags = np.zeros((len(powers), settings['sweep_points']))
phases = np.zeros((len(powers), settings['sweep_points']))

t0 = time.time()
for powerind in range(len(powers)):
    power = powers[powerind]
    print('Power: {}, final power: {}'.format(power, powers[-1]))
    settings['RFpower'] = power
    data = vna.spec_meas(settings)
    
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
