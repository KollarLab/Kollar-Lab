# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 14:07:57 2020

@author: Kollarlab
"""
import time
import os
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
import VNAplottingTools as VNAplots

project_dir = r'Z:\Data\HouckQuadTransmon'
meas_type = 'Trans'
save_Dir = userfuncs.saveDir(project_dir, meas_type)

name = 'powersweep'
stamp = userfuncs.timestamp()
scanname = name + '_' + stamp

CAV_Attenuation = -30

start_power = -40
stop_power = 0
power_points = 21

settings = vna.trans_default_settings()

settings['channel'] = 1
settings['avg_time'] = 10
settings['measurement'] = 'S21'
settings['start'] = 8.05e9
settings['stop'] = 8.15e9
settings['sweep_points'] = 2001
settings['RFpower'] = -30
settings['ifBW'] = 1e3

powers = np.linspace(start_power, stop_power, power_points)

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

VNAplots.power_plot(freqs, mags, phases, powers, scanname, CAV_Attenuation)

userfuncs.SaveFull(save_Dir, scanname, ['mags', 'phases', 'freqs', 'powers', 'CAV_Attenuation'], locals(), expsettings=settings)
plt.savefig(os.path.join(save_Dir, scanname+'.png'), dpi = 150)