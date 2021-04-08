# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 15:47:20 2020

@author: Kollarlab
"""
import userfuncs
import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

saveDir = r'Z:\Data\HouckTaTransmon\20201109'

stamp = userfuncs.timestamp()

name = 'cavityscan'

scanname = name + '_' + stamp

settings = vna.trans_default_settings()

#powers = np.linspace(-35,10,16)
numPowers= 20
powers = np.linspace(-45,10,numPowers)
averages = [int(x) for x in np.linspace(25,3,numPowers)]

settings['start'] = 7.165e9
settings['stop'] = 7.174e9
settings['sweep_points'] = 2001

HWattenuation = 0

mags = []
phases = []

for power, average in zip(powers, averages):
    print('Power: {}, final power: {}'.format(power, powers[-1]))
    settings['RFpower'] = power
    settings['averages'] = average
    
    mag, phase, freqs = vna.trans_meas(settings)
    
    mags.append(mag)
    phases.append(phase)
    
fig2 = plt.figure(8)
plt.clf()
ax = plt.subplot(1,1,1)
mags2 = numpy.asarray(mags)
freqs2 = np.zeros(len(freqs) + 1)
fdiff = freqs[1]-freqs[0]
powers2 = np.zeros(len(powers)+1)
pdiff = powers[1] - powers[0]

freqs2[0:-1] = freqs-fdiff/2
freqs2[-1] = freqs[-1] + fdiff/2

powers2[0:-1] = powers-pdiff/2
powers2[-1] = powers[-1] + pdiff/2

XX2,YY2 = np.meshgrid(freqs2,powers2+HWattenuation)
plt.pcolormesh(XX2/1e9,YY2,mags2, shading = 'nearest')



#ax.set_yaxis_ticks(powers, powers)

plt.xlabel("Frequency (GHz)")
plt.ylabel(r"Power (dBm)")
plt.title('RF power scan, {}'.format(scanname))

power_plot(freqs, mags, phases, powers, scanname, HWattenuation)

plt.show()

userfuncs.SaveFull(saveDir, scanname, ['mags', 'phases', 'freqs', 'powers', 'HWattenuation'], locals(), expsettings=settings)
plt.savefig(os.path.join(saveDir, scanname+'.png'), dpi = 150)







