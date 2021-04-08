# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 15:51:31 2021

@author: Kollarlab
"""

import matplotlib.pyplot as plt

import userfuncs
#import VNAplottingTools as VNAplots
import  plotting_tools as plots




#dat = userfuncs.LoadFull(r'Z:\Data\HouckDualHangerFluxonium\Spec\20210119\flxon_powersweep_20210119_150244.pkl')
dat = userfuncs.LoadFull(r'Z:\Data\HouckDualHangerFluxonium\Spec\20210119\flxon_powersweep_avgingTest_20210119_171748.pkl')

data = dat[0]
expsettings = dat[1]

freqs = data['freqs']
mags = data['mags'] 
powers = data['powers']- expsettings['Qbit_Attenuation']

plt.figure(14)
plt.clf()


ax = plt.subplot(1,2,1)
plots.general_colormap_subplot(ax,freqs, powers, mags, ['f', 'mag'], 'Hello World')
#plt.imshow(mags)
ax.set_aspect('auto')

ax = plt.subplot(1,2,2)
#powerList = [0,2,4,6]
powerList = [4]
waterfall = 12
for ind in range(0, len(powerList)):
    pind = powerList[ind]
    power = powers[pind]
    print(power)
    plt.plot(freqs/1e9, mags[pind,:] + pind*waterfall, label = str(np.round(power,3)) + ' 20s VNA', alpha = 0.7)
    
    
avgLine = np.mean(mags, 0)
plt.plot(freqs/1e9, avgLine + (pind+1)*waterfall, label = '160s hand average', color = 'firebrick', zorder = 5, linewidth = 2)


partialAvg = np.mean(mags[0:3,:], 0)
plt.plot(freqs/1e9, partialAvg + (pind+2)*waterfall, label = '60s hand average', color = 'mediumblue', zorder = 5, linewidth = 2)

hack = True
if hack:
    dat = userfuncs.LoadFull(r'Z:\Data\HouckDualHangerFluxonium\Spec\20210119\flxon_powersweep_zoom2_20210119_170547.pkl')
    data = dat[0]
    expsettings = dat[1]
    freqs = data['freqs']
    mags = data['mags'] 
    powers = data['powers']- expsettings['Qbit_Attenuation']
    
    power = powers[-2]
    print(power)
    plt.plot(freqs/1e9, mags[-2,:] + (pind+3)*waterfall, label = str(np.round(power,3)) + ' 60s VNA', alpha = 1, color = 'orange')



#ax.set_xlim([5.88, 5.93])
#ax.set_xlim([5.95, 6.05])
plt.ylabel('spec mag')
plt.xlabel('freq (GHz)')
ax.legend(loc = 'lower right')
plt.title('1 kHz thorlabs')

plt.show()