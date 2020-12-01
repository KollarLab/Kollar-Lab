# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 16:51:32 2020

'Takes in a list of hanger dips and goes to those to measure 'S21' signal 
for a center frequency and span'

Make sure to correct for the electrical delay to get the best phase results
@author: Kollarlab
"""
from VNAplottingTools import get_peaks
from Qcalcs import fit_Q
import userfuncs
import os
import numpy

saveDir = r'Z:\Data\PhosphorusDopedSilicon\SH2006A\20201125'
stamp = userfuncs.timestamp()
name = 'HangerDips'
temp = '227mK'
scanname = name + temp + '_' + stamp

broadscan_amp = vna.get_trace(1,'mag')
broadscan_f = vna.get_channel_axis(1)

hanger_dips_f, hanger_dips_amp = get_peaks(broadscan_f, broadscan_amp, window=1, polyorder=0, height=15, width=1, show_plot=True)
plt.savefig(os.path.join(saveDir, 'broadscan'+stamp+'.png'), dpi = 150)

#input('Make sure that the code has picked out all the peaks you were looking for, press any key to continue')

spans = [1e6]*len(hanger_dips_f)

settings = {}         
settings['ifBW'] = 100
settings['RFpower'] = -30
settings['channel'] = 1
settings['avg_time'] = 60
settings['sweep_points'] = 1001

mags, phases, freqs = vna.meas_lines(hanger_dips_f, spans, settings)

centers = numpy.zeros(len(hanger_dips_f))
Q_ints = numpy.zeros(len(hanger_dips_f))
Q_exts = numpy.zeros(len(hanger_dips_f))

for index in range(len(hanger_dips_f)):
    c, qi, qe = fit_Q(freqs[index], mags[index], hanger_dips_f[index])
    print('')
    plt.savefig(os.path.join(saveDir, 'F{}{}.png'.format(index, stamp)), dpi = 150)
    
    centers[index] = c
    Q_ints[index] = qi
    Q_exts[index] = qe

varsTosave = ['mags', 'phases', 'freqs', 'hanger_dips_f', 'hanger_dips_amp', 'broadscan_f', 'broadscan_amp', 'centers', 'Q_ints', 'Q_exts']

plt.figure(figsize=(12,7))
ax = plt.subplot(1,2,1)
ax.plot(centers/(2*numpy.pi), Q_ints,'x')
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Q_int')
ax.grid(True)
 
ax = plt.subplot(1,2,2)
ax.plot(centers/(2*numpy.pi), Q_exts,'x')
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Q_ext')
ax.grid(True)

plt.suptitle('SH2006A Temp:{}'.format(temp))
plt.savefig(os.path.join(saveDir, 'Qtrends'+stamp+'.png'), dpi = 150)
userfuncs.SaveFull(saveDir, scanname, varsTosave, locals(), expsettings=settings)