# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 11:47:12 2021

@author: Kollarlab
"""
import userfuncs
import plotting_tools as plots

from scipy.signal import savgol_filter
import os
import numpy as np
import matplotlib.pyplot as plt

projectfolder = r'Z:\Data\Fluxonium_Raman\CRF01_C3\Trans\20210127'
filename = 'initial_trans_sweep_m250mV_20210127_000346.pkl'
fullname = os.path.join(projectfolder, filename)


window = 11
polyorder = 5

data = dat[0]
#voltages = dat[1]['voltages']
powers = data['powers']

spec_mags = data['mags']
spec_freqs = data['freqs']

filtmag = np.zeros(spec_mags.shape)
spec_plt = np.zeros(spec_mags.shape)
#plt.figure()
#plt.plot(spec_freqs, spec_mags[0])
#for i in range(len(window)):
#    filtmag = savgol_filter(spec_mags[0], window[i], polyorder)
#    plt.plot(spec_freqs, filtmag+5*i)
#legend = [str(w) for w in window]
#plt.legend(legend)
#plt.xlabel('Freq (GHz)')
#plt.ylabel('Mag')

mins = np.zeros(len(spec_mags))
maxs = np.zeros(len(spec_mags))
for i in range(len(spec_mags)):
    baseline = np.mean(spec_mags[i][0:50])
    filtmag[i] = savgol_filter(spec_mags[i]/baseline, window, polyorder)
    spec_plt[i] = spec_mags[i]/baseline
#    filtmag[i] = spec_mags[i]-baseline
    mins[i] = np.mean(np.partition(filtmag[i],4)[:1])
    maxs[i] = np.mean(np.partition(filtmag[i],-4)[-1:])

#labels = ['Freq(GHz)', 'Voltage(V)']
labels = ['Freq(GHz)', 'Power(dBm)']

cmap_min = np.mean(mins)
cmap_max = np.mean(maxs)

cmap_range = [abs(cmap_min), abs(cmap_max)]
plt.figure()
ax = plt.subplot(211)
plots.general_colormap_subplot(ax,spec_freqs, powers, spec_plt, labels, 'Unfiltered', 'Blues')
ax = plt.subplot(212)
plots.general_colormap_subplot(ax,spec_freqs, powers, filtmag, labels, 'Filtered (baseline removed), window:{} polyorder:{}'.format(window, polyorder),'Blues')
#plots.general_colormap_subplot(ax,spec_freqs, voltages, spec_freqs, labels, 'Baseline removed', 'Blues')