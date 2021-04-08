# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 15:19:57 2021

@author: Kollarlab
"""

import userfuncs
import plotting_tools as plots

import os
import numpy as np
import matplotlib.pyplot as plt

projectfolder = r'Z:\Data\HouckDualHangerFluxonium\spec_flux_scan\20210203'
filename = 'tracking_fluxon_below_cav_20210203_124112.pkl'
fullname = os.path.join(projectfolder, filename)

raw_data = userfuncs.LoadFull(fullname)
variables = raw_data[0]
settings  = raw_data[1]

voltages = variables['voltages']
spec_data  = variables['specdata']
trans_data = variables['transdata']

spec_mags  = spec_data['mags']
spec_freqs = spec_data['xaxis']

filtmag = np.zeros(spec_mags.shape)
spec_plt = np.zeros(spec_mags.shape)


mins = np.zeros(len(spec_mags))
maxs = np.zeros(len(spec_mags))
for i in range(len(spec_mags)):
    baseline = np.mean(spec_mags[i][0:50])
    spec_plt[i] = spec_mags[i]-baseline
    filtmag[i] = spec_mags[i]
    mins[i] = np.mean(np.partition(spec_plt[i],10)[:10])
    maxs[i] = np.mean(np.partition(spec_plt[i],-10)[-10:])

labels = ['Freq(GHz)', 'Voltage(V)']

cmap_min = np.mean(mins)
cmap_max = np.mean(maxs)

cmap_range = [abs(cmap_min), abs(cmap_max)]
plt.figure()
ax = plt.subplot(211)
plots.general_colormap_subplot(ax,spec_freqs, voltages, spec_plt, labels, 'Unfiltered', 'viridis')
ax = plt.subplot(212)
plots.general_colormap_subplot(ax,spec_freqs, voltages, filtmag, labels, 'Improved', 'afmhot_r')