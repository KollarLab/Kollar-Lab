# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:35:43 2021

@author: Kollarlab
"""

from overlay_helpers import stitch_colormaps, overlay_GUI
import matplotlib.pyplot as plt

saveDir = r'C:\Users\Kollarlab\Downloads\FluxoniumPlotting\fluxonium_raman_20210305'

wide_scans = ['above_cavity_20210304_135121.pkl', 
              'fluxon_below_20210304_142416.pkl',
              'fluxon_below_track_20210304_153824.pkl',
              'cavity_assisted_search_below_20210308_231507.pkl', 
              'fine_around_cavity_20210304_221046.pkl',
              'fluxon_tracking_above_20210304_185318.pkl', 
              'cavity_assisted_search_above_20210309_015119.pkl',
              'cavity_assisted_search_half_flux_20210309_033536.pkl',
              'spec_fine_check_ultrafine_higher_power_20210310_014931.pkl',
              'Raman_assisted_search_5_spec_-65_cavity_20210325_172754.pkl',
              'Raman_assisted_search_5_spec_-65_cavity_20210326_010030.pkl',
              'Raman_assisted_search_-5_spec_-65_cavity_20210326_042649.pkl',
              'Raman_assisted_search_5_spec_-65_cavity_20210326_055342.pkl',
              'High_freq_sweep_5_spec_-65_cavity_20210326_100553.pkl',
              'Integer_flux_below_cavity_5_spec_-65_cavity_20210326_113257.pkl',
              ]
fig = plt.figure()
ax = plt.subplot(111)

volts_per_flux = 0.750116
offset = 0.2423

EC = 0.74
EJ = 9.06
EL = 0.54
cav_freq = 6.56
raman_freq = 6
min_flux = -0.6
max_flux = 0.6
flux_points = 51
min_freq = 3
max_freq = 14

twoPhoton = True
cavAssist = True
g2Transitions = True
Raman = True

ax.set_ylim([min_freq, max_freq])
ax.set_xlim([min_flux, max_flux])

stitch_colormaps(ax, wide_scans, saveDir, 'Blues', limits=[-5,5], convert_to_flux=False, offset=offset, volts_per_flux=volts_per_flux)

UI = []

ax.set_ylim([min_freq, max_freq])
ax.set_xlim([min_flux, max_flux])
overlay_GUI(fig, ax, UI, EC, EJ, EL, cav_freq, raman_freq, min_flux, max_flux, 51)