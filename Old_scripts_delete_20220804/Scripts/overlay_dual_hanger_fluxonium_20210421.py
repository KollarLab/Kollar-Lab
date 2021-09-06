# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 17:36:24 2021

@author: Kollarlab
"""

from overlay_helpers import stitch_colormaps, overlay_GUI, fluxonium_energies, fluxonium_plot_levels
import matplotlib.pyplot as plt

saveDir = r'C:\Users\Kollarlab\Downloads\FluxoniumPlotting\dual_hanger_fluxonium\low_cav'

wide_scans = ['wide_scan_20210202_192204.pkl', 
              'tracking_fluxon_below_cav_m15dBm_20210203_153809.pkl',
              'higher_power_above_cav_m15dBm_20210203_163109.pkl',
              'higher_power_above_cav_m10dBm_20210204_003416.pkl',
              'high_power_sweep_below_cavity_20210216_200242.pkl',
              'high_power_sweep_above_cavity_20210217_011205.pkl'
              ]
fig = plt.figure()
ax = plt.subplot(111)

volts_per_flux = 3.12629
offset = 1.215

EC = 1.18
EJ = 10.30
EL = 0.87
cav_freq = 7.57696
raman_freq = 6
min_flux = -0.4
max_flux = 0.25
flux_points = 51
min_freq = 4
max_freq = 14

twoPhoton = True
cavAssist = True
g2Transitions = True
Raman = True

ax.set_ylim([min_freq, max_freq])
#ax.set_xlim([min_flux, max_flux])

stitch_colormaps(ax, wide_scans, saveDir, 'Blues', limits=[-5,5], convert_to_flux=True, offset=offset, volts_per_flux=volts_per_flux)

UI = []

ax.set_ylim([min_freq, max_freq])
ax.set_xlim([min_flux, max_flux])

lines = []
energies, flux = fluxonium_energies(EC, EJ, EL, min_flux, max_flux, flux_points)
lines = fluxonium_plot_levels(ax, lines, energies, flux, cav_freq, raman_freq, cavAssist=True, g2Transitions=True)
#overlay_GUI(fig, ax, UI, EC, EJ, EL, cav_freq, raman_freq, min_flux, max_flux, 51)