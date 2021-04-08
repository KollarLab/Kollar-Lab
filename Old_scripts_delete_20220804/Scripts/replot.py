# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 12:17:47 2021

@author: Kollarlab
"""
import userfuncs
import os
import plotting_tools as plots

saveDir = r'Z:\Data\HouckDualHangerFluxonium\spec_flux_scan\20210212'
filename = 'wide_sweep_cavity1_20210212_174851.pkl'
fullpath = os.path.join(saveDir, filename)

raw_data = userfuncs.LoadFull(fullpath)
data = raw_data[0]

#transdata = data['transdata']
#specdata = data['specdata']
#singledata = data['singledata']
#voltages = data['voltages']
#trans_labels = data['trans_labels']
#spec_labels = data['spec_labels']
#plots.autoscan_plot(transdata, specdata, singledata, voltages, trans_labels, spec_labels, 'replot')

transdata = {}
transdata['xaxis'] = data['trans_freqs']
transdata['mags'] = data['trans_mags']
transdata['phases'] = data['trans_phases']

specdata = {}
specdata['xaxis'] = data['freqs']
specdata['mags'] = data['mags']
specdata['phases'] = data['phases']

singledata = {}
singledata['xaxis'] = data['freqs']
singledata['mag'] = data['mags'][-1]
singledata['phase'] = data['phases'][-1]

trans_labels = ['Freq (GHz)','Voltage (V)']
spec_labels = ['Freq (GHz)','Voltage (V)']
voltages = raw_data[1]['spec']['voltages']

plots.autoscan_plot(transdata, specdata, singledata, voltages[0:len(transdata)], filename, trans_labels, spec_labels, 'replot')