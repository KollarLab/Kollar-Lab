# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 10:46:45 2020

@author: Kollarlab
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.signal import find_peaks, savgol_filter

def general_VNAplot(xaxis, mags, phases, yaxis, scanname, HWattenuation = 0, 
                 labels=['Freq (GHz)', 'Power (dBm)'], identifier = '', 
                 fig_num = ''):
    if fig_num == '':
        fig = plt.figure(figsize=(13,8))
    else:
        fig = plt.figure(fig_num, figsize=(13,8))
    plt.clf()
    
    ax = plt.subplot(1,2,1)
    general_colormap_subplot(ax,xaxis, yaxis-HWattenuation, mags, labels, 'S21 mag')
    
    ax = plt.subplot(1,2,2)
    general_colormap_subplot(ax,xaxis, yaxis-HWattenuation, phases, labels, 'S21 phase')
    
    plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return

def simplescan_plot(full_data, singledata, 
                    yaxis, 
                    scanname, 
                    labels, 
                    identifier='', 
                    fig_num='', 
                    cmap='hot', 
                    vmin=np.NaN, vmax=np.NaN,
                    IQdata = False):
    if fig_num == '':
        fig = plt.figure(figsize=(13,8))
    else:
        fig = plt.figure(fig_num, figsize=(13,8))
    plt.clf()
    
    #modifying this function so that it can plot mag/phase data or I/Q data
    if not IQdata:
        key1 = 'mags'
        title1 = 'Mag'
        
        key2 = 'phases'
        title2 = 'Phase'
    else:
        key1 = 'Is'
        title1 = 'I'
        
        key2 = 'Qs'
        title2 = 'Q'
    
    ax = plt.subplot(2,2,1)
#    general_colormap_subplot(ax,full_data['xaxis'], yaxis, full_data['mags'], labels, 'Mag', cmap, vmin, vmax)
    general_colormap_subplot(ax,full_data['xaxis'], yaxis, full_data[key1], labels, title1, cmap, vmin, vmax)
    
    ax = plt.subplot(2,2,2)
#    general_colormap_subplot(ax,full_data['xaxis'], yaxis, full_data['phases'], labels, 'Phase', cmap)
    general_colormap_subplot(ax,full_data['xaxis'], yaxis, full_data[key2], labels, title2, cmap)
    
    ax = plt.subplot(2,2,3)
    key = key1[0:-1]
#    plt.plot(singledata['xaxis'], singledata['mag'])
    plt.plot(singledata['xaxis'], singledata[key])
    plt.xlabel(labels[0])
#    plt.title('Single shot mag')
    plt.title('Single shot ' + key)
    
    ax = plt.subplot(2,2,4)
    key = key2[0:-1]
#    plt.plot(singledata['xaxis'], singledata['phase'])
    plt.plot(singledata['xaxis'], singledata[key])
    plt.xlabel(labels[0])
#    plt.title('Single shot phase')
    plt.title('Single shot ' + key)
    
    plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return

def simplescan_plot_update(full_data, 
                           singledata, 
                           yaxis, 
                           scanname, 
                           labels, 
                           datanames,
                           titles,
                           identifier='',
                           fig_num='', 
                           cmap='hot', 
                           vmin=np.NaN, 
                           vmax=np.NaN):
    if fig_num == '':
        fig = plt.figure(figsize=(13,8))
    else:
        fig = plt.figure(fig_num, figsize=(13,8))
    plt.clf()
    
    key1 = datanames[0]
    key2 = datanames[1]
    title1 = titles[0]
    title2 = titles[1]
    ax = plt.subplot(2,2,1)
    general_colormap_subplot(ax,full_data['xaxis'], yaxis, full_data[key1], labels, title1, cmap, vmin, vmax)
    
    ax = plt.subplot(2,2,2)
    general_colormap_subplot(ax,full_data['xaxis'], yaxis, full_data[key2], labels, title2, cmap)
    
    ax = plt.subplot(2,2,3)
    key = key1[0:-1]
    plt.plot(singledata['xaxis'], singledata[key])
    plt.xlabel(labels[0])
    plt.title('Single shot ' + key)
    
    ax = plt.subplot(2,2,4)
    key = key2[0:-1]
    plt.plot(singledata['xaxis'], singledata[key])
    plt.xlabel(labels[0])
    plt.title('Single shot ' + key)
    
    plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return

def autoscan_plot(transdata, specdata, singledata, yaxis, scanname, trans_labels, spec_labels, identifier, fig_num = ''):
    if fig_num == '':
        fig = plt.figure(figsize=(13,8))
    else:
        fig = plt.figure(fig_num, figsize=(13,8))
    plt.clf()
    
    ax = plt.subplot(3,2,1)
    general_colormap_subplot(ax,transdata['xaxis'], yaxis, transdata['mags'], trans_labels, 'Trans mag')
    
    ax = plt.subplot(3,2,2)
    general_colormap_subplot(ax,transdata['xaxis'], yaxis, transdata['phases'], trans_labels, 'Trans phase')
    
    ax = plt.subplot(3,2,3)
    general_colormap_subplot(ax,specdata['xaxis'], yaxis, specdata['mags'], spec_labels, 'Spec mag')
    
    ax = plt.subplot(3,2,4)
    general_colormap_subplot(ax,specdata['xaxis'], yaxis, specdata['phases'], spec_labels, 'Spec phase')
    
    ax = plt.subplot(3,2,5)
    plt.plot(singledata['xaxis'], singledata['mag'])
    plt.xlabel(spec_labels[0])
    plt.title('Single shot spec mag')
    
    ax = plt.subplot(3,2,6)
    plt.plot(singledata['xaxis'], singledata['phase'])
    plt.xlabel(spec_labels[0])
    plt.title('Single shot spec phase')
    
    plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return

def general_colormap_subplot(ax, xaxis, yaxis, data, labels, 
                             title, cmap = 'hot', vmin = np.NaN, vmax = np.NaN, cbar=None):
    '''given a subplot object, it should make a decent imshow plot of the data.
    
    You have to handle subtracting out any attenuation and do the labeling after.
    
    Trying to make this the most generic function possible.'''
    
    plt.sca(ax)
    
    #this version will handle the color plot if there is only one row.
    if (type(xaxis) == list) or (type(xaxis) == np.ndarray):
        if (len(xaxis) == 1): 
            xlimits = [xaxis[0]-1e-3, xaxis[0]+1e-3]
        else:
            xlimits = [xaxis[0], xaxis[-1]]
    else:
        xlimits =  [xaxis-1e-3, xaxis+1e-3]
        
    if (type(yaxis) == list) or (type(yaxis) == np.ndarray):
        if (len(yaxis)) == 1:
            ylimits = [yaxis[0]-1e-3, yaxis[0]+1e-3]
        else:
            ylimits = [yaxis[0], yaxis[-1]]  
    else:
        ylimits = [yaxis-1e-3, yaxis+1e-3]
        
    limits = np.concatenate((xlimits, ylimits))
    
    if (not np.isnan(vmax)) and (not np.isnan(vmin)):
#        print('version with color scale limits')
        im = plt.imshow(data, extent = limits, origin='lower', aspect='auto', cmap=cmap, vmin = vmin, vmax = vmax)
    else:
        im = plt.imshow(data, extent = limits, origin='lower', aspect='auto', cmap=cmap)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.title(title)
    if not cbar:
        cbar = plt.colorbar()
        return cbar
    else:
        cbar.update_normal(im)
        return
    
def base_power_plot_imshow(fig, ax, xdata, ydata, zdata, labels, attenuation=0):
    #limits = [xdata[0], xdata[-1], ydata[0] + attenuation, ydata[-1] + attenuation]
    
    #this version will handle the color plot if there is only one row.
    if (type(xdata) == list) or (type(xdata) == np.ndarray):
        if (len(xdata) == 1): 
            xlimits = [xdata[0]-1e-3, xdata[0]+1e-3]
        else:
            xlimits = [xdata[0], xdata[-1]]
    else:
        xlimits =  [xdata-1e-3, xdata+1e-3]
        
    if (type(ydata) == list) or (type(ydata) == np.ndarray):
        if (len(ydata)) == 1:
            ylimits = [ydata[0]+ attenuation-1e-3, ydata[0]+ attenuation+1e-3]
        else:
            ylimits = [ydata[0] +attenuation, ydata[-1]+ attenuation]  
    else:
        ylimits = [ydata+attenuation-1e-3, ydata+attenuation+1e-3]
        
        
    limits = np.concatenate((xlimits, ylimits))
    
    pos = ax.imshow(zdata, extent = limits, origin='lower', aspect='auto', cmap='viridis')
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_title(labels[2])
    
    fig.colorbar(pos, ax = ax)
    
def pulsed_debug(fig, freqs, powers, mags, phases, amp, phase, scanname, power):
    
    ax = plt.subplot(2,2,1)
    labels = ['Frequency (GHz)', 'Power (dBm)', 'Voltage (mag)']
    base_power_plot_imshow(fig, ax, freqs/1e9, powers, mags, labels)
    ax = plt.subplot(2,2,2)
    labels = ['Frequency (GHz)', 'Power (dBm)', 'Voltage (phase)']
    base_power_plot_imshow(fig, ax, freqs/1e9, powers, phases, labels)
    ax = plt.subplot(2,2,3)
    plt.plot(freqs/1e9, amp)
    plt.xlabel('Freq (GHz)')
    plt.ylabel('Amp (V)')
    ax = plt.subplot(2,2,4)
    plt.plot(freqs/1e9, phase)
    plt.xlabel('Freq (GHz)')
    plt.ylabel('Phase (deg)')
    
    plt.suptitle('Filename: {}, Power: {}dB'.format(scanname, power))
    
def get_peaks(freqs, mags, window, polyorder, height, width, show_plot = False):
    
    filtmag = savgol_filter(mags, window, polyorder)
    
    peaks,_ = find_peaks(-filtmag, height = height, width = width)
    
    if show_plot:
        fig = plt.figure()
        plt.plot(freqs, mags)
        plt.plot(freqs, filtmag)
        plt.plot(freqs[peaks], filtmag[peaks], 'x')
        plt.title('Peaks in signal')
    
    return [freqs[peaks], mags[peaks]]