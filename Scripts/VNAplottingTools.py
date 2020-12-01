# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 10:46:45 2020

@author: Kollarlab
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.signal import find_peaks, savgol_filter

def base_power_plot(fig, ax, freqs, ydata, powers, scanname, scanformat, HWattenuation):
    mags2 = np.asarray(ydata)
    freqs2 = np.zeros(len(freqs) + 1)
    fdiff = freqs[1]-freqs[0]
    powers2 = np.zeros(len(powers)+1)
    
    if len(powers) == 1:
        pdiff = 0.01
    else:
        pdiff = powers[1] - powers[0]
    
    freqs2[0:-1] = freqs-fdiff/2
    freqs2[-1] = freqs[-1] + fdiff/2
    
    powers2[0:-1] = powers-pdiff/2
    powers2[-1] = powers[-1] + pdiff/2
    
    XX2,YY2 = np.meshgrid(freqs2,powers2+HWattenuation)
    im = plt.pcolormesh(XX2/1e9,YY2,mags2, shading = 'nearest')
        
    plt.xlabel("Frequency (GHz)")
    plt.ylabel(r"Power (dBm)")
    plt.title('S21 {}'.format(scanformat, scanname))  
    
    plt.show()

def base_power_plot_imshow(fig, ax, xdata, ydata, zdata, labels, attenuation=0):
    limits = [xdata[0], xdata[-1], ydata[0] + attenuation, ydata[-1] + attenuation]
    pos = ax.imshow(zdata, extent = limits, origin='lower', aspect='auto', cmap='viridis')
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_title(labels[2])
    
    fig.colorbar(pos, ax = ax)
    
def base_power_plot_spec(fig, ax, freqs, ydata, powers, scanname, scanformat, HWattenuation):
    mags2 = np.asarray(ydata)
    freqs2 = np.zeros(len(freqs) + 1)
    fdiff = freqs[1]-freqs[0]
    powers2 = np.zeros(len(powers)+1)
    pdiff = powers[1] - powers[0]
    
    freqs2[0:-1] = freqs-fdiff/2
    freqs2[-1] = freqs[-1] + fdiff/2
    
    powers2[0:-1] = powers-pdiff/2
    powers2[-1] = powers[-1] + pdiff/2
    
    XX2,YY2 = np.meshgrid(freqs2*1e6,powers2/1e9)
    im = plt.pcolormesh(XX2/1e9,YY2,mags2, shading = 'nearest')
        
    plt.xlabel("Time (us)")
    plt.ylabel(r"Freq (GHz)")
    plt.title('Voltage {}, {}'.format(scanformat, scanname))  
    
#    ax = plt.gca()
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    axins = inset_axes(ax,
                        width="20%",  # width = 50% of parent_bbox width
                        height="2%",  # height : 5%
                        loc='upper right',
                        bbox_to_anchor=(-0.02, -0.12, 1.0, 1.05),
                        bbox_transform=ax.transAxes,
                        borderpad=0,)
    cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
    cbar.set_label(r"Volts",labelpad = -10,y=1,x=-0.35)
    axins.xaxis.set_ticks_position("top")
    axins.tick_params(labelsize=9)
    
    plt.show()
    
def base_raw_time_plot_spec(fig, ax, times, ydata, ys, scanname, scanformat):
    ''' 
    times is the values of the time axis of raw daa plots
    please give this already in us so that the plot and the axis have the same units
    
    ys is the control parameter values for the y axis (usually frequency)
    because this could either be frequency or power, I will not handle the units nicely.
    '''
    
    mags2 = np.asarray(ydata)
    times2 = np.zeros(len(times) + 1)
    tdiff = times[1]-times[0]
    ys2 = np.zeros(len(ys)+1)
    ydiff = ys[1] - ys[0]
    
    times2[0:-1] = times-tdiff/2
    times2[-1] = times[-1] + tdiff/2
    
    ys2[0:-1] = ys-ydiff/2
    ys2[-1] = ys[-1] + ydiff/2
    
    XX2,YY2 = np.meshgrid(times2,ys2)
    im = plt.pcolormesh(XX2,YY2,mags2, shading = 'nearest')
        
    plt.xlabel("Time (us?)")
    plt.ylabel(r"control param")
    plt.title('Voltage {}, {}'.format(scanformat, scanname))  
    
#    ax = plt.gca()
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    axins = inset_axes(ax,
                        width="20%",  # width = 50% of parent_bbox width
                        height="2%",  # height : 5%
                        loc='upper right',
                        bbox_to_anchor=(-0.02, -0.12, 1.0, 1.05),
                        bbox_transform=ax.transAxes,
                        borderpad=0,)
    cbar = fig.colorbar(im,cax=axins,orientation="horizontal",ticks=plt.MaxNLocator(2))
    cbar.set_label(r"Volts",labelpad = -10,y=1,x=-0.35)
    axins.xaxis.set_ticks_position("top")
    axins.tick_params(labelsize=9)
    
    plt.show()
    
    
def power_plot(freqs, mags, phases, powers, scanname, HWattenuation = 0):
    fig = plt.figure(figsize=(13,8))
    plt.clf()
    
    ax = plt.subplot(1,2,1)
    labels = ['Frequency (GHz)','Power (dBm)', 'S21 mag']
    base_power_plot_imshow(fig, ax, freqs, powers, mags, labels, HWattenuation)
    ax = plt.subplot(1,2,2)
    labels = ['Frequency (GHz)','Power (dBm)', 'S21 phase']
    base_power_plot_imshow(fig, ax, freqs, powers, phases, labels, HWattenuation)
    
    plt.suptitle('Filename: {}, HWattenuation: {}dB'.format(scanname, HWattenuation))

def pulsed_debug(freqs, mags, phases, time, amp, phase, scanname, power):
    fig = plt.figure(figsize=(13,8))
    plt.clf()
    
    ax = plt.subplot(2,2,1)
    labels = ['Time (us)', 'Frequency (GHz)', 'Voltage (mag)']
    base_power_plot_imshow(fig, ax, time, freqs, mags, labels)
    ax = plt.subplot(2,2,2)
    labels = ['Time (us)', 'Frequency (GHz)', 'Voltage (phase)']
    base_power_plot_imshow(fig, ax, time, freqs, phases, labels)
    ax = plt.subplot(2,2,3)
    plt.plot(time, amp)
    plt.xlabel('Time (us)')
    plt.ylabel('Amp (V)')
    ax = plt.subplot(2,2,4)
    plt.plot(time, amp)
    plt.xlabel('Time (us)')
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
    