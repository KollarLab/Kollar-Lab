# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 10:46:45 2020

@author: Kollarlab
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.signal import find_peaks, savgol_filter



#def general_VNAplot(xaxis, mags, phases, yaxis, scanname, HWattenuation = 0, 
#                 xlabel = 'Frequency (GHz)', ylabel = 'Power (dBm)', identifier = '', 
#                 fig_num = ''):
#    if fig_num == '':
#        fig = plt.figure(figsize=(13,8))
#    else:
#        fig = plt.figure(fig_num, figsize=(13,8))
#    plt.clf()
#    
#    ax = plt.subplot(1,2,1)
#    labels = [xlabel,ylabel, 'S21 mag']
#    base_power_plot_imshow(fig, ax, xaxis, yaxis, mags, labels, HWattenuation)
#    ax = plt.subplot(1,2,2)
#    labels = [xlabel,ylabel, 'S21 phase']
#    base_power_plot_imshow(fig, ax, xaxis, yaxis, phases, labels, HWattenuation)
#    
#    plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
#    
#    fig.canvas.draw()
#    fig.canvas.flush_events()
#    return


def general_VNAplot(xaxis, mags, phases, yaxis, scanname, HWattenuation = 0, 
                 xlabel = 'Frequency (GHz)', ylabel = 'Power (dBm)', identifier = '', 
                 fig_num = ''):
    if fig_num == '':
        fig = plt.figure(figsize=(13,8))
    else:
        fig = plt.figure(fig_num, figsize=(13,8))
    plt.clf()
    
    ax = plt.subplot(1,2,1)
    general_colormap_subplot(ax,xaxis, yaxis-HWattenuation, mags)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('S21 mag')
    
    ax = plt.subplot(1,2,2)
    general_colormap_subplot(ax,xaxis, yaxis-HWattenuation, phases)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('S21 phase')
    
    
    plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return

def spec_fluxscanplot(transdata, specdata, singledata, yaxis, scanname, trans_labels, spec_labels, identifier, fig_num = ''):
    if fig_num == '':
        fig = plt.figure(figsize=(13,8))
    else:
        fig = plt.figure(fig_num, figsize=(13,8))
    plt.clf()
    
    ax = plt.subplot(3,2,1)
    general_colormap_subplot(ax,transdata['xaxis'], yaxis, transdata['mags'])
    plt.xlabel(trans_labels[0])
    plt.ylabel(trans_labels[1])
    plt.title('Trans mag')
    
    ax = plt.subplot(3,2,2)
    general_colormap_subplot(ax,transdata['xaxis'], yaxis, transdata['phases'])
    plt.xlabel(trans_labels[0])
    plt.ylabel(trans_labels[1])
    plt.title('Trans phase')
    
    ax = plt.subplot(3,2,3)
    general_colormap_subplot(ax,specdata['xaxis'], yaxis, specdata['mags'])
    plt.xlabel(spec_labels[0])
    plt.ylabel(spec_labels[1])
    plt.title('Spec mag')
    
    ax = plt.subplot(3,2,4)
    general_colormap_subplot(ax,specdata['xaxis'], yaxis, specdata['phases'])
    plt.xlabel(spec_labels[0])
    plt.ylabel(spec_labels[1])
    plt.title('Spec phase')
    
    ax = plt.subplot(3,2,5)
    plt.plot(singledata['xaxis'], singledata['mag'])
    plt.xlabel(spec_labels[0])
    plt.title('Single shot spec mag')
    
    ax = plt.subplot(3,2,6)
    plt.plot(singledata['xaxis'], singledata['phase'])
    plt.xlabel(spec_labels[0])
    plt.title('Single shot spec mag')
    
    
    plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
    
    fig.canvas.draw()
    fig.canvas.flush_events()
    return

def general_colormap_subplot(ax,xaxis, yaxis, data):
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
    
    plt.imshow(data, extent = limits, origin='lower', aspect='auto', cmap='jet_r')
    plt.colorbar()
        
    return
    


def base_power_plot_imshow(fig, ax, xdata, ydata, zdata, labels, attenuation=0):
#    limits = [xdata[0], xdata[-1], ydata[0] + attenuation, ydata[-1] + attenuation]
    
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
    plt.title('S21 {}'.format(scanformat))  
    
    plt.show()


    
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
    
    #ax = plt.gca()
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
    
def base_raw_time_plot_spec(fig, ax, times, ydata, ys, ylabel, zlabel, scanname, scanformat):
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
        
    plt.xlabel("Time (us)")
    plt.ylabel(ylabel)
    plt.title('Voltage {}, {}'.format(scanformat, scanname))  
    
    #ax = plt.gca()
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
    cbar.set_label(zlabel,labelpad = -10,y=1,x=-0.35)
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