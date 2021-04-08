# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:50:09 2021

@author: Kollarlab
"""

import os
import numpy as np

from matplotlib.widgets import Slider, Button, CheckButtons

from FluxoniumFunc import flux_sweep
import userfuncs

def stitch_colormaps(ax, files, saveDir, cmap, limits=[-5,5], convert_to_flux=False, offset=0, volts_per_flux=1):
    '''
    Takes a list of files and combines them into on big colormap plot
    The code will subtract the mean from every row of the data to smooth out 
    variations in the transmission plot
    The xaxis can be volts or flux quanta
    Params:
        ax: axis of the figure where plots should be made
        files: list of file names (ordered by the user to get best results) of the colormaps we want to plot
        saveDir: directory containing all the files
        cmap: colormap color scheme ('Blues' seems to be quite nice)
        limits: colormap limits, tighter limits can highlight faint features but also makes the plots noisier
        convert_to_flux: boolean to toggle the x axis
        offset: offset voltage corresponding to integer flux
        volts_per_flux: conversion between voltage and flux quanta
    '''
    for file in files:
        print(file)
        dat = userfuncs.LoadFull(os.path.join(saveDir, file))
        
        data = dat[0]
        
        voltages = data['voltages']
        specdata = data['specdata']
        mags  = specdata['mags']
        freqs = specdata['xaxis']
        
        for ind, mag in enumerate(mags):
            mags[ind] = mag - np.mean(mag)
        
        mags = np.transpose(mags)
        if convert_to_flux:
            print('flux')
            xaxis = volt_to_flux(voltages, offset, volts_per_flux)
        else:
            xaxis = voltages
        im=ax.pcolormesh(xaxis, freqs/1e9, mags, cmap=cmap, vmin=limits[0], vmax=limits[1])
        ax.set_ylabel('Frequency (GHz)')
        if convert_to_flux:
            ax.set_xlabel('Flux (phi_ext/phi_0)')
        else:
            ax.set_xlabel('Voltage (V)')
        
#    plt.colorbar(im, ax=ax)

def bare_transitions(ax, lines, full_energies, flux):
    '''
    Plot the n lowest energy levels of the fluxonium with the groundstate subtracted
    '''
    labels = ['g2', 'e1', 'e2', 'e3']
    colors = ['navy', 'green', 'red', 'purple']
    marker = 'o'
    energies = full_energies[1:] - full_energies[0]
    for energy, label, color in zip(energies, labels, colors):
        line = plot_line(ax, flux, energy, color=color, label=label, marker=marker)
        lines.append(line)
    ax.legend()
        
def two_photon(ax, lines, full_energies, flux):
    '''
    Plot two photon transitions (essentially higher energy/2)
    '''
    labels = ['e2/2', 'e3/2', 'e4/2']
    colors = ['red', 'purple', 'black']
    marker = 'd'
    energies = (full_energies[3:] - full_energies[0])/2
    for energy, label, color in zip(energies, labels, colors):
        line = plot_line(ax, flux, energy, color=color, label=label, marker=marker)
        lines.append(line)
    ax.legend()
    
def cavity_assisted(ax, lines, full_energies, flux, cavity_freq):
    '''
    Plot the cavity assisted transitions
    '''
    labels = ['e2 cav assist', 'e3 cav assist', 'e4 cav assist']
    colors = ['red', 'purple', 'black']
    marker = 'd'
    linestyle = '-.'
    energies = (full_energies[3:] - full_energies[0]) - cavity_freq
    for energy, label, color in zip(energies, labels, colors):
        line = plot_line(ax, flux, energy, color=color, label=label, marker=marker, linestyle=linestyle)
        lines.append(line)
    ax.legend()
    
def g2_transitions(ax, lines, full_energies, flux):
    '''
    Plot transitions starting from the g2 level
    '''
    labels = ['g2->e1', 'g2->e2', 'g2->e3']
    colors = ['orange', 'brown', 'limegreen']
    marker = 'x'
    linestyle = ':'
    energies = full_energies[2:] - full_energies[1]
    for energy, label, color in zip(energies, labels, colors):
        line = plot_line(ax, flux, energy, color=color, label=label, marker=marker, linestyle=linestyle)
        lines.append(line)
    ax.legend()
    
def g2_Raman(ax, lines, full_energies, flux, cavity_freq):
    '''
    Plot the g2 raman transitions
    '''
    labels = ['cav+g2', 'cav-g2']
    colors = ['magenta', 'deepskyblue']
    marker = 'x'
    linestyle = ':'
    energies = [full_energies[1] + cavity_freq- full_energies[0], cavity_freq - (full_energies[1] - full_energies[0])]
    for energy, label, color in zip(energies, labels, colors):
        line = plot_line(ax, flux, energy, color=color, label=label, marker=marker, linestyle=linestyle)
        lines.append(line)
    ax.legend()

def overlay_GUI(fig, ax, UI, EC, EJ, EL, cav_freq, raman_freq=6, min_flux=-0.5, max_flux=0.5, flux_points=51, bare=True, twoPhoton=False, cavAssist=False, g2Transitions=False, g2Raman=False):
    '''
    GUI allowing the user to tweak the fluxonium values to get the best data fit by eye
    There are also options to change which transitions are plotted to minimize 
    clutter in the early fitting stages
    Params:
        fig: reference figure
        ax: reference axes
        UI: list containing pointers to all the widget objects (important to keep them active)
        EC (GHz)
        EJ (GHz)
        EL (GHz)
        cav_freq (GHz)
        raman_freq (GHz)
        min_flux (phi/2pi)
        max_flux (phi/2pi)
        flux_points (int)
        bare, twoPhoton, cavAssist, g2Transitions , g2Raman: booleans controlling 
        which transitions are plotted
    '''
    # Get the starting energy spectrum and plot the lines
    # lines, energies and flux are variables used throughout the code so need to be global (I think?)
    global lines
    global energies
    global flux
    
    lines = []
    energies, flux = fluxonium_energies(EC, EJ, EL, min_flux, max_flux, flux_points)
    lines = fluxonium_plot_levels(ax, lines, energies, flux, cav_freq, raman_freq, bare, twoPhoton, cavAssist, g2Transitions, g2Raman)
    
    # Adjust the subplots region to leave some space for the sliders and buttons
    fig.subplots_adjust(left=0.1, bottom=0.35)
    axis_color = 'lightgoldenrodyellow'
    
    #Sliders to perform fit of fluxonium spectrum (all units are GHz)
    slider_length = 0.4
    slider_height = 0.03
    slider_xpos = 0.2
    slider_ypos = [0.2, 0.15, 0.1, 0.05]
    
    cav_slider_range = 0.2
    EC_slider_range  = 1.
    EJ_slider_range  = 4.
    EL_slider_range  = 1.

    Cav_slider = custom_slider(fig, [slider_xpos, slider_ypos[0], slider_length, slider_height], 
                               axis_color, cav_freq, cav_slider_range, 'Cav')
    EC_slider = custom_slider(fig, [slider_xpos, slider_ypos[1], slider_length, slider_height], 
                               axis_color, EC, EC_slider_range, 'EC')
    EJ_slider = custom_slider(fig, [slider_xpos, slider_ypos[2], slider_length, slider_height], 
                               axis_color, EJ, EJ_slider_range, 'EJ')
    EL_slider = custom_slider(fig, [slider_xpos, slider_ypos[3], slider_length, slider_height], 
                               axis_color, EL, EL_slider_range, 'EL')    
    
    #Reset button, will return all the sliders to their middle position and reset the checkboxes
    reset_button_ax = fig.add_axes([0.85, 0.025, 0.1, 0.04])
    reset_button = Button(reset_button_ax, 'Reset', color=axis_color,
                          hovercolor='0.975')
    
    #Checkboxes controlling which transitions we want to plot
    activated = [True, False, False, False, False]
    labels = ['Bare transitions', 'Two photon', 'Cav assist', 'g2 trans', 'Raman']
    ax_check = fig.add_axes([0.7, 0.05, 0.08, 0.1])
    plot_button = CheckButtons(ax_check,labels, activated)
    
    #Update plot everytime that the user changes the qubit parameters
    def sliders_on_changed(val):
        
        global lines
        global energies
        global flux
        
        EC = EC_slider.val
        EJ = EJ_slider.val
        EL = EL_slider.val
        cav_freq = Cav_slider.val
        
        status = plot_button.get_status()
        bare, twoPhoton, cavAssist, g2Transitions, g2Raman = status
        
        energies, flux = fluxonium_energies(EC, EJ, EL, min_flux, max_flux, flux_points)
        lines = fluxonium_plot_levels(ax, lines, energies, flux, cav_freq, raman_freq, bare, twoPhoton, cavAssist, g2Transitions, g2Raman)
        
        fig.canvas.draw_idle()

    #Reset all values to starting values
    def reset(mouse_event):
        global lines
        
        EJ_slider.reset()
        EC_slider.reset()
        EL_slider.reset()
        Cav_slider.reset()
        print('resetting')
        
        bare, twoPhoton, cavAssist, g2Transitions, g2Raman = [True, False, False, False, False]
        plot_button.set_active(1)
        plot_button.set_active(2)
        plot_button.set_active(3)
        plot_button.set_active(4)
        energies, flux = fluxonium_energies(EC, EJ, EL, min_flux, max_flux, flux_points)
        lines = fluxonium_plot_levels(ax, lines, energies, flux, cav_freq, raman_freq, bare, twoPhoton, cavAssist, g2Transitions, g2Raman)

    #Plot transitions/ spectral lines specified by the user (does not recompute spectrum)
    def select_plot(label):
        global lines
        global energies
        global flux
        
        status = plot_button.get_status()
        bare, twoPhoton, cavAssist, g2Transitions, g2Raman = status
        lines = fluxonium_plot_levels(ax, lines, energies, flux, cav_freq, raman_freq, bare, twoPhoton, cavAssist, g2Transitions, g2Raman)
        fig.canvas.draw_idle()
    
    #Activate all the buttons and sliders
    Cav_slider.on_changed(sliders_on_changed)
    EC_slider.on_changed(sliders_on_changed)
    EJ_slider.on_changed(sliders_on_changed)
    EL_slider.on_changed(sliders_on_changed)
    plot_button.on_clicked(select_plot)
    reset_button.on_clicked(reset)
    
    #CRITICAL: widgets need to have a reference in the 'main' script (UI) here otherwise
    #they become unresponsive 
    sliders = [EL_slider, EJ_slider, EC_slider, Cav_slider]
    buttons = [plot_button, reset_button]
    UI.append(sliders + buttons)
    
#### 
#HELPER FUNCTIONS
####

def custom_slider(fig, position, color, center, slider_range, name):
    '''
    Makes a custom slider at a specified position/ size (x,y,length, height)
    center/ initial value at center and specified slider range
    '''
    slider_ax = fig.add_axes(position, facecolor=color)
    
    start, stop = np.round(center-slider_range/2,3), np.round(center+slider_range/2,3)
    slider = Slider(slider_ax, name+' [{}:{}]'.format(start, stop), start, stop, valinit=center)
    
    return slider

def fluxonium_energies(EC, EJ, EL, min_flux, max_flux, flux_points):
    '''
    Computes the spectrum of the fluxonium specified by EC, EJ, EL over the 
    flux range [min_flux, max_flux] with flux_points points
    Params:
        EC (GHz)
        EJ (GHz)
        EL (GHz)
        min_flux: min ext flux/2pi 
        max_flux: max ext flux/2pi
        flux_points: positive int (should be odd to avoid weird aliasing)
    Returns:
        flux: np array of fluxes
        energies: matrix holding the computed eigen vals over flux range
    '''
    flux = np.linspace(min_flux, max_flux, flux_points)
    (evals_s, _) = flux_sweep(flux*2*np.pi, EC, EJ, EL)
    energies = evals_s.T
    
    return energies, flux
    
def fluxonium_plot_levels(ax, lines, energies, flux, cav_freq, raman_freq=6, bare=True, twoPhoton=False, cavAssist=False, g2Transitions=False, Raman=False):
    '''
    Plot all the user specified lines from provided spectral data:
        lines: list of lines to keep track of so we can remove them before plotting new ones
        energies: matrix of the n lowest energies that are being considered
        flux: np array of flux values at which the energies were computed
        cav_freq: frequency of the cavity in GHz
        raman_freq: frequency of the raman cavity in GHz
        bare: boolean indicating if the bare transistions should be plotted
        twoPhoton: boolean indicating if the two photon transistions should be plotted
        cavAssist: boolean indicating if the cavity assisted transistions should be plotted
        g2Transitions: boolean indicating if the transistions starting from g2 should be plotted
        g2Raman: boolean indicating if the Raman transitions from g2 should be plotted
    '''
    for line in lines:
        try:
            line.remove()
        except:
            continue
    lines = []
    
    if bare:
        bare_transitions(ax, lines, energies, flux)
    if twoPhoton:
        two_photon(ax, lines, energies, flux)
    if cavAssist:
        cavity_assisted(ax, lines, energies, flux, cav_freq)
    if g2Transitions:
        g2_transitions(ax, lines, energies, flux)
    if Raman:
        g2_Raman(ax, lines, energies, flux, raman_freq)
    
    ax.legend(
        bbox_to_anchor = (1.1,1),
        loc = 1,
        borderaxespad = 0.0
    )
        
    return lines

def plot_line(ax, flux, energy, linewidth=2, alpha=0.5, color='red', linestyle='-', marker='o', label='default'):
    '''
    Standard line plotting function to keep everything consistent and uncluttered late in the code
    '''
    [line, ] = ax.plot(
        flux,
        energy,
        linewidth=linewidth,
        alpha=alpha,
        color=color,
        linestyle=linestyle, 
        marker=marker,
        label = label
        )
    return line

def volt_to_flux(volts, offset, volts_per_flux):
    '''
    Helper function to convert volts to flux
    Params:
        volts: applied external voltage (double or array of doubles)
        offset: voltage where we have integer flux (double)
        volts_per_flux: conversion from volts to flux quantum (double)
    Returns:
        flux or flux array rounded to 4 places
    '''
        
    return np.round((volts-offset)/volts_per_flux, 4)

def flux_to_volt(flux, offset, volts_per_flux):
    '''
    Helper function to convert volts to flux
    Params:
        flux: external flux (double or array of doubles)
        offset: voltage where we have integer flux (double)
        volts_per_flux: conversion from volts to flux quantum (double)
    Returns:
        volts or volts array rounded to 4 places
    '''
    return np.round(flux*volts_per_flux + offset, 4)