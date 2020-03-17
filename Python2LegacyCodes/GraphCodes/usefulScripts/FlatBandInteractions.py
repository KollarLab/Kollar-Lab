#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 15:30:44 2017

@author: kollar2
"""

import re
import scipy
import pylab
import numpy
import time

import pickle
import datetime
import os
import sys

import scipy.linalg

#from LayoutGenerator3 import PlanarLayout
from LayoutGenerator4 import PlanarLayout



#def get_index(num, itt, az = True):
#    '''
#    get the index location of a semidual point. Point spcified by
#    
#    number within the itteration points/radials
#    itteration number
#    az (+True for azimuthal or false for radial)
#    '''
#    currInd = 0
#    for i in range(0, itt):
#        currInd = currInd + len(test.SDpoints[i]) + test.SDradials[i].shape[0]
#        
#    if az == False:
#        currInd = currInd + len(test.SDpoints[itt])
#
#    currInd = currInd + num
#    
#    return currInd

#def V_int(ind1, ind2, states):
#    '''
#    Calculate total interaction enegery of two particles at lattice sites
#    indexed by index 1 and index 2
#    
#    states is the set of eigenvectors that you want to include
#    '''
#    psis_1 = test.Psis[ind1,states]
#    psis_2 = test.Psis[ind2,states]
#    
#    return numpy.dot(numpy.conj(psis_2), psis_1)

#def V_int_map(source_ind):
#    '''
#    calculate a map of the interaction energy for a given location of the first
#    qubit.
#    Lattice sites specified by index in semidual points array
#    '''
#    
#    int_vals = numpy.zeros(len(sdx))
#    for ind2 in range(0, len(sdx)):
#        int_vals[ind2] = V_int(source_ind, ind2)
#    
#    return int_vals

#def plot_layout_state(state_vect, ax, title = 'state weight', colorbar = False, plot_links = False, cmap = 'Wistia'):
#    '''
#    plot a state on the graph
#    '''
#    Amps = state_vect
#    Probs = numpy.abs(Amps)**2
#    mSizes = Probs * len(Probs)*30
#    mColors = numpy.angle(Amps)/numpy.pi
#    
#    cm = pylab.cm.get_cmap(cmap)
#    
#    sdx, sdy = test.get_all_semidual_points()
#    
#    pylab.sca(ax)
#    pylab.scatter(sdx, sdy,c =  mColors, s = mSizes, marker = 'o', edgecolors = 'k', cmap = cm, vmin = -0.5, vmax = 1.5)
#    if colorbar:
#        cbar = pylab.colorbar(fraction=0.046, pad=0.04)
#        cbar.set_label('phase (pi radians)', rotation=270)
#          
#    if plot_links:
#        if test.itter>3:
#            test.draw_SDlinks(ax, 3, linewidth = 0.5, color = 'firebrick')
#        else:
#            test.draw_SDlinks(ax, test.itter,linewidth = 0.5, color = 'firebrick')
#    
#    pylab.title(title, fontsize=8)
#    ax.xaxis.set_visible(False)
#    ax.yaxis.set_visible(False)
#    ax.set_aspect('equal')
#    return
    
#def plot_map_state(map_vect, ax, title = 'ineraction weight', colorbar = False, plot_links = False, cmap = 'winter'):
#    '''plot an interaction map on the graph
#    '''
#    Amps = map_vect
#    
#    mSizes = 10
#    mColors = Amps
#    
##    cm = pylab.cm.get_cmap('seismic')
#    cm = pylab.cm.get_cmap(cmap)
##    cm = pylab.cm.get_cmap('RdBu')
#    
#    sdx, sdy = test.get_all_semidual_points()
#    
#    pylab.sca(ax)
#    pylab.scatter(sdx, sdy,c =  mColors, s = mSizes, marker = 'o', edgecolors = 'k', cmap = cm, vmin = -0.5, vmax = 1.5)
#    if colorbar:
#        cbar = pylab.colorbar(fraction=0.046, pad=0.04)
#        cbar.set_label('phase (pi radians)', rotation=270)
#          
#    if plot_links:
#        if test.itter>3:
#            test.draw_SDlinks(ax, 3, linewidth = 0.5, color = 'firebrick')
#        else:
#            test.draw_SDlinks(ax, test.itter,linewidth = 0.5, color = 'firebrick')
#    
#    pylab.title(title, fontsize=8)
#    ax.xaxis.set_visible(False)
#    ax.yaxis.set_visible(False)
#    ax.set_aspect('equal')
#    return

    


t0 = time.time()
source = '7gon_3vertex_ 2.pkl'
#source = '7gon_3vertex_ 3.pkl'
#source = '7gon_3vertex_ 4.pkl'

test = PlanarLayout(file_path = source)
t_elapsed = time.time()- t0

print 'load time = ' + str(t_elapsed)


#find the flat band
flat_band = numpy.where(test.Es < -1.95*test.t)[0]

fb_states = test.Psis[:,flat_band]

#get all the semidual points
sdx, sdy = test.get_all_semidual_points()



#######
itt1 = 0
#az1 = False
az1 = True
aznum1 = 1
radnum1 = 4
##
if az1:
    num1 = aznum1
else:
    num1 = radnum1
######
#######
#itt2 = 0
#az2 = True
#aznum2 = 4
#radnum2 = 0
#
#if az2:
#    num2 = aznum2
#else:
#    num2 = radnum2
    
    
state0 = numpy.zeros(len(sdx))
source_ind = test.get_SDindex(num1, itt1, az1)
#target_ind = test.get_SDindex(num2, itt2, az2)


#print test.V_int(source_ind, target_ind, flat_band)

#map1 = test.V_int_map(source_ind)

#map1 = test.V_int_map(source_ind, flat_band)

###more_states = scipy.arange(0, 150, 1)
#more_states = scipy.arange(0, 50, 1)
#map1 = test.V_int_map(source_ind, more_states)

map1 = test.V_int_map(source_ind, [-4])

pylab.figure(1)
pylab.clf()
ax1 = pylab.subplot(1,1,1)
#pylab.plot(numpy.abs(map1))
test.plot_map_state(map1, ax1, title = 'interaction', colorbar = True, plot_links = True)
#test.plot_map_state(map1, ax1, title = 'interaction', colorbar = True, plot_links = True, autoscale = True)
pylab.scatter([test.SDx[source_ind]], [test.SDy[source_ind]], c =  'gold', s = 150, edgecolors = 'k')
pylab.show()







