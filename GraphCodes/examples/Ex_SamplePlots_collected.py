#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 15:54:53 2017

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

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


def single_site(xgrid, ygrid, x0, y0, sigma, amp):
    vals = -amp*numpy.exp(-((xgrid-x0)**2 + (ygrid-y0)**2)/sigma**2 )
    return vals

#planningFolderPath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/'
#if not planningFolderPath in sys.path:
#    sys.path.append(planningFolderPath)
##from LayoutGenerator4 import PlanarLayout  
from LayoutGenerator5 import PlanarLayout 

planningFolderPath = r'/volumes/ourphoton/Alicia/Layouts/HyperbolicPickles/'




numPoints = 200
#numSites = 5
#amp = 0.1
#sigma = 0.2

numSites = 5
amp = 0.1
sigma = 0.13

xs = numpy.linspace(-1.3, 1.3, numPoints)
ys = numpy.linspace(-1.3, 1.3, numPoints)

Xgrid, Ygrid = numpy.meshgrid(xs, ys)


regularLattice = numpy.zeros((numPoints, numPoints))
siteXs = numpy.linspace(-1,1, numSites)
siteYs = siteXs
siteGridX, siteGridY = numpy.meshgrid(siteXs, siteYs)



for indn in range(0, numSites):
    for indm in range(0, numSites):
        tempSite = single_site(Xgrid, Ygrid, siteXs[indn], siteYs[indm], sigma, amp)
        regularLattice = regularLattice+tempSite

#fig = plt.figure(1)
#fig.clf()
#ax = fig.gca(projection='3d')
#
#surf = ax.plot_surface(Xgrid, Ygrid, regularLattice, cmap=cm.afmhot,
#                       linewidth=0, antialiased=False)
#
#plt.title('regularLattice')
#plt.show()




noiseSigma = 0.20*2/numSites
#xnoise = numpy.random.normal(0,noiseSigma, numSites**2 )
#ynoise = numpy.random.normal(0,noiseSigma, numSites**2 )
#xnoise = (numpy.random.rand(numSites**2 )*2-1)*noiseSigma
#ynoise = (numpy.random.rand(numSites**2 )*2-1)*noiseSigma
xnoise = (2*numpy.floor(2*numpy.random.rand(numSites**2))   -1  )*noiseSigma
ynoise = (2*numpy.floor(2*numpy.random.rand(numSites**2))   -1  )*noiseSigma
siteGridX_irr = siteGridX + numpy.reshape(xnoise, [numSites, numSites])
siteGridY_irr = siteGridY + numpy.reshape(ynoise, [numSites, numSites])

irregularLattice = numpy.zeros((numPoints, numPoints))
for indn in range(0, numSites):
    for indm in range(0, numSites):
#        tempSite = single_site(Xgrid, Ygrid, siteXs[indn] + xnoise[ind], siteYs[indm]+ynoise[ind], sigma, amp)
        tempSite = single_site(Xgrid, Ygrid, siteGridX_irr[indn, indm], siteGridY_irr[indn, indm], sigma, amp)
        irregularLattice = irregularLattice+tempSite
        
#fig = plt.figure(2)
#fig.clf()
#ax = fig.gca(projection='3d')
#
##cmap=cm.YlGnBu_r
#surf = ax.plot_surface(Xgrid, Ygrid, irregularLattice, cmap=cm.YlGnBu_r,
#                       linewidth=0, antialiased=True)
#
#plt.title('irregular Lattice')
#plt.show()




def drawgraph(axis, mode = 'regular', linecolor = 'k'):
    if mode == 'regular':
        gridx = siteGridX
        gridy = siteGridY
    else:
        gridx = siteGridX_irr
        gridy = siteGridY_irr
        
    for indn in range(0, numSites):
        for indm in range(0, numSites):
            
            #draw vertical links
            if indn< numSites-1:
                xlims = [gridx[indn, indm], gridx[indn+1, indm]]
                ylims = [gridy[indn, indm], gridy[indn+1, indm]]
                axis.plot(xlims, ylims, color = linecolor)
            
            #draw horizontal links
            if indm < numSites-1:
                xlims = [gridx[indn, indm], gridx[indn, indm+1]]
                ylims = [gridy[indn, indm], gridy[indn, indm+1]]
                axis.plot(xlims, ylims, color = linecolor)
                
    for indn in range(0, numSites):
        for indm in range(0, numSites):
            axis.plot(gridx[indn, indm], gridy[indn, indm], marker = '.', color = 'r', markersize = 14)
    return






#mapp = cm.viridis_r
#mapp = cm.YlGnBu_r
mapp = cm.jet_r
fig = pylab.figure(3)
pylab.clf()
ax1 = pylab.subplot(2,2,1)
#pylab.pcolor(regularLattice, cmap = mapp)
pylab.imshow(numpy.flipud(regularLattice), cmap = mapp, interpolation = 'hanning')
ax1.set_aspect('equal')
pylab.title('regular lattice potential')
#ax1.axes.get_xaxis().set_visible(False)
#ax1.axes.get_yaxis().set_visible(False)
ax1.axis('off')

ax3 = pylab.subplot(2,2,3)
#pylab.pcolor(irregularLattice, cmap = mapp)
pylab.imshow(numpy.flipud(irregularLattice), cmap = mapp, interpolation = 'hanning')
#pylab.imshow(numpy.diag(numpy.linspace(1,10,10)), cmap = cm.viridis)
ax3.set_aspect('equal')
pylab.title('disordered lattice potential')
#ax3.axes.get_xaxis().set_visible(False)
#ax3.axes.get_yaxis().set_visible(False)
ax3.axis('off')


ax2 = pylab.subplot(2,2,2)
drawgraph(ax2, 'regular')
ax2.set_aspect('equal')
pylab.title('TB graph of regular lattice')
#ax2.axes.get_xaxis().set_visible(False)
#ax2.axes.get_yaxis().set_visible(False)
ax2.axis('off')

ax4 = pylab.subplot(2,2,4)
drawgraph(ax4, 'irregular')
ax4.set_aspect('equal')
pylab.title('TB graph of regular lattice')
#ax4.axes.get_xaxis().set_visible(False)
#ax4.axes.get_yaxis().set_visible(False)
ax4.axis('off')

pylab.show()









######
#simple hyperbolic layout lattice
######
#test = PlanarLayout(file_path = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/7gon_3vertex_ 1.pkl')
test = PlanarLayout(gon = 7, vertex = 3,side = 1, radius_method = 'exp')
test.populate(maxItter = 1)

pylab.figure(4)
pylab.clf()
ax = pylab.subplot(1,1,1)
ax.set_aspect('equal')
#test.draw_all_radials(ax, color = 'k', linewidth = 2)
#test.draw_all_azimuthals(ax, color = 'k', linewidth = 2)
test.draw_resonator_lattice(ax, color = 'k', linewidth = 2)
test.draw_resonator_end_points(ax, color = 'b', edgecolor = 'c', size = 150)



ax.axis('off')
pylab.show()








####
#plot eigenstates
####
pylab.figure(2)
pylab.clf()
ax = pylab.subplot(1,1,1)
#plot state
test.plot_layout_state(test.Psis[:,41], ax, title = 'state 41', colorbar = False, plot_links = False, cmap = 'Wistia')
#plot graph
test.draw_SDlinks(ax, 1, color = 'firebrick', linewidth = 0.5)

ax.axis('off')
pylab.show()











########
#spectrum versus itteration number
########
#test = PlanarLayout(file_path = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/7gon_3vertex_ 5.pkl')
#test = PlanarLayout(file_path = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/7gon_3vertex_ 4.pkl')
test = PlanarLayout(file_path = os.path.join(planningFolderPath, '7gon_3vertex_ 4.pkl'))
#test = PlanarLayout(file_path = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/7gon_3vertex_ 4HW.pkl')
#test = PlanarLayout(file_path = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/8gon_3vertex_ 3.pkl')
#test = PlanarLayout(file_path = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/4gon_4vertex_ 4.pkl')

#test = PlanarLayout(gon = 8, vertex = 3, modeType = 'FW', radius_method = 'lin')
#test.populate(maxItter = 4)

maxItter = test.itter



colors = {1:'r', 2:'b', 3:'deepskyblue', 4:'goldenrod', 5:'indigo'}

pylab.figure(5)
pylab.clf()
for itt in range(1, maxItter+1):
    [Evals, Psis, Eorder] = test.get_eigs(itt)
    xaxis = numpy.linspace(0,1, len(Eorder))
    ax = pylab.subplot(1,maxItter,itt)
    
#    pylab.plot(Evals, 'b.')
    pylab.plot(xaxis, Evals[Eorder], color = colors[itt], marker = '.', linestyle = '', label = str(itt), alpha = 1)
    pylab.xlabel('Eigenvalue number')
    pylab.ylabel('Energy (J)')
    pylab.title(str(itt))
  
pylab.suptitle('Eigenergy Spectrum v. System Size')
pylab.show()


pylab.figure(6)
pylab.clf()


[Evals, Psis, Eorder] = test.get_eigs()
EValsMaxItt = Evals[Eorder]
MaxXAxis = numpy.linspace(0,1, len(Eorder))
ax = pylab.subplot(1,1,1)
for itt in scipy.arange(maxItter, 0, -1):
    [Evals, Psis, Eorder] = test.get_eigs(itt)
    
    ####very simple way to plot curves together. Alignment is not qutite right
    xaxis = numpy.linspace(0,1, len(Eorder))
    
#    ###trying to align curves more carefully by changing domain. Smaller system sizes plotted
#    ###over smaller domain because tail gets clipped.
#    ###This lines up better, but not quite sure how to justify this exact alignment
#    ###I'm approximately compensating for some particle in a box cutoff at slow phase winding, but I'm not
#    ###sure that I've actually don it correctly.
#    Emax = Evals[Eorder[-1]]
#    if itt<maxItter:
#        endInd = numpy.where(Emax < EValsMaxItt)[0][0] #index in final data where the smaller data set ends
#        endVal = MaxXAxis[endInd]
#        xaxis = numpy.linspace(0,endVal, len(Eorder))
#    else:
#        xaxis = MaxXAxis
    
#    pylab.plot(xaxis, Evals[Eorder], color = colors[itt], marker = '.', linestyle = '', label = str(itt), alpha = 0.6)
    pylab.plot(xaxis, Evals[Eorder], color = colors[itt], marker = '.', linestyle = '', label = str(itt), alpha = 1)
#    pylab.plot(xaxis, Evals[Eorder], color = colors[itt], marker = '.', linestyle = '-', linewidth = 0.3, label = str(itt), alpha = 0.6)
    
pylab.xlabel('Normalized Eigenvalue Number')
pylab.ylabel('Energy (|t|)')
pylab.title('Eigenergy Spectrum v. System Size')
ax.legend(loc = 'upper left', title = 'System Size')

pylab.show()




######
#trying DOS comparison versus sytem size instead
######

#set up frequency sweep
freq_range = 4.01
freq_res = 0.04
#freq_res = 0.12

freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)

[fullDOS, bins_out] = numpy.histogram(test.Es, freq_bins)

pylab.figure(7)
pylab.clf()

ax1 = pylab.subplot(1,1,1)
for itt in scipy.arange(maxItter, 0, -1):
    [Evals, Psis, Eorder] = test.get_eigs(itt)

#    pylab.hist(Evals[Eorder], freq_bins, histtype='stepfilled', orientation='vertical', facecolor = colors[itt], color = 'red', alpha = 0.6, label = str(itt))
    
    [DOS, bins_out] = numpy.histogram(Evals, freq_bins)

    bins_centers = (bins_out[0:-1] + bins_out[1:])/2
    binWidth = bins_out[1] - bins_out[0]

#    print numpy.max(DOS)
#    pylab.bar(bins_centers, 1.*DOS, width = binWidth, color =  colors[itt], label = str(itt), alpha = 0.6)
#    pylab.bar(bins_centers, 1.*DOS, width = binWidth, color =  colors[itt], label = str(itt), alpha = 1)
    
    pylab.bar(bins_centers, 1.*DOS/len(Evals), width = binWidth, color =  colors[itt], label = str(itt), alpha = 1)
    
    
pylab.xlabel('Energy (J)')
pylab.ylabel('DOS')
ax1.legend(loc = 'upper right')
pylab.title('DOs v system size')
ax1.set_xlim([-2.5,4])
ax1.set_ylim([0,40])
ax1.set_ylim([0,0.05])

pylab.show()



#############
#eigenstate array plots
###########
##try to build the double plaquette state
#def build_double_plaquette_state(maxItter = -1):
#    if maxItter > test.itter:
#            raise ValueError, 'dont have this many itterations'
#    elif maxItter <0:
#        maxItter = test.itter
#    [xs,ys] = test.get_all_semidual_points()
#    
#    state = numpy.zeros(len(xs))*(0+0j)
#    
#    currentInd = 0
#    for itteration in range(0, maxItter+1):
#        numPoints = len(test.points[itteration])
#        numRadials = test.radials[itteration].shape[0]
#        
#        
#        #azimuthal points
##        print 'azimuthal ' + str(itteration)
#        if itteration == 0:
#            #zeroth ring
#            ring_state = scipy.arange(0,numPoints,1)
#            ring_state[1:]= numpy.mod(ring_state[1:],2)*2-1
#            
#            state[0:numPoints] = ring_state
#        if itteration == 1:
#            ring_state = scipy.arange(0,test.gon-3,1)
#            ring_state= 1*(numpy.mod(ring_state,2)*2-1)
#            state[currentInd:currentInd + len(ring_state)] = ring_state
#        currentInd = currentInd + numPoints
#        
#        #radial points
#        if itteration !=0:
##            print 'radial ' + str(itteration)
#            if itteration == 1:
#                #radials over to next plaquette
#                rad_state = numpy.asarray([1,-1])
#                state[currentInd: currentInd +2] = rad_state
#                
#            currentInd = currentInd + numRadials
#                
#    #normalize
#    state = state/numpy.linalg.norm(state)
#    return state

def build_double_plaquette_state(layout, maxItter = -1):
    if maxItter > layout.itter:
            raise ValueError, 'dont have this many itterations'
    elif maxItter <0:
        maxItter = layout.itter
    [xs,ys] = layout.get_all_semidual_points()
    
    state = numpy.zeros(len(xs))*(0+0j)
    
    currentInd = 0
    for itteration in range(0, maxItter+1):
        numPoints = len(layout.points[itteration])
        numRadials = layout.radials[itteration].shape[0]
        
        
        #azimuthal points
#        print 'azimuthal ' + str(itteration)
        if itteration == 0:
            #zeroth ring
            ring_state = scipy.arange(0,numPoints,1)
            ring_state[1:]= numpy.mod(ring_state[1:],2)*2-1
            
            state[0:numPoints] = ring_state
        if itteration == 1:
            ring_state = scipy.arange(0,test.gon-3,1)
            ring_state= 1*(numpy.mod(ring_state,2)*2-1)
            state[currentInd:currentInd + len(ring_state)] = ring_state
        currentInd = currentInd + numPoints
        
        #radial points
        if itteration !=0:
#            print 'radial ' + str(itteration)
            if itteration == 1:
                #radials over to next plaquette
                rad_state = numpy.asarray([1,-1])
                state[currentInd: currentInd +2] = rad_state
                
            currentInd = currentInd + numRadials
                
    #normalize
    state = state/numpy.linalg.norm(state)
    return state


test = PlanarLayout(file_path = planningFolderPath + '7gon_3vertex_ 3.pkl')


FBstate = build_double_plaquette_state(test)



#Something is way wierd with this and the saved file doesn't quite look exactly like the one in the figure. !!!!
#!!!!!!!
#Hence the wacky scale factor

#latticeItt = 3 #this for nice pictures, but will take a long time to plot
latticeItt = 1  #use this one for quick checks

scaleFactor = 1.2

fig = pylab.figure(8)
pylab.clf()
ax = pylab.subplot(3,3,1)
test.plot_layout_state(scaleFactor*test.Psis[:,398], ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')
test.draw_SDlinks(ax, latticeItt , color = 'firebrick', linewidth = 0.5)
ax.axis('off')

ax = pylab.subplot(3,3,2)
test.plot_layout_state(scaleFactor*test.Psis[:,397], ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')
test.draw_SDlinks(ax, latticeItt , color = 'firebrick', linewidth = 0.5)
ax.axis('off')

ax = pylab.subplot(3,3,3)
test.plot_layout_state(scaleFactor*test.Psis[:,394], ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')
test.draw_SDlinks(ax, latticeItt , color = 'firebrick', linewidth = 0.5)
ax.axis('off')



ax = pylab.subplot(3,3,4)
test.plot_layout_state(scaleFactor*test.Psis[:,142], ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')
test.draw_SDlinks(ax, latticeItt, color = 'firebrick', linewidth = 0.5)
ax.axis('off')

ax = pylab.subplot(3,3,5)
test.plot_layout_state(scaleFactor*test.Psis[:,114], ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')
test.draw_SDlinks(ax, latticeItt, color = 'firebrick', linewidth = 0.5)
ax.axis('off')

ax = pylab.subplot(3,3,6)
test.plot_layout_state(scaleFactor*test.Psis[:,85], ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')
test.draw_SDlinks(ax, latticeItt, color = 'firebrick', linewidth = 0.5)
ax.axis('off')


ax = pylab.subplot(3,3,8)
test.plot_layout_state(scaleFactor*FBstate*0.35, ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')
test.draw_SDlinks(ax, latticeItt, color = 'firebrick', linewidth = 0.5)
ax.axis('off')
cbar = pylab.colorbar(fraction=0.03, pad=0.04) 
cbar.set_label('phase ($\pi$ radians)', rotation=270, labelpad= 10)



fig.set_size_inches(9, 9)

#pylab.tight_layout()
pylab.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
pylab.show()

#save figure
#fig.savefig('try1.png' , dpi = 300)



########
#plot eigenspectra
########

test7 = PlanarLayout(file_path = planningFolderPath + '7gon_3vertex_ 3.pkl')

test6 = PlanarLayout(file_path = planningFolderPath + '6gon_3vertex_ 3.pkl')


fig1 = pylab.figure(9)
pylab.clf()
pylab.plot(test7.Es[test7.Eorder], 'b.')
pylab.ylabel('Energy (t)')
pylab.xlabel('Eigenvector')

fig1.set_size_inches(6, 7)
#pylab.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=1, hspace=1)
pylab.show()

#fig1.savefig('heptagonSpectrum.png' , dpi = 300)



fig1 = pylab.figure(10)
pylab.clf()
pylab.plot(test6.Es[test6.Eorder], 'b.')
pylab.ylabel('Energy (t)')
pylab.xlabel('Eigenvector')

fig1.set_size_inches(6, 7)
#pylab.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=1, hspace=1)
pylab.show()

#fig1.savefig('hexagonSpectrum.png' , dpi = 300)



#####
#show how itterations of the lattice are made and a nice plot of layout v effective
#######
test2 = PlanarLayout(gon = 7, vertex = 3,side = 1, radius_method = 'Mattias', file_path = '', modeType = 'FW')
#test2 = PlanarLayout(gon = 7, vertex = 3,side = 1, radius_method = 'lin', file_path = '', modeType = 'FW')
test2.populate(2)


#test3 = PlanarLayout(gon = 7, vertex = 3,side = 1, radius_method = 'lin', file_path = '', modeType = 'FW')
#test3.populate(3)
####Shell consruction
#fig3 = pylab.figure(11)
#pylab.clf()
#for ind in range(1,4):
#    ax = pylab.subplot(1,3,ind)
#    test3.draw_resonator_lattice(ax, mode = 'line', maxItter = ind, color = 'mediumblue', alpha = 1 , linewidth = 1.5)
#    test3.draw_SDlinks(ax, ind, color = 'deepskyblue', linewidth = 2.5, minus_links = False, minus_color = 'goldenrod')
#    sdx, sdy = test3.get_all_semidual_points(maxItter = ind)
##    pylab.plot(sdx, sdy,color =  'goldenrod', marker = 'o', linestyle = '')
#    pylab.scatter(sdx, sdy,c =  'goldenrod', marker = 'o', edgecolors = 'k', s = 5,  zorder=5)
#    
#    ax.set_aspect('equal')
#    ax.axis('off')
#    
#pylab.tight_layout()
#fig3.set_size_inches(13, 5)
#####fig3.savefig('ShellConstruction.png', dpi = 200)


#
#nice overlay of layout and effective lattices
fig4 = pylab.figure(12)
pylab.clf()
ind = 2
ax = pylab.subplot(1,1,1)
test2.draw_resonator_lattice(ax, mode = 'line', maxItter = ind, color = 'mediumblue', alpha = 1 , linewidth = 1.5)
test2.draw_SDlinks(ax, ind, color = 'deepskyblue', linewidth = 2.5, minus_links = False, minus_color = 'goldenrod')
sdx, sdy = test2.get_all_semidual_points(maxItter = ind)
pylab.scatter(sdx, sdy,c =  'goldenrod', marker = 'o', edgecolors = 'k', s = 5,  zorder=5)

ax.set_aspect('equal')
ax.axis('off')
pylab.tight_layout()
pylab.show()
fig4.set_size_inches(5, 5)
#fig4.savefig(r'/Users/kollar2/Desktop/Overlay_Both.png', dpi = 200, transparent=True)
#fig4.savefig(r'/Users/kollar2/Desktop/Mattias_Overlay_Both.png', dpi = 200, transparent=True)


#nice overlay of layout
fig5 = pylab.figure(13)
pylab.clf()
ind = 2
ax = pylab.subplot(1,1,1)
#test2.draw_resonator_lattice(ax, mode = 'line', maxItter = ind, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
#test2.draw_resonator_lattice(ax, mode = 'line', maxItter = ind, color = 'mediumblue', alpha = 1 , linewidth = 5)
#test2.draw_resonator_lattice(ax, mode = 'line', maxItter = ind, color = 'paleturquoise', alpha = 1 , linewidth = 5)
#test2.draw_resonator_lattice(ax, mode = 'line', maxItter = ind, color = 'white', alpha = 1 , linewidth = 5)
test2.draw_resonator_lattice(ax, mode = 'line', maxItter = ind, color = 'white', alpha = 1 , linewidth = 7)

xs = test2.coords[:,0]
ys = test2.coords[:,1]
pylab.sca(ax)
#pylab.scatter(xs, ys ,c =  'goldenrod', s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
#pylab.scatter(xs, ys ,c =  'royalblue', s = 40, marker = 'o', edgecolors = 'k', zorder = 5)
pylab.scatter(xs, ys ,c =  'royalblue', s = 60, marker = 'o', edgecolors = 'k', zorder = 5)

ax.set_aspect('equal')
ax.axis('off')
pylab.tight_layout()
pylab.show()
fig5.set_size_inches(5, 5)
#fig5.savefig(r'/Users/kollar2/Desktop/Overlay_Layout.png', dpi = 200, transparent=True)
#fig5.savefig(r'/Users/kollar2/Desktop/Mattias_Overlay_Layout.png', dpi = 200, transparent=True)
#fig5.savefig(r'/Users/kollar2/Desktop/Mattias_Overlay_Layout2.png', dpi = 200, transparent=True)


#nice overlay of effective
fig6 = pylab.figure(14)
pylab.clf()
ind = 2
ax = pylab.subplot(1,1,1)
#test2.draw_SDlinks(ax, ind, color = 'deepskyblue', linewidth = 2.5, minus_links = False, minus_color = 'goldenrod')
#test2.draw_SDlinks(ax, ind, color = 'dodgerblue', linewidth = 5, minus_links = False, minus_color = 'goldenrod')
test2.draw_SDlinks(ax, ind, color = 'dodgerblue', linewidth = 7, minus_links = False, minus_color = 'goldenrod')

sdx, sdy = test2.get_all_semidual_points(maxItter = ind)
#pylab.scatter(sdx, sdy,c =  'goldenrod', marker = 'o', edgecolors = 'k', s = 20,  zorder=5)
#pylab.scatter(sdx, sdy,c =  'lavender', marker = 'o', edgecolors = 'k', s = 40,  zorder=5)
pylab.scatter(sdx, sdy,c =  'lavender', marker = 'o', edgecolors = 'k', s = 60,  zorder=5)

ax.set_aspect('equal')
ax.axis('off')
pylab.tight_layout()
pylab.show()
fig6.set_size_inches(5, 5)
#fig6.savefig(r'/Users/kollar2/Desktop/Overlay_Effective.png', dpi = 200, transparent=True)
#fig6.savefig(r'/Users/kollar2/Desktop/Mattias_Overlay_Effective.png', dpi = 200, transparent=True)
#fig6.savefig(r'/Users/kollar2/Desktop/Mattias_Overlay_Effective2.png', dpi = 200, transparent=True)

#nice overlay of layout and effective lattices
fig7 = pylab.figure(15)
pylab.clf()
ind = 2
ax = pylab.subplot(1,1,1)
test2.draw_SDlinks(ax, ind, color = 'deepskyblue', linewidth = 1.5, minus_links = False, minus_color = 'goldenrod')
test2.draw_resonator_lattice(ax, mode = 'line', maxItter = ind, color = 'mediumblue', alpha = 1 , linewidth = 2.5)

#sdx, sdy = test2.get_all_semidual_points(maxItter = ind)
#pylab.scatter(sdx, sdy,c =  'goldenrod', marker = 'o', edgecolors = 'k', s = 5,  zorder=5)

xs = test2.coords[:,0]
ys = test2.coords[:,1]
pylab.sca(ax)
#pylab.scatter(xs, ys ,c =  'goldenrod', s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)
#pylab.scatter(xs, ys ,c =  'goldenrod', s = 40, marker = 'o', edgecolors = 'k', zorder = 5)

ax.set_aspect('equal')
ax.axis('off')
pylab.tight_layout()
pylab.show()
fig7.set_size_inches(5, 5)
#fig7.savefig(r'/Users/kollar2/Desktop/Overlay_Both2.png', dpi = 200, transparent=True)
#fig7.savefig(r'/Users/kollar2/Desktop/Overlay_Both2_3.png', dpi = 200, transparent=True)
#fig7.savefig(r'/Users/kollar2/Desktop/Mattias_Overlay_Both2.png', dpi = 200, transparent=True)






#####
#plot falt band state
#####
test4 = PlanarLayout(gon = 7, vertex = 3,side = 1, radius_method = 'lin', file_path = '', modeType = 'FW')
test4.populate(2)
#test4 = PlanarLayout(file_path = planningFolderPath + '7gon_3vertex_ 2.pkl') #don't use this old pickle. Cetner Az ring positions are a bit off.

plotFBstate = build_double_plaquette_state(test4)
#scaleFactor4 = 2.4
scaleFactor4 = 4.2


FBfig = pylab.figure(16)
pylab.clf()
ax = pylab.subplot(1,1,1)
test4.plot_layout_state(scaleFactor4*plotFBstate*0.35, ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')
test4.draw_SDlinks(ax, 2, color = 'firebrick', linewidth = 0.5)
ax.axis('off')

pylab.tight_layout()

FBfig.set_size_inches(4, 4)
#FBfig.savefig(r'/Users/kollar2/Desktop/FB_2.png', dpi = 200, transparent=True)


