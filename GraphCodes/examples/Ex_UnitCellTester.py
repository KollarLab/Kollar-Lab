#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 12:09:58 2018

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

import scipy.io as sio
from scipy import signal

pkgDir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

if not pkgDir in sys.path:
    sys.path.append(pkgDir)

import context

from GraphCodes.EuclideanLayoutGenerator2 import EuclideanLayout
from GraphCodes.EuclideanLayoutGenerator2 import UnitCell

#generalFolderPath = r'/Users/kollar2/Documents/HouckLab/GeneralLayoutCode/'
#if not generalFolderPath in sys.path:
#    sys.path.append(generalFolderPath)
    
from GraphCodes.GeneralLayoutGenerator import *


#pylab.rcParams.update({'font.size': 14})
#pylab.rcParams.update({'font.size': 8})

modeType = 'FW'
#modeType = 'HW'

#testCell = UnitCell('Huse')
#testCell = UnitCell('74Huse')
#testCell = UnitCell('84Huse')
#testCell = UnitCell('123Huse')
#testCell = UnitCell('kagome')

#testCell = UnitCell('Huse1_1')
#testCell = UnitCell('Huse2_1')
#testCell = UnitCell('Huse2_2')
#testCell = UnitCell('Huse3_1')
#testCell = UnitCell('Huse3_3')

#testCell = UnitCell('84Huse2_1')

#testCell = UnitCell('PeterChain')
#testCell = UnitCell('square')

#cell0 = UnitCell('kagome')
#res0 = cell0.resonators
#res1 = split_resonators(res0)
#testCell = UnitCell('split_kagome', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)


#cell0 = UnitCell('square')
#res0 = cell0.resonators
#res1 = split_resonators(res0)
#testCell = UnitCell('split_square', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)



########split Euclidean, hoffman attemps
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#resonators0 = test1.resonators #graphene layout
#splitGraph = split_resonators(resonators0) #split graphene layout
#resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
#resonators2 = split_resonators(resonators1) #split further
#resonators = resonators2
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
##sites = [80,81,82, 83,79,78,77,76, 71, 84,85,86,87,114,115, 116, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
#sites = [80,81,83,79,78,77,76, 71, 50, 49, 48, 47, 46, 44, 55,54, 53,52] #DO NOT MODIFY !!!!!
#state = numpy.zeros(len(testLattice.SDx))
#state[sites] = 0.1
#newCell = testLattice.resonators[sites, :]
##fig1 = pylab.figure(114)
##pylab.clf()
##ax = pylab.subplot(1,1,1)
##testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
##xs = testLattice.coords[:,0]
##ys = testLattice.coords[:,1]
##pylab.sca(ax)
##pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)
##pylab.scatter(testLattice.SDx[sites], testLattice.SDy[sites], c =  'firebrick', s = 75, marker = 'o', edgecolors = 'maroon', zorder = 5,linewidth = 1)
##ax.set_aspect('equal')
##ax.axis('off')
##pylab.title(testLattice.name)
##pylab.tight_layout()
##pylab.show()
#testCell = UnitCell('extremal_hofmann_attempt', resonators = newCell, a1 = test1.unitcell.a1, a2 = test1.unitcell.a2)



###failed Hofman unit cell attempt. I think it's not ablt to add the necessary resonators
#cell0 = UnitCell('kagome')
#res0 = cell0.resonators
#res1 = split_resonators(res0) #these are now a split graphene layout
#res2 = generate_line_graph(res1) #now this is McLaughlin-like kagome
#res3 = split_resonators(res2)#split further
#res = res3
#testCell = UnitCell('extremal_hofmann_attempt', resonators = res3, a1 = cell0.a1, a2 = cell0.a2)

cell0 = UnitCell('kagome')
decCell = UnitCell('PeterChain_tail', side = 1)
res0 = cell0.resonators
res1 = decorate_layout(res0, decCell.resonators)
testCell = UnitCell('Peter_graphene', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)


#cell0 = UnitCell('Huse')
#res0 = cell0.resonators
#res1 = split_resonators(res0)
#testCell = UnitCell('split_HPG', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)

#cell0 = UnitCell('Huse')
#decCell = UnitCell('PeterChain_tail', side = 1)
#res0 = cell0.resonators
#res1 = decorate_layout(res0, decCell.resonators)
#testCell = UnitCell('Peter_HPG', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)


#cell0 = UnitCell('kagome')
#res0 = cell0.resonators
#res1 = generate_line_graph(res0,)
#testCell = UnitCell('LG_kagome', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)

#testCell = UnitCell('square')


plotAllBands = False

######
#test the unit cell
#######
pylab.figure(1)
pylab.clf()
ax = pylab.subplot(1,2,1)
testCell.draw_sites(ax)
pylab.title('Sites of Huse Cell')

ax = pylab.subplot(1,2,2)
testCell.draw_sites(ax,color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 20)
testCell.draw_resonators(ax, color = 'cornflowerblue', linewidth = 1)
testCell.draw_SDlinks(ax, color = 'firebrick', linewidth = 1)
testCell.draw_resonator_end_points(ax, color = 'deepskyblue', edgecolor = 'k',  marker = 'o' , size = 20)
pylab.title('Links of Unit Cell')
pylab.show()


######
#show the orientations
######
#alternate version
fig = pylab.figure(2)
pylab.clf()
ax = pylab.subplot(1,1,1)
testCell.draw_resonators(ax, color = 'cornflowerblue', linewidth = 1)
testCell.draw_resonator_end_points(ax, color = 'indigo', edgecolor = 'indigo',  marker = '+' , size = 20)
testCell.draw_site_orientations(ax, title = 'unit cell convention', colorbar = False, plot_links = False, cmap = 'jet', scaleFactor = 0.5)
testCell.draw_SDlinks(ax, linewidth = 1.5, HW = True , minus_color = 'goldenrod')
pylab.title('site orientations : ' + testCell.type)
#ax.set_aspect('auto')
ax.set_aspect('equal')
#    fig.savefig('HW.png', dpi = 200)

pylab.show()







#####
#testing bloch theory
####

Hmat = testCell.generate_Bloch_matrix(0,0,  modeType = modeType)
pylab.figure(3)
pylab.clf()
ax = pylab.subplot(1,2,1)
pylab.imshow(numpy.abs(Hmat))
pylab.colorbar()
pylab.title('|H|')

ax = pylab.subplot(1,2,2)
pylab.imshow(numpy.real(Hmat - numpy.transpose(numpy.conj(Hmat))))
pylab.title('H - Hdagger')

pylab.show()



#kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = 100, modeType = modeType)
#kx_y, ky_y, cuty = testCell.compute_band_structure(0, -8./3*numpy.pi, 0, 8./3*numpy.pi, numsteps = 100, modeType = modeType)
kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = 100, modeType = modeType)
kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = 100, modeType = modeType)
#kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = 100, modeType = modeType)
#kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2*numpy.pi, 0, 2*numpy.pi, numsteps = 100, modeType = modeType)


fig2 = pylab.figure(4)
pylab.clf()
ax = pylab.subplot(1,2,1)
testCell.plot_band_cut(ax, cutx)
pylab.title('xcut')
#ax.set_ylim([-2.1,4.05])

ax = pylab.subplot(1,2,2)
testCell.plot_band_cut(ax, cuty)
pylab.title('ycut')
#ax.set_ylim([-2.1,4.05])

titleStr = testCell.type + ', modeType: ' + modeType + ' (Made with UnitCell class)' 
pylab.suptitle(titleStr)

pylab.show()



#########
#surface plots
########
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

numSurfPoints = 50
xs = numpy.linspace(-2*numpy.pi, 2*numpy.pi, numSurfPoints)
ys = numpy.linspace(-2.5*numpy.pi, 2.5*numpy.pi, numSurfPoints)
#ys = numpy.linspace(-0.5*numpy.pi, 0.5*numpy.pi, numSurfPoints)

Xgrid, Ygrid = numpy.meshgrid(xs, ys)

bands = numpy.zeros((testCell.numSites, len(xs), len(ys)))

#####compute surfaces
for xind in range(0,len(xs)):
    for yind in range(0,len(ys)):
        xval = xs[xind]
        yval = ys[yind]
        kx,yk, Es = testCell.compute_band_structure(xval, yval, xval, yval, numsteps = 1, modeType = 'FW')
        bands[:, xind, yind] = numpy.transpose(Es)

####plot all surfaces together
fig = plt.figure(17)
fig.clf()
ax = fig.gca(projection='3d')
for ind in range(0,testCell.numSites):
    surf = ax.plot_surface(Xgrid, Ygrid, bands[ind,:,:], cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.set_zlim(-4, 4)
#    ax.set_zlim(-4*J, 4*J)

plt.title('dispersion curves')
plt.show()

###plot each surface
if plotAllBands:
    for ind in range(0,testCell.numSites):
        fig = plt.figure(40+ind)
        fig.clf()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(Xgrid, Ygrid, bands[ind,:,:], cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(-4, 4)
        plt.title(str(ind))
        plt.show()







######
##look at above gap state at k= 0
######
#Es, Psis = scipy.linalg.eigh(Hmat)
#
#stateInd = 0
#aboveGap = Psis[:,stateInd]
#print Es[stateInd]
#print aboveGap
#
#pylab.figure(5)
#pylab.clf()
#
#ax = pylab.subplot(1,1,1)
##testCell.draw_sites(ax,color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 20)
#testCell.draw_SDlinks(ax, color = 'firebrick', linewidth = 1)
#testCell.draw_resonators(ax, color = 'cornflowerblue', linewidth = 1)
#testCell.draw_resonator_end_points(ax, color = 'deepskyblue', edgecolor = 'k',  marker = 'o' , size = 20)
##testCell.plot_bloch_wave(aboveGap*2, ax, title = 'state weight', colorbar = False, plot_links = False, cmap = 'Wistia')
#temp = testCell.plot_bloch_wave_end_state(aboveGap*2, ax,modeType = modeType,  title = modeType + '_' + str(stateInd), colorbar = False, plot_links = False, cmap = 'Wistia')
#ax.set_aspect('equal')
#pylab.show()
#
#
#####try to plot all the unit cell wave functions. Doesn't work very well. You can't see anything
##pylab.figure(6)
##pylab.clf()
##for ind in range(0, testCell.numSites):
##    ax = pylab.subplot(1,testCell.numSites,ind+1)
##    testCell.draw_SDlinks(ax, color = 'firebrick', linewidth = 1)
##    testCell.draw_resonators(ax, color = 'cornflowerblue', linewidth = 1)
##    testCell.draw_resonator_end_points(ax, color = 'deepskyblue', edgecolor = 'k',  marker = 'o' , size = 20)
###    testCell.plot_bloch_wave(Psis[:,ind], ax, title = 'state weight', colorbar = False, plot_links = False, cmap = 'Wistia')
##    testCell.plot_bloch_wave_end_state(Psis[:,ind], ax,modeType = modeType,  title = str(ind), colorbar = False, plot_links = False, cmap = 'Wistia')
##    ax.set_aspect('equal')
##pylab.show()





######
#craxy flux lattice analogy
######
def get_FLMat(theta, kx= 0 , ky = 0):
    t = -1#flip the sign of the tunnel couplings
    startMat = fluxLatticeCell.generate_Bloch_matrix(kx, ky, modeType = 'FW', t = t)+ (0 +0j)
    
    phaseFactor = numpy.exp(1j*theta)
    
    #put in the fluxes
    startMat[1,7] = phaseFactor*startMat[1,7]
    startMat[7,1] = numpy.conj(phaseFactor)*startMat[7,1]
    
    startMat[7,9] = phaseFactor*startMat[7,9]
    startMat[9,7] = numpy.conj(phaseFactor)*startMat[9,7]
    
    startMat[8,3] = phaseFactor*startMat[8,3]
    startMat[3,8] = numpy.conj(phaseFactor)*startMat[3,8]
    
    startMat[9,11] = phaseFactor*startMat[9,11]
    startMat[11,9] = numpy.conj(phaseFactor)*startMat[11,9]
    
    startMat[3,10] = phaseFactor*startMat[3,10]
    startMat[10,3] = numpy.conj(phaseFactor)*startMat[10,3]
    
    startMat[11,5] = phaseFactor*startMat[11,5]
    startMat[5,11] = numpy.conj(phaseFactor)*startMat[5,11]
    
    startMat[5,6] = phaseFactor*startMat[5,6]
    startMat[6,5] = numpy.conj(phaseFactor)*startMat[6,5]
    
    startMat[1,6] = phaseFactor*startMat[1,6]
    startMat[6,1] = numpy.conj(phaseFactor)*startMat[6,1]
    
    FLMat = startMat
    
    return FLMat

fluxLatticeCell = UnitCell('Huse')

kx = 0
ky = 0

#kx = 1.
#ky = 0

#kx = 0.
#ky = 1.



FWMat = fluxLatticeCell.generate_Bloch_matrix(kx,ky, modeType = 'FW', t = 1)
HWMat = fluxLatticeCell.generate_Bloch_matrix(kx,ky, modeType = 'HW', t = 1)
#FLMat = get_FLMat(numpy.pi)

EsFW,PsisFW = scipy.linalg.eigh(FWMat)
EsHW,PsisHW = scipy.linalg.eigh(HWMat)

#EsFL,PsisFL = scipy.linalg.eigh(FLMat)

pylab.figure(7)
pylab.clf()
xs = numpy.zeros(len(EsFW))
pylab.scatter(xs, EsFW, c =  'firebrick', s = 75, edgecolors = 'k')
pylab.scatter(xs +1, EsHW, c =  'deepskyblue', s = 75, edgecolors = 'k')

numSteps = 100
steps = numpy.linspace(0.0, numpy.pi + numpy.pi/numSteps, numSteps)
for ind in range(0,len(steps)):
    step = steps[ind]
    offset = ind*1./(len(steps)-1)
    FLMat = get_FLMat(step, kx, ky)
    EsFL,PsisFL = scipy.linalg.eigh(FLMat)
    pylab.plot(xs +offset, EsFL, color = 'cornflowerblue' , marker = '.', markersize = '5', linestyle = '')
    
pylab.xlabel('flux (pi)')
pylab.ylabel('eigen energies')
pylab.title('Evolution of flux lattice: kx = ' + str(kx) + ' , ky = ' +str(ky))
pylab.show()


#pylab.figure(8)
#pylab.clf()
#numSubPlots = 5
#indStep = int(len(steps)/(numSubPlots-1))
#spind = 1
#ind = 0
#while ind< len(steps):
#    step = steps[ind]
#    offset = ind*1./(len(steps)-1)
#    FLMat = get_FLMat(step, kx, ky)
#    EsFL,PsisFL = scipy.linalg.eigh(FLMat)
#    
#    ax = pylab.subplot(1, numSubPlots, spind)
#    titleStr = 'angle: '+  str(numpy.round(step,2))
#    fluxLatticeCell.plot_bloch_wave(PsisFL[:,stateInd]*2, ax, title =titleStr, colorbar = False, plot_links = True, cmap = 'Wistia')    
#    ind = ind+ indStep
#    spind = spind +1
#    
#ax = pylab.subplot(1, numSubPlots, numSubPlots)
#titleStr = 'angle: '+  str(numpy.round(numpy.pi,2))
#fluxLatticeCell.plot_bloch_wave(PsisHW[:,stateInd]*2, ax, title =titleStr, colorbar = False, plot_links = True, cmap = 'Wistia')  
#
#pylab.suptitle('FL state evolution: ' +str(stateInd) + ' : kx = ' + str(kx)+ ' , ky = ' +str(ky))
#pylab.show()










###########
###plot for poster
##########
#pylab.rcParams.update({'font.size': 14})
#
#fig2 = pylab.figure(14)
#pylab.clf()
#ax = pylab.subplot(1,2,1)
##plot_band_cut_mc(ax, cutx)
#testCell.plot_band_cut(ax, cutx)
#pylab.title('xcut')
#pylab.ylabel('Energy (t)')
#pylab.xlabel('$k_x$ ($\pi$/a)')
#pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2,0,2], rotation='horizontal')
#
#ax = pylab.subplot(1,2,2)
##plot_band_cut_mc(ax, cuty)
#testCell.plot_band_cut(ax, cuty)
#pylab.title('ycut')
#pylab.ylabel('Energy (t)')
#pylab.xlabel('$k_y$ ($\pi$/a)')
#pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
#
#titleStr = testCell.type + ', modeType: ' + modeType + ' (Made with UnitCell class)' 
#pylab.suptitle(titleStr)
#
#pylab.show()
#
####fig2.savefig('HL_BS_1.png', dpi = 200)


