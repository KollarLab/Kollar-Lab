#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 15:53:30 2018

@author: kollar2
"""

import pylab
import numpy
from scipy import signal
import scipy
import pickle
import sys
import os.path
import matplotlib.gridspec as gridspec

#hyperbolicFolderPath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/'
#euclideanFolderPath = r'/Users/kollar2/Documents/HouckLab/EuclideanLatticePlanning/'
#generalFolderPath = r'/Users/kollar2/Documents/HouckLab/GeneralLayoutCode/'
#if not hyperbolicFolderPath in sys.path:
#    sys.path.append(hyperbolicFolderPath)
#if not euclideanFolderPath in sys.path:
#    sys.path.append(euclideanFolderPath)
#if not generalFolderPath in sys.path:
#    sys.path.append(generalFolderPath)
    
    
#FunctionFolderPath = r'/home/pranav/PhotonicLattices'
#DataPickleFolderPath = r'/volumes/ourphoton/Alicia/Layouts/HyperbolicPickles'
#if not FunctionFolderPath in sys.path:
#    sys.path.append(FunctionFolderPath)

pkgDir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

if not pkgDir in sys.path:
    sys.path.append(pkgDir)

import context

from GraphCodes.GeneralLayoutGenerator import GeneralLayout
from GraphCodes.GeneralLayoutGenerator import TreeResonators

from GraphCodes.EuclideanLayoutGenerator2 import UnitCell
from GraphCodes.EuclideanLayoutGenerator2 import EuclideanLayout

from GraphCodes.LayoutGenerator5 import PlanarLayout


from GraphCodes.GeneralLayoutGenerator import split_resonators
from GraphCodes.GeneralLayoutGenerator import generate_line_graph
from GraphCodes.GeneralLayoutGenerator import decorate_layout



#############
#defaults
##########
    
    
bigCdefault = 110
smallCdefault = 30

layoutLineColor = 'mediumblue'
layoutCapColor = 'goldenrod'

FWlinkAlpha = 0.7
FWsiteAlpha = 0.6

HWlinkAlpha = 0.8
HWsiteAlpha = 0.6

FWlinkColor = 'dodgerblue'
FWsiteColor = 'lightsteelblue'
FWsiteEdgeColor = 'mediumblue'

HWlinkColor = 'lightsteelblue'
HWminusLinkColor = 'b'
HWsiteColor = 'lightgrey'
HWsiteEdgeColor = 'midnightblue'

stateColor1 = 'gold'
stateEdgeColor1 = 'darkgoldenrod'
stateEdgeWidth1 = 1.25

stateColor2 = 'firebrick'
stateEdgeColor2 = 'maroon'
stateEdgeWidth2 = 1

DefaultFig = True
if DefaultFig:
    fig444 = pylab.figure(444)
    pylab.clf()
    
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    kagomeLattice = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'HW')
    kagomeLatticeHW = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    
    red = [3,5,13,18,20,21,24,32,33,35]
    yellow = [4,8,9,19,22,25,31,36,34,37]
    
    red2 = [9,12,15,0,18,27,34]
    yellow2 = [10,13,16,5,14,23]
    
    ####layout graph
    ax1 = pylab.subplot(1,5,1, adjustable='box', aspect=1)
    pylab.cla()
    kagomeLattice.draw_resonator_lattice(ax1, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    #testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
    xs = kagomeLattice.coords[:,0]
    ys = kagomeLattice.coords[:,1]
    pylab.sca(ax1)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    #ax.set_aspect('equal')
    ax1.axis('off')
    pylab.title('layout')
    
    ax2 = pylab.subplot(1,5,2, adjustable='box', aspect=1)
    pylab.cla()
    kagomeLattice.draw_SDlinks(ax2, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=1)
    pylab.sca(ax2)
    pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
    pylab.title('FW lattice')
    ax2.set_aspect('equal')
    ax2.axis('off')
    
    ax3 = pylab.subplot(1,5,3, adjustable='box', aspect=1)
    pylab.cla()
    kagomeLatticeHW.draw_SDlinks(ax3, color = HWlinkColor, linewidth = 3, minus_links = True, minus_color = HWminusLinkColor, alpha=1)
    pylab.sca(ax3)
    pylab.scatter(kagomeLatticeHW.SDx, kagomeLatticeHW.SDy ,c =  HWsiteColor, s = smallCdefault, marker = 'o', edgecolors = HWsiteEdgeColor, zorder = 5, alpha=1)
    pylab.title('HW lattice')
    ax3.set_aspect('equal')
    ax3.axis('off')
    
    ax4 = pylab.subplot(1,5,4, adjustable='box', aspect=1)
    pylab.cla()
    kagomeLattice.draw_SDlinks(ax4, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=FWlinkAlpha)
    pylab.sca(ax4)
    pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=FWsiteAlpha)
    pylab.scatter(kagomeLattice.SDx[yellow], kagomeLattice.SDy[yellow], c =  stateColor1, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(kagomeLattice.SDx[red], kagomeLattice.SDy[red], c =  stateColor2, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    pylab.title('FW state')
    ax4.set_aspect('equal')
    ax4.axis('off')
    
    ax5 = pylab.subplot(1,5,5, adjustable='box', aspect=1)
    pylab.cla()
    kagomeLatticeHW.draw_SDlinks(ax5, color = HWlinkColor, linewidth = 3, minus_links = True, minus_color = HWminusLinkColor, alpha=HWlinkAlpha)
    pylab.sca(ax5)
    pylab.scatter(kagomeLatticeHW.SDx, kagomeLatticeHW.SDy ,c =  HWsiteColor, s = smallCdefault, marker = 'o', edgecolors = HWsiteEdgeColor, zorder = 5, alpha=HWsiteAlpha)
    pylab.scatter(kagomeLatticeHW.SDx[yellow], kagomeLatticeHW.SDy[yellow], c =  stateColor1, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(kagomeLatticeHW.SDx[red], kagomeLatticeHW.SDy[red], c =  stateColor2, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    pylab.title('HW state')
    ax5.set_aspect('equal')
    ax5.axis('off')
    
#    pylab.suptitle('Defaults')
    fig444.set_size_inches(14.5, 3.5)
    pylab.tight_layout()
else:
    pass



###############
#make a unit cell
cell = UnitCell('PeterChain')

cell.a1
cell.a2
cell.resonators #(x0,y0)  (x1,y1)

pylab.figure(1)
pylab.clf()

ax = pylab.subplot(1,3,1)
cell.draw_resonators(ax, linewidth = 2, color = 'mediumblue')
cell.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 80)
#cell.draw_site_orientations(ax)
ax.set_aspect('equal')
ax.axis('off')
pylab.title('resonators')

ax = pylab.subplot(1,3,2)
cell.draw_SDlinks(ax, color = 'firebrick', linewidth = 2.5, minus_color = 'goldenrod')
ax.set_aspect('equal')
ax.axis('off')
pylab.title('links')


ax = pylab.subplot(1,3,3)
cell.draw_resonators(ax, linewidth = 2, color = 'mediumblue')
cell.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 80)
cell.draw_SDlinks(ax, color = 'firebrick', linewidth = 2.5, minus_color = 'goldenrod')
#cell.draw_site_orientations(ax)
ax.set_aspect('equal')
ax.axis('off')
pylab.title('links')

pylab.show()



###################
#make a 1D lattice (not periodic boundary conditions)
testLattice = EuclideanLayout(4,1,lattice_type = 'PeterChain', modeType = 'FW')
#testLattice = EuclideanLayout(3,1,lattice_type = 'PeterChain', modeType = 'HW')
genlattice = GeneralLayout(testLattice.resonators, modeType = testLattice.modeType)




pylab.figure(2)
pylab.clf()

ax = pylab.subplot(1,3,1)
testLattice.draw_resonator_lattice(ax, linewidth = 2, color = 'mediumblue')
testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 80)
#testLattice.draw_site_orientations(ax)
#ax.set_aspect('equal')
ax.axis('off')
pylab.title('resonators')


ax = pylab.subplot(1,3,2)
testLattice.draw_resonator_lattice(ax, linewidth = 2, color = 'mediumblue')
#testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 80)
testLattice.draw_SDlinks(ax, color = 'firebrick', linewidth = 2.5, minus_color = 'goldenrod', minus_links = True)
testLattice.plot_end_layout_state(0.35*numpy.ones(len(testLattice.SDx)), ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia', scaleFactor = 0.5)
ax.axis('off')
pylab.title('links')


ax = pylab.subplot(1,3,3)
testLattice.draw_SDlinks(ax, color = 'firebrick', linewidth = 0.5, minus_color = 'goldenrod', minus_links = True)
#testLattice.plot_layout_state(0.5*numpy.ones(len(testLattice.SDx)), ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')

state = testLattice.Psis[:,-2]
testLattice.plot_layout_state(state, ax, title = '', colorbar = False, plot_links = False, cmap = 'Wistia')


ax.axis('off')

pylab.tight_layout()
pylab.show()







##############
#compute a simple band structure
cellModeType = 'FW'

cell2  = UnitCell('PeterChain_tail', side = 1)


#band structure
numSurfPoints = 300
kx_x, ky_y, cutx = cell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
#kx_y, ky_y, cuty = cell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
kx_x, ky_y, cutx2 = cell2.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
#kx_y, ky_y, cuty2 = cell2.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)



pylab.figure(3)
pylab.clf()

ax = pylab.subplot(1,4,1)
cell.draw_resonators(ax, linewidth = 2, color = 'mediumblue')
cell.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 80)
#cell.draw_site_orientations(ax)
ax.set_aspect('equal')
ax.axis('off')
pylab.title('resonators')


ax = pylab.subplot(1,4,2)
pylab.cla()
cell.plot_band_cut(ax, cutx)
pylab.title('')
pylab.ylabel('Energy (|t|)')
pylab.xlabel('$k_x$ ($\pi$/a)')
pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')





ax = pylab.subplot(1,4,3)
cell2.draw_resonators(ax, linewidth = 2, color = 'mediumblue')
cell2.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 80)
#cell.draw_site_orientations(ax)
ax.set_aspect('equal')
ax.axis('off')
pylab.title('resonators')



ax = pylab.subplot(1,4,4)
pylab.cla()
cell.plot_band_cut(ax, cutx2)
pylab.title('')
pylab.ylabel('Energy (|t|)')
pylab.xlabel('$k_x$ ($\pi$/a)')
pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')


pylab.suptitle('Basic Band Structures')
pylab.show()







###############
#decorate and compute
##########

squareLattice = EuclideanLayout(4,4,lattice_type = 'square', modeType = 'FW')
#squareLattice = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'FW')

squareResonators = squareLattice.resonators
peterResonators = cell.resonators
peterResonators2 = cell2.resonators


decorated_resonators = decorate_layout(squareResonators, peterResonators2)
resonators = decorated_resonators
decoratedLattice2 = GeneralLayout(resonators , modeType = squareLattice.modeType, name =  'Peter2_square')

decorated_resonators = decorate_layout(squareResonators, peterResonators)
resonators = decorated_resonators
decoratedLattice = GeneralLayout(resonators , modeType = squareLattice.modeType, name =  'Peter_square')


testLattice = decoratedLattice2




#need a new cell
decorated_cell_resonators = decorate_layout(squareLattice.unitcell.resonators, peterResonators2)
#peterSquareCell = UnitCell('peterSquare',resonators = decorated_cell_resonators, a1 = cell2.a1, a2 = cell2.a2)
peterSquareCell = UnitCell(resonators = decorated_cell_resonators, a1 = squareLattice.unitcell.a1, a2 = squareLattice.unitcell.a2, lattice_type = 'peterSquare')




##cell0 = UnitCell('kagome')
##cell0 = UnitCell('square')
#cell0 = squareLattice.unitcell
#decCell = UnitCell('PeterChain_tail', side = 1)
#res0 = cell0.resonators
##res1 = decorate_layout(res0, decCell.resonators)
##res1 = decorate_layout(res0, peterResonators2)
#res1 = decorate_layout(squareLattice.unitcell.resonators, peterResonators2)
#testCell = UnitCell('Peter_graphene', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
#peterSquareCell = testCell



#band structure
numSurfPoints = 300
#kx_x, ky_y, cutx = peterSquareCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
#kx_y, ky_y, cuty = peterSquareCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)

#kx_x, ky_y, cutx = peterSquareCell.compute_band_structure(-2*numpy.pi, numpy.pi, 2*numpy.pi, numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
#kx_y, ky_y, cuty = peterSquareCell.compute_band_structure(-numpy.pi, -2.5*numpy.pi, numpy.pi, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)

kx_x, ky_y, cutx = peterSquareCell.compute_band_structure(-2*numpy.pi, -2*numpy.pi, 2*numpy.pi, 2*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
kx_y, ky_y, cuty = peterSquareCell.compute_band_structure(-2*numpy.pi, 2*numpy.pi, 2*numpy.pi, -2*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)



pylab.figure(4)
pylab.clf()

ax = pylab.subplot(1,4,1)
pylab.cla()
squareLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
xs = squareLattice.coords[:,0]
ys = squareLattice.coords[:,1]
pylab.sca(ax)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#ax.set_aspect('equal')
ax.axis('off')
pylab.title('square')


ax = pylab.subplot(1,4,2)
pylab.cla()
testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#peterSquareCell.draw_resonators(ax)
#testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
xs = testLattice.coords[:,0]
ys = testLattice.coords[:,1]
pylab.sca(ax)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#ax.set_aspect('equal')
ax.axis('off')
pylab.title('decorated square')



ax = pylab.subplot(1,4,3)
pylab.cla()
peterSquareCell.plot_band_cut(ax, cutx)
pylab.title('')
pylab.ylabel('Energy (|t|)')
pylab.xlabel('$k_x$ ($\pi$/a)')
pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
#ax.set_ylim([-2.5, 4])


ax = pylab.subplot(1,4,4)
pylab.cla()
peterSquareCell.plot_band_cut(ax, cuty)
pylab.title('')
pylab.ylabel('Energy (|t|)')
pylab.xlabel('$k_y$ ($\pi$/a)')
pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
#ax.set_ylim([-2.5, 4])



pylab.show()





