#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 16:48:15 2018

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

hyperbolicFolderPath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/'
euclideanFolderPath = r'/Users/kollar2/Documents/HouckLab/EuclideanLatticePlanning/'
generalFolderPath = r'/Users/kollar2/Documents/HouckLab/GeneralLayoutCode/'
if not hyperbolicFolderPath in sys.path:
    sys.path.append(hyperbolicFolderPath)
if not euclideanFolderPath in sys.path:
    sys.path.append(euclideanFolderPath)
if not generalFolderPath in sys.path:
    sys.path.append(generalFolderPath)
    
    
#FunctionFolderPath = r'/home/pranav/PhotonicLattices'
DataPickleFolderPath = r'/volumes/ourphoton/Alicia/Layouts/HyperbolicPickles'
#if not FunctionFolderPath in sys.path:
#    sys.path.append(FunctionFolderPath)
   
from GeneralLayoutGenerator import GeneralLayout
from GeneralLayoutGenerator import TreeResonators

from EuclideanLayoutGenerator2 import UnitCell
from EuclideanLayoutGenerator2 import EuclideanLayout

from LayoutGenerator5 import PlanarLayout


from GeneralLayoutGenerator import split_resonators
from GeneralLayoutGenerator import generate_line_graph
#from GeneralLayoutGenerator import decorate_layout



cell0 = UnitCell('kagome')
res0 = cell0.resonators
res1 = split_resonators(res0)
#testCell = UnitCell('split_kagome', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
testCell = cell0.split_cell(name = 'kaphene', splitIn = 2)

kapheneCell = cell0.split_cell(name = 'kaphene', splitIn = 2)
newRes = kapheneCell.line_graph_cell(name = 'Charon', resonatorsOnly = True)
testCell_auto = kapheneCell.line_graph_cell(name = 'Charon', resonatorsOnly = False)

tempCell = UnitCell(lattice_type = 'tempCell', side = 1, resonators = newRes, a1 =kapheneCell.a1, a2 = kapheneCell.a2)
#tempLattice = EuclideanLayout(xcells = 4, ycells = 4, modeType = 'FW', resonatorsOnly=False, initialCell = tempCell )
tempLattice = EuclideanLayout(xcells = 4, ycells = 4, modeType = 'FW', resonatorsOnly=False, initialCell = testCell_auto)



#1) extremal hofman
########split Euclidean, hoffman attemps
##!!!!!!!!needs to be this size for the cell to come out right. DO NOT CHANGE!!!
test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
resonators0 = test1.resonators #graphene layout
splitGraph = split_resonators(resonators0) #split graphene layout
resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome, kaphene
resonators2 = split_resonators(resonators1) #split further
resonators = resonators2
testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
#sites = [80,81,82, 83,79,78,77,76, 71, 84,85,86,87,114,115, 116, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
sites = [80,81,83,79,78,77,76, 71, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
state = numpy.zeros(len(testLattice.SDx))
state[sites] = 0.1
newCell = testLattice.resonators[sites, :]
testCell = UnitCell('extremal_hofmann', resonators = newCell, a1 = test1.unitcell.a1, a2 = test1.unitcell.a2)


#make all the cells that I want
kagomeCell = UnitCell('kagome') #graphene layout, kagome effective
kapheneCell = kagomeCell.split_cell(name = 'kaphene', splitIn = 2) #split graphene layout, kaphene effective
lineKapheneCell = kapheneCell.line_graph_cell(name = 'lineKaphene', resonatorsOnly = False) #kaphene layout, line kaffene effective

hpkCell = UnitCell('Huse') #HPG layout, HPK effective
hpKapheneCell = hpkCell.split_cell(name = 'HPKaphene', splitIn = 2) #split hpg layout, heptagon-pentagon kaphene effective
lineHPKapheneCell = hpKapheneCell.line_graph_cell(name = 'lineHPKaphene') #HPKaphene layout, line HP kapheen effective

CharonCell = lineKapheneCell.split_cell(name = 'Charon', splitIn = 2) #split kaphene layout, Charon effective
AliciumCell = CharonCell.line_graph_cell(name = 'Alicium', resonatorsOnly = False) #charonLayout, Alicium effective

lineKagomeCell = kagomeCell.line_graph_cell(name = 'lineKagome', resonatorsOnly = False) #kagome layout, soemthing effective

print 'computing absurd cell'
absurdity = 1 #line kaphene
#absurdity = 2 #Alicium
#absurdity = 3
#absurdity = 4 #who the hell knows
#absurdity = 5 #who the hell knows
startCell = kagomeCell
for ind in range(0, absurdity):
    cell1 = startCell.split_cell(name = 'TBD', splitIn = 2)
    cell2 = cell1.line_graph_cell(name = 'TBD')
    startCell = cell2
    
absurdCell = startCell
absurdCell.type = 'absurdCell_' + str(absurdity) 



#plotCell = kagomeCell
plotCell = kapheneCell
#plotCell = lineKapheneCell

#plotCell = hpkCell
#plotCell = hpKapheneCell
#plotCell = lineHPKapheneCell

#plotCell = hpkCell = CharonCell
#plotCell = hpkCell = AliciumCell

#plotCell = absurdCell

#plotLattice = EuclideanLayout(xcells = 4, ycells = 4, modeType = 'FW', resonatorsOnly=False, initialCell = plotCell)
#plotLattice = EuclideanLayout(xcells = 2, ycells = 2, modeType = 'FW', resonatorsOnly=False, initialCell = plotCell)
plotLattice = EuclideanLayout(xcells = 1, ycells = 1, modeType = 'FW', resonatorsOnly=False, initialCell = plotCell)






#pylab.figure(1)
#pylab.clf()
#ax = pylab.subplot(1,2,1)
#testCell.draw_resonators(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
#testCell.draw_resonator_end_points(ax, color = 'darkgoldenrod', edgecolor = 'k',  marker = 'o' , size = 35, zorder = 5)
##testCell.draw_SDlinks(ax)
#
#pylab.title('original manual construction')
#
#ax = pylab.subplot(1,2,2)
#kapheneCell.draw_resonators(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
#kapheneCell.draw_resonator_end_points(ax, color = 'darkgoldenrod', edgecolor = 'k',  marker = 'o' , size = 35, zorder = 2)
##kapheneCell.draw_SDlinks(ax)
#for ind in range(0, newRes.shape[0]):
#    res = newRes[ind,:]
#    xs = [res[0], res[2]]
#    ys = [res[1], res[3]]
#    pylab.plot(xs, ys)
#pylab.title('new Auto construction')
#
#pylab.show()
#
#
#pylab.figure(2)
#pylab.clf()
#ax = pylab.subplot(1,1,1)
#tempLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
#tempLattice.draw_resonator_end_points(ax, color = 'darkgoldenrod', edgecolor = 'k',  marker = 'o' , size = 35, zorder = 2)
#pylab.show()














print 'computing band cut'
modeType = 'FW'
#kx_x, ky_y, cutx = plotCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = 100, modeType = modeType)
#kx_y, ky_y, cuty = plotCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = 100, modeType = modeType)
kx_x, ky_y, cutx = plotCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = 25, modeType = modeType)
kx_y, ky_y, cuty = plotCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = 25, modeType = modeType)
print 'plotting'


fig2 = pylab.figure(4)
pylab.clf()

ax = pylab.subplot(1,4,1)
plotLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
plotLattice.draw_resonator_end_points(ax, color = 'darkgoldenrod', edgecolor = 'k',  marker = 'o' , size = 35, zorder = 2)
ax.axis('off')

ax = pylab.subplot(1,4,2)
plotLattice.draw_SDlinks(ax, color = 'dodgerblue', linewidth = 1.5)
pylab.scatter(plotLattice.SDx, plotLattice.SDy ,c =  'lightsteelblue', s = 25, marker = 'o', edgecolors = 'mediumblue', zorder = 5, alpha=1)
ax.axis('off')

ax = pylab.subplot(1,4,3)
plotCell.plot_band_cut(ax, cutx)
pylab.title('xcut')
#ax.set_ylim([-2.1,4.05])

ax = pylab.subplot(1,4,4)
plotCell.plot_band_cut(ax, cuty)
pylab.title('ycut')
#ax.set_ylim([-2.1,4.05])

titleStr = plotCell.type + ', modeType: ' + modeType + ' (Made with UnitCell class)' 
pylab.suptitle(titleStr)

pylab.tight_layout()

pylab.show()





