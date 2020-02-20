#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 16:09:15 2018

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

pkgDir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

if not pkgDir in sys.path:
    sys.path.append(pkgDir)

from GeneralLayoutGenerator import GeneralLayout
from GeneralLayoutGenerator import TreeResonators

#hyperbolicFolderPath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/'
#euclideanFolderPath = r'/Users/kollar2/Documents/HouckLab/EuclideanLatticePlanning/'
#if not hyperbolicFolderPath in sys.path:
#    sys.path.append(hyperbolicFolderPath)
#if not euclideanFolderPath in sys.path:
#    sys.path.append(euclideanFolderPath)

from EuclideanLayoutGenerator2 import UnitCell
from EuclideanLayoutGenerator2 import EuclideanLayout

from LayoutGenerator5 import PlanarLayout


from GeneralLayoutGenerator import split_resonators
from GeneralLayoutGenerator import generate_line_graph
from GeneralLayoutGenerator import decorate_layout


#######hyperbolic
#test1 = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
#test1.populate(2, resonatorsOnly=False)
#resonators = test1.get_all_resonators()
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  '7gon_3vertex_2')

#######single ring
#test1 = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
#test1.populate(2, resonatorsOnly=False)
#resonators = test1.get_all_resonators()
#resonators = resonators[0:test1.gon, :]
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'ring')

#######split hyperbolic
#test1 = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin')
#test1.populate(2, resonatorsOnly=False)
#resonators = test1.get_all_resonators()
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'split7gon_3vertex_2')
#
#######Euclidean
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#resonators = test1.resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'kagome')

#######Euclidean
#test1 = EuclideanLayout(4,4,lattice_type = 'square', modeType = 'FW')
#resonators = test1.resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'square')

#######Euclidean
##test1 = EuclideanLayout(4,3,lattice_type = 'Huse', modeType = 'FW')
#test1 = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'FW')
#resonators = test1.resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Huse')


#######tree
#test1 = TreeResonators(degree = 3, iterations = 5, side = 1, file_path = '', modeType = 'FW')
#resonators = test1.get_all_resonators()
###Tree2 = TreeResonators(file_path = '3regularTree_ 3_.pkl')
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'TREEEEE')


#######split tree
#test1 = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
#resonators = test1.get_all_resonators()
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'McLaughlinTree')

#######line graph of split tree
#test1 = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
#resonators = test1.get_all_resonators()
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_McLaughlinTree')

    
########split Euclidean
#test1 = EuclideanLayout(3,2,lattice_type = 'kagome', modeType = 'FW')
#resonators = test1.resonators
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Splitkagome')

########split Euclidean
#test1 = EuclideanLayout(3,3,lattice_type = 'square', modeType = 'FW')
#resonators = test1.resonators
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'SplitSquare')

#########split Euclidean, hoffman attemps
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#resonators0 = test1.resonators #graphene layout
#splitGraph = split_resonators(resonators0) #split graphene layout
#resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
#resonators2 = split_resonators(resonators1) #split further
#resonators = resonators2
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')

#########split Euclidean, hoffman attemps
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#resonators0 = test1.resonators #graphene layout
#splitGraph = split_resonators(resonators0) #split graphene layout
#resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
#resonators2 = split_resonators(resonators1) #split further
#resonators3 = generate_line_graph(resonators2)
#resonators = resonators3
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'L(hofmannAttempt)')



########split Euclidean 1.5 (n-way split)
#splitIn = 3
#test1 = EuclideanLayout(3,2,lattice_type = 'kagome', modeType = 'FW')
#resonators = test1.resonators
#splitGraph = split_resonators(resonators, splitIn = splitIn)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  str(splitIn) + 'SplitGraphene')

########split Euclidean2
#test1 = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'FW')
#resonators = test1.resonators
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'split_HPG')


########line graph of Euyclidean
#test1 = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'FW')
#resonators = test1.resonators
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_kagome')

########line graph of Euyclidean2
#test1 = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'FW')
#resonators = test1.resonators
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_Huse')

########line graph of line graph of Euyclidean
#test1 = EuclideanLayout(2,2,lattice_type = 'kagome', modeType = 'FW')
#resonators = test1.resonators
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_LG_kagome')


######non-trivial tree
#test1 = TreeResonators(cell ='Peter', degree = 3, iterations = 3, side = 1, file_path = '', modeType = 'FW')
#resonators = test1.get_all_resonators()
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'PeterTREEEEE')
#
#######non-trivial tree, 2
#test1 = TreeResonators(cell ='', degree = 3, iterations = 3, side = 1, file_path = '', modeType = 'FW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.get_all_resonators()
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'PeterTREEEEE2')
#
#########decorated Euyclidean
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.resonators
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Peter-graphene')
#
#########decorated Euyclidean
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.resonators
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Peter-kagome')
#
#########decorated Euyclidean
#test1 = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'FW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.resonators
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Peter-HPG')
#
#
######decorated hyperbolic
test1 = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin')
test1.populate(2, resonatorsOnly=False)
#resonators = test1.get_all_resonators()
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
ucell = UnitCell('PeterChain_tail', side = 1)
resonators = test1.get_all_resonators()
decorated_resonators = decorate_layout(resonators, ucell.resonators)
resonators = decorated_resonators
testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'decorated_hyperbolic')


    


######
#plot results
######





fig1 = pylab.figure(1)
pylab.clf()
ax = pylab.subplot(1,1,1)
testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
testLattice.draw_SDlinks(ax, color = 'deepskyblue', linewidth = 1.5, minus_links = True, minus_color = 'firebrick')
xs = testLattice.coords[:,0]
ys = testLattice.coords[:,1]
pylab.sca(ax)
#pylab.scatter(xs, ys ,c =  'goldenrod', s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)
#pylab.scatter(xs, ys ,c =  'goldenrod', s = 40, marker = 'o', edgecolors = 'k', zorder = 5)
ax.set_aspect('equal')
ax.axis('off')
pylab.title(testLattice.name)
fig1.set_size_inches(5, 5)
pylab.tight_layout()
pylab.show()


eigNum = 168
eigNum = 167
eigNum = 0

pylab.figure(2)
pylab.clf()
ax = pylab.subplot(1,2,1)
pylab.imshow(testLattice.H,cmap = 'winter')
pylab.title('Hamiltonian')
ax = pylab.subplot(1,2,2)
pylab.imshow(testLattice.H - numpy.transpose(testLattice.H),cmap = 'winter')
pylab.title('H - Htranspose')
pylab.suptitle(testLattice.name)
pylab.show()



xs = scipy.arange(0,len(testLattice.Es),1)
eigAmps = testLattice.Psis[:,testLattice.Eorder[eigNum]]

pylab.figure(3)
pylab.clf()
ax1 = pylab.subplot(1,2,1)
pylab.plot(testLattice.Es, 'b.')
pylab.plot(xs[eigNum],testLattice.Es[testLattice.Eorder[eigNum]], color = 'firebrick' , marker = '.', markersize = '10' )
pylab.title('eigen spectrum')
pylab.ylabel('Energy (t)')
pylab.xlabel('eigenvalue number')

ax2 = pylab.subplot(1,2,2)
titleStr = 'eigenvector weight : ' + str(eigNum)
testLattice.plot_layout_state(eigAmps, ax2, title = titleStr, colorbar = True, plot_links = True, cmap = 'Wistia')

pylab.suptitle(testLattice.name)
pylab.show()




########teting vertex finding codes
##vertexInd = 11
#vertexInd = int(numpy.floor(testLattice.coords.shape[0]/2.))
#
#fig4 = pylab.figure(4)
#pylab.clf()
#ax = pylab.subplot(1,1,1)
#testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
#testLattice.draw_SDlinks(ax, color = 'deepskyblue', linewidth = 1.5, minus_links = True, minus_color = 'goldenrod')
#
##newCoords = get_coords(testLattice.resonators)
#vertexx = testLattice.coords[vertexInd,0]
#vertexy = testLattice.coords[vertexInd,1]
#pylab.scatter(vertexx, vertexy ,c =  'firebrick', s = 200, marker = 'o', edgecolors = 'k', zorder = 0, alpha = 1)
#
##testLattice.generate_vertex_dict()
#connectedResonators = testLattice.vertexDict[vertexInd]
#for rind in connectedResonators:
#    [x0,y0,x1,y1] = testLattice.resonators[rind]
#    pylab.plot([x0, x1],[y0, y1] , color = 'gold', linewidth = 4, alpha = 1, zorder = 10)
#
#xs = testLattice.coords[:,0]
#ys = testLattice.coords[:,1]
#pylab.sca(ax)
#pylab.scatter(xs, ys ,c =  'goldenrod', s = 34, marker = 'o', edgecolors = 'k', zorder = 5)
#
#
#ax.set_aspect('equal')
#ax.axis('off')
#pylab.title('checking vertex dictionary : ' + testLattice.name)
#fig1.set_size_inches(5, 5)
#pylab.tight_layout()
#pylab.show()






fig5 = pylab.figure(5)
pylab.clf()
ax = pylab.subplot(1,2,1)
pylab.cla()
testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
xs = testLattice.coords[:,0]
ys = testLattice.coords[:,1]
pylab.sca(ax)
pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)
ax.set_aspect('equal')
ax.axis('off')
pylab.title('layout graph')


ax2 = pylab.subplot(1,2,2)
testLattice.draw_SDlinks(ax2, color = 'dodgerblue', linewidth = 2.5, minus_links = True, minus_color = 'gold')
#testLattice.draw_SD_points(ax, color = 'mediumblue', edgecolor = 'goldenrod',  marker = 'o' , size = 30,  extra = False)
pylab.sca(ax2)
pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  'lightsteelblue', s = 30, marker = 'o', edgecolors = 'mediumblue', zorder = 5)
ax2.set_aspect('equal')
ax2.axis('off')
pylab.title('effective (line) graph')

pylab.suptitle(testLattice.name)
fig5.set_size_inches(10, 5)
pylab.tight_layout()
pylab.show()





fig6 = pylab.figure(6)
pylab.clf()
ax = pylab.subplot(1,1,1)
#testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 0.5)
testLattice.draw_SDlinks(ax, color = 'firebrick', linewidth = 1., minus_links = True, minus_color = 'goldenrod')
xs = testLattice.coords[:,0]
ys = testLattice.coords[:,1]
pylab.sca(ax)
##pylab.scatter(xs, ys ,c =  'goldenrod', s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
#pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)
##pylab.scatter(xs, ys ,c =  'goldenrod', s = 40, marker = 'o', edgecolors = 'k', zorder = 5)
ax.set_aspect('equal')
ax.axis('off')
#pylab.title(testLattice.name)
pylab.tight_layout()
pylab.show()




