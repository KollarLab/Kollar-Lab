#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 15:49:59 2018

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

from GraphCodes.GeneralLayoutGenerator import GeneralLayout     #the most general layout that we have
from GraphCodes.GeneralLayoutGenerator import TreeResonators     # special case for making trees

from GraphCodes.EuclideanLayoutGenerator2 import UnitCell    #fundamental domain object
from GraphCodes.EuclideanLayoutGenerator2 import EuclideanLayout   #normal Euclidean lattices

from GraphCodes.LayoutGenerator5 import PlanarLayout   #original hyperbolic

from GraphCodes.GeneralLayoutGenerator import split_resonators     #functions to modify existing layouts
from GraphCodes.GeneralLayoutGenerator import generate_line_graph
from GraphCodes.GeneralLayoutGenerator import decorate_layout






#######
#hyperbolic lattices
#######
test = PlanarLayout(gon = 7, vertex = 3,side = 1, modeType = 'FW')
test.populate(maxItter = 2)


pylab.figure(1)
pylab.clf()
ax = pylab.subplot(1,1,1)
test.draw_resonator_lattice(ax, color = 'mediumblue', linewidth = 2)
test.draw_resonator_end_points(ax, color = 'cornflowerblue', edgecolor = 'c', size = 150)

ax.axis('off') #reomves the outside box
ax.set_aspect('equal')

pylab.tight_layout() #don't waste space
pylab.show()



#######
#Euclidean lattice
#######

testCell = UnitCell('kagome') #fundamental domain object
testLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW') #The full lattice

pylab.figure(2)
pylab.clf()
ax = pylab.subplot(1,2,1)
testLattice.draw_resonator_lattice(ax, color = 'mediumblue', linewidth = 2)
testLattice.draw_resonator_lattice(ax, color = 'firebrick', linewidth = 2, extras = True)
testLattice.draw_resonator_end_points(ax, color = 'cornflowerblue', edgecolor = 'c', size = 150)

ax.axis('off') #reomves the outside box
ax.set_aspect('equal')
pylab.title('lattice')



ax = pylab.subplot(1,2,2)
testLattice.unitcell.draw_resonators(ax, color = 'mediumblue', linewidth = 2)
testLattice.unitcell.draw_resonator_end_points(ax, color = 'cornflowerblue', edgecolor = 'c', size = 150)
ax.axis('off') #reomves the outside box
ax.set_aspect('equal')
pylab.title('cell')


pylab.show()



######
#General layout
######


#######hyperbolic
#test1 = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
#test1.populate(2, resonatorsOnly=False)
#resonators = test1.get_all_resonators()
#genLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  '7gon_3vertex_2')


#######euclidean
#resonators = testLattice.get_all_resonators()
#genLattice = GeneralLayout(resonators , modeType = testLattice.modeType, name =  'kagome33')


#######modified euclidean
#resonators = testLattice.get_all_resonators()
#resonators = split_resonators(resonators)
#genLattice = GeneralLayout(resonators , modeType = testLattice.modeType, name =  'splitkagome33')



#
#######tree
#test1 = TreeResonators(degree = 3, iterations = 5, side = 1, file_path = '', modeType = 'FW')
#resonators = test1.get_all_resonators()
#genLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'TREEEEE')


######split tree
#test1 = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
#resonators = test1.get_all_resonators()
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#genLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'McLaughlinTree')


#######line graph of split tree
test1 = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
resonators = test1.get_all_resonators()
splitGraph = split_resonators(resonators)
resonators = splitGraph
LGresonators = generate_line_graph(resonators)
resonators = LGresonators
genLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_McLaughlinTree')





pylab.figure(3)
pylab.clf()
ax = pylab.subplot(1,1,1)
genLattice.draw_resonator_lattice(ax, color = 'mediumblue', linewidth = 2)
genLattice.draw_resonator_end_points(ax, color = 'cornflowerblue', edgecolor = 'c', size = 150)

ax.axis('off') #reomves the outside box
ax.set_aspect('equal')
pylab.title('general')


pylab.show()







########
#line graph a.k.a. medial lattice
#######
#####split tree
test1 = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
resonators = test1.get_all_resonators()
splitGraph = split_resonators(resonators)
resonators = splitGraph
genLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'McLaughlinTree')


#######line graph of split tree
test1 = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
resonators = test1.get_all_resonators()
splitGraph = split_resonators(resonators)
resonators = splitGraph
LGresonators = generate_line_graph(resonators)
resonators = LGresonators
genLattice2 = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_McLaughlinTree')

pylab.figure(4)
pylab.clf()
ax = pylab.subplot(1,2,1)
genLattice.draw_resonator_lattice(ax, color = 'mediumblue', linewidth = 2)
genLattice.draw_resonator_end_points(ax, color = 'cornflowerblue', edgecolor = 'c', size = 150)

ax.axis('off') #reomves the outside box
ax.set_aspect('equal')
pylab.title('hardware layout')


ax = pylab.subplot(1,2,2)
genLattice.draw_resonator_lattice(ax, color = 'mediumblue', linewidth = 2, alpha = 0.5)
genLattice.draw_resonator_end_points(ax, color = 'cornflowerblue', edgecolor = 'c', size = 75 )

genLattice2.draw_resonator_lattice(ax, color = 'k', linewidth = 2)
genLattice2.draw_resonator_end_points(ax, color = 'firebrick', edgecolor = 'c', size = 75)

ax.axis('off') #reomves the outside box
ax.set_aspect('equal')
pylab.title('both')

pylab.show()










