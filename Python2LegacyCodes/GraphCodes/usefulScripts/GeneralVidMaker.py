#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  8 15:53:03 2018

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



generalFolderPath = r'/Users/kollar2/Documents/HouckLab/GeneralLayoutCode/'
if not generalFolderPath in sys.path:
    sys.path.append(generalFolderPath)
    
hyperbolicFolderPath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/'
euclideanFolderPath = r'/Users/kollar2/Documents/HouckLab/EuclideanLatticePlanning/'
if not hyperbolicFolderPath in sys.path:
    sys.path.append(hyperbolicFolderPath)
if not euclideanFolderPath in sys.path:
    sys.path.append(euclideanFolderPath)
    
from GeneralLayoutGenerator import*

from EuclideanLayoutGenerator2 import EuclideanLayout
from EuclideanLayoutGenerator2 import UnitCell

from LayoutGenerator5 import PlanarLayout







def make_layout_vid(GenLayout, startInd =0, stopInd = -1, figNum = 8, xsize = 14, ysize = 4):
    
    
    if stopInd <0:
        #set a break point that will not be reached
        breakPoint = len(GenLayout.Es)+1
    else:
        breakPoint = stopInd
        
    if startInd > len(GenLayout.Es):
        raise ValueError , 'dont have this many eigenvectors'
        
    
    folder = GenLayout.name + '_' + GenLayout.modeType

    
    currDir = os.getcwd()
    saveDir = os.path.join(currDir,folder)
    if os.path.isdir(saveDir):
        pass
    else:
        os.mkdir(saveDir)
    
    
    print folder + '\n'
        
    fig = pylab.figure(figNum)
    #######################################################
    for eigNum in range(startInd, len(GenLayout.Es)):
        print eigNum
        
        eigAmps = GenLayout.Psis[:,GenLayout.Eorder[eigNum]]
        
        
        
        fig.clf()
#        fig.figsize = (100,100)
        ax1 = pylab.subplot(1,2,1)
        xs = scipy.arange(0,len(GenLayout.Es),1)
        pylab.plot(xs, GenLayout.Es[GenLayout.Eorder], 'b.')
        
        pylab.plot(xs[eigNum],GenLayout.Es[GenLayout.Eorder[eigNum]], color = 'firebrick' , marker = '.', markersize = '10' )
        pylab.title(' full tight-binding Spectrum')
        ax1.set_ylim([-2.2, 4.1])
        
#        ax1 = pylab.subplot(1,2,2)
#        xs = scipy.arange(0,len(GenLayout.Es),1)
#        pylab.plot(xs, GenLayout.Es[GenLayout.Eorder], 'b.')
#        
#        pylab.plot(xs[eigNum],GenLayout.Es[GenLayout.Eorder[eigNum]], color = 'firebrick' , marker = '.', markersize = '10' )
#        pylab.title(' full tight-binding Spectrum')
#        ax1.set_ylim([-2.2, 4.1])
        
        ax1 = pylab.subplot(1,2,2)
        titleStr = 'eigenvector weight : ' + str(eigNum) + ' ; E = ' + str(GenLayout.Es[GenLayout.Eorder[eigNum]])
##        GenLayout.plot_layout_state(eigAmps, ax1, title = titleStr, colorbar = True, plot_links = True, cmap = 'Wistia')
        GenLayout.plot_layout_state(eigAmps, ax1, title = titleStr, colorbar = True, plot_links = False, cmap = 'Wistia')
        GenLayout.draw_SDlinks(ax1, color = 'firebrick', linewidth = 0.5, minus_links = True, minus_color = 'dodgerblue', NaNs = True)
        
        
        ax1.set_aspect('auto')
        

        fig.set_size_inches(xsize, ysize)
        
        pylab.tight_layout()
    
        fig_name = GenLayout.name + '_' + GenLayout.modeType +  '_eig_' + str(eigNum) + '.png' 
        fig_path = os.path.join(saveDir, fig_name)
        pylab.savefig(fig_path)
        
        if eigNum == breakPoint:
            break
        
    pylab.show()
    
    return


#######hyperbolic
#test1 = PlanarLayout(gon = 10, vertex = 3, side =1, radius_method = 'lin')
#test1.populate(2, resonatorsOnly=False)
#resonators = test1.get_all_resonators()
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  '10gon_3vertex_2')
#xsize = 12.7
#ysize = 6.3

#######hyperbolic, HW
#test1 = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'HW')
#test1.populate(3, resonatorsOnly=False)
#resonators = test1.get_all_resonators()
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  '7gon_3vertex_3')
#xsize = 12.7
#ysize = 6.3
    
#######hyperbolic
#test1 = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin')
#test1.populate(2, resonatorsOnly=False)
#resonators = test1.get_all_resonators()
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  '7gon_3vertex_2')
#xsize = 14
#ysize = 4

#######Euclidean
#test1 = EuclideanLayout(4,3,lattice_type = 'Huse', modeType = 'FW')
#resonators = test1.resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Huse_43')
#xsize = 10
#ysize = 10
    
#######Euclidean
#test1 = EuclideanLayout(6,6,lattice_type = 'kagome', modeType = 'FW')
#resonators = test1.resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'kagome_66')
#xsize = 10.4
#ysize = 6.3
    
#######Euclidean2
#test1 = EuclideanLayout(4,3,lattice_type = '123Huse', modeType = 'FW')
#resonators = test1.resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  '123Huse_43')
#xsize = 9
#ysize = 7.5


########tree
#test1 = TreeResonators(degree = 3, iterations = 5, side = 1, file_path = '', modeType = 'FW')
#resonators = test1.get_all_resonators()
###Tree2 = TreeResonators(file_path = '3regularTree_ 3_.pkl')
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  '3regTree_5')
#xsize = 11.5
#ysize = 5.4


#######split tree (effective= McLaughlin)
#test1 = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
#resonators = test1.get_all_resonators()
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'McLaughlinTree_4')
#xsize = 10
#ysize = 5

########line graph of split tree (layout = McLaughlin)
#test1 = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
#resonators = test1.get_all_resonators()
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_McLaughlinTree_4')
#xsize = 10
#ysize = 5

    
########split Euclidean
#test1 = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'FW')
#resonators = test1.resonators
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'split_kagome_44')
#xsize = 9
#ysize = 6

########split Euclidean2
#test1 = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'FW')
#resonators = test1.resonators
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'SplitHPG')
#xsize = 10
#ysize = 10


########line graph of Euyclidean
#test1 = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'FW')
#resonators = test1.resonators
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_kagome')
#xsize = 14
#ysize = 4

########line graph of Euyclidean2
#test1 = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'FW')
#resonators = test1.resonators
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_Huse')
#xsize = 10
#ysize = 10



#######non-trivial tree, 2
#test1 = TreeResonators(cell ='', degree = 3, iterations = 3, side = 1, file_path = '', modeType = 'FW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.get_all_resonators()
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'PeterTree_3')
#xsize = 11.5
#ysize = 6

#########decorated Euyclidean
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.resonators
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Peter_graphene_33')
#xsize = 10
#ysize = 8

#########decorated Euyclidean
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.resonators
#LGresonators = generate_line_graph(resonators)
#resonators = LGresonators
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Peter_kagome_33')
#xsize = 14
#ysize = 4

#########decorated Euyclidean
#test1 = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'FW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.resonators
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Peter_HPG_32')
#xsize = 10
#ysize = 9

#########decorated Euyclidean, HW
#test1 = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'HW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.resonators
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Peter_HPG_32')
#xsize = 10
#ysize = 9

#########decorated Euyclidean
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'HW')
#ucell = UnitCell('PeterChain_tail', side = 1)
#resonators = test1.resonators
#decorated_resonators = decorate_layout(resonators, ucell.resonators)
#resonators = decorated_resonators
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'Peter_graphene_33')
#xsize = 10
#ysize = 8


########split Euclidean with defect
#test1 = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'FW')
#resonators = test1.resonators
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'split_kagome_44_defect')
#testLattice.H[0,0] = 100
##for ind in range(0, testLattice.H.shape[0]):
##    testLattice.H[ind, ind] = 1.*ind
#testLattice.Es, testLattice.Psis = scipy.linalg.eigh(testLattice.H)
#testLattice.Eorder = numpy.argsort(testLattice.Es)
#xsize = 9
#ysize = 6
    



#########interface of Euyclidean
#xwidth = 3
#test1 = EuclideanLayout(xwidth,6,lattice_type = 'kagome', modeType = 'FW')
#resonators1 = test1.resonators
####testLattice = GeneralLayout(resonators1 , modeType = test1.modeType, name =  'developing')
#test2 = EuclideanLayout(xwidth,2,lattice_type = 'Huse', modeType = 'FW')
#resonators2 = test2.resonators
####testLattice = GeneralLayout(resonators2, modeType = test1.modeType, name =  'developing')
##shift second set of resonators to be next to first set
#resonators2 = shift_resonators(resonators2, xwidth * test1.unitcell.a1[0], xwidth*test1.unitcell.a1[1])
#
#
##remove the resonators that are double counted
#removals = []
#for rind in range(0, resonators2.shape[0]):
#    newres = resonators2[rind,:]
#    
#    temp = resonators1-newres
#    hitsVec = numpy.all(temp == 0, axis=1)
#    hits = numpy.where(hitsVec == True)[0] #resonators in resonators 1 that are double counted from resonators 2
#    
#    if not(hits.shape[0] == 0):
#        print 'found an edge resonator'
#        print hits
#        removals.append(hits[0])        
#resonators1[removals,:] = 0 #zero out the interface resoantors
#resonators1 = resonators1[~numpy.all(resonators1 == 0, axis=1)] #cut them out
#    
#resonators = numpy.concatenate((resonators1, resonators2)) #merge the two types of lattice
#    
#testLattice = GeneralLayout(resonators, modeType = test1.modeType, name =  'Huse_kagome_interface')
#xsize = 12
#ysize = 6.5


#########interface of Euyclidean2
#xwidth = 3
#test1 = EuclideanLayout(xwidth,3,lattice_type = '123Huse', modeType = 'FW')
#resonators1 = test1.resonators
####testLattice = GeneralLayout(resonators1 , modeType = test1.modeType, name =  'developing')
#test2 = EuclideanLayout(xwidth,3,lattice_type = 'Huse', modeType = 'FW')
#resonators2 = test2.resonators
####testLattice = GeneralLayout(resonators2, modeType = test1.modeType, name =  'developing')
##shift second set of resonators to be next to first set
#resonators2 = shift_resonators(resonators2, xwidth * test1.unitcell.a1[0], xwidth*test1.unitcell.a1[1])
#
#
##remove the resonators that are double counted
#removals = []
#for rind in range(0, resonators2.shape[0]):
#    newres = resonators2[rind,:]
#    
#    temp = resonators1-newres
#    hitsVec = numpy.all(temp == 0, axis=1)
#    hits = numpy.where(hitsVec == True)[0] #resonators in resonators 1 that are double counted from resonators 2
#    
#    if not(hits.shape[0] == 0):
#        print 'found an edge resonator'
#        print hits
#        removals.append(hits[0])        
#resonators1[removals,:] = 0 #zero out the interface resoantors
#resonators1 = resonators1[~numpy.all(resonators1 == 0, axis=1)] #cut them out
#    
#resonators = numpy.concatenate((resonators1, resonators2)) #merge the two types of lattice
#    
#testLattice = GeneralLayout(resonators, modeType = test1.modeType, name =  'Huse_123Huse_interface')
#xsize = 9.5
#ysize = 7.5




########split HPG, with hair
#cell0 = UnitCell('Huse')
#res0 = cell0.resonators
#res1 = split_resonators(res0)
#testCell = UnitCell('split_HPG', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
#
#ysize = 2
#test1 = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'FW')
##test1 = EuclideanLayout(5,6,lattice_type = 'Huse', modeType = modeType)
#resonators = test1.resonators
#splitGraph = split_resonators(resonators)
#resonators = splitGraph
#
#newResonators = numpy.zeros((3*ysize, 4))
#offset = numpy.asarray([testCell.coords[0,0], testCell.coords[0,1], testCell.coords[0,0], testCell.coords[0,1]])
#for rind in range(0, newResonators.shape[0]):
#    temp = numpy.asarray([-0.5,rind*testCell.a2[1]/3., 0,rind*testCell.a2[1]/3.])
#    newResonators[rind,:] =temp + offset
#resonators = numpy.concatenate((resonators, newResonators))
#testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'split_HPG_withhair_32')
#xsize = 10.5
#ysize = 7.



########split HPG, with hair everywhere
#cell0 = UnitCell('Huse')
#res0 = cell0.resonators
#res1 = split_resonators(res0)
#testCell = UnitCell('split_HPG', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
#
#testLattice = GeneralLayout(file_path = 'split_HPG_WHE_32.pkl')
##testLattice = GeneralLayout(file_path = 'split_HPG_WHE_32_HW.pkl')
#resonators = testLattice.resonators
#
#xsize = 10.5
#ysize = 7.


########Huse (HPG/HPK), with hair everywhere
#testCell = UnitCell('Huse')
#
#testLattice = GeneralLayout(file_path = 'Huse_WHE_32.pkl')
##testLattice = GeneralLayout(file_path = 'Huse_WHE_32_HW.pkl')
#resonators = testLattice.resonators
#
#xsize = 10.5
#ysize = 7.


########123Huse (DTG/DTK), with hair everywhere
#testCell = UnitCell('123Huse')
#
#testLattice = GeneralLayout(file_path = '123Huse_WHE_32.pkl')
##testLattice = GeneralLayout(file_path = '123Huse_WHE_32_HW.pkl')
#resonators = testLattice.resonators
#
#xsize = 10.5
#ysize = 7.



#########hairy hyperbolic
#hairy0 = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin')
#hairy0.populate(3)
#
#
#tempResonators =hairy0.get_all_resonators()
#totalres = tempResonators.shape[0]
#numOuterAz = hairy0.points[hairy0.itter].shape[0]
#numOuterRad = hairy0.radials[hairy0.itter].shape[0]
#resonators0 = tempResonators[0:(totalres  - numOuterAz-numOuterRad),:]
#resonators1 = tempResonators[(totalres-numOuterRad):,:]
#resonators = numpy.concatenate((resonators0, resonators1))
#
#testLattice = GeneralLayout(resonators, modeType = hairy0.modeType, name =  '7gon_3vertex_3_WHE')
#xsize = 10.5
#ysize = 5.















#######
#actually make videos
#######



#make_layout_vid(testLattice, startInd = 0, stopInd = 3, figNum = 1, xsize = xsize, ysize = ysize)
#make_layout_vid(testLattice, startInd = 0, stopInd = -1, figNum = 1, xsize = xsize, ysize = ysize)





######
#view lattices
######



pylab.figure(3)
pylab.clf()
ax = pylab.subplot(1,1,1)
testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
testLattice.draw_SDlinks(ax, color = 'deepskyblue', linewidth = 1.5, minus_links = True, minus_color = 'firebrick')
xs = testLattice.coords[:,0]
ys = testLattice.coords[:,1]
pylab.sca(ax)
#pylab.scatter(xs, ys ,c =  'goldenrod', s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 3)
#pylab.scatter(xs, ys ,c =  'goldenrod', s = 40, marker = 'o', edgecolors = 'k', zorder = 5)

siteNum = 0
pylab.scatter(testLattice.SDx[siteNum], testLattice.SDy[siteNum] ,c =  'firebrick', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)


ax.set_aspect('equal')
ax.axis('off')
pylab.title(testLattice.name )
pylab.tight_layout()
pylab.show()




