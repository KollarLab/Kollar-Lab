#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:23:05 2018

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
import scipy.io as sio

hyperbolicFolderPath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/'
euclideanFolderPath = r'/Users/kollar2/Documents/HouckLab/EuclideanLatticePlanning/'
generalFolderPath = r'/Users/kollar2/Documents/HouckLab/GeneralLayoutCode/'
if not hyperbolicFolderPath in sys.path:
    sys.path.append(hyperbolicFolderPath)
if not euclideanFolderPath in sys.path:
    sys.path.append(euclideanFolderPath)
if not generalFolderPath in sys.path:
    sys.path.append(generalFolderPath)

from GeneralLayoutGenerator import GeneralLayout
from GeneralLayoutGenerator import TreeResonators

from EuclideanLayoutGenerator2 import UnitCell
from EuclideanLayoutGenerator2 import EuclideanLayout

from LayoutGenerator5 import PlanarLayout


from GeneralLayoutGenerator import split_resonators
from GeneralLayoutGenerator import generate_line_graph
from GeneralLayoutGenerator import decorate_layout

#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#from matplotlib import cm



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













def plot_any_band_cut(ax, cut, connected = True):
    
    colorlist = ['firebrick', 'dodgerblue', 'blueviolet', 'mediumblue', 'goldenrod', 'cornflowerblue']
    
    pylab.sca(ax)
    
    shape = cut.shape
    
    for ind in range(0,shape[0]):
        colorInd = numpy.mod(ind, len(colorlist))
        if connected:
            pylab.plot(cut[ind,:], color = colorlist[colorInd] , marker = '', markersize = '5', linestyle = '-', linewidth = 3.5)
        else:
            pylab.plot(cut[ind,:], color = colorlist[colorInd] , marker = '.', markersize = '5', linestyle = '')
        
    pylab.title('some momentum cut')
    pylab.ylabel('Energy')
    pylab.xlabel('k_something')
        
    return



forceLims = True
#    forceLims = False
if forceLims:
    axmin = -2.75
    axmax = 3.25
    
    
    
    
    
    
    
    
    

##############1)
##kagome
#testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#resoantors1 = split_resonators(testEuclidLattice.resonators)
#splitLattice = GeneralLayout(resoantors1 , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
#
#newCell = split_resonators(testEuclidLattice.unitcell.resonators)
#testCell =  UnitCell('split_graphene', resonators = newCell, a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)
#
#parentCell = testEuclidLattice.unitcell
#
#grapheneCell = UnitCell('kagome', a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)


#############2)
##hpg
#testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'Huse', modeType = 'FW')
#resoantors1 = split_resonators(testEuclidLattice.resonators)
#splitLattice = GeneralLayout(resoantors1 , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
#
#newCell = split_resonators(testEuclidLattice.unitcell.resonators)
#testCell =  UnitCell('split_hpg', resonators = newCell, a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)
#
#parentCell = testEuclidLattice.unitcell
#
#grapheneCell = UnitCell('kagome', a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)



##############3)
##extremal hofman
#########split Euclidean, hoffman attemps
###!!!!!!!!needs to be this size for the cell to come out right. DO NOT CHANGE!!!
#test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#resonators0 = test1.resonators #graphene layout
#splitGraph = split_resonators(resonators0) #split graphene layout
#resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
#resonators2 = split_resonators(resonators1) #split further
#resonators = resonators2
#splitLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
##sites = [80,81,82, 83,79,78,77,76, 71, 84,85,86,87,114,115, 116, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
#sites = [80,81,83,79,78,77,76, 71, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
#state = numpy.zeros(len(splitLattice.SDx))
#state[sites] = 0.1
#newCell = splitLattice.resonators[sites, :]
##    
##    fig1 = pylab.figure(114)
##    pylab.clf()
##    ax = pylab.subplot(1,1,1)
##    testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
##    xs = testLattice.coords[:,0]
##    ys = testLattice.coords[:,1]
##    pylab.sca(ax)
##    pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)
##    pylab.scatter(testLattice.SDx[sites], testLattice.SDy[sites], c =  'firebrick', s = 75, marker = 'o', edgecolors = 'maroon', zorder = 5,linewidth = 1)
##    ax.set_aspect('equal')
##    ax.axis('off')
##    pylab.title(testLattice.name)
##    pylab.tight_layout()
##    pylab.show()
##    fig49 = pylab.figure(49)
#
#testCell = UnitCell('S(L(S(graphene)))', resonators = newCell, a1 = test1.unitcell.a1, a2 = test1.unitcell.a2)
#
#
##cheating way of getting the split graph band structure
#parentLattice = GeneralLayout(resonators1 , modeType = test1.modeType, name =  'hofmannParent')
#
##54, 55,56,57
##    sites = [50, 49, 48, 47, 46, 36,37, 51,52,53, 44,   67,64]
#sites = [50, 49, 48, 47, 46, 36,37,44,51]
#state = numpy.zeros(len(parentLattice.SDx))
#state[sites] = 0.1
#newCell = parentLattice.resonators[sites, :]
#
#parentCell = UnitCell('L(S(graphene))', resonators = newCell, a1 = test1.unitcell.a1, a2 = test1.unitcell.a2)



##############4)
##square
#testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'square', modeType = 'FW')
#resoantors1 = split_resonators(testEuclidLattice.resonators)
#splitLattice = GeneralLayout(resoantors1 , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
#
#newCell = split_resonators(testEuclidLattice.unitcell.resonators)
#testCell =  UnitCell('Lieb', resonators = newCell, a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)
#
#parentCell = testEuclidLattice.unitcell
#
#grapheneCell = UnitCell('kagome', a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)





#############5)
##mess
#testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#resoantors1 = split_resonators(testEuclidLattice.resonators)
#splitLattice = GeneralLayout(resoantors1 , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
#
#newCell = split_resonators(testEuclidLattice.unitcell.resonators)
##make all the cells that I want
#kagomeCell = UnitCell('kagome') #graphene layout, kagome effective
#kapheneCell = kagomeCell.split_cell(name = 'kaphene', splitIn = 2) #split graphene layout, kaphene effective
#lineKapheneCell = kapheneCell.line_graph_cell(name = 'lineKaphene', resonatorsOnly = False) #kaphene layout, line kaffene effective
#
#newCell = split_resonators(testEuclidLattice.unitcell.resonators)
#testCell =  lineKapheneCell.split_cell(name = 'split kaphene layout')
#
#parentCell = testEuclidLattice.unitcell
#parentCell = lineKapheneCell
#
#grapheneCell = UnitCell('kagome', a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)


##############6)
##mess
#testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#resoantors1 = split_resonators(testEuclidLattice.resonators)
#splitLattice = GeneralLayout(resoantors1 , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
#
#newCell = split_resonators(testEuclidLattice.unitcell.resonators)
##testCell =  UnitCell('split_graphene', resonators = newCell, a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)
#testCell =  UnitCell('graphene_layout', resonators = testEuclidLattice.unitcell.resonators, a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)
#
#parentCell = testEuclidLattice.unitcell
#
#grapheneCell = UnitCell('kagome', a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)
#



############5)
#split kaphene cell
testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
resoantors1 = split_resonators(testEuclidLattice.resonators)
splitLattice = GeneralLayout(resoantors1 , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)

newCell = split_resonators(testEuclidLattice.unitcell.resonators)
#make all the cells that I want
kagomeCell = UnitCell('kagome') #graphene layout, kagome effective
kapheneCell = kagomeCell.split_cell(name = 'kaphene', splitIn = 2) #split graphene layout, kaphene effective
lineKapheneCell = kapheneCell.line_graph_cell(name = 'lineKaphene', resonatorsOnly = False) #kaphene layout, line kaffene effective

newCell = split_resonators(testEuclidLattice.unitcell.resonators)
splitKapheneLayoutCell =  lineKapheneCell.split_cell(name = 'split kaphene layout')

testCell = splitKapheneLayoutCell
parentCell = splitKapheneLayoutCell

grapheneCell = UnitCell('kagome', a1 = testEuclidLattice.unitcell.a1, a2 = testEuclidLattice.unitcell.a2)

#kx,ky,temp = compute_layout_band_structure(splitKapheneLayoutCell, 0, 0, 0, 0, numsteps = 1, modeType = 'FW')







#band structure
numSurfPoints = 300
kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = 'FW')
kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = 'HW')







######plot structures
pylab.figure(1)
pylab.clf()
ax = pylab.subplot(1,3,1, adjustable='box', aspect=80)
pylab.cla()
#plot_band_cut_mc(ax, cuty)
testCell.plot_band_cut(ax, cuty)
pylab.title('effective band structure')
pylab.ylabel('Energy (|t|)')
pylab.xlabel('$k_y$ ($\pi$/a)')
pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
if forceLims:
    ax.set_ylim([axmin, axmax])




####hack the layout band strucutre
numSurfPoints = 300
#kx_x, ky_y, cutx = grapheneCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = 'FW')
#kx_y, ky_y, cuty = grapheneCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = 'FW')
kx_x, ky_y, cutx = parentCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = 'FW')
kx_y, ky_y, cuty = parentCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = 'FW')
layoutBandStructureX = cutx[1:,:] -1 #remove the flat bands and shift down
layoutBandStructureY = cuty[1:,:] -1 #remove the flat bands and shift down

stretchX1 = numpy.sqrt(layoutBandStructureX + 3)
stretchY1 = numpy.sqrt(layoutBandStructureY + 3)

FB = numpy.zeros((1, numSurfPoints))

splitCutX = numpy.concatenate((-stretchX1,FB,  stretchX1))
splitCutY = numpy.concatenate((-stretchY1, FB, stretchY1))
#    splitCutX = numpy.concatenate((-stretchX1, stretchX1))
#    splitCutY = numpy.concatenate((-stretchY1, stretchY1))



ax = pylab.subplot(1,3,2, adjustable='box', aspect=80)
pylab.cla()
#plot_band_cut_mc(ax, cuty)
plot_any_band_cut(ax, splitCutY)
pylab.title('hacked split band structure')
pylab.ylabel('Energy (|t|)')
pylab.xlabel('$k_y$ ($\pi$/a)')
pylab.xticks([0, splitCutY.shape[1]/2, splitCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
if forceLims:
    ax.set_ylim([axmin, axmax])
    
    
    
    
    
####actually compute the split band structure

#H_layout = numpy.zeros((testCell.coords.shape[0],testCell.coords.shape[0]))    

#def make_H(cell, kx, ky):
    


def find_layout_cell(cell):
    #determine the unit cell for the layout. The ends of the resonators is too big a set.
    allCoords = numpy.round(cell.coords[:,:], 3)
    svec_all = allCoords[:,0] + 1j*allCoords[:,1]
    
    
    
    def check_redundnacy(cell, site, svec_all, shift1, shift2):
        vec1 = numpy.round(cell.a1[0] + 1j*cell.a1[1], 3)
        vec2 = numpy.round(cell.a2[0] + 1j*cell.a2[1], 3)
        shiftedCoords = svec_all + shift1*(vec1) + shift2*(vec2)
    #    print shiftedCoords
    #    print -0.289-0.5j == shiftedCoords
        redundancies = numpy.where(numpy.round(site,3) == numpy.round(shiftedCoords,3))[0]
        return redundancies
        
    
    #determine coordinate equivalences
    redundancyDict = {}
    for cind in range(0, allCoords.shape[0]):
        site = svec_all[cind] #the site to compare
    #    print site
    #    print cind
        
        redundancyDict[cind] = []
        for shift1 in (-1,0,1):
            for shift2 in (-1,0,1):
                redundancies = check_redundnacy(cell, site, svec_all, shift1, shift2)
    #            print redundancies
                
                if len(redundancies) > 0:
                    if not (shift1 ==0 and shift2 == 0):
                        #found an actual redundancy
                        redundancyDict[cind] = numpy.concatenate((redundancyDict[cind], redundancies))
    #                    print redundancies
    #    print '   '
        
        
    #find the minimum cell
    minCellInds = [0.]
    for cind in range(1, allCoords.shape[0]):
        equivalentInds = redundancyDict[cind] #all the site that are the same as the one we are looking at
    #    print equivalentInds
    #    print minCellInds
    #    print '  '
        if len(equivalentInds)>0:
            for cind2 in range(0, len(equivalentInds)):
                currInd = equivalentInds[cind2]
                if currInd in minCellInds:
                    break
                if cind2 == len(equivalentInds)-1:
                    #no matches found for the site cind
                    minCellInds = numpy.concatenate((minCellInds, [cind]))
        else:
            #no redundant sites
            minCellInds = numpy.concatenate((minCellInds, [cind]))
            
    minCellInds = minCellInds.astype('int')
    #print minCellInds
    
    return minCellInds, redundancyDict
   

def generate_layout_Bloch_matrix(cell, kx, ky, modeType = 'FW', t = 1, phase = 0):
    
    allCoords = numpy.round(cell.coords[:,:], 3)
    svec_all = allCoords[:,0] + 1j*allCoords[:,1]
    
    minCellInds, redundancyDict = find_layout_cell(cell)
        
        
    #get the coordinates of the minimum unit cell
    #coords = numpy.round(cell.coords[[1,3],:], 3)
    coords = numpy.round(cell.coords[minCellInds,:], 3)
    svec = numpy.zeros((coords.shape[0]))*(1 + 1j)
    svec[:] = coords[:,0] + 1j*coords[:,1]
#    print svec
#    print '   '
    
    
    BlochMat = numpy.zeros((coords.shape[0], coords.shape[0]))*(0 + 0j)
    
    #store away the resonators, which tell me about all possible links
    resonators = numpy.round(cell.resonators, 3)
    zmat = numpy.zeros((resonators.shape[0],2))*(1 + 1j)
    zmat[:,0] = resonators[:,0] + 1j*resonators[:,1]
    zmat[:,1] = resonators[:,2] + 1j*resonators[:,3]
#    print zmat
    
    #convert the resonator matrix to links, basically fold it back to the unit cell
    for rind in range(0, resonators.shape[0]):
#        print rind
        source = zmat[rind,0]
        target = zmat[rind,1]
        
        sourceInd = numpy.where(numpy.round(source,3) == numpy.round(svec_all,3))[0][0]
        targetInd = numpy.where(numpy.round(target,3) == numpy.round(svec_all,3))[0][0]
        
        #figure out which types of points in the unit cell we are talking about
        #will call these variables the source class and target class
        if sourceInd in minCellInds:
            #this guy is in the basic unit cell
            sourceClass = sourceInd
        else:
            for cind in minCellInds:
                if sourceInd in redundancyDict[cind]:
                    sourceClass = cind
                    
        if targetInd in minCellInds:
            #this guy is in the basic unit cell
            targetClass = targetInd
        else:
            for cind in minCellInds:
                if targetInd in redundancyDict[cind]:
                    targetClass = cind
                    
#        print [sourceClass, targetClass]
                    
#        #update the source and target #ACTUALY. This appears to be bad. I want to direction of the 
        #actual bond, and the knowledge of which two classes of points I'm going between
#        source = svec_all[sourceClass]
#        target = svec_all[targetClass]
        
        #absolute corrdiates of origin site
        x0 = numpy.real(source)
        y0 = numpy.imag(source)
        
        #absolute coordinates of target site
        x1 = numpy.real(target)
        y1 = numpy.imag(target)
        
        deltaX = x1-x0
        deltaY = y1-y0
        
#        print deltaX
#        print deltaY
        
        #convert from minCellInds which tells which entries of the
        #total coords form a unit cell to
        #indices that label the entires for the matrix
        sourceMatInd = numpy.where(sourceClass == minCellInds)[0][0]
        targetMatInd = numpy.where(targetClass == minCellInds)[0][0]
        
        phaseFactor = numpy.exp(1j*kx*deltaX)*numpy.exp(1j*ky*deltaY)
        BlochMat[sourceMatInd, targetMatInd] = BlochMat[sourceMatInd, targetMatInd]+ t*phaseFactor
        BlochMat[targetMatInd, sourceMatInd] = BlochMat[targetMatInd, sourceMatInd]+ t*numpy.conj(phaseFactor)
        
#        print BlochMat
    return BlochMat

def compute_layout_band_structure(cell, kx_0, ky_0, kx_1, ky_1, numsteps = 100, modeType = 'FW', returnStates = False, phase  = 0):
        '''
        from scipy.linalg.eigh:
        The normalized selected eigenvector corresponding to the eigenvalue w[i] is the column v[:,i].
        
        This returns same format with two additional kx, ky indices
        '''
        
        kxs = numpy.linspace(kx_0, kx_1,numsteps)
        kys = numpy.linspace(ky_0, ky_1,numsteps)
        
        minCellInds, redundancyDict = find_layout_cell(cell)
        numLayoutSites = len(minCellInds)
        
        bandCut = numpy.zeros((numLayoutSites, numsteps))
        
        stateCut = numpy.zeros((numLayoutSites, numLayoutSites, numsteps)).astype('complex')
        
        for ind in range(0, numsteps):
            kvec = [kxs[ind],kys[ind]]
            
            H = generate_layout_Bloch_matrix(cell, kvec[0], kvec[1], modeType = modeType, phase  = phase)
        
            #Psis = numpy.zeros((self.numSites, self.numSites)).astype('complex')
            Es, Psis = scipy.linalg.eigh(H)
            
            bandCut[:,ind] = Es
            stateCut[:,:,ind] = Psis
        if returnStates:
            return kxs, kys, bandCut, stateCut
        else:
            return kxs, kys, bandCut        


#deltaA1 = 0
#deltaA2 = 0
#for sind in range(0, coords.shape[0]):         
#    source = svec[sind]
#
#    
#    hits0 = numpy.where(source == zmat[:,0])[0] #find all resonators with the source as the first point
##            print hits0
#    
#    hits1 = numpy.where(source == zmat[:,1])[0] #find all resonators with the source as the second point
##            print hits1
#
#    hits = numpy.concatenate((hits0, hits1))
#    for hind in range(0, len(hits)):
#        tind_res = hits[hind] #index of the target point in the resonator matrix
#        
#        #get coordinaes of the target point
#        if hind < len(hits0):
#            target = zmat[tind_res,0]
#        else:
#            target = zmat[tind_res,1]
#        print target
#        
#        #check if target is in the unit cell. Might be redundant
#        if target in svec:
#            print 'yes this is a link'
#            #this target state is in the unit cell, not one of the extras
#            #so this is a valid link that needs to be stored
#            
#            tind = numpy.where(target == svec)[0][0]
#            if sind == tind:
#                #found yourself
#                pass
#            else:
#            
#                #corrdiates of origin site
#                x0 = numpy.real(source)
#                y0 = numpy.imag(source)
#                
#                #coordinates of target site
#                x1 = cell.coords[tind,0] + deltaA1*cell.a1[0] + deltaA2*cell.a2[0]
#                y1 = cell.coords[tind,0] + deltaA1*cell.a1[1] + deltaA2*cell.a2[1]
#                
#                deltaX = x1-x0
#                deltaY = y1-y0
#                
#                phaseFactor = numpy.exp(1j*kx*deltaX)*numpy.exp(1j*ky*deltaY)
#                BlochMat[sind, tind] = BlochMat[sind, tind]+ t*phaseFactor

            




numSurfPoints = 300
#kx_x, ky_y, cutx = compute_layout_band_structure(grapheneCell, -2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = 'FW')
#kx_y, ky_y, cuty = compute_layout_band_structure(grapheneCell, 0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = 'FW')
kx_x, ky_y, cutx = compute_layout_band_structure(testCell, -2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = 'FW')
kx_y, ky_y, cuty = compute_layout_band_structure(testCell, 0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = 'FW')



ax = pylab.subplot(1,3,3, adjustable='box', aspect=80)
pylab.cla()
#plot_band_cut_mc(ax, cuty)
plot_any_band_cut(ax, cuty)
pylab.title('layout band structure')
pylab.ylabel('Energy (|t|)')
pylab.xlabel('$k_y$ ($\pi$/a)')
pylab.xticks([0, splitCutY.shape[1]/2, splitCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
if forceLims:
    ax.set_ylim([axmin, axmax])

        



pylab.suptitle('layout = ' + testCell.type)
pylab.show()


MAT = generate_layout_Bloch_matrix(grapheneCell, 0,0)
#print MAT
#print BlochMat





pylab.figure(2)
pylab.clf()
####layout graph
ax1 = pylab.subplot(1,2,1, adjustable='box', aspect=1)
pylab.cla()
splitLattice.draw_resonator_lattice(ax1, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
xs = splitLattice.coords[:,0]
ys = splitLattice.coords[:,1]
pylab.sca(ax1)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#ax.set_aspect('equal')
ax1.axis('off')
pylab.title('layout')

ax2 = pylab.subplot(1,2,2, adjustable='box', aspect=1)
pylab.cla()
splitLattice.draw_SDlinks(ax2, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=1)
pylab.sca(ax2)
pylab.scatter(splitLattice.SDx, splitLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
pylab.title('FW lattice')
ax2.set_aspect('equal')
ax2.axis('off')

pylab.suptitle('layout = ' + testCell.type)
pylab.show()



pylab.figure(3)
pylab.clf()

numSurfPoints = 300
#kx_x, ky_y, cutx = compute_layout_band_structure(grapheneCell, -2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = 'FW')
#kx_y, ky_y, cuty = compute_layout_band_structure(grapheneCell, 0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = 'FW')
kx_x, ky_y, cutx = compute_layout_band_structure(testCell, -2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = 'FW')
kx_y, ky_y, cuty = compute_layout_band_structure(testCell, 0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = 'FW')



ax = pylab.subplot(1,1,1, adjustable='box', aspect=80)
pylab.cla()
#plot_band_cut_mc(ax, cuty)
plot_any_band_cut(ax, cuty)
pylab.title('layout band structure')
pylab.ylabel('Energy (|t|)')
pylab.xlabel('$k_y$ ($\pi$/a)')
pylab.xticks([0, splitCutY.shape[1]/2, splitCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
if forceLims:
    ax.set_ylim([-3.2, 3.2])
    
pylab.tight_layout()
pylab.show()





