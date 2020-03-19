#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 21 10:25:51 2018

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





recompute = True
#recompute = False










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



def make_layout_vid(GenLayout, startInd =0, stopInd = -1, figNum = 8, xsize = 14, ysize = 4):
    
    
    if stopInd <0:
        #set a break point that will not be reached
        breakPoint = len(GenLayout.Es)+1
    else:
        breakPoint = stopInd
        
    if startInd > len(GenLayout.Es):
        raise ValueError('dont have this many eigenvectors')
        
    
    folder = GenLayout.name + '_' + GenLayout.modeType

    
    currDir = os.getcwd()
    saveDir = os.path.join(currDir,folder)
    if os.path.isdir(saveDir):
        pass
    else:
        os.mkdir(saveDir)
    
    
    print(folder + '\n')
        
    fig = pylab.figure(figNum)
    #######################################################
    for eigNum in range(startInd, len(GenLayout.Es)):
        print(eigNum)
        
        eigAmps = GenLayout.Psis[:,GenLayout.Eorder[eigNum]]
        
        
        
        fig.clf()
#        fig.figsize = (100,100)
        ax1 = pylab.subplot(1,2,1)
        xs = scipy.arange(0,len(GenLayout.Es),1)
        pylab.plot(xs, GenLayout.Es[GenLayout.Eorder], 'b.')
        
        pylab.plot(xs[eigNum],GenLayout.Es[GenLayout.Eorder[eigNum]], color = 'firebrick' , marker = '.', markersize = '10' )
        pylab.title(' full tight-binding Spectrum')
        ax1.set_ylim([-2.2, 4.1])
        
        ax1 = pylab.subplot(1,2,2)
        titleStr = 'eigenvector weight : ' + str(eigNum) + ' ; E = ' + str(GenLayout.Es[GenLayout.Eorder[eigNum]])
##        GenLayout.plot_layout_state(eigAmps, ax1, title = titleStr, colorbar = True, plot_links = True, cmap = 'Wistia')
#        GenLayout.plot_layout_state(eigAmps, ax1, title = titleStr, colorbar = True, plot_links = False, cmap = 'Wistia')
#        GenLayout.draw_SDlinks(ax1, color = 'firebrick', linewidth = 0.5, minus_links = True, minus_color = 'dodgerblue', NaNs = True)
        ax1.set_aspect('auto')
        
        
        if GenLayout.modeType == 'FW':
            pylab.cla()
#            GenLayout.draw_SDlinks(ax1, color = FWlinkColor, linewidth = 2, minus_links = True, minus_color = 'gold', alpha=FWlinkAlpha)
#            pylab.scatter(GenLayout.SDx, GenLayout.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 3, alpha=FWsiteAlpha)
            GenLayout.draw_SDlinks(ax1, color = FWlinkColor, linewidth = 2, minus_links = True, minus_color = 'gold', alpha=1)
            pylab.scatter(GenLayout.SDx, GenLayout.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 3, alpha=1)
        else:
            pylab.cla()
            GenLayout.draw_SDlinks(ax1, color = HWlinkColor, linewidth = 2, minus_links = True, minus_color = HWminusLinkColor, alpha=HWlinkAlpha)
            pylab.scatter(GenLayout.SDx, GenLayout.SDy ,c =  HWsiteColor, s = smallCdefault, marker = 'o', edgecolors = HWsiteEdgeColor, zorder = 3, alpha=HWsiteAlpha)
#            GenLayout.draw_SDlinks(ax1, color = HWlinkColor, linewidth = 2, minus_links = True, minus_color = HWminusLinkColor, alpha=1)
#            pylab.scatter(GenLayout.SDx, GenLayout.SDy ,c =  HWsiteColor, s = smallCdefault, marker = 'o', edgecolors = HWsiteEdgeColor, zorder = 3, alpha=1)
        GenLayout.plot_layout_state(eigAmps, ax1, title = titleStr, colorbar = True, plot_links = False, cmap = 'Wistia', zorder = 5)
        
        

        fig.set_size_inches(xsize, ysize)
        
        pylab.tight_layout()
    
        fig_name = GenLayout.name + '_' + GenLayout.modeType +  '_eig_' + str(eigNum) + '.png' 
        fig_path = os.path.join(saveDir, fig_name)
        pylab.savefig(fig_path)
        
        if eigNum == breakPoint:
            break
        
    pylab.show()
    
    return




def itter_generate_C60():
        '''
        full general (I hope) itteration function to generate the lattice
        by propagation of the tiling rule.
        
        Should be abe to handle and flat or hyperbolic tiling.
        Spherical tilings tend to crash eventudally because the tiling rule ends
        '''
        if test.itter == 0:
            #set value for the 1st itteration, whichis about to be done
            test.gon = 6
        elif test.itter > 1:
            #set it to 6 so it allocates more space. Will be messing with this
            test.gon = 6
            
        #make the last ring look a little better
        if test.itter == 3:
            test.radius = 1.25
        
        print('finished itteration equals = ' + str(test.itter))
        print('attempting itteration ' + str(test.itter+1))
        numVertices = len(test.points[test.itter])
        numFaces = len(test.points[test.itter])
        numRadials = test.radials[test.itter].shape[0]
        newItter = test.itter+1
        
        #find the new radial lines
        outgoing_radials_at_each_vertex = numpy.zeros(numVertices) #array to tell how manyrials emerge from each point
        for ind in range(0, numVertices):
            #see if this point is the end of a previous itterations radial
            currentAngle = numpy.mod(test.points[test.itter][ind], 2*numpy.pi)
            incoming = len(numpy.where(numpy.mod(test.radials[test.itter][:,1],2*numpy.pi)==currentAngle)[0])
            outgoing_radials_at_each_vertex[ind] = test.vertex-2-incoming
                
        print(outgoing_radials_at_each_vertex[0:7])
        print('\n')
        
        
        #find the number of new points
        currentIndexP = 0
        currentIndexR = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
                #do each gon that only has a corner touching the previous itteration (mostly)
                #move to the next site
                currentIndexP = currentIndexP + (test.gon-2)
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
            currentIndexP = currentIndexP - 1
            
        if test.itter < 1:
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = int(currentIndexP)
        elif test.itter == 1:
            numNewRadials = 10
            numNewPoints = 20
        elif test.itter == 2:
            numNewRadials = 10
            numNewPoints = 15
        elif test.itter == 3:
            numNewRadials = 5
            numNewPoints = 5
#        print numNewRadials
#        print numNewPoints
        
        
        #find the l ocations of the new azimuthal vertices
#        self.points[newItter ] = 0 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) #don't want 2pi included
#        self.points[newItter ] = -self.alpha/2 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints)
        if test.gon == 7:
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 5:
            if test.itter == 0:
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 0.
            else:
#                self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 1.0*self.alpha/2.
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 3:  
            test.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - test.alpha/2/1.5 #old as of 5-23-18
        else:
#            self.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - self.alpha/2/1.5 #old as of 5-23-18
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        
        #make the angles prettier
        if test.itter == 3:
            test.points[newItter ] = test.points[newItter ] + test.alpha*0.29
        
        if test.radius_method == 'lin':
            test.radii[newItter ] = test.radii[test.itter] + test.radius
        if test.radius_method == 'exp':
            test.radii[newItter ] = test.radius*numpy.exp(test.itter+1)
        if test.radius_method == 'Mattias':
            MattiasRadii = [1.0, 3.2, 5.1, 7.,8.]
            test.radii[newItter ] = test.radius*MattiasRadii[newItter]
        
        #store the cartesian locations of the new vertices
        xs = test.radii[newItter]*numpy.cos(test.points[newItter])
        ys = test.radii[newItter]*numpy.sin(test.points[newItter])
        newCoords = numpy.zeros((len(xs),2))
        newCoords[:,0] = xs
        newCoords[:,1] = ys
        test.coords = numpy.concatenate((test.coords, newCoords))
        
        
        test.radials[newItter] = numpy.zeros((numNewRadials,2))
        currentIndexP = 0
        currentIndexR = 0
        gonind = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
#                print ind
#                print currentIndexR
##                print currentIndexP
#                print '\n'
                #do each gon that only has a corner touching the previous itteration (mostly)
                test.radials[newItter][currentIndexR, :] = numpy.asarray([test.points[test.itter][ind], test.points[newItter][currentIndexP]])
                
                gonind = gonind+1
                if test.itter ==1:
                    if numpy.mod(gonind,2) == 0:
                        test.gon = 5
                    else:
                        test.gon = 6
                if test.itter ==2:
                    if numpy.mod(gonind,2) == 0:
                        test.gon = 6
                    else:
                        test.gon = 5
#                if test.itter ==2:
#                    if numpy.mod(gonind,2) == 0:
#                        test.gon = 5
#                    else:
#                        test.gon = 6
#                
                #ove to the next site
                currentIndexP = numpy.mod(currentIndexP + (test.gon-2), numNewPoints)
#                currentIndexP = currentIndexP + (self.gon-2)
                currentIndexR = currentIndexR + 1
            
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
#            currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)
            if ind == 0 and int(outgoing_radials_at_each_vertex[ind]) ==0:
                #haven't actually stepped, so I don't want to remove anything
                pass
            else:
                currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)

                

        
        test.itter = test.itter+1
        return


def itter_generate_C72():
        '''
        full general (I hope) itteration function to generate the lattice
        by propagation of the tiling rule.
        
        Should be abe to handle and flat or hyperbolic tiling.
        Spherical tilings tend to crash eventudally because the tiling rule ends
        '''
        if test.itter == 0:
            #set value for the 1st itteration, whichis about to be done
            test.gon = 6
        elif test.itter > 1:
            #set it to 6 so it allocates more space. Will be messing with this
            test.gon = 6
        
        print('finished itteration equals = ' + str(test.itter))
        print('attempting itteration ' + str(test.itter+1))
        numVertices = len(test.points[test.itter])
        numFaces = len(test.points[test.itter])
        numRadials = test.radials[test.itter].shape[0]
        newItter = test.itter+1
        
        #find the new radial lines
        outgoing_radials_at_each_vertex = numpy.zeros(numVertices) #array to tell how manyrials emerge from each point
        for ind in range(0, numVertices):
            #see if this point is the end of a previous itterations radial
            currentAngle = numpy.mod(test.points[test.itter][ind], 2*numpy.pi)
            incoming = len(numpy.where(numpy.mod(test.radials[test.itter][:,1],2*numpy.pi)==currentAngle)[0])
            outgoing_radials_at_each_vertex[ind] = test.vertex-2-incoming
                
        print(outgoing_radials_at_each_vertex[0:7])
        print('\n')
        
        
        #find the number of new points
        currentIndexP = 0
        currentIndexR = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
                #do each gon that only has a corner touching the previous itteration (mostly)
                #move to the next site
                currentIndexP = currentIndexP + (test.gon-2)
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
            currentIndexP = currentIndexP - 1
            
        if test.itter < 1:
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = int(currentIndexP)
        elif test.itter == 1:
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = 24
        elif test.itter == 2:
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = 18
        elif test.itter == 3:
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = int(currentIndexP)
#        print numNewRadials
#        print numNewPoints
        
        
        #find the l ocations of the new azimuthal vertices
#        self.points[newItter ] = 0 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) #don't want 2pi included
#        self.points[newItter ] = -self.alpha/2 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints)
        if test.gon == 7:
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 5:
            if test.itter == 0:
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 0.
            else:
#                self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 1.0*self.alpha/2.
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 3:  
            test.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - test.alpha/2/1.5 #old as of 5-23-18
        else:
#            self.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - self.alpha/2/1.5 #old as of 5-23-18
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        
        if test.itter == 3:
            test.points[newItter ] = test.points[newItter ] + test.alpha*0.3
        
        if test.radius_method == 'lin':
            test.radii[newItter ] = test.radii[test.itter] + test.radius
        if test.radius_method == 'exp':
            test.radii[newItter ] = test.radius*numpy.exp(test.itter+1)
        if test.radius_method == 'Mattias':
#            MattiasRadii = [1.0, 3, 4.5, 7.,8.]
#            MattiasRadii = [1.0, 3.1, 4.8, 7.,8.]
            MattiasRadii = [1.0, 3.2, 5.1, 7.,8.]
            test.radii[newItter ] = test.radius*MattiasRadii[newItter]
        
        #store the cartesian locations of the new vertices
        xs = test.radii[newItter]*numpy.cos(test.points[newItter])
        ys = test.radii[newItter]*numpy.sin(test.points[newItter])
        newCoords = numpy.zeros((len(xs),2))
        newCoords[:,0] = xs
        newCoords[:,1] = ys
        test.coords = numpy.concatenate((test.coords, newCoords))
        
        
        test.radials[newItter] = numpy.zeros((numNewRadials,2))
        currentIndexP = 0
        currentIndexR = 0
        gonind = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
#                print ind
#                print currentIndexR
##                print currentIndexP
#                print '\n'
                #do each gon that only has a corner touching the previous itteration (mostly)
                test.radials[newItter][currentIndexR, :] = numpy.asarray([test.points[test.itter][ind], test.points[newItter][currentIndexP]])
                
                gonind = gonind+1
                if test.itter ==1:
                    if numpy.mod(gonind,2) == 0:
                        test.gon = 5
                    else:
                        test.gon = 6
                if test.itter ==2:
                    if numpy.mod(gonind,2) == 0:
                        test.gon = 6
                    else:
                        test.gon = 5
                if test.itter ==3:
                    test.gon = 6
#                
                #ove to the next site
                currentIndexP = numpy.mod(currentIndexP + (test.gon-2), numNewPoints)
#                currentIndexP = currentIndexP + (self.gon-2)
                currentIndexR = currentIndexR + 1
            
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
#            currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)
            if ind == 0 and int(outgoing_radials_at_each_vertex[ind]) ==0:
                #haven't actually stepped, so I don't want to remove anything
                pass
            else:
                currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)

                

        
        test.itter = test.itter+1
        return





def itter_generate_C84_D7():
        '''
        full general (I hope) itteration function to generate the lattice
        by propagation of the tiling rule.
        
        Should be abe to handle and flat or hyperbolic tiling.
        Spherical tilings tend to crash eventudally because the tiling rule ends
        '''
        if test.itter == 0:
            #set value for the 1st itteration, whichis about to be done
            test.gon = 6
        elif test.itter > 1:
            #set it to 6 so it allocates more space. Will be messing with this
            test.gon = 6
        
        print('finished itteration equals = ' + str(test.itter))
        print('attempting itteration ' + str(test.itter+1))
        numVertices = len(test.points[test.itter])
        numFaces = len(test.points[test.itter])
        numRadials = test.radials[test.itter].shape[0]
        newItter = test.itter+1
        
        #find the new radial lines
        outgoing_radials_at_each_vertex = numpy.zeros(numVertices) #array to tell how manyrials emerge from each point
        for ind in range(0, numVertices):
            #see if this point is the end of a previous itterations radial
            currentAngle = numpy.mod(test.points[test.itter][ind], 2*numpy.pi)
            incoming = len(numpy.where(numpy.mod(test.radials[test.itter][:,1],2*numpy.pi)==currentAngle)[0])
            outgoing_radials_at_each_vertex[ind] = test.vertex-2-incoming
                
        print(outgoing_radials_at_each_vertex[0:7])
        print('\n')
        
        
        #find the number of new points
        currentIndexP = 0
        currentIndexR = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
                #do each gon that only has a corner touching the previous itteration (mostly)
                #move to the next site
                currentIndexP = currentIndexP + (test.gon-2)
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
            currentIndexP = currentIndexP - 1
            
        if test.itter < 1:
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = int(currentIndexP)
        elif test.itter == 1:
#            numNewRadials = 10
#            numNewPoints = 20
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = 28
        elif test.itter == 2:
#            numNewRadials = 10
#            numNewPoints = 15
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = 21
        elif test.itter == 3:
#            numNewRadials = 14
#            numNewPoints = 10
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = int(currentIndexP)
        
        
        #find the l ocations of the new azimuthal vertices
#        self.points[newItter ] = 0 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) #don't want 2pi included
#        self.points[newItter ] = -self.alpha/2 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints)
        if test.gon == 7:
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 5:
            if test.itter == 0:
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 0.
            else:
#                self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 1.0*self.alpha/2.
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 3:  
            test.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - test.alpha/2/1.5 #old as of 5-23-18
        else:
#            self.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - self.alpha/2/1.5 #old as of 5-23-18
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        
#        if test.itter == 1:
#            test.points[newItter ] = test.points[newItter ] + test.alpha*0.1
#        if test.itter == 2:
#            test.points[newItter ] = test.points[newItter ] - test.alpha*0.05
        if test.itter == 3:
            test.points[newItter ] = test.points[newItter ] + test.alpha*0.3
        
        if test.radius_method == 'lin':
            test.radii[newItter ] = test.radii[test.itter] + test.radius
        if test.radius_method == 'exp':
            test.radii[newItter ] = test.radius*numpy.exp(test.itter+1)
        if test.radius_method == 'Mattias':
#            MattiasRadii = [1.0, 3, 4.5, 7.,8.]
#            MattiasRadii = [1.0, 3.1, 4.8, 7.,8.]
            MattiasRadii = [1.0, 3.2, 5.1, 7.,8.]
            test.radii[newItter ] = test.radius*MattiasRadii[newItter]
        
        #store the cartesian locations of the new vertices
        xs = test.radii[newItter]*numpy.cos(test.points[newItter])
        ys = test.radii[newItter]*numpy.sin(test.points[newItter])
        newCoords = numpy.zeros((len(xs),2))
        newCoords[:,0] = xs
        newCoords[:,1] = ys
        test.coords = numpy.concatenate((test.coords, newCoords))
        
        
        test.radials[newItter] = numpy.zeros((numNewRadials,2))
        currentIndexP = 0
        currentIndexR = 0
        gonind = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
#                print ind
#                print currentIndexR
##                print currentIndexP
#                print '\n'
                #do each gon that only has a corner touching the previous itteration (mostly)
                test.radials[newItter][currentIndexR, :] = numpy.asarray([test.points[test.itter][ind], test.points[newItter][currentIndexP]])
                
                gonind = gonind+1
                if test.itter ==1:
                    if numpy.mod(gonind,2) == 0:
                        test.gon = 5
                    else:
                        test.gon = 6
                if test.itter ==2:
                    if numpy.mod(gonind,2) == 0:
                        test.gon = 6
                    else:
                        test.gon = 5
                if test.itter ==3:
                    test.gon = 6
#                if test.itter == 3:
#                    test.gon == 6
#                if test.itter ==2:
#                    if numpy.mod(gonind,2) == 0:
#                        test.gon = 5
#                    else:
#                        test.gon = 6
#                
                #ove to the next site
                currentIndexP = numpy.mod(currentIndexP + (test.gon-2), numNewPoints)
#                currentIndexP = currentIndexP + (self.gon-2)
                currentIndexR = currentIndexR + 1
            
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
#            currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)
            if ind == 0 and int(outgoing_radials_at_each_vertex[ind]) ==0:
                #haven't actually stepped, so I don't want to remove anything
                pass
            else:
                currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)

                

        
        test.itter = test.itter+1
        return



def itter_generate_C84():
        '''
        full general (I hope) itteration function to generate the lattice
        by propagation of the tiling rule.
        
        Should be abe to handle and flat or hyperbolic tiling.
        Spherical tilings tend to crash eventudally because the tiling rule ends
        '''
        print('running generate_c60 for ' + str(test.itter))
        if test.itter == 0:
            #set value for the 1st itteration, whichis about to be done
            test.gon = 6
        elif test.itter > 1:
            #set it to 6 so it allocates more space. Will be messing with this
            test.gon = 6
            
        #make the last ring look a little better
        if test.itter == 3:
            test.radius = 1.
        
        print('finished itteration equals = ' + str(test.itter))
        print('attempting itteration ' + str(test.itter+1))
        numVertices = len(test.points[test.itter])
        numFaces = len(test.points[test.itter])
        numRadials = test.radials[test.itter].shape[0]
        newItter = test.itter+1
        
        #find the new radial lines
        outgoing_radials_at_each_vertex = numpy.zeros(numVertices) #array to tell how manyrials emerge from each point
        for ind in range(0, numVertices):
            #see if this point is the end of a previous itterations radial
            currentAngle = numpy.mod(test.points[test.itter][ind], 2*numpy.pi)
            incoming = len(numpy.where(numpy.mod(test.radials[test.itter][:,1],2*numpy.pi)==currentAngle)[0])
            outgoing_radials_at_each_vertex[ind] = test.vertex-2-incoming
                
        print(outgoing_radials_at_each_vertex[0:7])
        print('\n')
        
        
        #find the number of new points
        currentIndexP = 0
        currentIndexR = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
                #do each gon that only has a corner touching the previous itteration (mostly)
                #move to the next site
                currentIndexP = currentIndexP + (test.gon-2)
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
            currentIndexP = currentIndexP - 1
            
        if test.itter < 1:
#            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
#            numNewPoints = int(currentIndexP)
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = 15
        elif test.itter == 1:
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints = 21
        elif test.itter == 2:
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints =  24 #21
        elif test.itter == 3:
            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
            numNewPoints =  18  #int(currentIndexP)
        elif test.itter == 4:
            numNewRadials = 6
            numNewPoints = 6
#            numNewRadials = 30
#            numNewPoints = 24
#        print numNewRadials
#        print numNewPoints
        
        
        #find the l ocations of the new azimuthal vertices
#        self.points[newItter ] = 0 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) #don't want 2pi included
#        self.points[newItter ] = -self.alpha/2 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints)
        if test.gon == 7:
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 5:
            if test.itter == 0:
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 0.
            else:
#                self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 1.0*self.alpha/2.
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 3:  
            test.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - test.alpha/2/1.5 #old as of 5-23-18
        else:
#            self.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - self.alpha/2/1.5 #old as of 5-23-18
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        
        if test.itter == 0:
            test.points[newItter ] = test.points[newItter ] - test.alpha*0.1
        if test.itter == 1:
            test.points[newItter ] = test.points[newItter ] - test.alpha*0.05
        if test.itter == 2:
            test.points[newItter ] = test.points[newItter ] - test.alpha*0.05
        if test.itter == 3:
            test.points[newItter ] = test.points[newItter ] - test.alpha*0.2
        if test.itter == 4:
            test.points[newItter ] = test.points[newItter ] + test.alpha*0.3
        
        if test.radius_method == 'lin':
            test.radii[newItter ] = test.radii[test.itter] + test.radius
        if test.radius_method == 'exp':
            test.radii[newItter ] = test.radius*numpy.exp(test.itter+1)
        if test.radius_method == 'Mattias':
#            MattiasRadii = [1.0, 3, 4.5, 7.,8.]
#            MattiasRadii = [1.0, 3.1, 4.8, 7.,8.]
            MattiasRadii = [1.0, 3.2, 5.1, 7.,8.]
            test.radii[newItter ] = test.radius*MattiasRadii[newItter]
        
        #store the cartesian locations of the new vertices
        xs = test.radii[newItter]*numpy.cos(test.points[newItter])
        ys = test.radii[newItter]*numpy.sin(test.points[newItter])
        newCoords = numpy.zeros((len(xs),2))
        newCoords[:,0] = xs
        newCoords[:,1] = ys
        test.coords = numpy.concatenate((test.coords, newCoords))
        
        
        test.radials[newItter] = numpy.zeros((numNewRadials,2))
        currentIndexP = 0
        currentIndexR = 0
        gonind = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
#                print ind
#                print currentIndexR
##                print currentIndexP
#                print '\n'
                #do each gon that only has a corner touching the previous itteration (mostly)
                test.radials[newItter][currentIndexR, :] = numpy.asarray([test.points[test.itter][ind], test.points[newItter][currentIndexP]])
                
                gonind = gonind+1
                if test.itter ==0:
                    if numpy.mod(gonind,2) == 0:
                        test.gon = 5
                    else:
                        test.gon = 6
                if test.itter ==1:
                    test.gon = 6
                if test.itter ==2:
#                    test.gon = 6
                    if numpy.mod(gonind,4) == 3:
                        test.gon = 5
                    else:
                        test.gon = 6
                if test.itter ==3:
                    if numpy.mod(gonind,2) == 0:
                        test.gon = 5
                    else:
                        test.gon = 6
                if test.itter ==4:
                    test.gon = 6

#                
                #ove to the next site
                currentIndexP = numpy.mod(currentIndexP + (test.gon-2), numNewPoints)
#                currentIndexP = currentIndexP + (self.gon-2)
                currentIndexR = currentIndexR + 1
            
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
#            currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)
            if ind == 0 and int(outgoing_radials_at_each_vertex[ind]) ==0:
                #haven't actually stepped, so I don't want to remove anything
                pass
            else:
                currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)

                

        
        test.itter = test.itter+1
        return
    
    
def itter_generate_interface_tiling(gon2, transition):
        '''
        full general (I hope) itteration function to generate the lattice
        by propagation of the tiling rule.
        
        Should be abe to handle and flat or hyperbolic tiling.
        Spherical tilings tend to crash eventudally because the tiling rule ends
        '''
        if test.itter < transition:
            pass
            #use initial value
        else:
            test.gon = gon2
        
        print('finished itteration equals = ' + str(test.itter))
        print('attempting itteration ' + str(test.itter+1))
        numVertices = len(test.points[test.itter])
        numFaces = len(test.points[test.itter])
        numRadials = test.radials[test.itter].shape[0]
        newItter = test.itter+1
        
        #find the new radial lines
        outgoing_radials_at_each_vertex = numpy.zeros(numVertices) #array to tell how manyrials emerge from each point
        for ind in range(0, numVertices):
            #see if this point is the end of a previous itterations radial
            currentAngle = numpy.mod(test.points[test.itter][ind], 2*numpy.pi)
            incoming = len(numpy.where(numpy.mod(test.radials[test.itter][:,1],2*numpy.pi)==currentAngle)[0])
            outgoing_radials_at_each_vertex[ind] = test.vertex-2-incoming
                
        print(outgoing_radials_at_each_vertex[0:7])
        print('\n')
        
        
        #find the number of new points
        currentIndexP = 0
        currentIndexR = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
                #do each gon that only has a corner touching the previous itteration (mostly)
                #move to the next site
                currentIndexP = currentIndexP + (test.gon-2)
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
            currentIndexP = currentIndexP - 1
            
        numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
        numNewPoints = int(currentIndexP)
#        print numNewRadials
#        print numNewPoints
#        if test.itter < 1:
#            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
#            numNewPoints = int(currentIndexP)
#        elif test.itter == 1:
##            numNewRadials = 10
##            numNewPoints = 20
#            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
#            numNewPoints = 28
#        elif test.itter == 2:
##            numNewRadials = 10
##            numNewPoints = 15
#            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
#            numNewPoints = 21
#        elif test.itter == 3:
##            numNewRadials = 14
##            numNewPoints = 10
#            numNewRadials = int(numpy.sum(outgoing_radials_at_each_vertex))
#            numNewPoints = int(currentIndexP)
        
        
        #find the l ocations of the new azimuthal vertices
#        self.points[newItter ] = 0 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) #don't want 2pi included
#        self.points[newItter ] = -self.alpha/2 + scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints)
        if test.gon == 7:
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 5:
            if test.itter == 0:
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 0.
            else:
#                self.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) + 1.0*self.alpha/2.
                test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        elif test.gon == 3:  
            test.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - test.alpha/2/1.5 #old as of 5-23-18
        else:
#            self.points[newItter ] =  scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - self.alpha/2/1.5 #old as of 5-23-18
            test.points[newItter ] = scipy.arange(0,2*numpy.pi, 2*numpy.pi/numNewPoints) - 0*test.alpha/2 + 0.25*(test.itter)**0.5*test.alpha/2.
        
##        if test.itter == 1:
##            test.points[newItter ] = test.points[newItter ] + test.alpha*0.1
##        if test.itter == 2:
##            test.points[newItter ] = test.points[newItter ] - test.alpha*0.05
#        if test.itter == 3:
#            test.points[newItter ] = test.points[newItter ] + test.alpha*0.3
#        
        if test.radius_method == 'lin':
            test.radii[newItter ] = test.radii[test.itter] + test.radius
        if test.radius_method == 'exp':
            test.radii[newItter ] = test.radius*numpy.exp(test.itter+1)
        if test.radius_method == 'Mattias':
#            MattiasRadii = [1.0, 3, 4.5, 7.,8.]
#            MattiasRadii = [1.0, 3.1, 4.8, 7.,8.]
            MattiasRadii = [1.0, 3.2, 5.1, 7.,8.]
            test.radii[newItter ] = test.radius*MattiasRadii[newItter]
        
        #store the cartesian locations of the new vertices
        xs = test.radii[newItter]*numpy.cos(test.points[newItter])
        ys = test.radii[newItter]*numpy.sin(test.points[newItter])
        newCoords = numpy.zeros((len(xs),2))
        newCoords[:,0] = xs
        newCoords[:,1] = ys
        test.coords = numpy.concatenate((test.coords, newCoords))
        
        
        test.radials[newItter] = numpy.zeros((numNewRadials,2))
        currentIndexP = 0
        currentIndexR = 0
        gonind = 0
        for ind in range(0, numVertices):
            for rad in range(0, int(outgoing_radials_at_each_vertex[ind])):
#                print ind
#                print currentIndexR
##                print currentIndexP
#                print '\n'
                #do each gon that only has a corner touching the previous itteration (mostly)
                test.radials[newItter][currentIndexR, :] = numpy.asarray([test.points[test.itter][ind], test.points[newItter][currentIndexP]])
                
#                gonind = gonind+1
#                if test.itter ==1:
#                    if numpy.mod(gonind,2) == 0:
#                        test.gon = 5
#                    else:
#                        test.gon = 6
#                if test.itter ==2:
#                    if numpy.mod(gonind,2) == 0:
#                        test.gon = 6
#                    else:
#                        test.gon = 5
#                if test.itter ==3:
#                    test.gon = 6
##                if test.itter == 3:
##                    test.gon == 6
##                if test.itter ==2:
##                    if numpy.mod(gonind,2) == 0:
##                        test.gon = 5
##                    else:
##                        test.gon = 6
#                
                #ove to the next site
                currentIndexP = numpy.mod(currentIndexP + (test.gon-2), numNewPoints)
#                currentIndexP = currentIndexP + (self.gon-2)
                currentIndexR = currentIndexR + 1
            
            #back up one site, because the next gon touches the previous itteration at a whole face, so we need a smaller step
#            currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)
            if ind == 0 and int(outgoing_radials_at_each_vertex[ind]) ==0:
                #haven't actually stepped, so I don't want to remove anything
                pass
            else:
                currentIndexP = numpy.mod(currentIndexP - 1, numNewPoints)

                

        
        test.itter = test.itter+1
        return


def draw_arc(ax, resonator, linewidth = 3, color = 'firebrick', numPoints = 100, direction = 1):
    x1,y1,x2,y2 = resonator
    
    xCenter = (x1+x2)/2
    yCenter = (y1+y2)/2
    
    r = 0.5* numpy.sqrt((x1-x2)**2 + (y1-y2)**2)
    
    theta0 = numpy.arctan((y1-yCenter) / (x1-xCenter))
    if direction == 1:
        thetas = numpy.linspace(0, numpy.pi, numPoints) + theta0
    if direction == -1:
        thetas = numpy.linspace(0, -numpy.pi, numPoints) + theta0
    
    arcXs = numpy.cos(thetas)*r + xCenter
    arcYs = numpy.sin(thetas)*r + yCenter
    
    pylab.sca(ax)
    pylab.plot(arcXs, arcYs, linewidth = linewidth, color = color)
#    pylab.scatter(xCenter, yCenter)
    
    return

def draw_variable_arc(ax, resonator, d = 1, linewidth = 3, color = 'firebrick', numPoints = 100, direction = 1):
    x1,y1,x2,y2 = resonator
    
    
    
    X = 0.5* numpy.sqrt((x1-x2)**2 + (y1-y2)**2)
    
    theta = 2*numpy.arctan(d/X)
    
    r = X/(numpy.sin(theta))
    
    xCenter = 0
    yCenter = 0

    
    xMid = (x1+x2)/2
    yMid = (y1+y2)/2
    print(xMid)
    print(yMid)
    
    bondOrientation = numpy.arctan((y1-yMid) / (x1-xMid))
    print(180*bondOrientation/numpy.pi)
    
    if direction == 1:
        thetas = numpy.linspace(-theta, theta, numPoints) + bondOrientation + numpy.pi/2
    if direction == -1:
        thetas = numpy.linspace(theta, -theta, numPoints) + bondOrientation - numpy.pi/2
#        thetas = numpy.linspace(theta, -theta, numPoints) + numpy.pi/2
    
#    arcXs = numpy.cos(thetas)*r + xCenter
#    arcYs = numpy.sin(thetas)*r + yCenter
    arcXs = numpy.cos(thetas)*r
    arcYs = numpy.sin(thetas)*r
#    arcXs = numpy.cos(thetas)*r 
#    arcYs = numpy.sin(thetas)*r
    
    newxMid = (arcXs[0]+arcXs[-1])/2
    newyMid = (arcYs[0]+arcYs[-1])/2
    
    arcXs = arcXs + xMid - newxMid
    arcYs = arcYs + yMid - newyMid
    
    pylab.sca(ax)
    pylab.plot(arcXs, arcYs, linewidth = linewidth, color = color)
#    pylab.scatter(xCenter, yCenter)
    
    return



#test = PlanarLayout(gon = 5, vertex = 3, side =1, radius_method = 'lin')
#test.populate(2, resonatorsOnly=False)
        
xsize = 11.25
ysize = 6

if recompute:
    #####
    #bucky ball
    test = PlanarLayout(gon = 5, vertex = 3, side =1, radius_method = 'lin')
    itter_generate_C60()
    itter_generate_C60()
    itter_generate_C60()
    itter_generate_C60()
    test.generate_semiduals()
    C60 = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  'C60')
    C60_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  'C60')
    
    ######
    ##C72
    test = PlanarLayout(gon = 6, vertex = 3, side =1, radius_method = 'lin')
    itter_generate_C72()
    itter_generate_C72()
    itter_generate_C72()
    itter_generate_C72()
    test.generate_semiduals()
    C72 = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  'C72')
    C72_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  'C72')
    
    ######
    ##C84, 7 fold version
    test = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin')
    itter_generate_C84_D7()
    itter_generate_C84_D7()
    itter_generate_C84_D7()
    itter_generate_C84_D7()
    test.generate_semiduals()
    C84 = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  'C84')
    C84_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  'C84')
    
    
    ######
    ##C84, no heptagon version
    test = PlanarLayout(gon = 6, vertex = 3, side =1, radius_method = 'lin')
    itter_generate_C84()
    itter_generate_C84()
    itter_generate_C84()
    itter_generate_C84()
    #itter_generate_C84()
    test.generate_semiduals()
    baseLayout = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  'C84_6_intermediate')
    newResonators = numpy.zeros((3,4))
    #newResonators[0,:] = numpy.concatenate(( baseLayout.coords[1,:], baseLayout.coords[34,:]))
    #newResonators[1,:] = numpy.concatenate(( baseLayout.coords[5,:], baseLayout.coords[-7,:]))
    #newResonators[2,:] = numpy.concatenate(( baseLayout.coords[52,:], baseLayout.coords[-3,:]))
    ind1 = 1
    ind2 = 14
    newResonators[0,:] = test.radii[4]* numpy.asarray([numpy.cos(test.points[4][ind1]), numpy.sin(test.points[4][ind1]), numpy.cos(test.points[4][ind2]), numpy.sin(test.points[4][ind2])   ])
    ind1 = 2
    ind2 = 7
    newResonators[1,:] = test.radii[4]* numpy.asarray([numpy.cos(test.points[4][ind1]), numpy.sin(test.points[4][ind1]), numpy.cos(test.points[4][ind2]), numpy.sin(test.points[4][ind2])   ])
    ind1 = 8
    ind2 = 13
    newResonators[2,:] = test.radii[4]* numpy.asarray([numpy.cos(test.points[4][ind1]), numpy.sin(test.points[4][ind1]), numpy.cos(test.points[4][ind2]), numpy.sin(test.points[4][ind2])   ])
    
    
    allResonators = numpy.concatenate((baseLayout.resonators, newResonators))
    
    C84_6 = GeneralLayout(allResonators,  modeType = 'FW', name =  'C84_6')
    C84_6_HW = GeneralLayout(allResonators,  modeType = 'HW', name =  'C84_6')
    
    extras = GeneralLayout(newResonators,  modeType = 'FW', name =  'C84_6_extras')
    #
    


######
##interface
#gon1 = 7
#gon2 = 6
#transition = 2
#test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#test.generate_semiduals()
#name_str = 'hybrid_' + str(gon1) + '_' +str(gon2) + '_' + str(transition) 
#hybrid = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
#hybrid_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)

#gon1 = 9
#gon2 = 6
#transition = 2
#test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#test.generate_semiduals()
#name_str = 'hybrid_' + str(gon1) + '_' +str(gon2) + '_' + str(transition) 
#hybrid = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
#hybrid_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)


#gon1 = 6
#gon2 = 5
#transition = 2
#test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
##itter_generate_interface_tiling(gon2, transition)
##itter_generate_interface_tiling(gon2, transition)
#test.generate_semiduals()
#name_str = 'hybrid_' + str(gon1) + '_' +str(gon2) + '_' + str(transition) 
#hybrid = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
#hybrid_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)

#gon1 = 7
#gon2 = 5
#transition = 2
#test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
##itter_generate_interface_tiling(gon2, transition)
##itter_generate_interface_tiling(gon2, transition)
#test.generate_semiduals()
#name_str = 'hybrid_' + str(gon1) + '_' +str(gon2) + '_' + str(transition) 
#hybrid2 = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
#hybrid2_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)

#gon1 = 5
#gon2 = 6
#transition = 2
#test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
##itter_generate_interface_tiling(gon2, transition)
##itter_generate_interface_tiling(gon2, transition)
##itter_generate_interface_tiling(gon2, transition)
#test.generate_semiduals()
#name_str = 'hybrid_' + str(gon1) + '_' +str(gon2) + '_' + str(transition) 
#hybrid = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
#hybrid_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)


#gon1 = 6
#gon2 = 7
#transition = 2
#test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#test.generate_semiduals()
#name_str = 'hybrid_' + str(gon1) + '_' +str(gon2) + '_' + str(transition) 
#hybrid = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
#hybrid_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)

#gon1 = 7
#gon2 = 8
#transition = 2
#test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
##itter_generate_interface_tiling(gon2, transition)
#test.generate_semiduals()
#name_str = 'hybrid_' + str(gon1) + '_' +str(gon2) + '_' + str(transition) 
#hybrid = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
#hybrid_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)



gon1 = 7
gon2 = 6
gon3 = 5
transition = 1
transition2 = 3
test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
itter_generate_interface_tiling(gon2, transition) #0,1
itter_generate_interface_tiling(gon2, transition) #2
itter_generate_interface_tiling(gon3, transition2) #3
itter_generate_interface_tiling(gon3, transition2)
#itter_generate_interface_tiling(8, transition2)
#itter_generate_interface_tiling(6, transition2)
test.generate_semiduals()
name_str = 'hybrid_' + str(gon1) +str(gon2) + str(gon3) 
hybrid = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
hybrid_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)


#gon1 = 5
#gon2 = 6
#gon3 = 7
#transition = 1
#transition2 = 3
#test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
#itter_generate_interface_tiling(gon2, transition) #0,1
#itter_generate_interface_tiling(gon2, transition) #2
#itter_generate_interface_tiling(gon3, transition2) #3
#itter_generate_interface_tiling(gon3, transition2)
#itter_generate_interface_tiling(gon3, transition2)
##itter_generate_interface_tiling(8, transition2)
##itter_generate_interface_tiling(6, transition2)
#test.generate_semiduals()
#name_str = 'hybrid_' + str(gon1) +str(gon2) + str(gon3) 
#hybrid = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
#hybrid_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)



#gon1 = 5
#gon2 = 6
#transition = 2
#test = PlanarLayout(gon = gon1, vertex = 3, side =1, radius_method = 'lin')
#itter_generate_C60()
#itter_generate_C60()
#itter_generate_interface_tiling(gon2, transition)
#itter_generate_interface_tiling(gon2, transition)
#test.generate_semiduals()
#name_str = 'hybrid_' + 'soccerBall' + '_' +str(gon2) + '_' + str(transition) 
#hybrid = GeneralLayout(test.get_all_resonators(),  modeType = 'FW', name =  name_str)
#hybrid_HW = GeneralLayout(test.get_all_resonators(),  modeType = 'HW', name =  name_str)




fig1=pylab.figure(1)
pylab.clf()

####layout graph
ax1 = pylab.subplot(1,2,1, adjustable='box', aspect=1)
pylab.cla()
C60.draw_resonator_lattice(ax1, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#C60.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
xs = C60.coords[:,0]
ys = C60.coords[:,1]
pylab.sca(ax1)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#C60.draw_resonator_end_points(ax1, color = layoutCapColor, edgecolor = 'k',  marker = 'o' , size = smallCdefault, zorder = 5)
#ax.set_aspect('equal')
ax1.axis('off')
pylab.title('C60: Buckminsterfullerene')

ax2 = pylab.subplot(1,2,2)
pylab.cla()
pylab.plot(C60.Es, 'b.', label = 'H_s')
pylab.plot(C60_HW.Es,color =  'deepskyblue', marker = '.', linestyle = '', label = 'H_a')
pylab.xlabel('eigenvalue index')
pylab.ylabel('Energy (|t|)')
pylab.title('$\sigma$')
ax2.legend(loc = 'upper left')

pylab.tight_layout()
pylab.show()
    


fig2=pylab.figure(2)
pylab.clf()

####layout graph
ax1 = pylab.subplot(1,2,1, adjustable='box', aspect=1)
pylab.cla()
C72.draw_resonator_lattice(ax1, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#C60.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
xs = C72.coords[:,0]
ys = C72.coords[:,1]
pylab.sca(ax1)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#ax.set_aspect('equal')
ax1.axis('off')
pylab.title('C72: fullerene')

ax2 = pylab.subplot(1,2,2)
pylab.cla()
pylab.plot(C72.Es, 'b.', label = 'H_s')
pylab.plot(C72_HW.Es,color =  'deepskyblue', marker = '.', linestyle = '', label = 'H_a')
pylab.xlabel('eigenvalue index')
pylab.ylabel('Energy (|t|)')
pylab.title('$\sigma$')
ax2.legend(loc = 'upper left')

pylab.tight_layout()
pylab.show()



fig3=pylab.figure(3)
pylab.clf()

####layout graph
ax1 = pylab.subplot(1,2,1, adjustable='box', aspect=1)
pylab.cla()
C84.draw_resonator_lattice(ax1, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#C60.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
xs = C84.coords[:,0]
ys = C84.coords[:,1]
pylab.sca(ax1)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#ax.set_aspect('equal')
ax1.axis('off')
pylab.title('C84: fullerene')

ax2 = pylab.subplot(1,2,2)
pylab.cla()
pylab.plot(C84.Es, 'b.', label = 'H_s')
pylab.plot(C84_HW.Es,color =  'deepskyblue', marker = '.', linestyle = '', label = 'H_a')
pylab.xlabel('eigenvalue index')
pylab.ylabel('Energy (|t|)')
pylab.title('$\sigma$')
ax2.legend(loc = 'upper left')

pylab.tight_layout()
pylab.show()
    




fig4=pylab.figure(4)
pylab.clf()

####layout graph
ax1 = pylab.subplot(1,2,1, adjustable='box', aspect=1)
pylab.cla()
baseLayout.draw_resonator_lattice(ax1, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#extras.draw_resonator_lattice(ax1, color = 'firebrick', alpha = 1 , linewidth = 2.5)
draw_arc(ax1, newResonators[0,:], linewidth = 2.5, color = layoutLineColor, numPoints = 100, direction = -1)
draw_arc(ax1, newResonators[1,:], linewidth = 2.5, color = layoutLineColor, numPoints = 100, direction = 1)
draw_arc(ax1, newResonators[2,:], linewidth = 2.5, color = layoutLineColor, numPoints = 100, direction = -1)
xs = C84_6.coords[:,0]
ys = C84_6.coords[:,1]
pylab.sca(ax1)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#ax.set_aspect('equal')
ax1.axis('off')
pylab.title('C84: fullerene')

ax2 = pylab.subplot(1,2,2)
pylab.cla()
pylab.plot(C84_6.Es, 'b.', label = 'H_s')
pylab.plot(C84_6_HW.Es,color =  'deepskyblue', marker = '.', linestyle = '', label = 'H_a')
pylab.xlabel('eigenvalue index')
pylab.ylabel('Energy (|t|)')
pylab.title('$\sigma$')
ax2.legend(loc = 'upper left')

pylab.tight_layout()
pylab.show()
    



fig6=pylab.figure(6)
pylab.clf()

####layout graph
ax1 = pylab.subplot(1,2,1, adjustable='box', aspect=1)
pylab.cla()
hybrid.draw_resonator_lattice(ax1, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
xs = hybrid.coords[:,0]
ys = hybrid.coords[:,1]
pylab.sca(ax1)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#ax.set_aspect('equal')
ax1.axis('off')
pylab.title('hybrid')

ax2 = pylab.subplot(1,2,2)
pylab.cla()
pylab.plot(hybrid.Es, 'b.', label = 'H_s')
pylab.plot(hybrid_HW.Es,color =  'deepskyblue', marker = '.', linestyle = '', label = 'H_a')
pylab.xlabel('eigenvalue index')
pylab.ylabel('Energy (|t|)')
pylab.title('$\sigma$')
ax2.legend(loc = 'upper left')

pylab.tight_layout()
pylab.show()




fig1.set_size_inches(11.5,6.25)
fig2.set_size_inches(11.5,6.25)
fig3.set_size_inches(11.5,6.25)
fig4.set_size_inches(11.5,6.25)

#fig1.savefig('C60.png',transparent= False, dpi = 200)
#fig2.savefig('C72.png',transparent= False, dpi = 200)
#fig3.savefig('C84.png',transparent= False, dpi = 200)
#fig4.savefig('C84_Ramanujan.png',transparent= False, dpi = 200)





########
##actually make videos
########
#
#
#
##make_layout_vid(C60, startInd = 0, stopInd = 1, figNum = 1, xsize = xsize, ysize = ysize)
#make_layout_vid(C60, startInd = 0, stopInd = -1, figNum = 1, xsize = xsize, ysize = ysize)
#
##make_layout_vid(C60_HW, startInd = 0, stopInd = 1, figNum = 1, xsize = xsize, ysize = ysize)
#make_layout_vid(C60_HW, startInd = 0, stopInd = -1, figNum = 1, xsize = xsize, ysize = ysize)
#
#
#
##make_layout_vid(C84, startInd = 0, stopInd = 1, figNum = 1, xsize = xsize, ysize = ysize)
#make_layout_vid(C84, startInd = 0, stopInd = -1, figNum = 1, xsize = xsize, ysize = ysize)
#
##make_layout_vid(C84_HW, startInd = 0, stopInd = 1, figNum = 1, xsize = xsize, ysize = ysize)
#make_layout_vid(C84_HW, startInd = 0, stopInd = -1, figNum = 1, xsize = xsize, ysize = ysize)


##make_layout_vid(C84, startInd = 0, stopInd = 1, figNum = 5, xsize = xsize, ysize = ysize)
#make_layout_vid(C84_6, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)
#
##make_layout_vid(C84_6_HW, startInd = 0, stopInd = 1, figNum = 5, xsize = xsize, ysize = ysize)
#make_layout_vid(C84_6_HW, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)


#make_layout_vid(hybrid, startInd = 0, stopInd = 3, figNum = 5, xsize = xsize, ysize = ysize)
#make_layout_vid(hybrid2, startInd = 0, stopInd = 3, figNum = 5, xsize = xsize, ysize = ysize)

#make_layout_vid(hybrid, startInd = 200, stopInd = 210, figNum = 5, xsize = xsize, ysize = ysize)
#make_layout_vid(hybrid2, startInd = 200, stopInd = 210, figNum = 5, xsize = xsize, ysize = ysize)

#make_layout_vid(hybrid, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)
##make_layout_vid(hybrid2, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)
#
#make_layout_vid(hybrid_HW, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)
#make_layout_vid(hybrid2_HW, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)


#make_layout_vid(hybrid, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)
#make_layout_vid(hybrid_HW, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)
#
#make_layout_vid(hybrid2, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)
#make_layout_vid(hybrid2_HW, startInd = 0, stopInd = -1, figNum = 5, xsize = xsize, ysize = ysize)




markersize = 5
#gs = gridspec.GridSpec(1, 4,
#                           width_ratios=[1,0.5, 1, 0.5],
#                           height_ratios=[1]
#                           )
gs = gridspec.GridSpec(1, 4,
                           width_ratios=[1,0.5, 0.5, 1.0],
                           height_ratios=[1]
                           )
fig5 = pylab.figure(5)
pylab.clf()

ax = pylab.subplot(gs[0], adjustable='box', aspect=1)
pylab.cla()
C60.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
xs = C60.coords[:,0]
ys = C60.coords[:,1]
pylab.sca(ax)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#ax.set_aspect('equal')
ax.axis('off')

#ax = pylab.subplot(gs[1], adjustable='box', aspect=45*len(C60.SDx)/len(C84.SDx))
#ax = pylab.subplot(gs[1], adjustable='box', aspect=55*(len(C60.SDx)+30)/(len(C84.SDx)+45))
ax = pylab.subplot(gs[1], adjustable='box', aspect=60*(len(C60.SDx)+37)/(len(C84.SDx)+52))
pylab.cla()
xs = scipy.arange(0, len(C60.Es), 1)
pylab.plot(xs, C60.Es, 'b.', label = 'H_s',markersize = markersize)
#pylab.plot(xs+30, C60_HW.Es,color =  'deepskyblue', marker = '.', linestyle = '', label = 'H_a',markersize = markersize)
pylab.plot(xs+37, C60_HW.Es,color =  'deepskyblue', marker = '.', linestyle = '', label = 'H_a',markersize = markersize)
pylab.xlabel('Eigenvalue Index')
pylab.ylabel('Energy (|t|)')
#pylab.title('$\sigma$')
ax.legend(loc = 'upper left')


ax = pylab.subplot(gs[3], adjustable='box', aspect=1)
pylab.cla()
baseLayout.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
dd = 2.75
draw_variable_arc(ax, newResonators[0,:], d = dd, linewidth = 2.5, color = layoutLineColor, numPoints = 100, direction = -1)
draw_variable_arc(ax, newResonators[1,:], d = dd, linewidth = 2.5, color = layoutLineColor, numPoints = 100, direction = 1)
draw_variable_arc(ax, newResonators[2,:], d = dd, linewidth = 2.5, color = layoutLineColor, numPoints = 100, direction = -1)
#draw_arc(ax, newResonators[0,:], linewidth = 2.5, color = layoutLineColor, numPoints = 100, direction = -1)
#draw_arc(ax, newResonators[1,:], linewidth = 2.5, color = layoutLineColor, numPoints = 100, direction = 1)
#draw_arc(ax, newResonators[2,:], linewidth = 2.5, color = layoutLineColor, numPoints = 100, direction = -1)
xs = C84_6.coords[:,0]
ys = C84_6.coords[:,1]
pylab.sca(ax)
pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#ax.set_aspect('equal')
ax.axis('off')


#ax = pylab.subplot(gs[2], adjustable='box', aspect=45)
#ax = pylab.subplot(gs[2], adjustable='box', aspect=55)
ax = pylab.subplot(gs[2], adjustable='box', aspect=60)
#ax = pylab.subplot(gs[2])
pylab.cla()
xs = scipy.arange(0, len(C84_6.Es), 1)
pylab.plot(xs, C84_6.Es, 'b.', label = 'H_s', markersize = markersize)
#pylab.plot(xs+45, C84_6_HW.Es,color =  'deepskyblue', marker = '.', linestyle = '', label = 'H_a', markersize = markersize)
pylab.plot(xs+52, C84_6_HW.Es,color =  'deepskyblue', marker = '.', linestyle = '', label = 'H_a', markersize = markersize)
pylab.xlabel('Eigenvalue Index')
pylab.ylabel('Energy (|t|)')
#pylab.title('$\sigma$')
ax.legend(loc = 'upper left')




fig5.set_size_inches([12.2, 4.5])
pylab.tight_layout()
pylab.show()


#fig5.savefig('Fullerenes.png',transparent= False, dpi = 200)
#fig5.savefig('Fullerenes.svg',transparent= True)






##random debugging
#testLattice = PlanarLayout(vertex = 3, gon = 7)
#testLattice.populate(2, resonatorsOnly = True)
#testLattice.get_all_resonators()
#
##testLattice = EuclideanLayout(xcells = 4, ycells = 4, lattice_type = 'Huse',resonatorsOnly=True)
##        
##testLattice = GeneralLayout(testLattice.resonators,name =  'NameMe', resonatorsOnly = True)


