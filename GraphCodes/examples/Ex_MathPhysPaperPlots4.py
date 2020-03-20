#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 16:47:27 2018

@author: kollar2
Modified on May 7, 2018:
@author: Pranav Mundada

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
DataPickleFolderPath = r'/volumes/ourphoton/Alicia/Layouts/HyperbolicPickles'
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
#from GeneralLayoutGenerator import decorate_layout


#%%
###########
#Chose which figures to make
##########


DefaultFig = True
#DefaultFig = False

GraphLineGraphFig = True
#GraphLineGraphFig = False

LGBandStructureFig = True
#LGBandStructureFig = False

TsubKFig = True
#TsubKFig = False

McLaughlinFig = True
#McLaughlinFig = False

KagomeTopologyFig = True
#KagomeTopologyFig = False

HPKoskFig = True
#HPKoskFig = False

HPKoskHWFig = True
#HPKoskHWFig = False

FWHWFBFig = True
#FWHWFBFig = False

FiniteHyperbolicFig = True
#FiniteHyperbolicFig = False

FiniteHyperbolicFig2 = True
#FiniteHyperbolicFig2 = False

ExtraFBstatesFig = True
#ExtraFBstatesFig = False

KagomeUnitCellFig = True
#KagomeUnitCellFig = False

SplitGraphFig = True
#SplitGraphFig = False
if SplitGraphFig:
#    effectiveInstead = True
    effectiveInstead = False
    

SplitGraphLineGraphFig = True
#SplitGraphLineGraphFig = False

extraPlots = True
#extraPlots = False

LLSplots = True
#LLSplots = False

AliciumFig  = True
#AliciumFig = False





#%%
    
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


#####
#trial split graph defaults
#####

###old
#splitCenterColor = 'gold'
#splitEdgeColor = 'midnightblue'
#splitCenterColor2 = 'lightsteelblue'
#splitEdgeColor2 = 'midnightblue'
#LGsplitCenterColor = layoutCapColor
#splitStateColor1 = stateColor1 
#splitStateEdgeColor1 = stateEdgeColor1

####new 1
n1_splitCenterColor = 'gold'
n1_splitEdgeColor = 'midnightblue'
n1_splitCenterColor2 = 'lightsteelblue'
n1_splitEdgeColor2 = 'midnightblue'
n1_LGsplitCenterColor = layoutCapColor
n1_splitStateColor1 = 'b'
n1_splitStateEdgeColor1 = 'mediumblue'

###new 2
n2_splitCenterColor = 'dodgerblue'
n2_splitEdgeColor = 'navy'
n2_splitCenterColor2 = 'lightsteelblue'
n2_splitEdgeColor2 = 'midnightblue'
n2_LGsplitCenterColor = 'gainsboro'
n2_splitStateColor1 = stateColor1 
n2_splitStateEdgeColor1 = stateEdgeColor1







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
    kagomeLattice.draw_resonator_end_points(ax1, color = layoutCapColor, edgecolor = 'k',  marker = 'o' , size = smallCdefault, zorder = 5)
#    xs = kagomeLattice.coords[:,0]
#    ys = kagomeLattice.coords[:,1]
#    pylab.sca(ax1)
#    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    
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
    fig444 = pylab.figure(444)
    pylab.clf()
    pylab.close(fig444)













#%%


########
#figure 1: resoantors, and effective lattices
#######
if GraphLineGraphFig:
    bigC = 90
    
    smallGap = 0.24
    bigGap = 0.52+ 0.24
    
    plotStates = True
#    plotStates = False
    
    gs = gridspec.GridSpec(4, 4,
                           width_ratios=[1, 1, 1, 1],
                           height_ratios=[1, 1, 0.1, 0.1]
                           )
    fig1 = pylab.figure(1)
    pylab.clf()
    #1) kagome
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    kagomeLattice = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    ax = pylab.subplot(gs[0], adjustable='box', aspect=1)
    pylab.cla()
    kagomeLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    xs = kagomeLattice.coords[:,0]
    ys = kagomeLattice.coords[:,1]
    pylab.sca(ax)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    #ax.set_aspect('equal')
    ax.axis('off')
    #pylab.title('layout graph')
    
    
    ax5 = pylab.subplot(gs[4+0], adjustable='box', aspect=1)
    kagomeLattice.draw_SDlinks(ax5, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax5)
    pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
    red = [3,5,13]
    yellow = [4,8,9]
    shiftx, shifty = testEuclidLattice.unitcell.a1
    if plotStates:
        pylab.scatter(kagomeLattice.SDx[yellow]+shiftx, kagomeLattice.SDy[yellow]+shifty, c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
        pylab.scatter(kagomeLattice.SDx[red]+shiftx, kagomeLattice.SDy[red]+shifty, c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    #ax5.set_aspect('equal')
    ax5.axis('off')



    
    #2) square
    testEuclidLattice = EuclideanLayout(4,4,lattice_type = 'square', modeType = 'FW')
    #testEuclidLattice = EuclideanLayout(4,4,lattice_type = 'square', modeType = 'HW')
    res_list = testEuclidLattice.get_all_resonators()
    mask = numpy.logical_and((res_list[:,0]+res_list[:,2])>0 ,(res_list[:,1]+res_list[:,3])>0)
    squareLattice = GeneralLayout(res_list[mask,:] , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    ax2 = pylab.subplot(gs[1])
    pylab.cla()
    squareLattice.draw_resonator_lattice(ax2, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    xs = squareLattice.coords[:,0]
    ys = squareLattice.coords[:,1]
    pylab.sca(ax2)
    pylab.scatter(xs, ys ,c = layoutCapColor, s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
#    testEuclidLattice.draw_resonator_end_points(ax2, color = layoutCapColor, edgecolor = 'k',  marker = 'o' , size = 20, zorder = 5)
#    squareLattice.draw_resonator_end_points(ax2, color = layoutCapColor, edgecolor = 'k',  marker = 'o' , size = 20, zorder = 5)
    ax2.set_aspect('equal')
    ax2.axis('off')
    #pylab.title('layout graph')
    
    
    ax6 = pylab.subplot(gs[4+1])
    squareLattice.draw_SDlinks(ax6, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax6)
    pylab.scatter(squareLattice.SDx, squareLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
#    red = [3,5,13]
#    yellow = [4,8,9]
    red = [13]
    yellow = [9]
    shiftx = testEuclidLattice.unitcell.a1[0]
    shifty = testEuclidLattice.unitcell.a2[1]
    if plotStates:
        pylab.scatter(squareLattice.SDx[yellow], squareLattice.SDy[yellow]-shifty, c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
        pylab.scatter(squareLattice.SDx[red]-shiftx, squareLattice.SDy[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
        pylab.scatter(squareLattice.SDx[yellow], squareLattice.SDy[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
        pylab.scatter(squareLattice.SDx[red], squareLattice.SDy[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    ax6.set_aspect('equal')
    ax6.axis('off')
 
    
    
    
    #3) tree
    testTree = TreeResonators(degree = 3, iterations = 5, side = 1, file_path = '', modeType = 'FW')
    #testTree = TreeResonators(degree = 3, iterations = 5, side = 1, file_path = '', modeType = 'HW')
    resonators = testTree.get_all_resonators()
    treeLattice = GeneralLayout(resonators , modeType = testTree.modeType, name =  'TREEEEE')
    ax3 = pylab.subplot(gs[2])
    pylab.cla()
    treeLattice.draw_resonator_lattice(ax3, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    xs = treeLattice.coords[:,0]
    ys = treeLattice.coords[:,1]
    pylab.sca(ax3)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    ax3.set_aspect('equal')
    ax3.axis('off')
    #pylab.title('layout graph')
    
    
    ax7 = pylab.subplot(gs[4+2])
    treeLattice.draw_SDlinks(ax7, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax7)
    pylab.scatter(treeLattice.SDx, treeLattice.SDy ,c =  FWsiteColor, s = 20, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
    state_vect = treeLattice.Psis[:,0]
    stepSize = 0.11/2
    tempVect = state_vect/stepSize
    maxVal = numpy.max(state_vect)
    for ind in range(0, len(state_vect)):
        rawVal = state_vect[ind]
        if numpy.abs(rawVal/stepSize)>8:
            #this is the peak
            pass
        elif numpy.abs(rawVal/stepSize)>4:
            #second ring
            state_vect[ind] = maxVal*numpy.sign(rawVal)/2
        elif numpy.abs(rawVal/stepSize)>2:
            #third ring
            state_vect[ind] = maxVal*numpy.sign(rawVal)/4
        elif numpy.abs(rawVal/stepSize)>1:
            #fourth ring
            state_vect[ind] = maxVal*numpy.sign(rawVal)/8
        else:
            state_vect[ind] = 0
            
    reds = numpy.where(state_vect <0)[0]
    yellows = numpy.where(state_vect > 0)[0]
    
    scaleFactor = 8
    
    redState = numpy.zeros(len(state_vect))
    redState[reds] = state_vect[reds]
    redAmps = redState
    redProbs = numpy.abs(redAmps)**2
    redSizes = redProbs * len(redProbs)*scaleFactor
    
    yellowState = numpy.zeros(len(state_vect))
    yellowState[yellows] = state_vect[yellows]
    yellowAmps = yellowState
    yellowProbs = numpy.abs(yellowAmps)**2
    yellowSizes = yellowProbs * len(yellowProbs)*scaleFactor
    
    if plotStates:
        pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  stateColor2, s = redSizes, marker = 'o', edgecolors = stateEdgeColor2,  zorder = 6, linewidth = stateEdgeWidth2)
        pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  stateColor1, s = yellowSizes, marker = 'o', edgecolors = stateEdgeColor1,  zorder = 6, linewidth = stateEdgeWidth1)
    
    ax7.set_aspect('equal')
    ax7.axis('off')




    
    #4) hyperbolic
    testHyperbolic = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
    testHyperbolic.populate(2, resonatorsOnly=False)
    resonators = testHyperbolic.get_all_resonators()
    nameStr = str(testHyperbolic.gon) + 'gon_' + str(testHyperbolic.vertex) + 'vertex_' + str(testHyperbolic.itter) + '_' + testHyperbolic.modeType
    hyperbolicLattice = GeneralLayout(resonators , modeType = testHyperbolic.modeType, name =  nameStr)
    
    
    
    ax4 = pylab.subplot(gs[3])
    pylab.cla()
    hyperbolicLattice.draw_resonator_lattice(ax4, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    xs = hyperbolicLattice.coords[:,0]
    ys = hyperbolicLattice.coords[:,1]
    pylab.sca(ax4)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
#    testHyperbolic.draw_resonator_end_points(ax4, color = layoutCapColor, edgecolor = 'k',  marker = 'o' , size = 20, zorder = 5)
    ax4.set_aspect('equal')
    ax4.axis('off')
    #pylab.title('layout graph')
    
    
    ax8 = pylab.subplot(gs[4+3])
    hyperbolicLattice.draw_SDlinks(ax8, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax8)
    pylab.scatter(hyperbolicLattice.SDx, hyperbolicLattice.SDy ,c =  FWsiteColor, s = 20, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
    red = [2,4,6,7,9,36]
    yellow = [1,3,5,8,10,35]
    if plotStates:
        pylab.scatter(hyperbolicLattice.SDx[yellow], hyperbolicLattice.SDy[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
        pylab.scatter(hyperbolicLattice.SDx[red], hyperbolicLattice.SDy[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
#    ax1.set_aspect('equal')
#    ax1.axis('off')
    
    ax8.set_aspect('equal')
    ax8.axis('off')
    #pylab.title('effective (line) graph')
    
    #pylab.suptitle(hyperbolicLattice.name)
    #fig1.set_size_inches(15, 8)
    #pylab.tight_layout()
    #pylab.show()
    #fig1.savefig('Figure1.pdf',transparent= True)
    
    
    ####   sample DOS plot
    #set up frequency sweep
    ## 1) Kagome 
    freq_start = -4.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = []#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-3.0,3.0)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    Amp_imp = 0.3
    Amp_band = 0.1
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[8+0], adjustable='box', aspect=6)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)    
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-3,3], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,6.5])
    #ax1.set_xlim([-4,4])
    #ax1.legend(loc= 'upper right')
    
    ## 2) Square 
    freq_start = -4.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = []#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-4.0,4.0)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[8+1], adjustable='box', aspect=6)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-4,4], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,6.5])
    
    ## 3) Tree 
    freq_start = -4.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = []#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-3 + smallGap,3 - smallGap)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[8+2], adjustable='box', aspect=6)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-3,3], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,6.5])
    
    ## 4) Hyperbolics 
    freq_start = -4.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = []#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-3 + bigGap ,3 - smallGap)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[8+3], adjustable='box', aspect=6)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-3,3], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,6.5])
    
    ## 5) Kagome Line 
    freq_start = -2.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [-2]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2,4)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[12+0], adjustable='box', aspect=6)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,4], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,6.5])
    
    ## 6) Square Line 
    freq_start = -2.00
    freq_stop = 6.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [-2]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2,6)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[12+1], adjustable='box', aspect=6)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,6], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,6.5])
    
    ## 7) Tree Line 
    freq_start = -2.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [-2]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2 + smallGap, 4 - smallGap)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[12+2], adjustable='box', aspect=6)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,4], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,6.5])
    
    ## 8) Hyperbolic Line 
    freq_start = -2.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [-2]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2 + bigGap,4 - smallGap)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[12+3], adjustable='box', aspect=6)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1)#, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6)#, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 0.6, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,4], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,6.5])
    
    fig1.set_size_inches(11.5, 7.2)
    pylab.tight_layout()
    pylab.show()
#    fig1.savefig('Figure1.svg',transparent= True)
#    fig1.savefig('GraphsAndLineGraphs.png',transparent= False, dpi = 200)
#    fig1.savefig('GraphsAndLineGraphs.png',transparent= False, dpi = 400)
#    fig1.savefig('GraphsAndLineGraphs.svg',transparent= True)
#    fig1.savefig('GraphsAndLineGraphs_noStates.png',transparent= False, dpi = 200)
else:
    fig1 = pylab.figure(1)
    pylab.clf()
    pylab.close(fig1)
    
    
    
    
    
    
    
    
    
    
    
    
#%%

#######
#band structure figure
#########    
if LGBandStructureFig:
    fig2 = pylab.figure(2)
    pylab.clf()
else:
    fig2 = pylab.figure(2)
    pylab.clf()
    pylab.suptitle('band structure fig')
    pylab.close(fig2)
    
    




















#%% 
    
############
#Figure 5: Topology argument in kagome
########################################
if KagomeTopologyFig:
##    bigC = 110
##    smallC = 30
#    fig5 = pylab.figure(5)
#    pylab.clf()
#    #1) kagome
#    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
#    #testEuclidLattice = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'HW')
#    kagomeLattice = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
#    pylab.cla()
#    red = [3,5,13,18,20,21,24,32,33,35]
#    yellow = [4,8,9,19,22,25,31,36,34,37]
#    
#    ax5 = pylab.subplot(1,2,1, adjustable='box', aspect=1)
#    kagomeLattice.draw_SDlinks(ax5, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=FWlinkAlpha)
#    pylab.sca(ax5)
#    pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=FWsiteAlpha)
#    pylab.scatter(kagomeLattice.SDx[yellow], kagomeLattice.SDy[yellow], c =  stateColor1, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
#    pylab.scatter(kagomeLattice.SDx[red], kagomeLattice.SDy[red], c =  stateColor2, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
#    #ax5.set_aspect('equal')
#    ax5.axis('off')
#    
#    red2 = [9,12,15,0,18,27,34]
#    yellow2 = [10,13,16,5,14,23]
#    
#    ax5 = pylab.subplot(1,2,2, adjustable='box', aspect=1)
#    kagomeLattice.draw_SDlinks(ax5, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=FWlinkAlpha)
#    pylab.sca(ax5)
#    pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=FWsiteAlpha)
#    pylab.scatter(kagomeLattice.SDx[yellow2], kagomeLattice.SDy[yellow2], c =  stateColor1, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
#    pylab.scatter(kagomeLattice.SDx[red2], kagomeLattice.SDy[red2], c =  stateColor2, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
#    #ax5.set_aspect('equal')
#    ax5.axis('off')
#    fig1.set_size_inches(13, 5)
#    #pylab.tight_layout()
#    pylab.show()
##    fig5.savefig('Figure5.pdf',transparent= True)
##    fig5.savefig('kagomeTopology.png',transparent= False, dpi = 200)
##    fig5.savefig('kagomeTopology.png',transparent= False, dpi = 400)
##    fig5.savefig('kagomeTopology.svg',transparent= True)
    

    fig5 = pylab.figure(5)
    pylab.clf()
    #1) kagome
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    #testEuclidLattice = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'HW')
    kagomeLattice = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    pylab.cla()
    red = [3,5,13,18,20,21,24,32,33,35]
    yellow = [4,8,9,19,22,25,31,36,34,37]
    
    ax5 = pylab.subplot(1,3,1, adjustable='box', aspect=1)
    kagomeLattice.draw_SDlinks(ax5, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=FWlinkAlpha)
    pylab.sca(ax5)
    pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=FWsiteAlpha)
    pylab.scatter(kagomeLattice.SDx[yellow], kagomeLattice.SDy[yellow], c =  stateColor1, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(kagomeLattice.SDx[red], kagomeLattice.SDy[red], c =  stateColor2, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    #ax5.set_aspect('equal')
    ax5.axis('off')
    
    red2 = [9,12,15,0,18,27,34]
    yellow2 = [10,13,16,5,14,23]
    
    ax5 = pylab.subplot(1,3,2, adjustable='box', aspect=1)
    kagomeLattice.draw_SDlinks(ax5, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=FWlinkAlpha)
    pylab.sca(ax5)
    pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=FWsiteAlpha)
    pylab.scatter(kagomeLattice.SDx[yellow2], kagomeLattice.SDy[yellow2], c =  stateColor1, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(kagomeLattice.SDx[red2], kagomeLattice.SDy[red2], c =  stateColor2, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    #ax5.set_aspect('equal')
    ax5.axis('off')



    ax = pylab.subplot(1,3,3, adjustable='box', aspect=1)
    testEuclidLattice = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'FW')
    kagomeLattice2 = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    keepers = [12,8, 18, 19, 20, 17,16,15, 11, 23, 22, 25, 24, 27, 28, 29, 30, 31, 32, 40, 39, 41, 36, 37]
    truncationRes = kagomeLattice2.resonators[keepers, :]
    truncation = GeneralLayout(truncationRes , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    
    truncation.draw_SDlinks(ax, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=FWlinkAlpha)
    pylab.sca(ax)
    pylab.scatter(truncation.SDx, truncation.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=FWsiteAlpha)
    ax.set_aspect('equal')
    ax.axis('off')

    fig5.set_size_inches(9.6, 3.6)
    pylab.tight_layout()
    pylab.show()
    
    
#    fig5.savefig('Figure5.pdf',transparent= True)
#    fig5.savefig('kagomeTopology.png',transparent= False, dpi = 200)
#    fig5.savefig('kagomeTopology.png',transparent= False, dpi = 400)
#    fig5.savefig('kagomeTopology.svg',transparent= True)
    
    
else:
    fig5 = pylab.figure(5)
    pylab.clf()
    pylab.suptitle('KagomeTopologyFig')
    pylab.close(fig5)










#%%
########
#figure 6: HPK, OSK
#######
if HPKoskFig:
    cellModeType = 'FW'
    #cellModeType = 'HW'
    
    bigC = 70
    smallC = 25
    #1) HPK
    HPKCell = UnitCell('Huse')
    testEuclidLattice = EuclideanLayout(4,2,lattice_type = 'Huse', modeType = cellModeType)
    resonators = testEuclidLattice.resonators
    HPKLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'Huse')
    
    
    gs = gridspec.GridSpec(2, 5,
                           width_ratios=[1, 1, 1, 1,1],
                           height_ratios=[1, 1]
                           )
    
    red2 = [13,21,27,29,31,35,38,40,49,56,65,70,73,89,86,92,94,97,98,111]
    yellow2 = [12,19,23,26,28,30,33,37,39,48,53,55,64,71,72,73,74,87,91,95,110]
    
    latticeStretch = 0.7
    
    fig6 = pylab.figure(6)
    pylab.clf()
#    ax = pylab.subplot(2,5,1)
    ax = pylab.subplot(gs[0], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    HPKLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    #add transparent bubbles so fiugres are sized correctly.
    pylab.scatter(HPKLattice.SDx[yellow2], HPKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=0)
    pylab.scatter(HPKLattice.SDx[red2], HPKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=0)
#    ax.set_aspect('equal')
    ax.axis('off')
#    pylab.title('effective (line) graph')
    
    
#    ax2 = pylab.subplot(2,5,2)
    ax2 = pylab.subplot(gs[1], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    HPKLattice.draw_SDlinks(ax2, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
    pylab.sca(ax2)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
    pylab.scatter(HPKLattice.SDx[yellow2], HPKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(HPKLattice.SDx[red2], HPKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
#    ax2.set_aspect('equal')
    ax2.axis('off')
#    pylab.title('Contractible')
    
#    ax3 = pylab.subplot(2,5,3)
    ax3 = pylab.subplot(gs[2], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    mask = numpy.logical_and(HPKLattice.SDx>0.433 ,HPKLattice.SDx<0.453)
    x1 = HPKLattice.SDx[mask]
    y1 = HPKLattice.SDy[mask]
    red3 = [0,24,48,72,113]
    yellow3 = [18,42,66,90]
    HPKLattice.draw_SDlinks(ax3, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
    pylab.sca(ax3)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
    pylab.scatter(x1[::2], y1[::2] ,c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
    pylab.scatter(x1[1::2], y1[1::2] ,c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(HPKLattice.SDx[yellow3], HPKLattice.SDy[yellow3],c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(HPKLattice.SDx[red3], HPKLattice.SDy[red3], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
    
#    ax3.set_aspect('equal')
    ax3.axis('off')
#    pylab.title('Incontractible')
    
    ##1) osK
    OSKCell = UnitCell('84Huse')
    testEuclidLattice = EuclideanLayout(4,2,lattice_type = '84Huse', modeType = cellModeType)
    resonators = testEuclidLattice.resonators
    OSKLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  '84Huse')
    
    red2 = [13,20,29,34,35,38,57,63,67,76,98,97]#,40,49,56,65,70,73,89,86,92,94,97,98,111]
    yellow2 = [12,19,21,28,30,32,37,53,59,62,70,75]#,39,48,53,55,64,71,72,73,74,87,91,95,110]
    
#    ax = pylab.subplot(2,5,6)
    ax = pylab.subplot(gs[5], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    OSKLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = 1)
    pylab.sca(ax)
    pylab.scatter(OSKLattice.SDx, OSKLattice.SDy , c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = 1)
    # add in transparent circles so that the figures are autosized correctly
    pylab.scatter(HPKLattice.SDx[yellow2], HPKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=0, linewidth = stateEdgeWidth1)
    pylab.scatter(HPKLattice.SDx[red2], HPKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=0, linewidth = stateEdgeWidth2)
#    ax.set_aspect('equal')
    ax.axis('off')
    
#    ax2 = pylab.subplot(2,5,7)
    ax2 = pylab.subplot(gs[6], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    OSKLattice.draw_SDlinks(ax2, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
    pylab.sca(ax2)
    pylab.scatter(OSKLattice.SDx, OSKLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
    pylab.scatter(OSKLattice.SDx[yellow2], OSKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(OSKLattice.SDx[red2], OSKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
#    ax2.set_aspect('equal')
    ax2.axis('off')
    
#    ax3 = pylab.subplot(2,5,8)
    ax3 = pylab.subplot(gs[7], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    mask = numpy.logical_and(OSKLattice.SDx>0.433 ,OSKLattice.SDx<0.453)
    x1 = OSKLattice.SDx[mask]
    y1 = OSKLattice.SDy[mask]
    red3 = [0,24,48,72,113]
    yellow3 = [18,42,66,90]
    OSKLattice.draw_SDlinks(ax3, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
    pylab.sca(ax3)
    pylab.scatter(OSKLattice.SDx, OSKLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
    pylab.scatter(x1[::2], y1[::2] ,c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
    pylab.scatter(x1[1::2], y1[1::2] ,c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(OSKLattice.SDx[yellow3], OSKLattice.SDy[yellow3], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(OSKLattice.SDx[red3], OSKLattice.SDy[red3], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
    
#    ax3.set_aspect('equal')
    ax3.axis('off')
    
    #pylab.suptitle(testLattice.name)
    
    #band structure
    numSurfPoints = 300
    kx_x, ky_y, cutx = HPKCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = HPKCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
    
    
    
#    ax4 = pylab.subplot(2,5,4, adjustable='box', aspect=80)
    ax4 = pylab.subplot(gs[3], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cutx)
    HPKCell.plot_band_cut(ax4, cutx)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_x$ ($\pi$/a)')
    pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2,0,2], rotation='horizontal')
    
#    ax5 = pylab.subplot(2,5,5, adjustable='box', aspect=80)
    ax5 = pylab.subplot(gs[4], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
    HPKCell.plot_band_cut(ax5, cuty)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
    
    numSurfPoints = 300
    kx_x, ky_y, cutx = OSKCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = OSKCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
    
    
    
#    ax9 = pylab.subplot(2,5,9, adjustable='box', aspect=80)
    ax9 = pylab.subplot(gs[8], adjustable='box', aspect=80)
    #plot_band_cut_mc(ax, cutx)
    pylab.cla()
    HPKCell.plot_band_cut(ax9, cutx)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_x$ ($\pi$/a)')
    pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2,0,2], rotation='horizontal')
    
#    ax10 = pylab.subplot(2,5,10, adjustable='box', aspect=80)
    ax10 = pylab.subplot(gs[9], adjustable='box', aspect=80)
    #plot_band_cut_mc(ax, cuty)
    pylab.cla()
    HPKCell.plot_band_cut(ax10, cuty)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
    #titleStr = testCell.type + ', modeType: ' + cellModeType + ' (Made with UnitCell class)' 
    #pylab.suptitle(titleStr)
    
    fig6.set_size_inches(12.5, 7.5)
    pylab.tight_layout()
    pylab.show()
#    fig6.savefig('Figure6.pdf',transparent= True)
#    fig6.savefig('HPK_OSK.png',transparent= False, dpi = 200)
#    fig6.savefig('HPK_OSK.png',transparent= False, dpi = 400)
#    fig6.savefig('HPK_OSK.svg',transparent= True)
    ###fig2.savefig('HL_BS_1.png', dpi = 200)
    
    
    
    
    #    def plot_band_cut(self, ax, bandCut, colorlist = ''):
    #        if colorlist == '':
    #            colorlist = ['firebrick', 'dodgerblue', 'blueviolet', 'mediumblue', 'goldenrod', 'cornflowerblue']
    #        
    #        for ind in range(0,self.numSites):
    #            colorInd = numpy.mod(ind, len(colorlist))
    #            pylab.plot(bandCut[ind,:], color = colorlist[colorInd] , marker = '.', markersize = '5', linestyle = '')
    ##            pylab.plot(bandCut[ind,:], '.')
    #        pylab.title('some momentum cut')
    #        pylab.ylabel('Energy')
    #        pylab.xlabel('k_something')
    
    
    
    #pylab.rcParams.update({'font.size': 14})
else:
    fig6 = pylab.figure(6)    
    pylab.clf()
    pylab.suptitle('HPK, OSK Fig')
    pylab.close(fig6)
    
    
    
    
    
    
    
    
#%%
########
#figure 6: HPK, OSK, alternate version with HW HPK
#######
if HPKoskHWFig:
    cellModeType = 'FW'
    #cellModeType = 'HW'
    
    bigC = 70
    smallC = 25
    
    allThree = True
#    allThree = False
    
    
    ##########1) HPK
    HPKCell = UnitCell('Huse')
    testEuclidLattice = EuclideanLayout(4,2,lattice_type = 'Huse', modeType = cellModeType)
    resonators = testEuclidLattice.resonators
    HPKLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'Huse')
    
    
    if allThree:
        gs = gridspec.GridSpec(3, 5,
                               width_ratios=[1, 1, 1, 1,1],
                               height_ratios=[1, 1,1]
                               )
    else:
        gs = gridspec.GridSpec(2, 5,
                               width_ratios=[1, 1, 1, 1,1],
                               height_ratios=[1, 1]
                               )
    
    red2 = [13,21,27,29,31,35,38,40,49,56,65,70,73,89,86,92,94,97,98,111]
    yellow2 = [12,19,23,26,28,30,33,37,39,48,53,55,64,71,72,73,74,87,91,95,110]
    
    latticeStretch = 0.7
    
    fig6 = pylab.figure(6)
    pylab.clf()
#    ax = pylab.subplot(2,5,1)
    ax = pylab.subplot(gs[0], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    HPKLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    #add transparent bubbles so fiugres are sized correctly.
    pylab.scatter(HPKLattice.SDx[yellow2], HPKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=0)
    pylab.scatter(HPKLattice.SDx[red2], HPKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=0)
#    ax.set_aspect('equal')
    ax.axis('off')
#    pylab.title('effective (line) graph')
    
    ####hack to make matching layout for a talk
#    pylab.cla()
#    HPKLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#    xs = HPKLattice.coords[:,0]
#    ys = HPKLattice.coords[:,1]
#    pylab.sca(ax)
#    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#    #ax.set_aspect('equal')
#    ax.axis('off')
    
    
    
    
#    ax2 = pylab.subplot(2,5,2)
    ax2 = pylab.subplot(gs[1], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    HPKLattice.draw_SDlinks(ax2, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
    pylab.sca(ax2)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
    pylab.scatter(HPKLattice.SDx[yellow2], HPKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(HPKLattice.SDx[red2], HPKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
#    ax2.set_aspect('equal')
    ax2.axis('off')
#    pylab.title('Contractible')
    
#    ax3 = pylab.subplot(2,5,3)
    ax3 = pylab.subplot(gs[2], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    mask = numpy.logical_and(HPKLattice.SDx>0.433 ,HPKLattice.SDx<0.453)
    x1 = HPKLattice.SDx[mask]
    y1 = HPKLattice.SDy[mask]
    red3 = [0,24,48,72,113]
    yellow3 = [18,42,66,90]
    HPKLattice.draw_SDlinks(ax3, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
    pylab.sca(ax3)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
    pylab.scatter(x1[::2], y1[::2] ,c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
    pylab.scatter(x1[1::2], y1[1::2] ,c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(HPKLattice.SDx[yellow3], HPKLattice.SDy[yellow3],c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(HPKLattice.SDx[red3], HPKLattice.SDy[red3], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
    
#    ax3.set_aspect('equal')
    ax3.axis('off')
#    pylab.title('Incontractible')
    
    
    #band structure
    numSurfPoints = 300
    kx_x, ky_y, cutx = HPKCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = HPKCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
    
    
    #    ax4 = pylab.subplot(2,5,4, adjustable='box', aspect=80)
    ax4 = pylab.subplot(gs[3], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cutx)
    HPKCell.plot_band_cut(ax4, cutx)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_x$ ($\pi$/a)')
    pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2,0,2], rotation='horizontal')
    
#    ax5 = pylab.subplot(2,5,5, adjustable='box', aspect=80)
    ax5 = pylab.subplot(gs[4], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
    HPKCell.plot_band_cut(ax5, cuty)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
    
    
    
    
    
    ##########1) HW HPK
    HPKCell = UnitCell('Huse')
    testEuclidLattice = EuclideanLayout(4,2,lattice_type = 'Huse', modeType = 'HW')
    resonators = testEuclidLattice.resonators
    HPKLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'Huse')



#    #heptagon1
#    red2 = [48,42,65]
#    yellow2 = [40,41,47,46]
#    #pentagon1
#    red2 = [27,28,11]
#    yellow2 = [8,9]
#    #pentagon2
#    red2 = [55,51,50]
#    yellow2 = [57,58]
#    #alternate pentagon
#    red2 = [75,74,79]
#    yellow2 = [81,82]
#    #hetagon2
#    red2 = [79,80,104,103]
#    yellow2 = [73,72,90]
    
    red2 = [48,42,65] + [27,28,11]+ [55,51,50] + [79,80,104,103]
    yellow2 = [40,41,47,46] + [8,9] + [57,58] + [73,72,90]
    
    
#    ax = pylab.subplot(2,5,1)
    ax = pylab.subplot(gs[5], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    HPKLattice.draw_SDlinks(ax, color = HWlinkColor, linewidth = 2.5, minus_links = True, minus_color = HWminusLinkColor)
    pylab.sca(ax)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c =  HWsiteColor, s = smallCdefault, marker = 'o', edgecolors = HWsiteEdgeColor, zorder = 5)
    #add transparent bubbles so fiugres are sized correctly.
    pylab.scatter(HPKLattice.SDx[yellow2], HPKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=0)
    pylab.scatter(HPKLattice.SDx[red2], HPKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=0)
#    ax.set_aspect('equal')
    ax.axis('off')
#    pylab.title('effective (line) graph')
    

    
    
#    ax2 = pylab.subplot(2,5,2)
    ax2 = pylab.subplot(gs[6], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    HPKLattice.draw_SDlinks(ax2, color = HWlinkColor, linewidth = 2.5, minus_links = True, minus_color = HWminusLinkColor, alpha = HWlinkAlpha)
    pylab.sca(ax2)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c =  HWsiteColor, s = smallCdefault, marker = 'o', edgecolors = HWsiteEdgeColor, zorder = 5, alpha = HWsiteAlpha)
    pylab.scatter(HPKLattice.SDx[yellow2], HPKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
    pylab.scatter(HPKLattice.SDx[red2], HPKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
#    ax2.set_aspect('equal')
    ax2.axis('off')
#    pylab.title('Contractible')
    

    
    
#    ax3 = pylab.subplot(2,5,3)
    ax3 = pylab.subplot(gs[7], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    mask = numpy.logical_and(HPKLattice.SDx>0.433 ,HPKLattice.SDx<0.453)
    x1 = HPKLattice.SDx[mask]
    y1 = HPKLattice.SDy[mask]
    red3 = [0,24,48,72,113, 18,42,66,90]
    HPKLattice.draw_SDlinks(ax3, color = HWlinkColor, linewidth = 2.5, minus_links = True, minus_color = HWminusLinkColor, alpha = HWlinkAlpha)
    pylab.sca(ax3)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c =  HWsiteColor, s = smallCdefault, marker = 'o', edgecolors = HWsiteEdgeColor, zorder = 5, alpha = HWsiteAlpha)
    pylab.scatter(x1, y1 ,c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
    pylab.scatter(HPKLattice.SDx[red3], HPKLattice.SDy[red3], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
    
#    ax3.set_aspect('equal')
    ax3.axis('off')
#    pylab.title('Incontractible')
    
    
    #band structure
    numSurfPoints = 300
    kx_x, ky_y, cutx = HPKCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = 'HW')
    kx_y, ky_y, cuty = HPKCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = 'HW')
    
    
    def plot_band_cut_hack(ax, bandCut, colorlist = ''):
        if colorlist == '':
            colorlist = ['firebrick', 'dodgerblue', 'blueviolet', 'mediumblue',     'cornflowerblue', 'goldenrod', 'mediumblue', 'blueviolet', 'dodgerblue', 'firebrick', 'cornflowerblue', 'goldenrod']
        
        for ind in range(0,HPKCell.numSites):
            colorInd = numpy.mod(ind, len(colorlist))
            pylab.plot(bandCut[ind,:], color = colorlist[colorInd] , marker = '.', markersize = '5', linestyle = '')
    #            pylab.plot(bandCut[ind,:], '.')
        pylab.title('some momentum cut')
        pylab.ylabel('Energy')
        pylab.xlabel('k_something')
    
    #    ax4 = pylab.subplot(2,5,4, adjustable='box', aspect=80)
    ax4 = pylab.subplot(gs[8], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cutx)
#    HPKCell.plot_band_cut(ax4, cutx)
    plot_band_cut_hack(ax4, cutx, colorlist = '')
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_x$ ($\pi$/a)')
    pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2,0,2], rotation='horizontal')
    
#    ax5 = pylab.subplot(2,5,5, adjustable='box', aspect=80)
    ax5 = pylab.subplot(gs[9], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
#    HPKCell.plot_band_cut(ax5, cuty)
    plot_band_cut_hack(ax5, cuty, colorlist = '')
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
    
    
    #    #    def plot_band_cut(self, ax, bandCut, colorlist = ''):
#    #        if colorlist == '':
#    #            colorlist = ['firebrick', 'dodgerblue', 'blueviolet', 'mediumblue', 'goldenrod', 'cornflowerblue']
#    #        
#    #        for ind in range(0,self.numSites):
#    #            colorInd = numpy.mod(ind, len(colorlist))
#    #            pylab.plot(bandCut[ind,:], color = colorlist[colorInd] , marker = '.', markersize = '5', linestyle = '')
#    ##            pylab.plot(bandCut[ind,:], '.')
#    #        pylab.title('some momentum cut')
#    #        pylab.ylabel('Energy')
#    #        pylab.xlabel('k_something')
    
    
    
    
    if allThree:
        ###########3) osK
        OSKCell = UnitCell('84Huse')
        testEuclidLattice = EuclideanLayout(4,2,lattice_type = '84Huse', modeType = cellModeType)
        resonators = testEuclidLattice.resonators
        OSKLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  '84Huse')
        
        red2 = [13,20,29,34,35,38,57,63,67,76,98,97]#,40,49,56,65,70,73,89,86,92,94,97,98,111]
        yellow2 = [12,19,21,28,30,32,37,53,59,62,70,75]#,39,48,53,55,64,71,72,73,74,87,91,95,110]
        
    #    ax = pylab.subplot(2,5,6)
        ax = pylab.subplot(gs[10], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        OSKLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = 1)
        pylab.sca(ax)
        pylab.scatter(OSKLattice.SDx, OSKLattice.SDy , c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = 1)
        # add in transparent circles so that the figures are autosized correctly
        pylab.scatter(HPKLattice.SDx[yellow2], HPKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=0, linewidth = stateEdgeWidth1)
        pylab.scatter(HPKLattice.SDx[red2], HPKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=0, linewidth = stateEdgeWidth2)
    #    ax.set_aspect('equal')
        ax.axis('off')
        
    #    ax2 = pylab.subplot(2,5,7)
        ax2 = pylab.subplot(gs[11], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        OSKLattice.draw_SDlinks(ax2, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
        pylab.sca(ax2)
        pylab.scatter(OSKLattice.SDx, OSKLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
        pylab.scatter(OSKLattice.SDx[yellow2], OSKLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
        pylab.scatter(OSKLattice.SDx[red2], OSKLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
    #    ax2.set_aspect('equal')
        ax2.axis('off')
        
    #    ax3 = pylab.subplot(2,5,8)
        ax3 = pylab.subplot(gs[12], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        mask = numpy.logical_and(OSKLattice.SDx>0.433 ,OSKLattice.SDx<0.453)
        x1 = OSKLattice.SDx[mask]
        y1 = OSKLattice.SDy[mask]
        red3 = [0,24,48,72,113]
        yellow3 = [18,42,66,90]
        OSKLattice.draw_SDlinks(ax3, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
        pylab.sca(ax3)
        pylab.scatter(OSKLattice.SDx, OSKLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
        pylab.scatter(x1[::2], y1[::2] ,c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
        pylab.scatter(x1[1::2], y1[1::2] ,c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
        pylab.scatter(OSKLattice.SDx[yellow3], OSKLattice.SDy[yellow3], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=1, linewidth = stateEdgeWidth1)
        pylab.scatter(OSKLattice.SDx[red3], OSKLattice.SDy[red3], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=1, linewidth = stateEdgeWidth2)
        
    #    ax3.set_aspect('equal')
        ax3.axis('off')
        
        #pylab.suptitle(testLattice.name)
        
        
        
        numSurfPoints = 300
        kx_x, ky_y, cutx = OSKCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
        kx_y, ky_y, cuty = OSKCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
        
        
    #    ax9 = pylab.subplot(2,5,9, adjustable='box', aspect=80)
        ax9 = pylab.subplot(gs[13], adjustable='box', aspect=80)
        #plot_band_cut_mc(ax, cutx)
        pylab.cla()
        HPKCell.plot_band_cut(ax9, cutx)
        pylab.title('')
        pylab.ylabel('Energy (|t|)')
        pylab.xlabel('$k_x$ ($\pi$/a)')
        pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2,0,2], rotation='horizontal')
        
    #    ax10 = pylab.subplot(2,5,10, adjustable='box', aspect=80)
        ax10 = pylab.subplot(gs[14], adjustable='box', aspect=80)
        #plot_band_cut_mc(ax, cuty)
        pylab.cla()
        HPKCell.plot_band_cut(ax10, cuty)
        pylab.title('')
        pylab.ylabel('Energy (|t|)')
        pylab.xlabel('$k_y$ ($\pi$/a)')
        pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        #titleStr = testCell.type + ', modeType: ' + cellModeType + ' (Made with UnitCell class)' 
        #pylab.suptitle(titleStr)
#    
        
        
    if allThree:
#        fig6.set_size_inches(10.55, 7.4)
#        pylab.tight_layout()
#        fig6.set_size_inches(10.55, 11.25)
        fig6.set_size_inches(12.5, 11.25)
        pylab.tight_layout()
        pylab.show()
        
#        fig6.savefig('HPK_OSK_HW.png',transparent= False, dpi = 200)
#        fig6.savefig('HPK_OSK_HW.svg',transparent= True)
        
    else: 
        fig6.set_size_inches(12.5, 7.5)
        pylab.tight_layout()
        pylab.show()
    #    fig6.savefig('Figure6.pdf',transparent= True)
    #    fig6.savefig('HPK_OSK.png',transparent= False, dpi = 200)
    #    fig6.savefig('HPK_OSK.png',transparent= False, dpi = 400)
    #    fig6.savefig('HPK_OSK.svg',transparent= True)
    
    
    
    

    
    
    
    #pylab.rcParams.update({'font.size': 14})
else:
    fig6 = pylab.figure(6)    
    pylab.clf()
    pylab.suptitle('HPK, OSK Fig')
    pylab.close(fig6)
    
    
    
    
    
    
    
#%%
    
#############
#Figure 7: FW vs HW
###############
    
if FWHWFBFig:
#    gs = gridspec.GridSpec(2, 2,
#                           width_ratios=[1, 1],
#                           height_ratios=[1, 1.8]
#                           )
#    gs = gridspec.GridSpec(1, 4,
#                           width_ratios=[1.8,1.8,1,1],
#                           height_ratios=[1]
#                           )
    gs = gridspec.GridSpec(1, 4,
                           width_ratios=[2,2,1,1],
                           height_ratios=[1]
                           )
    
    cellModeType = 'FW'
    #cellModeType = 'HW'
    
    bigC = 95
#    smallC = 30
    
    bigC2 = 95*0.6
    
    #1) hyperbolic
    testHyperbolic = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
    testHyperbolic.populate(2, resonatorsOnly=False)
    resonators = testHyperbolic.get_all_resonators()
    nameStr = str(testHyperbolic.gon) + 'gon_' + str(testHyperbolic.vertex) + 'vertex_' + str(testHyperbolic.itter) + '_' + testHyperbolic.modeType
    hyperbolicLattice = GeneralLayout(resonators , modeType = testHyperbolic.modeType, name =  nameStr)
    
    fig7 = pylab.figure(7)
    pylab.clf()
    red = [2,4,6,7,9,36]
    yellow = [1,3,5,8,10,35]
    #ax1 = pylab.subplot(2,2,1)
    ax1 = pylab.subplot(gs[0])
#    ax1 = pylab.subplot(gs[0], adjustable='box', aspect=2)
    pylab.cla()
    hyperbolicLattice.draw_SDlinks(ax1, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
    xs = hyperbolicLattice.SDx
    ys = hyperbolicLattice.SDy
    pylab.sca(ax1)
    pylab.scatter(xs, ys, c = FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
    pylab.scatter(hyperbolicLattice.SDx[yellow], hyperbolicLattice.SDy[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(hyperbolicLattice.SDx[red], hyperbolicLattice.SDy[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    ax1.set_aspect('equal')
    ax1.axis('off')
    
    testHyperbolic = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'HW')
    testHyperbolic.populate(2, resonatorsOnly=False)
    resonators = testHyperbolic.get_all_resonators()
    nameStr = str(testHyperbolic.gon) + 'gon_' + str(testHyperbolic.vertex) + 'vertex_' + str(testHyperbolic.itter) + '_' + testHyperbolic.modeType
    hyperbolicLattice = GeneralLayout(resonators , modeType = testHyperbolic.modeType, name =  nameStr)
    
    ax2 = pylab.subplot(gs[1])
    #ax2 = pylab.subplot(2,2,2)
    pylab.cla()
    hyperbolicLattice.draw_SDlinks(ax2, color = HWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'b', alpha = HWlinkAlpha)
    xs = hyperbolicLattice.SDx
    ys = hyperbolicLattice.SDy
    pylab.sca(ax2)
    pylab.scatter(xs, ys ,c =  HWsiteColor, s = smallCdefault, marker = 'o', edgecolors = HWsiteEdgeColor, zorder = 5,  linewidth = 1, alpha = HWsiteAlpha)
    pylab.scatter(hyperbolicLattice.SDx[0:7], hyperbolicLattice.SDy[0:7], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    ax2.set_aspect('equal')
    ax2.axis('off')
    
    #2) HPK
    HPKCell = UnitCell('Huse')
    testEuclidLattice = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'FW')
    resonators = testEuclidLattice.resonators
    HPKLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'Huse')
    
    red2 = [25,32,41,46,49,65]
    yellow2 = [24,31,40,47,48,50]
    red2b = [48,42,65]
    yellow2b = [40,41,47,46]
#    ax3 = pylab.subplot(gs[2])
    ax3 = pylab.subplot(gs[2], adjustable='box', aspect=0.8)
    #ax3 = pylab.subplot(2,2,3)
    pylab.cla()
    pylab.sca(ax3)
    HPKLattice.draw_SDlinks(ax3, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy, c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = FWsiteAlpha)
    pylab.scatter(HPKLattice.SDx[yellow2], HPKLattice.SDy[yellow2], c =  stateColor1, s = bigC2, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(HPKLattice.SDx[red2], HPKLattice.SDy[red2], c =  stateColor2, s = bigC2, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
#    ax3.set_aspect('equal')
    ax3.axis('off')
    
    
    HPKCell = UnitCell('Huse')
    testEuclidLattice = EuclideanLayout(3,2,lattice_type = 'Huse', modeType = 'HW')
    resonators = testEuclidLattice.resonators
    HPKLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'Huse')
    
#    ax4 = pylab.subplot(gs[3])
    ax4 = pylab.subplot(gs[3], adjustable='box', aspect=0.8)
    #ax4 = pylab.subplot(2,2,4)
    pylab.cla()
    pylab.sca(ax4)
    HPKLattice.draw_SDlinks(ax4, color = HWlinkColor, linewidth = 2.5, minus_links = True, minus_color = HWminusLinkColor, alpha = HWlinkAlpha)
    pylab.scatter(HPKLattice.SDx, HPKLattice.SDy ,c = FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha = HWsiteAlpha)
    pylab.scatter(HPKLattice.SDx[red2b], HPKLattice.SDy[red2b], c =  stateColor2, s = bigC2, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    pylab.scatter(HPKLattice.SDx[yellow2b], HPKLattice.SDy[yellow2b], c =  stateColor1, s = bigC2, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
#    ax4.set_aspect('equal')
    ax4.axis('off')
#    fig7.set_size_inches(10, 15)
    fig7.set_size_inches(10, 3.8)
    pylab.tight_layout()
    pylab.show()
#    fig7.savefig('Figure7.pdf',transparent= True)
#    fig7.savefig('FWHW_FB.png',transparent= False, dpi = 200)
#    fig7.savefig('FWHW_FB.svg',transparent= True)
else:
    fig7 = pylab.figure(7)
    pylab.clf()
    pylab.suptitle('FW HW Flat Band States Fig')
    pylab.close(fig7)
















#%%

########
#figure 9: McLaughlin
#######

if McLaughlinFig:
    #use the first color convention
    splitCenterColor = n1_splitCenterColor
    splitEdgeColor = n1_splitEdgeColor
    splitCenterColor2 = n1_splitCenterColor2
    splitEdgeColor2 = n1_splitEdgeColor2
    LGsplitCenterColor = n1_LGsplitCenterColor
    splitStateColor1 = n1_splitStateColor1
    splitStateEdgeColor1 = n1_splitStateEdgeColor1
    
#    splitCenterColor = n2_splitCenterColor
#    splitEdgeColor = n2_splitEdgeColor
#    splitCenterColor2 = n2_splitCenterColor2
#    splitEdgeColor2 = n2_splitEdgeColor2
#    LGsplitCenterColor = n2_LGsplitCenterColor
#    splitStateColor1 = n2_splitStateColor1
#    splitStateEdgeColor1 = n2_splitStateEdgeColor1
    
    
    
    ######split tree
#    gs = gridspec.GridSpec(2, 2,
#                           width_ratios=[1, 1],
#                           height_ratios=[1, 0.1]
#                           )
    gs = gridspec.GridSpec(2, 3,
                           width_ratios=[1, 1, 1],
                           height_ratios=[1, 0.1]
                           )
    test1 = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
    resonators = test1.get_all_resonators()
    splitGraph = split_resonators(resonators, splitIn = 2)
    resonators = splitGraph
    testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'McLaughlinTree')
    
    LGresonators = generate_line_graph(resonators)
    resonators = LGresonators
    testLattice2 = GeneralLayout(resonators , modeType = test1.modeType, name =  'LG_McLaughlinTree')
    
    
    fig17 = pylab.figure(17)
    pylab.clf()
    
#    #split tree
#    ax = pylab.subplot(gs[0])
#    pylab.cla()
##    testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#    testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
#    xs = testLattice.coords[:,0]
#    ys = testLattice.coords[:,1]
#    pylab.sca(ax)
##    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#    pylab.scatter(xs, ys ,c =  'skyblue', s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#    ax.set_aspect('equal')
#    ax.axis('off')
    
    #split tree
    ax = pylab.subplot(gs[0])
    pylab.cla()
#    testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    testLattice.draw_resonator_lattice(ax, color = 'b', alpha = 1 , linewidth = 2.5)
    xs = testLattice.coords[:,0]
    ys = testLattice.coords[:,1]
    pylab.sca(ax)
#    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#    pylab.scatter(xs, ys ,c =  'gold', s = smallCdefault, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
    pylab.scatter(xs, ys ,c =  splitCenterColor, s = smallCdefault, marker = 'o', edgecolors = splitEdgeColor, zorder = 5)

    
    xs = test1.coords[:,0]
    ys = test1.coords[:,1]
#    pylab.scatter(xs, ys ,c =  'deepskyblue', s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
#    pylab.scatter(xs, ys ,c = 'lightsteelblue', s = smallCdefault, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
    pylab.scatter(xs, ys ,c =  splitCenterColor2, s = smallCdefault, marker = 'o', edgecolors = splitEdgeColor2, zorder = 5)
    ax.set_aspect('equal')
    ax.axis('off')

    
    #McLaughlin
    ax2 = pylab.subplot(gs[1])
#    testLattice.draw_SDlinks(ax2, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    testLattice.draw_SDlinks(ax2, color = layoutLineColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax2)
#    pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    ax2.set_aspect('equal')
    ax2.axis('off')

    
    #L(McLaughlin)
    ax3 = pylab.subplot(gs[2])
    testLattice2.draw_SDlinks(ax3, color = FWlinkColor, linewidth = 1.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax3)
    pylab.scatter(testLattice2.SDx, testLattice2.SDy ,c =  FWsiteColor, s = 20, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    ax3.set_aspect('equal')
    ax3.axis('off')
    
    
    #pylab.suptitle(testLattice.name)
    #0) split DOS
    freq_start = -3.00
    freq_stop = 3.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [0.] # energies of flat bands
    N = len(freqs)
    a = numpy.sqrt(3 -2*numpy.sqrt(2))
    aprime = 0.4
    b = numpy.sqrt(3 +2*numpy.sqrt(2))
    bprime = 2.4
    bandX = [(-bprime,-aprime), (aprime,bprime)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    Amp_imp = 0.3
    Amp_band = 0.1
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[3], adjustable='box', aspect=3)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)     
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-3,0,3], rotation='horizontal')
#    pylab.xticks([-2,-1.965, -0.149, 0,3], rotation='horizontal')
    #ax1.set_ylim([0,0.3])
    
    
    #pylab.suptitle(testLattice.name)
    #1) DOS
    freq_start = -2.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-1.92,-0.16), (1.16,2.92)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    Amp_imp = 0.3
    Amp_band = 0.1
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[4], adjustable='box', aspect=3)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)     
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,0,3], rotation='horizontal')
#    pylab.xticks([-2,-1.965, -0.149, 0,3], rotation='horizontal')
    #ax1.set_ylim([0,0.3])
    
    #2) DOS
    freq_start = -2.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
#    impX = [-2., -1., 0.4] # energies of flat bands
    impX = [-2., -1., 1] # energies of flat bands
    N = len(freqs)
#    bandX = [(-0.92,0.32), (1.16,2.92)] # start and end points for normal bands
    bandX = [(-0.92,-0.16+1), (2.16,3.92)] # start and end points for normal bands
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    Amp_imp = 0.3
    Amp_band = 0.1
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
        
    ax1 = pylab.subplot(gs[5], adjustable='box', aspect=3)
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)        
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,-1,1,4], rotation='horizontal')
    
    
    
    #ax1.set_ylim([0,0.3])
    
    fig17.set_size_inches(12, 5.5)
    pylab.tight_layout()
    pylab.show()
#    fig17.savefig('Figure9.pdf',transparent= True)
#    fig17.savefig('McLaughlin.png',transparent= False, dpi = 200)
#    fig17.savefig('McLaughlin.svg',transparent= True)
else:
    fig17 = pylab.figure(17)
    pylab.clf()
    pylab.suptitle('McLaughlinFig')
    pylab.close(fig17)
    
    
    
    
    
    
    
    
    
    
    
    
    

#%% 
###########
#Figure 8, finite hyperbolic and regularized hyperbolic
###########



if FiniteHyperbolicFig:

    showRegularized = True
    showHardWall = True
    minimal = True
    #minimal = False
    
    if showHardWall:
        if minimal:
    #        testDict = pickle.load( open('8gon_3vertex_4_FW_DOSminimal.pkl', "rb"))
    #        testDict = pickle.load( open('8gon_3vertex_5_FW_DOSminimal.pkl', "rb"))
            
            testDict = pickle.load( open(os.path.join(DataPickleFolderPath,'7gon_3vertex_6_FW_DOSminimal.pkl'), "rb"))
    #        testDict = pickle.load( open(os.path.join(DataPickleFolderPath,'7gon_3vertex_6_FW_DOSminimal.pkl'), "rb"), encoding = 'latin1')
        else:
        #    test = PlanarLayout(file_path = '7gon_3vertex_ 3.pkl')
        #    test = PlanarLayout(file_path = '7gon_3vertex_ 4.pkl')
        #    test = PlanarLayout(file_path = '7gon_3vertex_ 5.pkl')
            #test = PlanarLayout(file_path = '7gon_3vertex_ 6.pkl')
            
            
            test = PlanarLayout(file_path = os.path.join(DataPickleFolderPath,'8gon_3vertex_ 3.pkl'))
        #    test = PlanarLayout(file_path = '8gon_3vertex_ 4.pkl')
        
        
       
        #set up frequency sweep
        freq_range = 4.01
        freq_res = 0.04
    
        
        freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
        freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
        
    
    
        if minimal:
            Evals = testDict['Es']
            Eorder = testDict['Eorder']
            itt = testDict['itter']
        else:
            [Evals, Eorder] =  [test.Es, test.Eorder]
            itt = test.itter
    
        
        [DOS, bins_out] = numpy.histogram(Evals, freq_bins)
    
        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
        binWidth = bins_out[1] - bins_out[0]
    
        fig8 = pylab.figure(45)
        pylab.clf()
        ax1 = pylab.subplot(1,4,1)
        pylab.bar(bins_centers, 1.*DOS/len(Evals), width = binWidth, color =  'firebrick', label = str(itt), alpha = 1)
        
            
        pylab.xlabel('Energy (|t|)')
        pylab.ylabel('Density of States')
#        ax1.legend(loc = 'upper right')
#        pylab.title('Heptagon')
        ax1.set_xlim([-2.5,4])
        ax1.set_ylim([0,40])
        ax1.set_ylim([0,0.04])
        
    
        testDict = pickle.load( open(os.path.join(DataPickleFolderPath,'8gon_3vertex_5_FW_DOSminimal.pkl'), "rb"))
    #    testDict = pickle.load( open(os.path.join(DataPickleFolderPath,'8gon_3vertex_5_FW_DOSminimal.pkl'), "rb"), encoding = 'latin1')
        #set up frequency sweep
        freq_range = 4.01
        freq_res = 0.04
    
        
        freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
        freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
        
    
    
        if minimal:
            Evals = testDict['Es']
            Eorder = testDict['Eorder']
            itt = testDict['itter']
        else:
            [Evals, Eorder] =  [test.Es, test.Eorder]
            itt = test.itter
    
        
        [DOS, bins_out] = numpy.histogram(Evals, freq_bins)
    
        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
        binWidth = bins_out[1] - bins_out[0]
    
        ax1 = pylab.subplot(1,4,2)
        pylab.bar(bins_centers, 1.*DOS/len(Evals), width = binWidth, color =  'dodgerblue', label = str(itt), alpha = 1)
        
            
        pylab.xlabel('Energy (|t|)')
        pylab.ylabel('Density of States')
#        ax1.legend(loc = 'upper right')
#        pylab.title('Octagon')
        ax1.set_xlim([-2.5,4])
        ax1.set_ylim([0,40])
        ax1.set_ylim([0,0.04])
    
    
    
    
    
    if showRegularized:
        ###sample figure of FW v HW regu;arized
    #    loaded = pickle.load(open('7gon_3vertex_4_FW_reg.pkl', "rb" ) , encoding = 'latin1')
    #    loaded2 = pickle.load(open('7gon_3vertex_4_HW_reg.pkl', "rb" ) , encoding = 'latin1')
        loaded = pickle.load(open(os.path.join(DataPickleFolderPath,'7gon_3vertex_5_FW_reg.pkl'), "rb" ))
        loaded2 = pickle.load(open(os.path.join(DataPickleFolderPath,'7gon_3vertex_5_HW_reg.pkl'), "rb" ))
    #    loaded = pickle.load(open(os.path.join(DataPickleFolderPath,'7gon_3vertex_5_FW_reg.pkl'), "rb" ) , encoding = 'latin1')
    #    loaded2 = pickle.load(open(os.path.join(DataPickleFolderPath,'7gon_3vertex_5_HW_reg.pkl'), "rb" ) , encoding = 'latin1')
    #    pylab.figure(44)
    #    pylab.clf()
    #    ax = pylab.subplot(1,1,1)
    #    
    #    xs = numpy.linspace(0,1,len(loaded['regEs']))
    #    xs2 = numpy.linspace(0,1,len(loaded2['regEs']))
    #    pylab.plot(xs, loaded['regEs'], color =  'mediumblue', marker = '.', linestyle = '', label = loaded['name'])
    #    pylab.plot(xs2, loaded2['regEs'], color =  'deepskyblue', marker = '.', linestyle = '', label = loaded2['name'])
    #    
    #    ax.set_ylim([-2.5, 4.2])
    #    pylab.ylabel('Energy (|t|)')
    #    pylab.xlabel('normalized eigenvector index')
    #    pylab.title('regularized layout')
    #    
    #    ax.legend(loc = 'upper left')
    #    pylab.show()
            
            
        ####   sample DOS plot
        #set up frequency sweep
        freq_range = 4.01
        freq_res = 0.04
        
        freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
        freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
        
        #[fullDOS, bins_out] = numpy.histogram(loaded['regEs'], freq_bins)
        #[fullDOS2, bins_out] = numpy.histogram(loaded2['regEs'], freq_bins)
        
        ax1 = pylab.subplot(1,4,4)
        
        
        [DOS_unnorm, bins_out] = numpy.histogram(loaded['regEs'], freq_bins)
        DOS= DOS_unnorm*1.0/len(loaded['regEs'])
        
        [DOS2_unnorm, bins_out] = numpy.histogram(loaded2['regEs'], freq_bins)
        DOS2= DOS2_unnorm*1.0/len(loaded2['regEs'])
        
        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
        binWidth = bins_out[1] - bins_out[0]
        
#        pylab.bar(bins_centers, 1.*DOS, width = binWidth, color = 'mediumblue', label = loaded['name'], alpha = 0.6)
#        pylab.bar(bins_centers, 1.*DOS2, width = binWidth, color = 'deepskyblue', label = loaded2['name'], alpha = 0.6) 
        pylab.bar(bins_centers, 1.*DOS, width = binWidth, color = 'mediumblue', label = 'full-wave', alpha = 0.6)
        pylab.bar(bins_centers, 1.*DOS2, width = binWidth, color = 'deepskyblue', label = 'half-wave', alpha = 0.6)  
            
        pylab.xlabel('Energy (|t|)')
        pylab.ylabel('Density of States')
#        pylab.title('DOs for FW and HW regularized')
        ax1.set_xlim([-2.5,4.1])
        #ax1.set_xlim([-2.5,4.0])
        ax1.set_ylim([0,0.035])
        ax1.legend(loc= 'upper right')
        
        
        
        
            
    
    
    
    
    #%######
    #figure 8, in part
    ######
    
    #test = PlanarLayout(file_path = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/7gon_3vertex_ 3.pkl')
    test = PlanarLayout(file_path = os.path.join(DataPickleFolderPath,'8gon_3vertex_ 3.pkl'))
    
    
    #ax = pylab.subplot(1,1,1)
    #xaxis = numpy.linspace(0,1, len(test.Eorder))
    #pylab.plot(xaxis, test.Es[test.Eorder], color = 'mediumblue', marker = '.', linestyle = '', label = str(test.gon), alpha = 1)
    #pylab.xlabel('Eigenvalue Index')
    #pylab.ylabel('Energy (|t|)')
    #pylab.title('finite hyperbolic')
    #pylab.show()
    
    
    
    
    
    
    
    ##find properties versus polygon
    itterGons = numpy.asarray([5,6,7,8,9,10,11,12,13,14,15,16,17])  #DO NOT CHANGE THIS!!!!
    itterDepth = 2 #OR THIS!!!!
    
    itterModeType = 'FW'
    pickleName = 'propVpolygon_reg_' + str(itterModeType) + '_' + str(itterDepth) + '.pkl'
    pickleName_Wall = 'propVpolygon_' + str(itterModeType) + '_' + str(itterDepth) + '.pkl'
    
    
    #pull old pickled dictionary
    #storageDict = pickle.load(open(os.path.join(DataPickleFolderPath, pickleName), "rb" ), encoding = 'latin1' )
    storageDict = pickle.load(open(os.path.join(DataPickleFolderPath, pickleName), "rb" ) )
    FBSizeArray = storageDict['FBsize']
    FBfracArray = storageDict['FBfraction'] 
    gapSizeArray = storageDict['gap']
    sysSizeArray = storageDict['SysSize']
    
    #storageDict_Wall = pickle.load(open(os.path.join(DataPickleFolderPath, pickleName_Wall), "rb" ), encoding = 'latin1' )
    storageDict_Wall = pickle.load(open(os.path.join(DataPickleFolderPath, pickleName_Wall), "rb" ))
    FBSizeArray_Wall = storageDict_Wall['FBsize']
    FBfracArray_Wall = storageDict_Wall['FBfraction'] 
    gapSizeArray_Wall = storageDict_Wall['gap']
    sysSizeArray_Wall = storageDict_Wall['SysSize']
    
    
    
    
    
    ax = pylab.subplot(1,4,3)
    gapSize = numpy.exp(-(itterGons/2.-2.2)) + 0.18
    gapSize = numpy.exp(-(itterGons*1.0-5.)**2) + 0.18
    gapSize = numpy.exp(-(itterGons/2.0 - 2)**1.25) + 0.175
    gapSize = 0.6/(itterGons-4.)**1.5 + 0.165
    
    fit_guess = [4., 2., 1., 0.175] # exponential
    
    if itterGons[0] == 6:
        pylab.plot(itterGons[1::2],gapSizeArray_Wall[1::2] , color = 'firebrick', marker = '.', linestyle = '', label = 'hard wall', alpha = 1)
        pylab.plot(itterGons[1::2],gapSizeArray[1::2] , color = 'mediumblue', marker = '.', linestyle = '', label = 'regularized odd', alpha = 1)
        pylab.plot(itterGons[2::2],gapSizeArray[2::2] , color = 'dodgerblue', marker = '.', linestyle = '', label = 'regularized even', alpha = 1)
    elif itterGons[0] == 5:
        oddGons = itterGons[0::2]
        oddGaps = gapSizeArray[0::2]
        oddGaps_Wall = gapSizeArray_Wall[0::2]
    #    fit_out, pcov = curve_fit(expff, oddGons, oddGaps, p0 = fit_guess)
        pylab.plot(oddGons,oddGaps_Wall , color = 'firebrick', marker = 'D', linestyle = '', label = 'hard wall', alpha = 1, markersize = 5)
        pylab.plot(oddGons,oddGaps , color = 'mediumblue', marker = 'v', linestyle = '', label = 'regularized odd', alpha = 1, markersize = 5)
        pylab.plot(itterGons[3::2],gapSizeArray[3::2] , color = 'dodgerblue', marker = '^', linestyle = '', label = 'regularized even', alpha = 1, markersize = 5)
    
    pylab.plot(itterGons, 1./itterGons**2, 'dodgerblue', linewidth = 0.75, label = '1/$k^2$ bound' )
    
    ax.set_ylim([0,0.9])
    
#    pylab.title('gap size')
    pylab.xlabel('Polygon')
    pylab.ylabel('Gap Energy (|t|)')
    ax.legend(loc = 'upper right')
    
    #pylab.suptitle('properties of regularized layouts: itter depth = ' + str(itterDepth))
    fig8.set_size_inches(12,3.8)
    pylab.tight_layout()
    pylab.show()
#    fig8.savefig('Figure8.pdf', transparent= True)
#    fig8.savefig('FiniteHyperbolic.png', transparent= False, dpi = 200)
#    fig8.savefig('FiniteHyperbolic.svg', transparent= True)
else:
    fig8= pylab.figure(45)
    pylab.clf()
    pylab.suptitle('Finite and REgularized Hyperbolics Fig')
    pylab.close(fig8)













#%% 
###########
#Figure 8, finite hyperbolic and regularized hyperbolic, alternate version with HW too
###########



if FiniteHyperbolicFig2:

    showRegularized = True
    showHardWall = True
    minimal = True #use small pickles made for this purpose
    #minimal = False
    
    numCols = 3
    numRows = 2
    
    Scale = 100.
    
    DataPickleFolderPath2 = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/'
    
    pylab.rcParams.update({'font.size': 12})
    
    if showHardWall:
        #####hard wall heptagon kagome, FW
        if minimal:
            testDict = pickle.load( open(os.path.join(DataPickleFolderPath,'7gon_3vertex_6_FW_DOSminimal.pkl'), "rb"))
        else:
            test = PlanarLayout(file_path = os.path.join(DataPickleFolderPath,'8gon_3vertex_ 3.pkl'))

        
        
       
        #set up frequency sweep
        freq_range = 4.01
        freq_res = 0.04
    
        
        freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
        freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
        
    
    
        if minimal:
            Evals = testDict['Es']
            Eorder = testDict['Eorder']
            itt = testDict['itter']
        else:
            [Evals, Eorder] =  [test.Es, test.Eorder]
            itt = test.itter
    
        
        [DOS, bins_out] = numpy.histogram(Evals, freq_bins)
    
        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
        binWidth = bins_out[1] - bins_out[0]
    
        fig8 = pylab.figure(45)
        pylab.clf()
        ax1 = pylab.subplot(numRows,numCols,1)
        pylab.bar(bins_centers, Scale*DOS/len(Evals), width = binWidth, color =  'mediumblue', label = str(itt), alpha = 1)
        
            
        pylab.xlabel('Energy (|t|)')
        pylab.ylabel('Density of States (%)')
#        ax1.legend(loc = 'upper right')
#        pylab.title('Heptagon')
        ax1.set_xlim([-2.5,4])
        ax1.set_ylim([0,40])
        ax1.set_ylim([0,0.04*Scale])
        
        #####hard wall ocatgon kagome, FW
        testDict = pickle.load( open(os.path.join(DataPickleFolderPath,'8gon_3vertex_5_FW_DOSminimal.pkl'), "rb"))
        #set up frequency sweep
        freq_range = 4.01
        freq_res = 0.04
    
        
        freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
        freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
        
    
    
        if minimal:
            Evals = testDict['Es']
            Eorder = testDict['Eorder']
            itt = testDict['itter']
        else:
            [Evals, Eorder] =  [test.Es, test.Eorder]
            itt = test.itter
    
        
        [DOS, bins_out] = numpy.histogram(Evals, freq_bins)
    
        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
        binWidth = bins_out[1] - bins_out[0]
    
        ax1 = pylab.subplot(numRows,numCols,4)
        pylab.bar(bins_centers, Scale*DOS/len(Evals), width = binWidth, color =  'darkred', label = str(itt), alpha = 1)
        
            
        pylab.xlabel('Energy (|t|)')
        pylab.ylabel('Density of States (%)')
#        ax1.legend(loc = 'upper right')
#        pylab.title('Octagon')
        ax1.set_xlim([-2.5,4])
        ax1.set_ylim([0,40])
        ax1.set_ylim([0,0.04*Scale])
        
        
        
        #####hard wall heptagon kagome, HW
        testDict = pickle.load( open(os.path.join(DataPickleFolderPath,'7gon_3vertex_6_HW_DOSminimal.pkl'), "rb"))
        #set up frequency sweep
        freq_range = 4.01
        freq_res = 0.04
    
        
        freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
        freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
        
    
    
        if minimal:
            Evals = testDict['Es']
            Eorder = testDict['Eorder']
            itt = testDict['itter']
        else:
            [Evals, Eorder] =  [test.Es, test.Eorder]
            itt = test.itter
    
        
        [DOS, bins_out] = numpy.histogram(Evals, freq_bins)
    
        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
        binWidth = bins_out[1] - bins_out[0]
    
        ax1 = pylab.subplot(numRows,numCols,2)
        pylab.bar(bins_centers, Scale*DOS/len(Evals), width = binWidth, color =  'dodgerblue', label = str(itt), alpha = 1)
        
            
        pylab.xlabel('Energy (|t|)')
        pylab.ylabel('Density of States (%)')
#        ax1.legend(loc = 'upper right')
#        pylab.title('Octagon')
        ax1.set_xlim([-2.5,4])
        ax1.set_ylim([0,40])
        ax1.set_ylim([0,0.04*Scale])
        
        
#        #####hard wall octagon kagome, HW
#        testDict = pickle.load( open(os.path.join(DataPickleFolderPath,'8gon_3vertex_5_HW_DOSminimal.pkl'), "rb"))
#        #set up frequency sweep
#        freq_range = 4.01
#        freq_res = 0.04
#    
#        
#        freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
#        freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
#        
#    
#    
#        if minimal:
#            Evals = testDict['Es']
#            Eorder = testDict['Eorder']
#            itt = testDict['itter']
#        else:
#            [Evals, Eorder] =  [test.Es, test.Eorder]
#            itt = test.itter
#    
#        
#        [DOS, bins_out] = numpy.histogram(Evals, freq_bins)
#    
#        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
#        binWidth = bins_out[1] - bins_out[0]
#    
#        ax1 = pylab.subplot(numRows,numCols,4)
#        pylab.bar(bins_centers, 1.*DOS/len(Evals), width = binWidth, color =  'firebrick', label = str(itt), alpha = 1)
#        
#            
#        pylab.xlabel('Energy (|t|)')
#        pylab.ylabel('Density of States')
##        ax1.legend(loc = 'upper right')
##        pylab.title('Octagon')
#        ax1.set_xlim([-2.5,4])
#        ax1.set_ylim([0,40])
#        ax1.set_ylim([0,0.04])
    
    
    
    
    
    if showRegularized:
        ### FW v HW regularized, heptagon kagome
        loaded = pickle.load(open(os.path.join(DataPickleFolderPath,'7gon_3vertex_5_FW_reg.pkl'), "rb" ))
        loaded2 = pickle.load(open(os.path.join(DataPickleFolderPath,'7gon_3vertex_5_HW_reg.pkl'), "rb" ))
            
            
        ####   sample DOS plot
        #set up frequency sweep
        freq_range = 4.01
        freq_res = 0.04
        
        freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
        freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
        
        #[fullDOS, bins_out] = numpy.histogram(loaded['regEs'], freq_bins)
        #[fullDOS2, bins_out] = numpy.histogram(loaded2['regEs'], freq_bins)
        
        ax1 = pylab.subplot(numRows,numCols,3)
        
        
        [DOS_unnorm, bins_out] = numpy.histogram(loaded['regEs'], freq_bins)
        DOS= DOS_unnorm*1.0/len(loaded['regEs'])
        
        [DOS2_unnorm, bins_out] = numpy.histogram(loaded2['regEs'], freq_bins)
        DOS2= DOS2_unnorm*1.0/len(loaded2['regEs'])
        
        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
        binWidth = bins_out[1] - bins_out[0]
        
#        pylab.bar(bins_centers, 1.*DOS, width = binWidth, color = 'mediumblue', label = loaded['name'], alpha = 0.6)
#        pylab.bar(bins_centers, 1.*DOS2, width = binWidth, color = 'deepskyblue', label = loaded2['name'], alpha = 0.6) 
        pylab.bar(bins_centers, Scale*DOS, width = binWidth, color = 'mediumblue', label = 'full-wave', alpha = 0.8)
        pylab.bar(bins_centers, Scale*DOS2, width = binWidth, color = 'deepskyblue', label = 'half-wave', alpha = 0.6)  
            
        pylab.xlabel('Energy (|t|)')
        pylab.ylabel('Density of States (%)')
#        pylab.title('DOs for FW and HW regularized')
        ax1.set_xlim([-2.5,4.1])
        #ax1.set_xlim([-2.5,4.0])
        ax1.set_ylim([0,0.035*Scale])
        ax1.legend(loc= 'upper right')
        
        
        
        
        ######## FW v HW regularized, ocatagon kagome
        loaded = pickle.load(open(os.path.join(DataPickleFolderPath,'8gon_3vertex_4_FW_reg.pkl'), "rb" ))
            
            
        ####   sample DOS plot
        #set up frequency sweep
        freq_range = 4.01
        freq_res = 0.04
        
        freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
        freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
        
        #[fullDOS, bins_out] = numpy.histogram(loaded['regEs'], freq_bins)
        #[fullDOS2, bins_out] = numpy.histogram(loaded2['regEs'], freq_bins)
        
        ax1 = pylab.subplot(numRows,numCols,5)
#        ax1 = pylab.subplot(numRows,numCols-1,5)
        
        
        [DOS_unnorm, bins_out] = numpy.histogram(loaded['regEs'], freq_bins)
        DOS= DOS_unnorm*1.0/len(loaded['regEs'])
        
        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
        binWidth = bins_out[1] - bins_out[0]
         
        pylab.bar(bins_centers, Scale*DOS, width = binWidth, color = 'firebrick', label = 'full-wave', alpha = 1)
            
        pylab.xlabel('Energy (|t|)')
        pylab.ylabel('Density of States (%)')
#        pylab.title('DOs for FW and HW regularized')
        ax1.set_xlim([-2.5,4.1])
        #ax1.set_xlim([-2.5,4.0])
        ax1.set_ylim([0,0.035*Scale])
#        ax1.legend(loc= 'upper right')
        
        
        
        
            
    
    
    
#    
#    #%######
#    #hardwall octagon DOS
#    ######
#    
#    #test = PlanarLayout(file_path = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/7gon_3vertex_ 3.pkl')
#    test = PlanarLayout(file_path = os.path.join(DataPickleFolderPath,'8gon_3vertex_ 3.pkl'))
    
    

    
    
    ##find properties versus polygon
    itterGons = numpy.asarray([5,6,7,8,9,10,11,12,13,14,15,16,17])  #DO NOT CHANGE THIS!!!!
    itterDepth = 2 #OR THIS!!!!
    
    itterModeType = 'FW'
    pickleName = 'propVpolygon_reg_' + str(itterModeType) + '_' + str(itterDepth) + '.pkl'
    pickleName_Wall = 'propVpolygon_' + str(itterModeType) + '_' + str(itterDepth) + '.pkl'
    
    pickleName2 = 'propVpolygon_reg_' + 'HW' + '_' + str(itterDepth) + '.pkl'
    
    
    #pull old pickled dictionary
    storageDict = pickle.load(open(os.path.join(DataPickleFolderPath, pickleName), "rb" ) )
    FBSizeArray = storageDict['FBsize']
    FBfracArray = storageDict['FBfraction'] 
    gapSizeArray = storageDict['gap']
    sysSizeArray = storageDict['SysSize']
    
    storageDict_Wall = pickle.load(open(os.path.join(DataPickleFolderPath, pickleName_Wall), "rb" ))
    FBSizeArray_Wall = storageDict_Wall['FBsize']
    FBfracArray_Wall = storageDict_Wall['FBfraction'] 
    gapSizeArray_Wall = storageDict_Wall['gap']
    sysSizeArray_Wall = storageDict_Wall['SysSize']
    
    storageDict2 = pickle.load(open(os.path.join(DataPickleFolderPath, pickleName2), "rb" ) )
    FBSizeArray2 = storageDict2['FBsize']
    FBfracArray2 = storageDict2['FBfraction'] 
    gapSizeArray2 = storageDict2['gap']
    sysSizeArray2 = storageDict2['SysSize']
    
    
    
    
    
    ax = pylab.subplot(numRows,numCols,6)
    gapSize = numpy.exp(-(itterGons/2.-2.2)) + 0.18
    gapSize = numpy.exp(-(itterGons*1.0-5.)**2) + 0.18
    gapSize = numpy.exp(-(itterGons/2.0 - 2)**1.25) + 0.175
    gapSize = 0.6/(itterGons-4.)**1.5 + 0.165
    
    fit_guess = [4., 2., 1., 0.175] # exponential
    
    if itterGons[0] == 6:
        pylab.plot(itterGons[1::2],gapSizeArray_Wall[1::2] , color = 'firebrick', marker = '.', linestyle = '', label = 'hard wall', alpha = 1)
        pylab.plot(itterGons[1::2],gapSizeArray[1::2] , color = 'mediumblue', marker = '.', linestyle = '', label = 'regularized odd', alpha = 1)
        pylab.plot(itterGons[2::2],gapSizeArray[2::2] , color = 'dodgerblue', marker = '.', linestyle = '', label = 'regularized even', alpha = 1)
    elif itterGons[0] == 5:
        oddGons = itterGons[0::2]
        oddGaps = gapSizeArray[0::2]
        oddGaps_Wall = gapSizeArray_Wall[0::2]
    #    fit_out, pcov = curve_fit(expff, oddGons, oddGaps, p0 = fit_guess)
        pylab.plot(itterGons[3:],gapSizeArray2[3:] , color = 'darkgoldenrod', marker = 'o', linestyle = '', label = 'HW regularized', alpha = 1, markersize = 5)
        pylab.plot(oddGons,oddGaps_Wall , color = 'firebrick', marker = 'D', linestyle = '', label = 'FW hard wall', alpha = 1, markersize = 5)
        pylab.plot(oddGons,oddGaps , color = 'mediumblue', marker = 'v', linestyle = '', label = 'FW regularized odd', alpha = 1, markersize = 5)
        pylab.plot(itterGons[3::2],gapSizeArray[3::2] , color = 'dodgerblue', marker = '^', linestyle = '', label = 'FW regularized even', alpha = 1, markersize = 5)
        
    
    pylab.plot(itterGons, 1./itterGons**2, 'dodgerblue', linewidth = 0.75, label = '1/$k^2$ bound' )
    
    ax.set_ylim([0,0.9])
    
#    pylab.title('gap size')
    pylab.xlabel('Polygon')
    pylab.ylabel('Gap Energy (|t|)')
    ax.legend(loc = 'upper right')
    
    #pylab.suptitle('properties of regularized layouts: itter depth = ' + str(itterDepth))
    fig8.set_size_inches(12,6.4)
    pylab.tight_layout()
    pylab.show()
#    fig8.savefig('Figure8.pdf', transparent= True)
#    fig8.savefig('FiniteHyperbolic.png', transparent= False, dpi = 200)
#    fig8.savefig('FiniteHyperbolic.svg', transparent= True)
else:
    fig8= pylab.figure(45)
    pylab.clf()
    pylab.suptitle('Finite and REgularized Hyperbolics Fig')
    pylab.close(fig8)















#%%

############
#sample graphs for peter notes
######

if TsubKFig:
    ##1 tetrahedron
    gs = gridspec.GridSpec(2, 3,
                           width_ratios=[1, 1, 1],
                           height_ratios=[0.7, 1]
                           )
    a = 1
    c = 0.8*a
    b = 2*a
    lilGraph = numpy.zeros((6,4))
    lilGraph[0,:] = [-a,0,a,0]
    lilGraph[1,:] = [-a,0,0,c]
    lilGraph[2,:] = [a,0,0,c]
    lilGraph[3,:] = [0, c,0 ,b]
    lilGraph[4,:] = [-a, 0,0 ,b]
    lilGraph[5,:] = [a, 0,0 ,b]
    nameStr = 'tetrahedron'
    testLattice = GeneralLayout(lilGraph , modeType = 'FW', name =  nameStr)
    
    fig47 = pylab.figure(47)
    pylab.clf()
    ax = pylab.subplot(gs[0])
    pylab.cla()
    testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    #testLattice.save(name='test.pkl', protocol=2)
    #testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
    xs = testLattice.coords[:,0]
    ys = testLattice.coords[:,1]
    pylab.sca(ax)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    ##1 cube
    a = 0.8
    b = 2
    lilGraph = numpy.zeros((12,4))
    lilGraph[0,:] = [-a,a,a,a]
    lilGraph[1,:] = [-a,-a,a,-a]
    lilGraph[2,:] = [a,a,a,-a]
    lilGraph[3,:] = [-a, a,-a,-a]
    lilGraph[4,:] = [-b,b,b,b]
    lilGraph[5,:] = [-b,-b,b,-b]
    lilGraph[6,:] = [b,b,b,-b]
    lilGraph[7,:] = [-b, b,-b,-b]
    lilGraph[8,:] = [a,a,b,b]
    lilGraph[9,:] = [a,-a,b,-b]
    lilGraph[10,:] = [-a,a,-b,b]
    lilGraph[11,:] = [-a,-a,-b,-b]
    
    nameStr = 'cube'
    testLattice = GeneralLayout(lilGraph , modeType = 'FW', name =  nameStr)
    ax = pylab.subplot(gs[1])
    pylab.cla()
    testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    #testLattice.save(name='test.pkl', protocol=2)
    #testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
    xs = testLattice.coords[:,0]
    ys = testLattice.coords[:,1]
    pylab.sca(ax)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    ##4 dodecahedron
    testLattice = GeneralLayout(file_path = os.path.join(DataPickleFolderPath,'dodecahedron.pkl'))
    
    ax = pylab.subplot(gs[2])
    pylab.cla()
    testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    #testLattice.save(name='test.pkl', protocol=2)
    #testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
    xs = testLattice.coords[:,0]
    ys = testLattice.coords[:,1]
    pylab.sca(ax)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    ax.set_aspect('equal')
    ax.axis('off')
    #pylab.title('layout graph')
    
    ##5 Graphene
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    #testEuclidLattice = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'HW')
    kagomeLattice = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    ax = pylab.subplot(gs[3])
    pylab.cla()
    kagomeLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    #testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
    xs = kagomeLattice.coords[:,0]
    ys = kagomeLattice.coords[:,1]
    pylab.sca(ax)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    ax.set_aspect('equal')
    ax.axis('off')
    #pylab.title('layout graph')
    
    
    fig47.set_size_inches(10, 5)
    pylab.tight_layout()
    pylab.show()
#    fig47.savefig('Figure10_Peter.pdf', transparent= True)
#    fig47.savefig('TsubK.png', transparent= False, dpi = 200)
#    fig47.savefig('TsubK.svg', transparent= True)
else:
    fig47 = pylab.figure(47)
    pylab.clf()
    pylab.suptitle('PeterNotesFig') 
    pylab.close(fig47)
    
    
    
    
    
#%%
    
    
########
#extra FB states for Mattias
#####
def build_double_plaquette_state(layout, maxItter = -1):
    if maxItter > layout.itter:
            raise ValueError('dont have this many itterations')
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
            ring_state = scipy.arange(0,layout.gon-3,1)
            print(len(ring_state))
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

def seperate_pm(state):
    plusInds = numpy.where(state>0)[0]
    minusInds = numpy.where(state<0)[0]
    
    return plusInds, minusInds
    
    
    
if ExtraFBstatesFig:
    fig74 = pylab.figure(74)
    pylab.clf()
    
    bigC = 90
    
    gs = gridspec.GridSpec(1, 4,
                           width_ratios=[0.8, 1, 1, 1],
                           height_ratios=[1]
                           )
    
    #1) kagome
#    ax = pylab.subplot(1,4,1)
    ax = pylab.subplot(gs[0])
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    kagomeLattice = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)

    kagomeLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax)
    pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
    red = [3,5,13]
    yellow = [4,8,9]
    shiftx, shifty = testEuclidLattice.unitcell.a1
    pylab.scatter(kagomeLattice.SDx[yellow]+shiftx, kagomeLattice.SDy[yellow]+shifty, c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(kagomeLattice.SDx[red]+shiftx, kagomeLattice.SDy[red]+shifty, c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    ax.set_aspect('equal')
    ax.axis('off')
    
    
    #1) heptagon kagome
#    ax = pylab.subplot(1,4,2)
    ax = pylab.subplot(gs[1])
    
    testHyperbolic = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
    testHyperbolic.populate(2)
    resonators = testHyperbolic.get_all_resonators()
    nameStr = str(testHyperbolic.gon) + 'gon_' + str(testHyperbolic.vertex) + 'vertex_' + str(testHyperbolic.itter) + '_' + testHyperbolic.modeType
    hyperbolicLattice = GeneralLayout(resonators , modeType = testHyperbolic.modeType, name =  nameStr)
    
    hyperbolicLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax)
    pylab.scatter(hyperbolicLattice.SDx, hyperbolicLattice.SDy ,c =  FWsiteColor, s = 20, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
#    red = [2,4,6,7,9,36]
#    yellow = [1,3,5,8,10,35]
    FBstate =  build_double_plaquette_state(testHyperbolic)
    red, yellow = seperate_pm(FBstate)
    pylab.scatter(hyperbolicLattice.SDx[yellow], hyperbolicLattice.SDy[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(hyperbolicLattice.SDx[red], hyperbolicLattice.SDy[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    ax.set_aspect('equal')
    ax.axis('off')
    
    
    
    #1) octagon kagome
#    ax = pylab.subplot(1,4,2)
    ax = pylab.subplot(gs[2])
    
    testHyperbolic = PlanarLayout(gon = 8, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
    testHyperbolic.populate(1, resonatorsOnly=False)
    resonators = testHyperbolic.get_all_resonators()
    nameStr = str(testHyperbolic.gon) + 'gon_' + str(testHyperbolic.vertex) + 'vertex_' + str(testHyperbolic.itter) + '_' + testHyperbolic.modeType
    hyperbolicLattice = GeneralLayout(resonators , modeType = testHyperbolic.modeType, name =  nameStr)
    
    hyperbolicLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax)
    pylab.scatter(hyperbolicLattice.SDx, hyperbolicLattice.SDy ,c =  FWsiteColor, s = 20, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
#    FBstate =  build_double_plaquette_state(testHyperbolic)
#    red, yellow = seperate_pm(FBstate)
    red = [0,2,4,6]
    yellow = [1,3,5,7]
    pylab.scatter(hyperbolicLattice.SDx[yellow], hyperbolicLattice.SDy[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(hyperbolicLattice.SDx[red], hyperbolicLattice.SDy[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    ax.set_aspect('equal')
    ax.axis('off')
    
    
    
        #1) nonagon kagome
#    ax = pylab.subplot(1,4,2)
    ax = pylab.subplot(gs[3])
    
    testHyperbolic = PlanarLayout(gon = 9, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
    testHyperbolic.populate(1, resonatorsOnly=False)
    resonators = testHyperbolic.get_all_resonators()
    nameStr = str(testHyperbolic.gon) + 'gon_' + str(testHyperbolic.vertex) + 'vertex_' + str(testHyperbolic.itter) + '_' + testHyperbolic.modeType
    hyperbolicLattice = GeneralLayout(resonators , modeType = testHyperbolic.modeType, name =  nameStr)
    
    hyperbolicLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
    pylab.sca(ax)
    pylab.scatter(hyperbolicLattice.SDx, hyperbolicLattice.SDy ,c =  FWsiteColor, s = 20, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
    FBstate =  build_double_plaquette_state(testHyperbolic)
    red, yellow = seperate_pm(FBstate)
#    hyperbolicLattice.plot_layout_state(FBstate, ax)
    pylab.scatter(hyperbolicLattice.SDx[yellow], hyperbolicLattice.SDy[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(hyperbolicLattice.SDx[red], hyperbolicLattice.SDy[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    ax.set_aspect('equal')
    ax.axis('off')
    
    
    
    
    
    fig74.set_size_inches([11.2, 3.2])
    pylab.tight_layout()
    
#    fig74.savefig('ManyFBstates_1.png', transparent= False, dpi = 200)
#    fig74.savefig('ManyFBstates_2.png', transparent= False, dpi = 400)
#    fig74.savefig('ManyFBstates_2.svg', transparent= True)
    
else:
    fig74 = pylab.figure(74)
    pylab.clf()
    pylab.close(fig74)
    
    
    
    
    
    
    
    
#%%
    
############
#Part of kaomge topology fig with unit cell
########################################
if KagomeUnitCellFig:

    fig55 = pylab.figure(55)
    pylab.clf()
    #1) kagome
#    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    testEuclidLattice = EuclideanLayout(4,4,lattice_type = 'kagome', modeType = 'FW')
    kagomeLattice = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    pylab.cla()
    red = [3,5,13,18,20,21,24,32,33,35]
    yellow = [4,8,9,19,22,25,31,36,34,37]
    
#    ax5 = pylab.subplot(1,2,1, adjustable='box', aspect=1)
#    kagomeLattice.draw_SDlinks(ax5, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=FWlinkAlpha)
#    pylab.sca(ax5)
#    pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=FWsiteAlpha)
#    pylab.scatter(kagomeLattice.SDx[yellow], kagomeLattice.SDy[yellow], c =  stateColor1, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
#    pylab.scatter(kagomeLattice.SDx[red], kagomeLattice.SDy[red], c =  stateColor2, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
#    #ax5.set_aspect('equal')
#    ax5.axis('off')
    
    
#    keepers = [12,13, 14, 9, 8]
    keepers = [12,8, 18, 19, 20, 17,16,15, 11, 23, 22, 25, 24, 27, 28, 29, 30, 31, 32, 40, 39, 41, 36, 37]
    
    ax = pylab.subplot(1,2,2, adjustable='box', aspect=1)
    pylab.cla()
    kagomeLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
    #testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
    xs = kagomeLattice.coords[:,0]
    ys = kagomeLattice.coords[:,1]
    pylab.sca(ax)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault, marker = 'o', edgecolors = 'k', zorder = 5)
    pylab.scatter(kagomeLattice.SDx[keepers], kagomeLattice.SDy[keepers], c =  stateColor2, s = bigCdefault, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    ax.set_aspect('equal')
    ax.axis('off')
    
    
    
    ax = pylab.subplot(1,2,1, adjustable='box', aspect=1)
    truncationRes = kagomeLattice.resonators[keepers, :]
    truncation = GeneralLayout(truncationRes , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    
    truncation.draw_SDlinks(ax, color = FWlinkColor, linewidth = 3, minus_links = True, minus_color = 'gold', alpha=FWlinkAlpha)
    pylab.sca(ax)
    pylab.scatter(truncation.SDx, truncation.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=FWsiteAlpha)
    ax.set_aspect('equal')
    ax.axis('off')
    

    pylab.tight_layout()
    pylab.show()
#    fig55.savefig('kagomeUnitCellFig.png',transparent= False, dpi = 400)
#    fig55.savefig('kagomeUnitCellFig.svg',transparent= True)
else:
    fig55 = pylab.figure(55)
    pylab.clf()
    pylab.suptitle('KagomeUnitCellFig')
    pylab.close(fig55)

    

 
    
    
    
    
#%%
########
#figure infinity: split graphene, split hpg, split thingy
#######
def plot_any_band_cut(ax, cut):
    
    colorlist = ['firebrick', 'dodgerblue', 'blueviolet', 'mediumblue', 'goldenrod', 'cornflowerblue']
    
    pylab.sca(ax)
    
    shape = cut.shape
    
    for ind in range(0,shape[0]):
        colorInd = numpy.mod(ind, len(colorlist))
        pylab.plot(cut[ind,:], color = colorlist[colorInd] , marker = '.', markersize = '5', linestyle = '')
        
    pylab.title('some momentum cut')
    pylab.ylabel('Energy')
    pylab.xlabel('k_something')
        
    return
    



if SplitGraphFig:
    cellModeType = 'FW'
    #cellModeType = 'HW'
    
    McLaughlinColors = True
#    McLaughlinColors= False
    
    allTogether = True
#    allTogether = False
    
    forceLims = True
#    forceLims = False
    if forceLims:
        axmin = -2.75
        axmax = 3.25
    
    
    bigC = 70
    smallC = 25
    smallC2 = 20
    
    splitCenterColor = n1_splitCenterColor
    splitEdgeColor = n1_splitEdgeColor
    splitCenterColor2 = n1_splitCenterColor2
    splitEdgeColor2 = n1_splitEdgeColor2
    LGsplitCenterColor = n1_LGsplitCenterColor
    splitStateColor1 = n1_splitStateColor1
    splitStateEdgeColor1 = n1_splitStateEdgeColor1
    
 
#    gs = gridspec.GridSpec(2, 6,
#                           width_ratios=[1, 1, 1, 1,1,1],
#                           height_ratios=[1, 1]
#                           )
#    gs = gridspec.GridSpec(3, 3,
#                           width_ratios=[1, 1, 1],
#                           height_ratios=[1, 1,1]
#                           )
#    gs = gridspec.GridSpec(2, 3,
#                           width_ratios=[ 1,1, 1],
#                           height_ratios=[1, 1]
#                           )
#    gs2 = gridspec.GridSpec(1, 3,
#                           width_ratios=[ 1,1, 1],
#                           height_ratios=[ 1]
#                           )
    
    gs = gridspec.GridSpec(2, 4,
                           width_ratios=[1, 1,1, 1],
                           height_ratios=[1, 1]
                           )
    gs2 = gridspec.GridSpec(1, 4,
                           width_ratios=[1, 1,1, 1],
                           height_ratios=[ 1]
                           )
    gs_t = gridspec.GridSpec(3, 4,
                           width_ratios=[1, 1,1, 1],
                           height_ratios=[1, 1, 1]
                           )
    
    
    
    fig48 = pylab.figure(48)
    pylab.clf()
    
    
    #1) split graphene
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = cellModeType)
    resonators = testEuclidLattice.resonators
    resonators = split_resonators(resonators)
    testLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'splitkagome')
    
    cell0 = UnitCell('kagome')
    res0 = cell0.resonators
    res1 = split_resonators(res0)
    testCell = UnitCell('split_kagome', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
    
    
    #latticeStretch = 0.7
    latticeStretch = 1.0
    
    
    
    
    #cheating way of getting the split graph band structure
    originalCell = testEuclidLattice.unitcell
    numSurfPoints = 300
    kx_x, ky_y, cutx = originalCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = originalCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
    layoutBandStructureX = cutx[1:,:] -1 #remove the flat bands and shift down
    layoutBandStructureY = cuty[1:,:] -1 #remove the flat bands and shift down
    
    stretchX1 = numpy.sqrt(layoutBandStructureX + 3)
    stretchY1 = numpy.sqrt(layoutBandStructureY + 3)
    
    FB = numpy.zeros((1, numSurfPoints))
    
    splitCutX = numpy.concatenate((-stretchX1,FB,  stretchX1))
    splitCutY = numpy.concatenate((-stretchY1, FB, stretchY1))
#    splitCutX = numpy.concatenate((-stretchX1, stretchX1))
#    splitCutY = numpy.concatenate((-stretchY1, stretchY1))
    
    if allTogether:
        ax = pylab.subplot(gs_t[0], adjustable='box', aspect=80)
    else:
        ax = pylab.subplot(gs[0], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
#    testCell.plot_band_cut(ax, splitCutY)
    plot_any_band_cut(ax, splitCutY)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, splitCutX.shape[1]/2, splitCutX.shape[1]], [-2.5,0,2.5], rotation='horizontal')
    if forceLims:
        ax.set_ylim([axmin, axmax])
    
    
    ####graphs
    if allTogether:
        ax = pylab.subplot(gs_t[1], adjustable='box', aspect=latticeStretch)
    else:
        ax = pylab.subplot(gs[1], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    if McLaughlinColors:
        #splitGraphPrecursor
        testLattice.draw_resonator_lattice(ax, color = 'b', alpha = 1 , linewidth = 2.5)
        xs = testLattice.coords[:,0]
        ys = testLattice.coords[:,1]
        pylab.sca(ax)
#        pylab.scatter(xs, ys ,c =  'gold', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor, s = smallC, marker = 'o', edgecolors = splitEdgeColor, zorder = 5)
        
        xs = testEuclidLattice.coords[:,0]
        ys = testEuclidLattice.coords[:,1]
#        pylab.scatter(xs, ys ,c = 'lightsteelblue', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor2, s = smallC, marker = 'o', edgecolors = splitEdgeColor2, zorder = 5)
#        ax.set_aspect('equal')
        ax.axis('off')
        
        
    else:
        testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
        xs = testLattice.coords[:,0]
        ys = testLattice.coords[:,1]
        pylab.sca(ax)
        pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
        ax.axis('off')
        
    
    if effectiveInstead:
        newgraph = generate_line_graph(testLattice.resonators)
        newLattice = GeneralLayout(newgraph , modeType = 'FW', name =  'kagome_McLaughlin')
        
        if allTogether:
            ax = pylab.subplot(gs_t[2], adjustable='box', aspect=latticeStretch)
        else:
            ax = pylab.subplot(gs[2], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        newLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha=1)
        pylab.sca(ax)
        pylab.scatter(newLattice.SDx, newLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
        ax.axis('off')
        
        #band structure
        numSurfPoints = 300
        kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
        kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)

        #hack the effective band structure
        effectiveBandStructureX = cutx[:,:] +1 #remove the flat bands and shift down
        effectiveBandStructureY = cuty[:,:] +1 #remove the flat bands and shift down
        FB = numpy.zeros((1, numSurfPoints))-2
        effectiveCutX = numpy.concatenate((FB, effectiveBandStructureX))
        effectiveCutY = numpy.concatenate((FB, effectiveBandStructureY))
        
        if allTogether:
            ax = pylab.subplot(gs_t[3], adjustable='box', aspect=80)
        else:
            ax = pylab.subplot(gs[3], adjustable='box', aspect=80)
        pylab.cla()
        #plot_band_cut_mc(ax, cuty)
        #testCell.plot_band_cut(ax, cuty)
        plot_any_band_cut(ax, effectiveCutY)
        pylab.title('')
        pylab.ylabel('Energy (|t|)')
        pylab.xlabel('$k_y$ ($\pi$/a)')
        pylab.xticks([0, effectiveCutY.shape[1]/2, effectiveCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        
    else:
        if allTogether:
            ax = pylab.subplot(gs_t[2], adjustable='box', aspect=latticeStretch)
        else:
            ax = pylab.subplot(gs[2], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        if McLaughlinColors:
            testLattice.draw_SDlinks(ax, color = layoutLineColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
            pylab.sca(ax)
            pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
        #    ax.set_aspect('equal')
            ax.axis('off')
        else:
            testLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
            pylab.sca(ax)
            pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
        #    ax.set_aspect('equal')
            ax.axis('off')
    
    
        #band structure
        numSurfPoints = 300
        kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
        kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
        
        if allTogether:
            ax = pylab.subplot(gs_t[3], adjustable='box', aspect=80)
        else:
            ax = pylab.subplot(gs[3], adjustable='box', aspect=80)
        pylab.cla()
        #plot_band_cut_mc(ax, cuty)
        testCell.plot_band_cut(ax, cuty)
        pylab.title('')
        pylab.ylabel('Energy (|t|)')
        pylab.xlabel('$k_y$ ($\pi$/a)')
        pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        if forceLims:
            ax.set_ylim([axmin, axmax])
    
    





    #2) split hpg
    testEuclidLattice = EuclideanLayout(2,1,lattice_type = 'Huse', modeType = cellModeType)
    resonators = testEuclidLattice.resonators
    resonators = split_resonators(resonators)
    testLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'splitHPG')
    
    cell0 = UnitCell('Huse')
    res0 = cell0.resonators
    res1 = split_resonators(res0)
    testCell = UnitCell('split_hpg', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
    
    #latticeStretch = 0.7
    latticeStretch = 1.0
    
    
    #cheating way of getting the split graph band structure
    originalCell = testEuclidLattice.unitcell
    numSurfPoints = 300
    kx_x, ky_y, cutx = originalCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = originalCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
    layoutBandStructureX = cutx[4:,:] -1 #remove the flat bands and shift down
    layoutBandStructureY = cuty[4:,:] -1 #remove the flat bands and shift down
    
    stretchX1 = numpy.sqrt(layoutBandStructureX + 3)
    stretchY1 = numpy.sqrt(layoutBandStructureY + 3)
    
    FB = numpy.zeros((1, numSurfPoints))
    
    splitCutX = numpy.concatenate((-stretchX1,FB,  stretchX1))
    splitCutY = numpy.concatenate((-stretchY1, FB, stretchY1))
#    splitCutX = numpy.concatenate((-stretchX1, stretchX1))
#    splitCutY = numpy.concatenate((-stretchY1, stretchY1))
    
    if allTogether:
        ax = pylab.subplot(gs_t[4], adjustable='box', aspect=80)
    else:
        ax = pylab.subplot(gs[4], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
#    testCell.plot_band_cut(ax, splitCutY)
    plot_any_band_cut(ax, splitCutY)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, splitCutX.shape[1]/2, splitCutX.shape[1]], [-2.5,0,2.5], rotation='horizontal')
    if forceLims:
        ax.set_ylim([axmin, axmax])
    
    
    
    ####graphs
    if allTogether:
        ax = pylab.subplot(gs_t[5], adjustable='box', aspect=latticeStretch)
    else:
        ax = pylab.subplot(gs[5], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    if McLaughlinColors:
        #splitGraphPrecursor
        testLattice.draw_resonator_lattice(ax, color = 'b', alpha = 1 , linewidth = 2.5)
        xs = testLattice.coords[:,0]
        ys = testLattice.coords[:,1]
        pylab.sca(ax)
#        pylab.scatter(xs, ys ,c =  'gold', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor, s = smallC, marker = 'o', edgecolors = splitEdgeColor, zorder = 5)
        
        xs = testEuclidLattice.coords[:,0]
        ys = testEuclidLattice.coords[:,1]
#        pylab.scatter(xs, ys ,c = 'lightsteelblue', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor2, s = smallC, marker = 'o', edgecolors = splitEdgeColor2, zorder = 5)
#        ax.set_aspect('equal')
        ax.axis('off')
    else:
        testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
        xs = testLattice.coords[:,0]
        ys = testLattice.coords[:,1]
        pylab.sca(ax)
        pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
        ax.axis('off')
    
    
    if effectiveInstead:
        newgraph = generate_line_graph(testLattice.resonators)
        newLattice = GeneralLayout(newgraph , modeType = 'FW', name =  'hpk_McLaughlin')
        
        if allTogether:
            ax = pylab.subplot(gs_t[6], adjustable='box', aspect=latticeStretch)
        else:
            ax = pylab.subplot(gs[6], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        newLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha=1)
        pylab.sca(ax)
        pylab.scatter(newLattice.SDx, newLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
        ax.axis('off')
        
        #band structure
        numSurfPoints = 300
        kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
        kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)

        #hack the effective band structure
        effectiveBandStructureX = cutx[:,:] +1 #remove the flat bands and shift down
        effectiveBandStructureY = cuty[:,:] +1 #remove the flat bands and shift down
        FB = numpy.zeros((1, numSurfPoints))-2
        effectiveCutX = numpy.concatenate((FB, effectiveBandStructureX))
        effectiveCutY = numpy.concatenate((FB, effectiveBandStructureY))
        
        if allTogether:
            ax = pylab.subplot(gs_t[7], adjustable='box', aspect=80)
        else:
            ax = pylab.subplot(gs[7], adjustable='box', aspect=80)
        pylab.cla()
        #plot_band_cut_mc(ax, cuty)
        #testCell.plot_band_cut(ax, cuty)
        plot_any_band_cut(ax, effectiveCutY)
        pylab.title('')
        pylab.ylabel('Energy (|t|)')
        pylab.xlabel('$k_y$ ($\pi$/a)')
        pylab.xticks([0, effectiveCutY.shape[1]/2, effectiveCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        
    else:
        if allTogether:
            ax = pylab.subplot(gs_t[6], adjustable='box', aspect=latticeStretch)
        else:
            ax = pylab.subplot(gs[6], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        if McLaughlinColors:
            testLattice.draw_SDlinks(ax, color = layoutLineColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
            pylab.sca(ax)
            pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
            ax.axis('off')        
        else:
            testLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold')
            pylab.sca(ax)
            pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
            ax.axis('off')
    
    
        #band structure
        numSurfPoints = 300
        kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
        kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
        
        
        if allTogether:
            ax = pylab.subplot(gs_t[7], adjustable='box', aspect=80)
        else:
            ax = pylab.subplot(gs[7], adjustable='box', aspect=80)
        pylab.cla()
        #plot_band_cut_mc(ax, cuty)
        testCell.plot_band_cut(ax, cuty)
    #    testCell.plot_band_cut(ax, cutx)
        pylab.title('')
        pylab.ylabel('Energy (|t|)')
        pylab.xlabel('$k_y$ ($\pi$/a)')
    #    pylab.xlabel('$k_x$ ($\pi$/a)')
        pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')   
        if forceLims:
            ax.set_ylim([axmin, axmax])
    
    
    
    if allTogether:
        pass
    else:
        fig49 = pylab.figure(49)
        pylab.clf()
    
    #1) extremal hofman
    ########split Euclidean, hoffman attemps
    ##!!!!!!!!needs to be this size for the cell to come out right. DO NOT CHANGE!!!
    test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    resonators0 = test1.resonators #graphene layout
    splitGraph = split_resonators(resonators0) #split graphene layout
    resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
    resonators2 = split_resonators(resonators1) #split further
    resonators = resonators2
    testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
    #sites = [80,81,82, 83,79,78,77,76, 71, 84,85,86,87,114,115, 116, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
    sites = [80,81,83,79,78,77,76, 71, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
    state = numpy.zeros(len(testLattice.SDx))
    state[sites] = 0.1
    newCell = testLattice.resonators[sites, :]
#    
#    fig1 = pylab.figure(114)
#    pylab.clf()
#    ax = pylab.subplot(1,1,1)
#    testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
#    xs = testLattice.coords[:,0]
#    ys = testLattice.coords[:,1]
#    pylab.sca(ax)
#    pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)
#    pylab.scatter(testLattice.SDx[sites], testLattice.SDy[sites], c =  'firebrick', s = 75, marker = 'o', edgecolors = 'maroon', zorder = 5,linewidth = 1)
#    ax.set_aspect('equal')
#    ax.axis('off')
#    pylab.title(testLattice.name)
#    pylab.tight_layout()
#    pylab.show()
#    fig49 = pylab.figure(49)
    
    testCell = UnitCell('extremal_hofmann', resonators = newCell, a1 = test1.unitcell.a1, a2 = test1.unitcell.a2)
    
    
    #latticeStretch = 0.7
    latticeStretch = 1.0
    
    

    #cheating way of getting the split graph band structure
    parentLattice = GeneralLayout(resonators1 , modeType = test1.modeType, name =  'hofmannParent')
    
    #54, 55,56,57
#    sites = [50, 49, 48, 47, 46, 36,37, 51,52,53, 44,   67,64]
    sites = [50, 49, 48, 47, 46, 36,37,44,51]
    state = numpy.zeros(len(parentLattice.SDx))
    state[sites] = 0.1
    newCell = parentLattice.resonators[sites, :]
    
    parentCell = UnitCell('extremal_hofmann_parent', resonators = newCell, a1 = test1.unitcell.a1, a2 = test1.unitcell.a2)
#    fig1 = pylab.figure(114)
#    pylab.clf()
#    ax = pylab.subplot(1,1,1)
#    parentLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
#    xs = parentLattice.coords[:,0]
#    ys = parentLattice.coords[:,1]
#    pylab.sca(ax)
#    pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)
#    pylab.scatter(parentLattice.SDx[sites], parentLattice.SDy[sites], c =  'firebrick', s = 75, marker = 'o', edgecolors = 'maroon', zorder = 5,linewidth = 1)
#    ax.set_aspect('equal')
#    ax.axis('off')
#    pylab.title(parentLattice.name)
#    pylab.tight_layout()
#    pylab.show()
#    fig49 = pylab.figure(49)
    
    originalCell = parentCell
    numSurfPoints = 300
    kx_x, ky_y, cutx = originalCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = originalCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
    layoutBandStructureX = cutx[3:,:] -1 #remove the flat bands and shift down
    layoutBandStructureY = cuty[3:,:] -1 #remove the flat bands and shift down
    
    stretchX1 = numpy.sqrt(layoutBandStructureX + 3)
    stretchY1 = numpy.sqrt(layoutBandStructureY + 3)
    
    FB = numpy.zeros((1, numSurfPoints))
    
    splitCutX = numpy.concatenate((-stretchX1,FB,  stretchX1))
    splitCutY = numpy.concatenate((-stretchY1, FB, stretchY1))
#    splitCutX = numpy.concatenate((-stretchX1, stretchX1))
#    splitCutY = numpy.concatenate((-stretchY1, stretchY1))
    
    if allTogether:
        ax = pylab.subplot(gs_t[8+0], adjustable='box', aspect=80)
    else:
        ax = pylab.subplot(gs2[0], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
#    testCell.plot_band_cut(ax, splitCutY)
    plot_any_band_cut(ax, splitCutY)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, splitCutX.shape[1]/2, splitCutX.shape[1]], [-2.5,0,2.5], rotation='horizontal')
    if forceLims:
        ax.set_ylim([axmin, axmax])
    
    
    
    ####redo the lattice for plotting
    test1 = EuclideanLayout(5,5,lattice_type = 'kagome', modeType = 'FW')
    resonators0 = test1.resonators #graphene layout
    splitGraph = split_resonators(resonators0) #split graphene layout
    resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
    resonators2 = split_resonators(resonators1) #split further
    resonators = resonators2
    testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
    
    
    if allTogether:
        ax = pylab.subplot(gs_t[8+1], adjustable='box', aspect=latticeStretch)
    else:
        ax = pylab.subplot(gs2[1], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    if McLaughlinColors:
        #splitGraphPrecursor
        testLattice.draw_resonator_lattice(ax, color = 'b', alpha = 1 , linewidth = 2.5)
        xs = testLattice.coords[:,0]
        ys = testLattice.coords[:,1]
        pylab.sca(ax)
#        pylab.scatter(xs, ys ,c =  'gold', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor, s = smallC, marker = 'o', edgecolors = splitEdgeColor, zorder = 5)
        
        tempGen = GeneralLayout(resonators1 , modeType = test1.modeType, name =  'temp')
        xs = tempGen.coords[:,0]
        ys = tempGen.coords[:,1]
#        pylab.scatter(xs, ys ,c = 'lightsteelblue', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor2, s = smallC, marker = 'o', edgecolors = splitEdgeColor2, zorder = 5)
#        ax.set_aspect('equal')
        ax.axis('off')
    else:
        testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
        xs = testLattice.coords[:,0]
        ys = testLattice.coords[:,1]
        pylab.sca(ax)
        pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallC2, marker = 'o', edgecolors = 'k', zorder = 5)
        ax.axis('off')
    width = 1.9 +.4
    xStart = -.28
    yStart = 0.35 + 1
    ax.set_xlim([xStart, xStart + width])
    ax.set_ylim([yStart, yStart + width])
    
    if effectiveInstead:
        newgraph = generate_line_graph(testLattice.resonators)
        newLattice = GeneralLayout(newgraph , modeType = test1.modeType, name =  'hofmannlayout')
        
        if allTogether:
            ax = pylab.subplot(gs_t[8+2], adjustable='box', aspect=latticeStretch)
        else:
            ax = pylab.subplot(gs2[2], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        newLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha=1)
        pylab.sca(ax)
        pylab.scatter(newLattice.SDx, newLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
        ax.axis('off')
        width = 1.7 +0.045
        xStart = -0.005
        yStart = 1.625
        ax.set_xlim([xStart, xStart + width])
        ax.set_ylim([yStart, yStart + width])
        
        #band structure
        numSurfPoints = 300
        kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
        kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)

        #hack the effective band structure
        effectiveBandStructureX = cutx[:,:] +1 #remove the flat bands and shift down
        effectiveBandStructureY = cuty[:,:] +1 #remove the flat bands and shift down
        FB = numpy.zeros((1, numSurfPoints))-2
        effectiveCutX = numpy.concatenate((FB, effectiveBandStructureX))
        effectiveCutY = numpy.concatenate((FB, effectiveBandStructureY))
        
        if allTogether:
            ax = pylab.subplot(gs_t[8+3], adjustable='box', aspect=80)
        else:
            ax = pylab.subplot(gs2[3], adjustable='box', aspect=80)
        pylab.cla()
        #plot_band_cut_mc(ax, cuty)
        #testCell.plot_band_cut(ax, cuty)
        plot_any_band_cut(ax, effectiveCutY)
        pylab.title('')
        pylab.ylabel('Energy (|t|)')
        pylab.xlabel('$k_y$ ($\pi$/a)')
        pylab.xticks([0, effectiveCutY.shape[1]/2, effectiveCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        
    else:
        if allTogether:
            ax = pylab.subplot(gs_t[8+2], adjustable='box', aspect=latticeStretch)
        else:
            ax = pylab.subplot(gs2[2], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        if McLaughlinColors:
            testLattice.draw_SDlinks(ax, color = layoutLineColor, linewidth = 1.5, minus_links = True, minus_color = 'gold')
            pylab.sca(ax)
            pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  layoutCapColor, s = smallC2, marker = 'o', edgecolors = 'k', zorder = 5)
            ax.axis('off')
        else:
            testLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 1.5, minus_links = True, minus_color = 'gold')
            pylab.sca(ax)
            pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  FWsiteColor, s = smallC2, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
            #add transparent bubbles so fiugres are sized correctly.
        #    pylab.scatter(testLattice.SDx[yellow2], testLattice.SDy[yellow2], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, alpha=0)
        #    pylab.scatter(testLattice.SDx[red2], testLattice.SDy[red2], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, alpha=0)
        #    ax.set_aspect('equal')
            ax.axis('off')
        width = 1.7 +0.045
        xStart = -0.005
        yStart = 1.625
        ax.set_xlim([xStart, xStart + width])
        ax.set_ylim([yStart, yStart + width])
    
    
        #band structure
        numSurfPoints = 300
        kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
        kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
        
        
        if allTogether:
            ax = pylab.subplot(gs_t[8+3], adjustable='box', aspect=80)
        else:
            ax = pylab.subplot(gs2[3], adjustable='box', aspect=80)
        pylab.cla()
        #plot_band_cut_mc(ax, cuty)
        testCell.plot_band_cut(ax, cuty)
        pylab.title('')
        pylab.ylabel('Energy (|t|)')
        pylab.xlabel('$k_y$ ($\pi$/a)')
        pylab.xticks([0, cutx.shape[1]/2, cutx.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        if forceLims:
            ax.set_ylim([axmin, axmax])
    
    

    if allTogether:
        fig48.set_size_inches(9.0*4/3, 7.5*3/2)
        pylab.figure(48)
        pylab.tight_layout()
        pylab.show()
        
#    fig48.savefig('SplitGraphs.png',transparent= False, dpi = 200)
#    fig48.savefig('SplitGraphs.svg',transparent= True)
        
        fig49 = pylab.figure(49)    
        pylab.clf()
        pylab.suptitle('Split Graphs2')
        pylab.close(fig49)
    else:
        fig48.set_size_inches(9.0*4/3, 7.5)
        pylab.figure(48)
        pylab.tight_layout()
        pylab.show()
        
        fig49.set_size_inches(9.0*4/3, 7.5/2)
    #    ig49.set_size_inches(13.6, 6.75)
        pylab.figure(49)
        pylab.tight_layout()
        pylab.show()
    
    
#    fig48.savefig('SplitGraphs1.png',transparent= False, dpi = 200)
#    fig48.savefig('SplitGraphs1.svg',transparent= True)
    
#    fig49.savefig('SplitGraphs2.png',transparent= False, dpi = 200)
#    fig49.savefig('SplitGraphs2.svg',transparent= True)
    
    
else:
    fig48 = pylab.figure(49)    
    pylab.clf()
    pylab.suptitle('Split Graphs')
    pylab.close(fig48)
    
    fig49 = pylab.figure(49)    
    pylab.clf()
    pylab.suptitle('Split Graphs2')
    pylab.close(fig49)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#%%   
#######
#figure 1 million: split graphs, and effective lattices, localized states and spectra
#######
if SplitGraphLineGraphFig:
#    bigC = 90
    bigC = 60
    smallC = 20
    
    smallGap = 0.24
    bigGap = 0.52+ 0.24 -0.24
    
    Amp_imp = 0.3
    Amp_band = 0.1

    
    
    McLaughlinColors = True
#    McLaughlinColors = False

    splitCenterColor = n2_splitCenterColor
    splitEdgeColor = n2_splitEdgeColor
    splitCenterColor2 = n2_splitCenterColor2
    splitEdgeColor2 = n2_splitEdgeColor2
    LGsplitCenterColor = n2_LGsplitCenterColor
    splitStateColor1 = n2_splitStateColor1
    splitStateEdgeColor1 = n2_splitStateEdgeColor1
    
    
    gs = gridspec.GridSpec(4, 4,
                           width_ratios=[1, 1, 1, 1],
                           height_ratios=[1, 1, 0.1, 0.1]
                           )
    fig111 = pylab.figure(111)
    pylab.clf()
    #1) kagome
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    tempRes = split_resonators(testEuclidLattice.get_all_resonators())
    kagomeLattice = GeneralLayout(tempRes, modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    ax = pylab.subplot(gs[0], adjustable='box', aspect=1)
    pylab.cla()
    if McLaughlinColors:
        #splitGraphPrecursor
        kagomeLattice.draw_resonator_lattice(ax, color = 'b', alpha = FWlinkAlpha , linewidth = 2.5)
        xs = kagomeLattice.coords[:,0]
        ys = kagomeLattice.coords[:,1]
        pylab.sca(ax)
#        pylab.scatter(xs, ys ,c =  'gold', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor, s = smallC, marker = 'o', edgecolors = splitEdgeColor, zorder = 5)
        
        xs = testEuclidLattice.coords[:,0]
        ys = testEuclidLattice.coords[:,1]
#        pylab.scatter(xs, ys ,c = 'lightsteelblue', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c = splitCenterColor2, s = smallC, marker = 'o', edgecolors = splitEdgeColor2, zorder = 5)
#        ax.set_aspect('equal')
        ax.axis('off')
    else:
        kagomeLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
        xs = kagomeLattice.coords[:,0]
        ys = kagomeLattice.coords[:,1]
        pylab.sca(ax)
        pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
        #ax.set_aspect('equal')
        ax.axis('off')

    red = [24, 14,6]
    yellow = [15, 23, 5]
    shiftx, shifty = testEuclidLattice.unitcell.a1
    coordsX = kagomeLattice.coords[:,0]
    coordsY = kagomeLattice.coords[:,1]
#    pylab.scatter(coordsX[yellow]+shiftx, coordsY[yellow]+shifty, c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(coordsX[yellow]+shiftx, coordsY[yellow]+shifty, c =  splitStateColor1, s = bigC, marker = 'o', edgecolors = splitStateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(coordsX[red]+shiftx, coordsY[red]+shifty, c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2) 
        
    
    
    
    ax5 = pylab.subplot(gs[4+0], adjustable='box', aspect=1)
    if McLaughlinColors:
        kagomeLattice.draw_SDlinks(ax5, color = layoutLineColor, linewidth = 2, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
        pylab.sca(ax5)
#        pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
        pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  LGsplitCenterColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
    else:
        kagomeLattice.draw_SDlinks(ax5, color = FWlinkColor, linewidth = 2, minus_links = True, minus_color = 'gold')
        pylab.sca(ax5)
        pylab.scatter(kagomeLattice.SDx, kagomeLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
        
#    red = [3,5,13]
#    yellow = [4,8,9]
    red = [10,11, 6,7, 26, 27]
    yellow = [8,9, 16, 17, 18, 19]
    shiftx, shifty = testEuclidLattice.unitcell.a1
#    pylab.scatter(kagomeLattice.SDx[yellow]+shiftx, kagomeLattice.SDy[yellow]+shifty, c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(kagomeLattice.SDx[yellow]+shiftx, kagomeLattice.SDy[yellow]+shifty, c =  splitStateColor1, s = bigC, marker = 'o', edgecolors = splitStateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(kagomeLattice.SDx[red]+shiftx, kagomeLattice.SDy[red]+shifty, c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    #ax5.set_aspect('equal')
    ax5.axis('off')



    
    #2) square
    testEuclidLattice = EuclideanLayout(4,4,lattice_type = 'square', modeType = 'FW')
    #testEuclidLattice = EuclideanLayout(4,4,lattice_type = 'square', modeType = 'HW')
    res_list = testEuclidLattice.get_all_resonators()
    mask = numpy.logical_and((res_list[:,0]+res_list[:,2])>0 ,(res_list[:,1]+res_list[:,3])>0)
    tempRes = split_resonators(res_list[mask,:])
    tempLattice = GeneralLayout(res_list[mask,:] , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    squareLattice = GeneralLayout(tempRes , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    ax = pylab.subplot(gs[1])
    pylab.cla()
    if McLaughlinColors:
        #splitGraphPrecursor
        squareLattice.draw_resonator_lattice(ax, color = 'b', alpha = FWlinkAlpha , linewidth = 2.5)
        xs = squareLattice.coords[:,0]
        ys = squareLattice.coords[:,1]
        pylab.sca(ax)
#        pylab.scatter(xs, ys ,c =  'gold', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor, s = smallC, marker = 'o', edgecolors = splitEdgeColor, zorder = 5)
        
        xs = tempLattice.coords[:,0]
        ys = tempLattice.coords[:,1]
#        pylab.scatter(xs, ys ,c = 'lightsteelblue', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor2, s = smallC, marker = 'o', edgecolors = splitEdgeColor2, zorder = 5)
#        ax.set_aspect('equal')
        ax.axis('off')
    else:
        squareLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
        xs = squareLattice.coords[:,0]
        ys = squareLattice.coords[:,1]
        pylab.sca(ax)
        pylab.scatter(xs, ys ,c = layoutCapColor, s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
        ax.set_aspect('equal')
        ax.axis('off')


    red = [11, 23]
    yellow = [16, 17]
    shiftx, shifty = 0,0
    coordsX = squareLattice.coords[:,0]
    coordsY = squareLattice.coords[:,1]
#    pylab.scatter(coordsX[yellow]+shiftx, coordsY[yellow]+shifty, c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(coordsX[yellow]+shiftx, coordsY[yellow]+shifty, c =  splitStateColor1, s = bigC, marker = 'o', edgecolors = splitStateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(coordsX[red]+shiftx, coordsY[red]+shifty, c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2) 
    
    
    
    
    ax6 = pylab.subplot(gs[4+1])
    if McLaughlinColors:
        squareLattice.draw_SDlinks(ax6, color = layoutLineColor, linewidth = 2, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
        pylab.sca(ax6)
#        pylab.scatter(squareLattice.SDx, squareLattice.SDy ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
        pylab.scatter(squareLattice.SDx, squareLattice.SDy ,c =  LGsplitCenterColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
    else:
        squareLattice.draw_SDlinks(ax6, color = FWlinkColor, linewidth = 2, minus_links = True, minus_color = 'gold')
        pylab.sca(ax6)
        pylab.scatter(squareLattice.SDx, squareLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
#    red = [3,5,13]
#    yellow = [4,8,9]
    red = [13,12, 26, 27]
    yellow = [15, 18, 19, 14]
    shiftx = testEuclidLattice.unitcell.a1[0]
    shifty = testEuclidLattice.unitcell.a2[1]
#    pylab.scatter(squareLattice.SDx[yellow], squareLattice.SDy[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(squareLattice.SDx[yellow], squareLattice.SDy[yellow], c =  splitStateColor1, s = bigC, marker = 'o', edgecolors = splitStateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(squareLattice.SDx[red], squareLattice.SDy[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    ax6.set_aspect('equal')
    ax6.axis('off')
 
    
    
    
    #3) tree
    testTree = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
    #testTree = TreeResonators(degree = 3, iterations = 5, side = 1, file_path = '', modeType = 'HW')
    resonators = testTree.get_all_resonators()
    tempRes = split_resonators(resonators)
    treeLattice = GeneralLayout(tempRes , modeType = testTree.modeType, name =  'TREEEEE')
    ax3 = pylab.subplot(gs[2])
    pylab.cla()
    if McLaughlinColors:
        #splitGraphPrecursor
        treeLattice.draw_resonator_lattice(ax3, color = 'b', alpha = FWlinkAlpha , linewidth = 2.5)
        xs = treeLattice.coords[:,0]
        ys = treeLattice.coords[:,1]
        pylab.sca(ax3)
#        pylab.scatter(xs, ys ,c =  'gold', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor, s = smallC, marker = 'o', edgecolors = splitEdgeColor, zorder = 5)
        
        xs = testTree.coords[:,0]
        ys = testTree.coords[:,1]
#        pylab.scatter(xs, ys ,c = 'lightsteelblue', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor2, s = smallC, marker = 'o', edgecolors = splitEdgeColor2, zorder = 5)
#        ax3.set_aspect('equal')
        ax3.axis('off')
    else:
        treeLattice.draw_resonator_lattice(ax3, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
        xs = treeLattice.coords[:,0]
        ys = treeLattice.coords[:,1]
        pylab.sca(ax3)
        pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
        ax3.set_aspect('equal')
        ax3.axis('off')
        #pylab.title('layout graph')
    
    
    ax7 = pylab.subplot(gs[4+2])
    if McLaughlinColors:
        treeLattice.draw_SDlinks(ax7, color = layoutLineColor, linewidth = 1.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
        pylab.sca(ax7)
#        pylab.scatter(treeLattice.SDx, treeLattice.SDy ,c =  layoutCapColor, s = 15, marker = 'o', edgecolors = 'k', zorder = 5)
        pylab.scatter(treeLattice.SDx, treeLattice.SDy ,c =  LGsplitCenterColor, s = 15, marker = 'o', edgecolors = 'k', zorder = 5)
    else:
        treeLattice.draw_SDlinks(ax7, color = FWlinkColor, linewidth = 1.5, minus_links = True, minus_color = 'gold')
        pylab.sca(ax7)
        pylab.scatter(treeLattice.SDx, treeLattice.SDy ,c =  FWsiteColor, s = 15, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
#    state_vect = treeLattice.Psis[:,44]
#    stepSize = 0.11/2
#    tempVect = state_vect/stepSize
#    maxVal = numpy.max(state_vect)
#    for ind in range(0, len(state_vect)):
#        rawVal = state_vect[ind]
#        if numpy.abs(rawVal/stepSize)>8:
#            #this is the peak
#            pass
#        elif numpy.abs(rawVal/stepSize)>4:
#            #second ring
#            state_vect[ind] = maxVal*numpy.sign(rawVal)/2
#        elif numpy.abs(rawVal/stepSize)>2:
#            #third ring
#            state_vect[ind] = maxVal*numpy.sign(rawVal)/4
#        elif numpy.abs(rawVal/stepSize)>1:
#            #fourth ring
#            state_vect[ind] = maxVal*numpy.sign(rawVal)/8
#        else:
#            state_vect[ind] = 0


    #make a flat band state
    state_vect = numpy.zeros(len(treeLattice.SDx))
    state_vect[0] = 1
    state_vect[1] = 1
    
    state_vect[2] = -1
    state_vect[3] = -1
    
    state_vect[6] = -0.5
    state_vect[7] = -0.5
    
    state_vect[8] = -0.5
    state_vect[9] = -0.5
    
    state_vect[10] = 0.5
    state_vect[11] = 0.5
    
    state_vect[12] = 0.5
    state_vect[13] = 0.5
    
    state_vect[18] = 0.25
    state_vect[19] = 0.25
    
    state_vect[20] = 0.25
    state_vect[21] = 0.25
    
    state_vect[22] = 0.25
    state_vect[23] = 0.25
    
    state_vect[24] = 0.25
    state_vect[25] = 0.25
    
    state_vect[26] = -0.25
    state_vect[27] = -0.25
    
    state_vect[28] = -0.25
    state_vect[29] = -0.25
    
    state_vect[30] = -0.25
    state_vect[31] = -0.25
    
    state_vect[32] = -0.25
    state_vect[33] = -0.25
    
    state_vect[42:58] = 0.5*numpy.concatenate((state_vect[26:34], state_vect[26:34]))

    state_vect[58:74] = - state_vect[42:58]
    
    state_vect = state_vect*0.45
            
    reds = numpy.where(state_vect <0)[0]
    yellows = numpy.where(state_vect > 0)[0]
    
    scaleFactor = 8
#    scaleFactor = 64
    
    redState = numpy.zeros(len(state_vect))
    redState[reds] = state_vect[reds]
    redAmps = redState
    redProbs = numpy.abs(redAmps)**2
    redSizes = redProbs * len(redProbs)*scaleFactor
    
    yellowState = numpy.zeros(len(state_vect))
    yellowState[yellows] = state_vect[yellows]
    yellowAmps = yellowState
    yellowProbs = numpy.abs(yellowAmps)**2
    yellowSizes = yellowProbs * len(yellowProbs)*scaleFactor
    
    pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  stateColor2, s = redSizes, marker = 'o', edgecolors = stateEdgeColor2,  zorder = 6, linewidth = stateEdgeWidth2)
#    pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  stateColor1, s = yellowSizes, marker = 'o', edgecolors = stateEdgeColor1,  zorder = 6, linewidth = stateEdgeWidth1)
    pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  splitStateColor1, s = yellowSizes, marker = 'o', edgecolors = splitStateEdgeColor1,  zorder = 6, linewidth = stateEdgeWidth1)
#    pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  stateColor2, s = redSizes, marker = 'o', edgecolors = stateEdgeColor2,  zorder = 6, linewidth = stateEdgeWidth2)
#    pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  'forestgreen', s = yellowSizes, marker = 'o', edgecolors = 'darkolivegreen',  zorder = 6, linewidth = stateEdgeWidth1)

    
    
    ax7.set_aspect('equal')
    ax7.axis('off')




    ###go back to the tree and build the flat band state, relies on the original being ina  sensible order
    pylab.sca(ax3)
    ends = numpy.where(state_vect !=0)[0]
    xcoords = numpy.zeros(len(ends)/2)
    ycoords = numpy.zeros(len(ends)/2)
    treeReds = numpy.zeros(len(xcoords))
    treeYellows = numpy.zeros(len(ycoords))
    
    for ind in range(0, len(ends)/2):
        site = ends[2*ind]
        site2 = ends[2*ind+ 1]
        x1 = treeLattice.SDx[site]
        x2 = treeLattice.SDx[site2]
        
        y1 = treeLattice.SDy[site]
        y2 = treeLattice.SDy[site2]
        
        xcoords[ind] = (x1 + x2)/2
        ycoords[ind] = (y1 + y2)/2
        
        if state_vect[site] > 0:
            #this is a yellow
            treeYellows[ind] = yellowSizes[site]
        else:
            treeReds[ind] = redSizes[site]
    pylab.scatter(xcoords, ycoords, c =  stateColor2, s = treeReds, marker = 'o', edgecolors = stateEdgeColor2,  zorder = 6, linewidth = stateEdgeWidth2)
#    pylab.scatter(xcoords, ycoords, c =  stateColor1, s = treeYellows, marker = 'o', edgecolors = stateEdgeColor1,  zorder = 6, linewidth = stateEdgeWidth1)
    pylab.scatter(xcoords, ycoords, c =  splitStateColor1, s = treeYellows, marker = 'o', edgecolors = splitStateEdgeColor1,  zorder = 6, linewidth = stateEdgeWidth1)
    ax3.set_aspect('equal')
    ax3.axis('off')


    
    #4) hyperbolic
    testHyperbolic = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
    testHyperbolic.populate(2, resonatorsOnly=False)
#    testHyperbolic.populate(1, resonatorsOnly=False)
    resonators = testHyperbolic.get_all_resonators()
    tempRes = split_resonators(resonators)
    nameStr = str(testHyperbolic.gon) + 'gon_' + str(testHyperbolic.vertex) + 'vertex_' + str(testHyperbolic.itter) + '_' + testHyperbolic.modeType
    hyperbolicLattice = GeneralLayout(tempRes , modeType = testHyperbolic.modeType, name =  nameStr)
    
    
    
    ax4 = pylab.subplot(gs[3])
    pylab.cla()
    if McLaughlinColors:
        #splitGraphPrecursor
        hyperbolicLattice.draw_resonator_lattice(ax4, color = 'b', alpha = FWlinkAlpha , linewidth = 2.5)
        xs = hyperbolicLattice.coords[:,0]
        ys = hyperbolicLattice.coords[:,1]
        pylab.sca(ax4)
#        pylab.scatter(xs, ys ,c =  'gold', s = 20, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor, s = smallC, marker = 'o', edgecolors = splitEdgeColor, zorder = 5)
        
        xs = testHyperbolic.coords[:,0]
        ys = testHyperbolic.coords[:,1]
#        pylab.scatter(xs, ys ,c = 'lightsteelblue', s = 20, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
        pylab.scatter(xs, ys ,c =  splitCenterColor2, s = smallC, marker = 'o', edgecolors = splitEdgeColor2, zorder = 5)
#        ax4.set_aspect('equal')
        ax4.axis('off')
    else:
        hyperbolicLattice.draw_resonator_lattice(ax4, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
        xs = hyperbolicLattice.coords[:,0]
        ys = hyperbolicLattice.coords[:,1]
        pylab.sca(ax4)
        pylab.scatter(xs, ys ,c =  layoutCapColor, s = 15, marker = 'o', edgecolors = 'k', zorder = 5)
        ax4.set_aspect('equal')
        ax4.axis('off')
        #pylab.title('layout graph')
     
    
    red = [150]
    yellow = [180]
    shiftx, shifty = 0,0
    coordsX = hyperbolicLattice.coords[:,0]
    coordsY = hyperbolicLattice.coords[:,1]
    abss = numpy.sqrt(coordsX**2 + coordsY**2)
    minAbs = numpy.min(abss)
    core = numpy.where(abss < 1.2 *minAbs)[0]
    red = [core[1], core[8], core[7], 178, 181, 198]
    yellow = [core[3], core[4], core[11], 161, 188,    204]
#    pylab.scatter(coordsX[yellow]+shiftx, coordsY[yellow]+shifty, c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(coordsX[yellow]+shiftx, coordsY[yellow]+shifty, c =  splitStateColor1, s = bigC, marker = 'o', edgecolors = splitStateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(coordsX[red]+shiftx, coordsY[red]+shifty, c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2) 
    
    
    
    ax8 = pylab.subplot(gs[4+3])
    if McLaughlinColors:
        hyperbolicLattice.draw_SDlinks(ax8, color = layoutLineColor, linewidth = 2.0, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
        pylab.sca(ax8)
#        pylab.scatter(hyperbolicLattice.SDx, hyperbolicLattice.SDy ,c =  layoutCapColor, s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
        pylab.scatter(hyperbolicLattice.SDx, hyperbolicLattice.SDy ,c =  LGsplitCenterColor, s = 20, marker = 'o', edgecolors = 'k', zorder = 5)
    else:
        hyperbolicLattice.draw_SDlinks(ax8, color = FWlinkColor, linewidth = 2.0, minus_links = True, minus_color = 'gold')
        pylab.sca(ax8)
        pylab.scatter(hyperbolicLattice.SDx, hyperbolicLattice.SDy ,c =  FWsiteColor, s = 20, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
    
#    red = [2,4,6,7,9,36]
#    yellow = [1,3,5,8,10,35]
    red = [2,3, 6, 7, 10, 11, 16, 17, 20, 21, 70,71]
    yellow = [4,5, 8,9, 12,13,14, 15, 18, 19, 72,73]
#    pylab.scatter(hyperbolicLattice.SDx[yellow], hyperbolicLattice.SDy[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(hyperbolicLattice.SDx[yellow], hyperbolicLattice.SDy[yellow], c =  splitStateColor1, s = bigC, marker = 'o', edgecolors = splitStateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(hyperbolicLattice.SDx[red], hyperbolicLattice.SDy[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5,linewidth = stateEdgeWidth2)
    ax8.set_aspect('equal')
    ax8.axis('off')
    
    ax8.set_aspect('equal')
    ax8.axis('off')
    #pylab.title('effective (line) graph')
    
    #pylab.suptitle(hyperbolicLattice.name)
    #fig1.set_size_inches(15, 8)
    #pylab.tight_layout()
    #pylab.show()
    #fig1.savefig('Figure1.pdf',transparent= True)
    
    
    
    
    
    
    ####   sample DOS plot
    #set up frequency sweep
    ## 1) Kagome 
    ax1 = pylab.subplot(gs[8+0], adjustable='box', aspect=6)
    
    freq_start = -3.00
    freq_stop = 3.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [0]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2.4,0)] # start and end points for normal bands
    bandX2 = [(0,2.4)]
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    yband2 = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    for vals in bandX2:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
#    pylab.bar(freqs+freq_res/8., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
       
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-3,0,3], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,4.5])

    
    ## 2) Square     
    ax1 = pylab.subplot(gs[8+1], adjustable='box', aspect=6)
    
    freq_start = -4.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [0]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2.68,0)] # start and end points for normal bands
    bandX2 = [(0,2.68)]
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    yband2 = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    for vals in bandX2:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
#    pylab.bar(freqs+freq_res/8., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
       
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-4,-2,0,2,4], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,4.5])
    
    
    ## 3) Tree 
    ax1 = pylab.subplot(gs[8+2], adjustable='box', aspect=6)

    freq_start = -3.00
    freq_stop = 3.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [0]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2.4 + smallGap,-smallGap)] # start and end points for normal bands
    bandX2 = [(+smallGap,2.4 - smallGap)]
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    yband2 = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    for vals in bandX2:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
#    pylab.bar(freqs+freq_res/8., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
        
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-3,0,3], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,4.5])
    
    
    
    ## 4) Hyperbolics    
    ax1 = pylab.subplot(gs[8+3], adjustable='box', aspect=6)

    freq_start = -3.00
    freq_stop = 3.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [0]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2.4 + smallGap, -bigGap)] # start and end points for normal bands
    bandX2 = [(bigGap,2.4 - smallGap)]
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    yband2 = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    for vals in bandX2:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
#    pylab.bar(freqs+freq_res/8., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-3,0,3], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,4.5])
    
    ####### 5) Kagome Line
    ax1 = pylab.subplot(gs[12+0], adjustable='box', aspect=6)
    
    freq_start = -2.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [-2,0]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2,0)] # start and end points for normal bands
    bandX2 = [(1,3)]
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    yband2 = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    for vals in bandX2:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
#    pylab.bar(freqs+freq_res/8., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,0,3], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,4.5])
    
    
    
    ## 6) Square Line   
    ax1 = pylab.subplot(gs[12+1], adjustable='box', aspect=6)
    
    freq_start = -2.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [-2,0]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2,0)] # start and end points for normal bands
    bandX2 = [(2,4)]
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    yband2 = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    for vals in bandX2:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
#    pylab.bar(freqs+freq_res/8., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,0,2,4], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,4.5])
    
    
    ## 7) Tree Line    
    ax1 = pylab.subplot(gs[12+2], adjustable='box', aspect=6)
    
    freq_start = -2.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [-2,0]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2+smallGap,0 - smallGap)] # start and end points for normal bands
    bandX2 = [(1+smallGap,3-smallGap)]
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    yband2 = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    for vals in bandX2:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
#    pylab.bar(freqs+freq_res/8., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,0,1,3], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,4.5])
    
    
    ## 8) Hyperbolic Line     
    ax1 = pylab.subplot(gs[12+3], adjustable='box', aspect=6)
    
    freq_start = -2.00
    freq_stop = 4.00
    freq_res = 0.04
    
    freqs = numpy.arange(freq_start, freq_stop+freq_res, freq_res).round(2)
    
    impX = [-2,0]#[-2., 0.] # energies of flat bands
    N = len(freqs)
    bandX = [(-2+smallGap,0 - bigGap)] # start and end points for normal bands
    bandX2 = [(1+bigGap,3-smallGap)]
    yimp = numpy.zeros(N)
    yband = numpy.zeros(N)
    yband2 = numpy.zeros(N)
    
    
    for val in impX:
        ind = numpy.where(freqs == val)
        yimp += Amp_imp*signal.unit_impulse(N,ind)
    #
    for vals in bandX:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    for vals in bandX2:
        indstart = numpy.where(freqs == vals[0])[0][0]
        indstop = numpy.where(freqs == vals[1])[0][0]
        yband[indstart :indstop] = Amp_band
    
#    pylab.bar(freqs+freq_res/8., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
#    pylab.bar(freqs+freq_res/8., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
#    pylab.bar(freqs+freq_res/8., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yimp, width = freq_res*2, color = 'mediumblue', label = 'flatband', alpha = 1, linewidth = 0)
    pylab.bar(freqs+freq_res/2., 1.*yband, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)  
    pylab.bar(freqs+freq_res/2., 1.*yband2, width = freq_res, color = 'deepskyblue', label = 'garden variety band', alpha = 0.6, linewidth = 0)
    
    
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_visible(False)
    pylab.xticks([-2,0,1,3], rotation='horizontal')
    ax1.set_ylim([0,0.3])
    ax1.set_xlim([-4.5,4.5])
    
    
    
    
    
    fig111.set_size_inches(11.5, 7.2)
    pylab.tight_layout()
    pylab.show()
    
#    fig111.savefig('SplitGraphsAndLineGraphs.png',transparent= False, dpi = 200)
#    fig111.savefig('SplitGraphsAndLineGraphs.png',transparent= False, dpi = 400)
#    fig111.savefig('SplitGraphsAndLineGraphs.svg',transparent= True)
else:
    fig111 = pylab.figure(111)
    pylab.clf()
    pylab.close(fig111)   
    
    
    
    
#%%
# 
    
    
###extra plots for talks and stuff
if extraPlots:
    plotStates = True
#    plotStates = False
    
    fig999 = pylab.figure(999)
    pylab.clf()
    
    bigC = bigCdefault*3
    smallC = smallCdefault*2.5
    linewidth = 4.5
    
    ax = pylab.subplot(1,2,1)
    pylab.cla()
    testHyperbolic = PlanarLayout(gon = 6, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
    testHyperbolic.populate(1, resonatorsOnly=False)
    resonators = testHyperbolic.get_all_resonators()
    resonators = resonators[0:6,:]
    nameStr = 'heptagon_ring'
    ring = GeneralLayout(resonators , modeType = testHyperbolic.modeType, name =  nameStr)
    
    
    

    ring.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = linewidth)
    #testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
    xs = ring.coords[:,0]
    ys = ring.coords[:,1]
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)

    red = [0,3,4]
    yellow = [2,1,5]
    if plotStates:
        pylab.scatter(xs[yellow], ys[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
        pylab.scatter(xs[red], ys[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    ax.set_aspect('equal')
    ax.axis('off')
    
    
    ax = pylab.subplot(1,2,2)
    pylab.cla()
    testHyperbolic = PlanarLayout(gon = 7, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
    testHyperbolic.populate(1, resonatorsOnly=False)
    resonators = testHyperbolic.get_all_resonators()
    resonators = resonators[0:7,:]
    nameStr = 'heptagon_ring'
    ring = GeneralLayout(resonators , modeType = testHyperbolic.modeType, name =  nameStr)
    
    ring.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = linewidth)
    #testLattice.draw_resonator_end_points(ax, color = 'goldenrod', edgecolor = 'k',  marker = 'o' , size = 30)
    xs = ring.coords[:,0]
    ys = ring.coords[:,1]
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)

    red = [0,3,4]
    yellow = [2,1,5]
    if plotStates:
        pylab.scatter(xs[yellow], ys[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
        pylab.scatter(xs[red], ys[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    ax.set_aspect('equal')
    ax.axis('off')
    
    fig999.set_size_inches(4.1, 2)
    pylab.tight_layout()
    pylab.show()
    
    #fig999.savefig('Cycles.png',transparent= False, dpi = 200)
    
    
    
    
    
    
    
    fig998 = pylab.figure(998)
    pylab.clf()
    
    bigC = bigCdefault*1.5
    
    #1) kagome
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    kagomeLattice = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
    ax = pylab.subplot(1,2,1)
    pylab.cla()
    kagomeLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 3.5)
    xs = kagomeLattice.coords[:,0]
    ys = kagomeLattice.coords[:,1]
    pylab.sca(ax)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault*2, marker = 'o', edgecolors = 'k', zorder = 5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    
    ax = pylab.subplot(1,2,2)
    pylab.cla()
    kagomeLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 3.5)
    xs = kagomeLattice.coords[:,0]
    ys = kagomeLattice.coords[:,1]
    pylab.sca(ax)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallCdefault*2, marker = 'o', edgecolors = 'k', zorder = 5)
    
    red = [0,3,4]
    yellow = [2,1,5]
    inds = scipy.arange(0,len(xs), 1)
    red = [0,1,2,7,8,9,10,15, 16, 17, 18,23,24,25,26]
    yellow = [3,4,5, 6,11,12,13,14,19,20,21,22,27,28,29]
    pylab.scatter(xs[yellow], ys[yellow], c =  stateColor1, s = bigC, marker = 'o', edgecolors = stateEdgeColor1, zorder = 5, linewidth = stateEdgeWidth1)
    pylab.scatter(xs[red], ys[red], c =  stateColor2, s = bigC, marker = 'o', edgecolors = stateEdgeColor2, zorder = 5, linewidth = stateEdgeWidth2)
    ax.set_aspect('equal')
    ax.axis('off')
    
    fig998.set_size_inches(6.4, 4.8)
    pylab.tight_layout()
    pylab.show()
    
    
    fig997 = pylab.figure(997)
    pylab.clf()
    
    bigC = bigCdefault*1
    smallC = smallCdefault*1.25
    linewidth = 2.75
    latticeStretch = 0.7
    
    #1) huse
    testEuclidLattice = EuclideanLayout(4,2,lattice_type = 'Huse', modeType = 'FW')
    testLattice = GeneralLayout(testEuclidLattice.get_all_resonators() , modeType = testEuclidLattice.modeType, name =  testEuclidLattice.lattice_type)
#    ax = pylab.subplot(1,2,1)
    ax = pylab.subplot(1,2,1, adjustable='box', aspect=latticeStretch)
    pylab.cla()
    testLattice.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = linewidth)
    xs = testLattice.coords[:,0]
    ys = testLattice.coords[:,1]
    pylab.sca(ax)
    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallC , marker = 'o', edgecolors = 'k', zorder = 5)
    ax.axis('off')
    
    
#    ax = pylab.subplot(1,2,2)
    ax = pylab.subplot(1,2,2, adjustable='box', aspect=latticeStretch)
    pylab.cla()
    testLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = linewidth, minus_links = True, minus_color = 'gold', alpha=1)
    pylab.sca(ax)
    pylab.scatter(testLattice.SDx, testLattice.SDy ,c =  FWsiteColor, s = smallC , marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
    ax.axis('off')
    
    fig997.set_size_inches(6.4, 4.8)
    pylab.tight_layout()
    pylab.show()
    
    #fig997.savefig('HPKlattice.png',transparent= False, dpi = 200)
    
    
    fig996 = pylab.figure(996)
    pylab.clf()
    
    test1 = PlanarLayout(gon =7, vertex = 3, modeType = 'FW')
    test1.populate(2)
    
    test1_HW = PlanarLayout(gon =7, vertex = 3, modeType = 'HW')
    test1_HW.populate(2)
    
    ax = pylab.subplot(1,2,1)
    pylab.cla()
    test1.draw_SDlinks(ax, color = FWlinkColor, linewidth = 3.5, minus_links = True, minus_color = 'gold', alpha=1)
    pylab.sca(ax)
    pylab.scatter(test1.SDx, test1.SDy ,c =  FWsiteColor, s = smallCdefault*1.5 , marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
    ax.axis('off')
    ax.set_aspect('equal')
    
    ax = pylab.subplot(1,2,2)
    pylab.cla()
    test1_HW.draw_SDlinks(ax, color = HWlinkColor, linewidth = 3.5, minus_links = True, minus_color = HWminusLinkColor, alpha=1)
    pylab.sca(ax)
    pylab.scatter(test1_HW.SDx, test1_HW.SDy ,c =  HWsiteColor, s = smallCdefault*1.5 , marker = 'o', edgecolors = HWsiteEdgeColor, zorder = 5, alpha=1)
    ax.axis('off')
    ax.set_aspect('equal')
    
    
    fig996.set_size_inches(11, 5.8)
    pylab.tight_layout()
    pylab.show()
    
    
    #fig996.savefig('FW_HW_hyperbolics.png',transparent= False, dpi = 200)
else:
    fig999 = pylab.figure(999)
    pylab.clf()
    pylab.close(fig999)     
    
    fig998 = pylab.figure(998)
    pylab.clf()
    pylab.close(fig998) 
    
    fig997 = pylab.figure(997)
    pylab.clf()
    pylab.close(fig997)
    
    fig996 = pylab.figure(996)
    pylab.clf()
    pylab.close(fig996)
    
    
   
    
    
    
    
    
    
    
#%%
########
#figure infinity: split graphene, split hpg, split thingy
#######
def plot_any_band_cut(ax, cut):
    
    colorlist = ['firebrick', 'dodgerblue', 'blueviolet', 'mediumblue', 'goldenrod', 'cornflowerblue']
    
    pylab.sca(ax)
    
    shape = cut.shape
    
    for ind in range(0,shape[0]):
        colorInd = numpy.mod(ind, len(colorlist))
        pylab.plot(cut[ind,:], color = colorlist[colorInd] , marker = '.', markersize = '5', linestyle = '')
        
    pylab.title('some momentum cut')
    pylab.ylabel('Energy')
    pylab.xlabel('k_something')
        
    return
    



if LLSplots:
    cellModeType = 'FW'
    #cellModeType = 'HW'
    
    

    bigC = 70
    smallC = 25
    smallC2 = 20
    
 

    
    gs = gridspec.GridSpec(1, 2,
                           width_ratios=[1, 0.8],
                           height_ratios=[1]
                           )
    gs2 = gridspec.GridSpec(1, 2,
                           width_ratios=[1, 0.8],
                           height_ratios=[1]
                           )
    gs3 = gridspec.GridSpec(1, 2,
                           width_ratios=[1, 0.8],
                           height_ratios=[1]
                           )
    
    
    fig4567 = pylab.figure(4567)
    pylab.clf()
    
    
    #1) split graphene
    testEuclidLattice = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = cellModeType)
    resonators = testEuclidLattice.resonators
    resonators = split_resonators(resonators)
    testLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'splitkagome')
    
    cell0 = UnitCell('kagome')
    res0 = cell0.resonators
    res1 = split_resonators(res0)
    testCell = UnitCell('split_kagome', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
    
    
    #latticeStretch = 0.7
    latticeStretch = 1.0
    
    
    
   
    
    ####graphs
    ax = pylab.subplot(gs[0], adjustable='box', aspect=latticeStretch)

    newgraph = generate_line_graph(testLattice.resonators)
    newLattice = GeneralLayout(newgraph , modeType = test1.modeType, name =  'kagome_McLaughlin')
    
    pylab.cla()
    newLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha=1)
    pylab.sca(ax)
    pylab.scatter(newLattice.SDx, newLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
    ax.axis('off')
    
    #band structure
    numSurfPoints = 300
    kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)

    #hack the effective band structure
    effectiveBandStructureX = cutx[:,:] +1 #remove the flat bands and shift down
    effectiveBandStructureY = cuty[:,:] +1 #remove the flat bands and shift down
    FB = numpy.zeros((1, numSurfPoints))-2
    effectiveCutX = numpy.concatenate((FB, effectiveBandStructureX))
    effectiveCutY = numpy.concatenate((FB, effectiveBandStructureY))
    
    ax = pylab.subplot(gs[1], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
    #testCell.plot_band_cut(ax, cuty)
    plot_any_band_cut(ax, effectiveCutY)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, effectiveCutY.shape[1]/2, effectiveCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        



    fig4568 = pylab.figure(4568)
    pylab.clf()

    #2) split hpg
    testEuclidLattice = EuclideanLayout(2,1,lattice_type = 'Huse', modeType = cellModeType)
    resonators = testEuclidLattice.resonators
    resonators = split_resonators(resonators)
    testLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'splitHPG')
    
    cell0 = UnitCell('Huse')
    res0 = cell0.resonators
    res1 = split_resonators(res0)
    testCell = UnitCell('split_hpg', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
    
    #latticeStretch = 0.7
    latticeStretch = 1.0
    
    
    
    
    
    ####graphs
    ax = pylab.subplot(gs2[0], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    

    newgraph = generate_line_graph(testLattice.resonators)
    newLattice = GeneralLayout(newgraph , modeType = test1.modeType, name =  'hpk_McLaughlin')
    
    pylab.cla()
    newLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha=1)
    pylab.sca(ax)
    pylab.scatter(newLattice.SDx, newLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
    ax.axis('off')
    
    #band structure
    numSurfPoints = 300
    kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)

    #hack the effective band structure
    effectiveBandStructureX = cutx[:,:] +1 #remove the flat bands and shift down
    effectiveBandStructureY = cuty[:,:] +1 #remove the flat bands and shift down
    FB = numpy.zeros((1, numSurfPoints))-2
    effectiveCutX = numpy.concatenate((FB, effectiveBandStructureX))
    effectiveCutY = numpy.concatenate((FB, effectiveBandStructureY))
    
    ax = pylab.subplot(gs2[1], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
    #testCell.plot_band_cut(ax, cuty)
    plot_any_band_cut(ax, effectiveCutY)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, effectiveCutY.shape[1]/2, effectiveCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
    
    
    
    
    
    fig4569 = pylab.figure(4569)
    pylab.clf()
    
    #1) extremal hofman
    ########split Euclidean, hoffman attemps
    ##!!!!!!!!needs to be this size for the cell to come out right. DO NOT CHANGE!!!
    test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    resonators0 = test1.resonators #graphene layout
    splitGraph = split_resonators(resonators0) #split graphene layout
    resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
    resonators2 = split_resonators(resonators1) #split further
    resonators = resonators2
    testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
    #sites = [80,81,82, 83,79,78,77,76, 71, 84,85,86,87,114,115, 116, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
    sites = [80,81,83,79,78,77,76, 71, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
    state = numpy.zeros(len(testLattice.SDx))
    state[sites] = 0.1
    newCell = testLattice.resonators[sites, :]
#    
#    fig1 = pylab.figure(114)
#    pylab.clf()
#    ax = pylab.subplot(1,1,1)
#    testLattice.draw_resonator_lattice(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
#    xs = testLattice.coords[:,0]
#    ys = testLattice.coords[:,1]
#    pylab.sca(ax)
#    pylab.scatter(xs, ys ,c =  'goldenrod', s = 30, marker = 'o', edgecolors = 'k', zorder = 5)
#    pylab.scatter(testLattice.SDx[sites], testLattice.SDy[sites], c =  'firebrick', s = 75, marker = 'o', edgecolors = 'maroon', zorder = 5,linewidth = 1)
#    ax.set_aspect('equal')
#    ax.axis('off')
#    pylab.title(testLattice.name)
#    pylab.tight_layout()
#    pylab.show()
#    fig49 = pylab.figure(49)
    
    testCell = UnitCell('extremal_hofmann', resonators = newCell, a1 = test1.unitcell.a1, a2 = test1.unitcell.a2)
    
    
    #latticeStretch = 0.7
    latticeStretch = 1.0
    
    


    
    
    
    ####redo the lattice for plotting
    test1 = EuclideanLayout(5,5,lattice_type = 'kagome', modeType = 'FW')
    resonators0 = test1.resonators #graphene layout
    splitGraph = split_resonators(resonators0) #split graphene layout
    resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
    resonators2 = split_resonators(resonators1) #split further
    resonators = resonators2
    testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
    
    
    
    ax = pylab.subplot(gs3[0], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    
    newgraph = generate_line_graph(testLattice.resonators)
    newLattice = GeneralLayout(newgraph , modeType = test1.modeType, name =  'hofmannlayout')
    
    pylab.cla()
    newLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = 2.5, minus_links = True, minus_color = 'gold', alpha=1)
    pylab.sca(ax)
    pylab.scatter(newLattice.SDx, newLattice.SDy ,c =  FWsiteColor, s = smallCdefault, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
    ax.axis('off')
    width = 1.7 +0.045
    xStart = -0.005
    yStart = 1.625
    ax.set_xlim([xStart, xStart + width])
    ax.set_ylim([yStart, yStart + width])
    
    #band structure
    numSurfPoints = 300
    kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)

    #hack the effective band structure
    effectiveBandStructureX = cutx[:,:] +1 #remove the flat bands and shift down
    effectiveBandStructureY = cuty[:,:] +1 #remove the flat bands and shift down
    FB = numpy.zeros((1, numSurfPoints))-2
    effectiveCutX = numpy.concatenate((FB, effectiveBandStructureX))
    effectiveCutY = numpy.concatenate((FB, effectiveBandStructureY))
    
    ax = pylab.subplot(gs3[1], adjustable='box', aspect=80)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
    #testCell.plot_band_cut(ax, cuty)
    plot_any_band_cut(ax, effectiveCutY)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, effectiveCutY.shape[1]/2, effectiveCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        

    

#    
    fig4567.set_size_inches(10.25, 7.5)
    pylab.figure(4567)
    pylab.tight_layout()
    pylab.show()
    
    fig4568.set_size_inches(10.25, 7.5)
    pylab.figure(4568)
    pylab.tight_layout()
    pylab.show()
    
    fig4569.set_size_inches(12.6, 7.45)
    pylab.figure(4569)
    pylab.tight_layout()
    pylab.show()
    
    
#    fig4567.savefig('graphene_McLaughlin.png',transparent= False, dpi = 200)
#    fig4567.savefig('graphene_McLaughlin.svg',transparent= True)
    
#    fig4568.savefig('hpg_McLaughlin.png',transparent= False, dpi = 200)
#    fig4568.savefig('hpg_McLaughlin.svg',transparent= True)

#    fig4569.savefig('hoffman_McLaughlin.png',transparent= False, dpi = 200)
#    fig4569.savefig('hoffman_McLaughlin.svg',transparent= True)    
    
else:
    fig4567 = pylab.figure(4567)    
    pylab.clf()
    pylab.suptitle('Split Graphs')
    pylab.close(fig4567)
    
    fig4568 = pylab.figure(4568)    
    pylab.clf()
    pylab.suptitle('Split Graphs2')
    pylab.close(fig4568)
    
    fig4569 = pylab.figure(4569)    
    pylab.clf()
    pylab.suptitle('Split Graphs3')
    pylab.close(fig4569)
    
    
    
    
    
    
    
    
    
    
    
#%%
########
#Alicum and maybe line graph of split graphene
#######
    
if AliciumFig:
    
#    Both = False
    Both = True
    
    plotSurfaces = True
    
    
    if Both:
        gs = gridspec.GridSpec(1, 4,
                       width_ratios=[0.7, 0.6, 1.5,0.7],
                       height_ratios=[1]
                           )
#        smallC = smallCdefault
#        linewidth = 2.5
        
#        smallC = 25
        smallC = 22
        linewidth = 1.5
    else:
        gs = gridspec.GridSpec(1, 2,
                       width_ratios=[1, 0.55],
                       height_ratios=[1]
                           ) 
        smallC = 25
        linewidth = 1.5
        
    
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
    #######
    
    
    
    
    fig994 = pylab.figure(994)
    pylab.clf()

    
    
    if Both:
        #1) split graphene
        testEuclidLattice = EuclideanLayout(3,5,lattice_type = 'kagome', modeType = cellModeType)
        resonators = testEuclidLattice.resonators
        resonators = split_resonators(resonators)
        testLattice = GeneralLayout(resonators , modeType = testEuclidLattice.modeType, name =  'splitkagome')
        
        cell0 = UnitCell('kagome')
        res0 = cell0.resonators
        res1 = split_resonators(res0)
        testCell = UnitCell('split_graphene', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
        
        
        #latticeStretch = 0.7
        latticeStretch = 1.0
            
        #take the line graph to get 3-regular layout
        newgraph = generate_line_graph(testLattice.resonators)
        newLattice = GeneralLayout(newgraph , modeType = 'FW', name =  'kagome_McLaughlin')
        
    
        ax = pylab.subplot(gs[0], adjustable='box', aspect=latticeStretch)
        pylab.cla()
        newLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = linewidth, minus_links = True, minus_color = 'gold', alpha=1)
        pylab.sca(ax)
        pylab.scatter(newLattice.SDx, newLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
        ax.axis('off')
        xwidth = 1.7 +0.175
        xStart = -0.07
#        ywidth = 1.7 +0.035
#        yStart = 0.635
#        ywidth = 2.7 +0.035
#        yStart = 0 + 0.635
        ywidth = 3.3 +0.035
        yStart = 0.2 + 0.635
        ax.set_xlim([xStart, xStart + xwidth])
        ax.set_ylim([yStart, yStart + ywidth])
        
        #band structure
        numSurfPoints = 300
        kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
        kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)
    
        #hack the effective band structure
        effectiveBandStructureX = cutx[:,:] +1 #remove the flat bands and shift down
        effectiveBandStructureY = cuty[:,:] +1 #remove the flat bands and shift down
        FB = numpy.zeros((1, numSurfPoints))-2
        effectiveCutX = numpy.concatenate((FB, effectiveBandStructureX))
        effectiveCutY = numpy.concatenate((FB, effectiveBandStructureY))
        
        ax = pylab.subplot(gs[1], adjustable='box', aspect=105)
        pylab.cla()
        #plot_band_cut_mc(ax, cuty)
        #testCell.plot_band_cut(ax, cuty)
        plot_any_band_cut(ax, effectiveCutY)
        pylab.title('')
        pylab.ylabel('Energy (|t|)')
        pylab.xlabel('$k_y$ ($\pi$/a)')
        pylab.xticks([0, effectiveCutY.shape[1]/2, effectiveCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        
        
        if plotSurfaces:
            fig = pylab.figure(479)
            pylab.clf()
        
            bands = numpy.zeros((testCell.numSites, len(xs), len(ys)))
            
            #####compute surfaces
            for xind in range(0,len(xs)):
                for yind in range(0,len(ys)):
                    xval = xs[xind]
                    yval = ys[yind]
                    kx,yk, Es = testCell.compute_band_structure(xval, yval, xval, yval, numsteps = 1, modeType = 'FW')
                    bands[:, xind, yind] = numpy.transpose(Es)
            
            #hack it so that I get the band strucutre of the line graph of the line graph
            #hack the effective band structure
            effectiveIsh = bands[:,:,:] +1 #remove the flat bands and shift down
            FB = numpy.zeros((1, len(xs), len(xs)))-2
            effectiveBands = numpy.concatenate((FB, effectiveIsh))
    
            
            ####plot all surfaces together
            ax = fig.gca(projection='3d')
            for ind in range(0,testCell.numSites+1):
#                surf = ax.plot_surface(Xgrid, Ygrid, effectiveBands[ind,:,:], cmap='autumn',
#                                       linewidth=0, antialiased=False)
#                surf = ax.plot_surface(Xgrid, Ygrid, effectiveBands[ind,:,:], cmap='Blues_r',
#                                       linewidth=0, antialiased=False)
                surf = ax.plot_surface(Xgrid, Ygrid, effectiveBands[ind,:,:], cmap='YlGnBu_r',
                                       linewidth=0, antialiased=False)
#                surf = ax.plot_surface(Xgrid, Ygrid, effectiveBands[ind,:,:], cmap='hot',
#                                       linewidth=0, antialiased=False)
#                surf = ax.plot_surface(Xgrid, Ygrid, effectiveBands[ind,:,:],color = 'mediumblue',
#                                       linewidth=0, antialiased=False)
#                surf = ax.plot_surface(Xgrid, Ygrid, effectiveBands[ind,:,:],color = 'mediumblue',
#                                       linewidth=0, antialiased=False)
                ax.set_zlim(-2.5, 4.5)
            #    ax.set_zlim(-4*J, 4*J)
            
            ax.view_init(4, 0)
            
            plt.title('dispersion curves')
            plt.show()
            
            #back to old figure
            pylab.figure(994)
        
#            fig = pylab.figure(479)
#            ax = fig.gca(projection='3d')
#            ax.view_init(4, 10)
            

    #1) extremal hofman
    ########split Euclidean, hoffman attemps
    ##!!!!!!!!needs to be this size for the cell to come out right. DO NOT CHANGE!!!
    test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    resonators0 = test1.resonators #graphene layout
    splitGraph = split_resonators(resonators0) #split graphene layout
    resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
    resonators2 = split_resonators(resonators1) #split further
    resonators = resonators2
    testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
    sites = [80,81,83,79,78,77,76, 71, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
    state = numpy.zeros(len(testLattice.SDx))
    state[sites] = 0.1
    newCell = testLattice.resonators[sites, :]
#    
    
    testCell = UnitCell('extremal_hofmann', resonators = newCell, a1 = test1.unitcell.a1, a2 = test1.unitcell.a2)
    
    
    #latticeStretch = 0.7
    latticeStretch = 1.0
    
    
    ####redo the lattice for plotting
    test1 = EuclideanLayout(5,5,lattice_type = 'kagome', modeType = 'FW')
    resonators0 = test1.resonators #graphene layout
    splitGraph = split_resonators(resonators0) #split graphene layout
    resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
    resonators2 = split_resonators(resonators1) #split further
    resonators = resonators2
    testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
    
    
    newgraph = generate_line_graph(testLattice.resonators)
    newLattice = GeneralLayout(newgraph , modeType = test1.modeType, name =  'hofmannlayout')
    
    if Both:
        ax = pylab.subplot(gs[2], adjustable='box', aspect=latticeStretch)
    else:
        ax = pylab.subplot(gs[0], adjustable='box', aspect=latticeStretch)
    pylab.cla()
    newLattice.draw_SDlinks(ax, color = FWlinkColor, linewidth = linewidth, minus_links = True, minus_color = 'gold', alpha=1)
    pylab.sca(ax)
    pylab.scatter(newLattice.SDx, newLattice.SDy ,c =  FWsiteColor, s = smallC, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5, alpha=1)
    ax.axis('off')
#    width = 1.7 +0.045
#    xStart = -0.005
#    yStart = 1.625
    xwidth = 1.7 +0.125
    ywidth = 1.7 -0.105
    xStart = -0.045
    yStart = 1.705
    ax.set_xlim([xStart, xStart + xwidth])
    ax.set_ylim([yStart, yStart + ywidth])
    
    #band structure
    numSurfPoints = 300
    kx_x, ky_y, cutx = testCell.compute_band_structure(-2*numpy.pi, 0, 2*numpy.pi, 0, numsteps = numSurfPoints, modeType = cellModeType)
    kx_y, ky_y, cuty = testCell.compute_band_structure(0, -2.5*numpy.pi, 0, 2.5*numpy.pi, numsteps = numSurfPoints, modeType = cellModeType)

    #hack the effective band structure
    effectiveBandStructureX = cutx[:,:] +1 #remove the flat bands and shift down
    effectiveBandStructureY = cuty[:,:] +1 #remove the flat bands and shift down
    FB = numpy.zeros((1, numSurfPoints))-2
    effectiveCutX = numpy.concatenate((FB, effectiveBandStructureX))
    effectiveCutY = numpy.concatenate((FB, effectiveBandStructureY))
    
    if Both:
        ax = pylab.subplot(gs[3], adjustable='box', aspect=90)
    else:
        ax = pylab.subplot(gs[1], adjustable='box', aspect=85)
    pylab.cla()
    #plot_band_cut_mc(ax, cuty)
    #testCell.plot_band_cut(ax, cuty)
    plot_any_band_cut(ax, effectiveCutY)
    pylab.title('')
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('$k_y$ ($\pi$/a)')
    pylab.xticks([0, effectiveCutY.shape[1]/2, effectiveCutY.shape[1]], [-2.5,0,2.5], rotation='horizontal')
        

    

    if plotSurfaces:
        fig = pylab.figure(478)
        pylab.clf()
    
        bands = numpy.zeros((testCell.numSites, len(xs), len(ys)))
        
        #####compute surfaces
        for xind in range(0,len(xs)):
            for yind in range(0,len(ys)):
                xval = xs[xind]
                yval = ys[yind]
                kx,yk, Es = testCell.compute_band_structure(xval, yval, xval, yval, numsteps = 1, modeType = 'FW')
                bands[:, xind, yind] = numpy.transpose(Es)
        
        #hack it so that I get the band strucutre of the line graph of the line graph
        #hack the effective band structure
        effectiveIsh = bands[:,:,:] +1 #remove the flat bands and shift down
        FB = numpy.zeros((1, len(xs), len(xs)))-2
        effectiveBands = numpy.concatenate((FB, effectiveIsh))
        
        ####plot all surfaces together
        ax = fig.gca(projection='3d')
        for ind in range(0,testCell.numSites+1):
#            surf = ax.plot_surface(Xgrid, Ygrid, effectiveBands[ind,:,:], cmap=cm.coolwarm,
#                                   linewidth=0, antialiased=False)
#            surf = ax.plot_surface(Xgrid, Ygrid, effectiveBands[ind,:,:], color = 'mediumblue',
#                                   linewidth=0, antialiased=False)
            surf = ax.plot_surface(Xgrid, Ygrid, effectiveBands[ind,:,:], cmap='YlGnBu_r',
                                       linewidth=0, antialiased=False)
            ax.set_zlim(-2.5, 4.5)
        #    ax.set_zlim(-4*J, 4*J)
        ax.view_init(0, 0)
#        ax.view_init(1, 0)
        
        plt.title('dispersion curves')
        plt.show()

#        fig = pylab.figure(478)
#        ax = fig.gca(projection='3d')
#        ax.view_init(4, 10)

    
    


    if Both:
#        fig994.set_size_inches(11.8, 3.9)
        fig994.set_size_inches(14.3, 5.0)
        pylab.figure(994)
        pylab.tight_layout()
        pylab.show()
        
#        fig994.savefig('OptimalGaps_both.png',transparent= False, dpi = 200)
#        fig994.savefig('OptimalGaps_both.svg',transparent= True)  
    else:
        fig994.set_size_inches(8.5, 5)
        pylab.figure(994)
        pylab.tight_layout()
        pylab.show()
        
        
    
else:
    fig994 = pylab.figure(994)    
    pylab.clf()
    pylab.suptitle('Alicium')
    pylab.close(fig994)    
    
    
    
    
    
    
    
#%%
#    
####doodling
    

    
#pylab.figure(1479)
#pylab.clf()
##ax = pylab.subplot(1,2,1)
# 
##3) tree
#testTree = TreeResonators(degree = 3, iterations = 4, side = 1, file_path = '', modeType = 'FW')
##testTree = TreeResonators(degree = 3, iterations = 5, side = 1, file_path = '', modeType = 'HW')
#resonators = testTree.get_all_resonators()
#tempRes = split_resonators(resonators)
#treeLattice = GeneralLayout(tempRes , modeType = testTree.modeType, name =  'TREEEEE')
##ax3 = pylab.subplot(gs[2])
#ax3 = pylab.subplot(1,2,1)
#pylab.cla()
#if McLaughlinColors:
#    #splitGraphPrecursor
#    treeLattice.draw_resonator_lattice(ax3, color = 'b', alpha = FWlinkAlpha , linewidth = 2.5)
#    xs = treeLattice.coords[:,0]
#    ys = treeLattice.coords[:,1]
#    pylab.sca(ax3)
##    pylab.scatter(xs, ys ,c =  'gold', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
##    pylab.scatter(xs, ys ,c =  'mediumblue', s = smallC, marker = 'o', edgecolors = 'navy', zorder = 5)
##    pylab.scatter(xs, ys ,c =  'b', s = smallC, marker = 'o', edgecolors = 'navy', zorder = 5)
#    pylab.scatter(xs, ys ,c =  'dodgerblue', s = smallC, marker = 'o', edgecolors = 'navy', zorder = 5)
#    
#    xs = testTree.coords[:,0]
#    ys = testTree.coords[:,1]
#    pylab.scatter(xs, ys ,c = 'lightsteelblue', s = smallC, marker = 'o', edgecolors = 'midnightblue', zorder = 5)
#    ax3.set_aspect('equal')
#    ax3.axis('off')
#else:
#    treeLattice.draw_resonator_lattice(ax3, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#    xs = treeLattice.coords[:,0]
#    ys = treeLattice.coords[:,1]
#    pylab.sca(ax3)
#    pylab.scatter(xs, ys ,c =  layoutCapColor, s = smallC, marker = 'o', edgecolors = 'k', zorder = 5)
#    ax3.set_aspect('equal')
#    ax3.axis('off')
#    #pylab.title('layout graph')
#
#
#ax7 = pylab.subplot(1,2,2)
##ax7 = pylab.subplot(gs[4+2])
#if McLaughlinColors:
#    treeLattice.draw_SDlinks(ax7, color = layoutLineColor, linewidth = 1.5, minus_links = True, minus_color = 'gold', alpha = FWlinkAlpha)
#    pylab.sca(ax7)
##    pylab.scatter(treeLattice.SDx, treeLattice.SDy ,c =  layoutCapColor, s = 15, marker = 'o', edgecolors = 'k', zorder = 5)
#    pylab.scatter(treeLattice.SDx, treeLattice.SDy ,c =  'gainsboro', s = 15, marker = 'o', edgecolors = 'k', zorder = 5)
#else:
#    treeLattice.draw_SDlinks(ax7, color = FWlinkColor, linewidth = 1.5, minus_links = True, minus_color = 'gold')
#    pylab.sca(ax7)
#    pylab.scatter(treeLattice.SDx, treeLattice.SDy ,c =  FWsiteColor, s = 15, marker = 'o', edgecolors = FWsiteEdgeColor, zorder = 5)
#
##    state_vect = treeLattice.Psis[:,44]
##    stepSize = 0.11/2
##    tempVect = state_vect/stepSize
##    maxVal = numpy.max(state_vect)
##    for ind in range(0, len(state_vect)):
##        rawVal = state_vect[ind]
##        if numpy.abs(rawVal/stepSize)>8:
##            #this is the peak
##            pass
##        elif numpy.abs(rawVal/stepSize)>4:
##            #second ring
##            state_vect[ind] = maxVal*numpy.sign(rawVal)/2
##        elif numpy.abs(rawVal/stepSize)>2:
##            #third ring
##            state_vect[ind] = maxVal*numpy.sign(rawVal)/4
##        elif numpy.abs(rawVal/stepSize)>1:
##            #fourth ring
##            state_vect[ind] = maxVal*numpy.sign(rawVal)/8
##        else:
##            state_vect[ind] = 0
#
#
##    #make a flat band state
#state_vect = numpy.zeros(len(treeLattice.SDx))
#state_vect[0] = 1
#state_vect[1] = 1
#
#state_vect[2] = -1
#state_vect[3] = -1
#
#state_vect[6] = -0.5
#state_vect[7] = -0.5
#
#state_vect[8] = -0.5
#state_vect[9] = -0.5
#
#state_vect[10] = 0.5
#state_vect[11] = 0.5
#
#state_vect[12] = 0.5
#state_vect[13] = 0.5
#
#state_vect[18] = 0.25
#state_vect[19] = 0.25
#
#state_vect[20] = 0.25
#state_vect[21] = 0.25
#
#state_vect[22] = 0.25
#state_vect[23] = 0.25
#
#state_vect[24] = 0.25
#state_vect[25] = 0.25
#
#state_vect[26] = -0.25
#state_vect[27] = -0.25
#
#state_vect[28] = -0.25
#state_vect[29] = -0.25
#
#state_vect[30] = -0.25
#state_vect[31] = -0.25
#
#state_vect[32] = -0.25
#state_vect[33] = -0.25
#
#state_vect[42:58] = 0.5*numpy.concatenate((state_vect[26:34], state_vect[26:34]))
#
#state_vect[58:74] = - state_vect[42:58]
#
#state_vect = state_vect*0.45
#        
#reds = numpy.where(state_vect <0)[0]
#yellows = numpy.where(state_vect > 0)[0]
#
#scaleFactor = 8
##    scaleFactor = 64
#
#redState = numpy.zeros(len(state_vect))
#redState[reds] = state_vect[reds]
#redAmps = redState
#redProbs = numpy.abs(redAmps)**2
#redSizes = redProbs * len(redProbs)*scaleFactor
#
#yellowState = numpy.zeros(len(state_vect))
#yellowState[yellows] = state_vect[yellows]
#yellowAmps = yellowState
#yellowProbs = numpy.abs(yellowAmps)**2
#yellowSizes = yellowProbs * len(yellowProbs)*scaleFactor
#
##    pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  stateColor2, s = redSizes, marker = 'o', edgecolors = stateEdgeColor2,  zorder = 6, linewidth = stateEdgeWidth2)
##    pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  stateColor1, s = yellowSizes, marker = 'o', edgecolors = stateEdgeColor1,  zorder = 6, linewidth = stateEdgeWidth1)
#
#####possibility 3
#fillColor = stateColor1
#edgeColor = 'goldenrod'
#
#####possibility 1
##fillColor = 'lightgrey'
##edgeColor = 'grey'
#
#####possibility 2
##fillColor = 'mediumblue'
##edgeColor = 'navy'
#
######possibility 3
##fillColor = 'b'
##edgeColor = 'mediumblue'
#
#
##fillColor = 'cornflowerblue'
##edgeColor = 'lightsteelblue'
##
##fillColor = 'lightskyblue'
##edgeColor = 'lightsteelblue'
##
##fillColor = 'deepskyblue'
##edgeColor = 'lightskyblue'
#
##fillColor = 'lightsteelblue'
##edgeColor = 'cornflowerblue'
#
##fillColor = 'whitesmoke'
##edgeColor = 'lightgrey'
##
##fillColor = 'gold'
##edgeColor = 'darkgoldenrod'
#
#pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  stateColor2, s = redSizes, marker = 'o', edgecolors = stateEdgeColor2,  zorder = 6, linewidth = stateEdgeWidth2)
#pylab.scatter(treeLattice.SDx, treeLattice.SDy, c =  fillColor, s = yellowSizes, marker = 'o', edgecolors = edgeColor,  zorder = 6, linewidth = stateEdgeWidth1)
#
#
#
#ax7.set_aspect('equal')
#ax7.axis('off')
#
#
#
####go back to the tree and build the flat band state, relies on the original being ina  sensible order
#pylab.sca(ax3)
#ends = numpy.where(state_vect !=0)[0]
#xcoords = numpy.zeros(len(ends)/2)
#ycoords = numpy.zeros(len(ends)/2)
#treeReds = numpy.zeros(len(xcoords))
#treeYellows = numpy.zeros(len(ycoords))
#
#for ind in range(0, len(ends)/2):
#    site = ends[2*ind]
#    site2 = ends[2*ind+ 1]
#    x1 = treeLattice.SDx[site]
#    x2 = treeLattice.SDx[site2]
#    
#    y1 = treeLattice.SDy[site]
#    y2 = treeLattice.SDy[site2]
#    
#    xcoords[ind] = (x1 + x2)/2
#    ycoords[ind] = (y1 + y2)/2
#    
#    if state_vect[site] > 0:
#        #this is a yellow
#        treeYellows[ind] = yellowSizes[site]
#    else:
#        treeReds[ind] = redSizes[site]
#pylab.scatter(xcoords, ycoords, c =  stateColor2, s = treeReds, marker = 'o', edgecolors = stateEdgeColor2,  zorder = 6, linewidth = stateEdgeWidth2)
#pylab.scatter(xcoords, ycoords, c =  fillColor, s = treeYellows, marker = 'o', edgecolors = edgeColor,  zorder = 6, linewidth = stateEdgeWidth1)
#ax3.set_aspect('equal')
#ax3.axis('off')