#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 14:59:11 2018

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
import time

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

#########
#surface plots
########
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm





def compute_DOS(cell, numSurfPoints = 300, modeType = 'FW', freq_res = 0.04, detectFlatBands = True, split_layout = False):
    '''
    split_layout means somput the DOS the comes from the split grph of the layout
    So, if you call this on effective = kagome, it will comput split graphene. Hopefully.
    It'll maybe get the number of flat bands worng though. Beware
    '''
    
    divisions1 = numpy.linspace(0, 1, numSurfPoints)
    divisions2 = numpy.linspace(0, 1, numSurfPoints)
    
    
    Xgrid, Ygrid = numpy.meshgrid(divisions1, divisions2)
    
    #lattice vectors
    a1 = numpy.zeros(3)
    a1[0:2] = cell.a1
    
    a2 = numpy.zeros(3)
    a2[0:2] = cell.a2
    
    a3 = numpy.zeros(3)
    a3[2] = 1
    
    #reciprocal lattice vectors
    denom = numpy.dot(a1, numpy.cross(a2,a3))
    b1 = 2*numpy.pi * numpy.cross(a2,a3) /denom
    b2 = 2*numpy.pi * numpy.cross(a3,a1) /denom
    b3 = 2*numpy.pi * numpy.cross(a1,a2) /denom
    
    #compute the grid in k space
    for n in range(0, numSurfPoints):
        for m in range(0, numSurfPoints):
            ind1 = divisions1[n]
            ind2 = divisions2[m]
            
            vec = ind1*b1 + ind2*b2
            
            Xgrid[n,m] = vec[0]
            Ygrid[n,m] = vec[1]
            
    #allocate space for the bands
    if split_layout:
        bands = numpy.zeros((2*cell.numSites, len(divisions1), len(divisions2)))
    else:
        bands = numpy.zeros((cell.numSites, len(divisions1), len(divisions2)))
    
    #####compute surfaces
    for xind in range(0,len(divisions1)):
        for yind in range(0,len(divisions2)):
            xval = Xgrid[xind, yind]
            yval = Ygrid[xind, yind]
            kx,yk, Es = cell.compute_band_structure(xval, yval, xval, yval, numsteps = 1, modeType = modeType)
            if split_layout:
                effectiveEs = Es
                layoutEs = effectiveEs -1
                splitEs1 = numpy.sqrt(layoutEs + 3 + 10**-7) #adding a tiny offset to try to stabilize the numerics
                splitEs2 = -numpy.sqrt(layoutEs + 3 + 10**-7)
                
                splitEs = numpy.concatenate((splitEs1, splitEs2))
                bands[:, xind, yind] = numpy.transpose(splitEs)
            else:
                bands[:, xind, yind] = numpy.transpose(Es)
    

    #bin and comput DOS
    freq_range = 4.0+ freq_res/2
    freqs = scipy.arange(-2-freq_res/2, 4+ freq_res/2,  freq_res) + freq_res/2.
    if split_layout:
        freq_bins = scipy.arange(-3-freq_res/2, 3+ 1.5*freq_res/2, freq_res)
    else:
        freq_bins = scipy.arange(-2-freq_res/2, 4+ 1.5*freq_res/2, freq_res)
    
#    print bands.shape
#    print freq_bins
    [DOS, bins_out] = numpy.histogram(bands, freq_bins)
    bins_centers = (bins_out[0:-1] + bins_out[1:])/2
    binWidth = bins_out[1] - bins_out[0]
    
    #normalize DOS
    DOS = 1.*DOS/bands.size
    
    #autodetect flat bands
    if detectFlatBands:
        FBs = numpy.zeros(len(DOS))
        FBEs = numpy.where(DOS > 0.05)[0]
        
        DOS[FBEs] = 0
#        for ind in FBEs:
#            if ind ==0:
#                DOS[ind] = DOS[ind+1]
#            elif ind == len(DOS)-1:
#                DOS[ind] = DOS[ind-1]
#            else:
#                DOS[ind] = numpy.max((DOS[ind-1], DOS[ind+1]))
        
        FBs[FBEs] = 1.
        
        return DOS, FBs, bins_centers, binWidth
    else:
        return DOS, bins_centers, binWidth
        









#########################################################
#main body
#########

#saveLots = True
saveLots = False


#testSingleCell = True
testSingleCell = False




if testSingleCell:
    #testCell = UnitCell('Huse')
    #testCell = UnitCell('74Huse')
    testCell = UnitCell('84Huse')
    #testCell = UnitCell('123Huse')
    #testCell = UnitCell('kagome')
    
    
    
    #cell0 = UnitCell('Huse')
    #res0 = cell0.resonators
    #res1 = split_resonators(res0)
    #testCell = UnitCell('split_HPG', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
    
    
    #cell0 = UnitCell('kagome')
    #res0 = cell0.resonators
    #res1 = split_resonators(res0)
    #testCell = UnitCell('split_graphene', resonators = res1, a1 = cell0.a1, a2 = cell0.a2)
    
    
    
    
    ##1) extremal hofman
    #########split Euclidean, hoffman attemps
    ###!!!!!!!!needs to be this size for the cell to come out right. DO NOT CHANGE!!!
    #test1 = EuclideanLayout(3,3,lattice_type = 'kagome', modeType = 'FW')
    #resonators0 = test1.resonators #graphene layout
    #splitGraph = split_resonators(resonators0) #split graphene layout
    #resonators1 = generate_line_graph(splitGraph) #McLaughlin-esque kagome
    #resonators2 = split_resonators(resonators1) #split further
    #resonators = resonators2
    #testLattice = GeneralLayout(resonators , modeType = test1.modeType, name =  'hofmannAttempt')
    ##sites = [80,81,82, 83,79,78,77,76, 71, 84,85,86,87,114,115, 116, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
    #sites = [80,81,83,79,78,77,76, 71, 50, 49, 48, 47, 46, 44, 55,54, 53,52]
    #state = numpy.zeros(len(testLattice.SDx))
    #state[sites] = 0.1
    #newCell = testLattice.resonators[sites, :]
    #testCell = UnitCell('extremal_hofmann', resonators = newCell, a1 = test1.unitcell.a1, a2 = test1.unitcell.a2)






    
    
    
    pylab.figure(4)
    pylab.clf()
    ax = pylab.subplot(1,1,1)
    testCell.draw_resonators(ax, color = 'mediumblue', alpha = 1 , linewidth = 2.5)
    testCell.draw_resonator_end_points(ax, color = 'darkgoldenrod', edgecolor = 'k',  marker = 'o' , size = 35, zorder = 5)
    pylab.show()
    
    
    
    
    
    DOS, FBs, binfreqs, res = compute_DOS(testCell, numSurfPoints = 5, modeType = 'FW', freq_res = 0.04, detectFlatBands = True)
    #DOS, FBs, binfreqs, res = compute_DOS(testCell, numSurfPoints = 50, modeType = 'FW', freq_res = 0.04, detectFlatBands = True)
    #DOS, FBs, binfreqs, res = compute_DOS(testCell, numSurfPoints = 300, modeType = 'FW', freq_res = 0.04, detectFlatBands = True)
    
    
    fig1 =  pylab.figure(1)
    pylab.clf()
    ax = pylab.subplot(1,1,1) 
    pylab.bar(binfreqs, 1.*DOS, width = res, color =  'deepskyblue', label = str('fullDOS'), alpha = 1, align = 'center')
    pylab.bar(binfreqs, 1*FBs, width = res, color =  'mediumblue', label = str('fullDOS'), alpha = 1, align = 'center')
    ax.set_ylim([0,0.03])
    pylab.title(testCell.type)
    fig1.set_size_inches(9.2,2.75)
    pylab.xlabel('Energy (|t|)')
    pylab.ylabel('DOS')
    pylab.tight_layout()
    pylab.show()
    
    fig2 =  pylab.figure(2)
    pylab.clf()
    ax = pylab.subplot(1,1,1) 
    #pylab.plot(binfreqs, 1.*DOS, color =  'deepskyblue', label = str('fullDOS'), alpha = 1, zorder = 1)
    ax.fill_between(binfreqs, 0, DOS, color = 'deepskyblue', zorder = 2)
    #unmaskinds = numpy.where(DOS > 0)[0]
    #mask = numpy.ones(len(DOS))*numpy.max(DOS)
    #mask[unmaskinds] = 0
    #ax.fill_between(binfreqs, 0, mask, color = 'w', zorder = 3)
    pylab.bar(binfreqs, 1*FBs, width = res, color =  'mediumblue', label = str('fullDOS'), alpha = 1, align = 'center', zorder = 4)
    ax.set_ylim([0,0.03])
    pylab.title(testCell.type)
    fig2.set_size_inches(9.2,2.75)
    pylab.xlabel('Energy (|t|)')
    pylab.ylabel('DOS')
    pylab.tight_layout()
    pylab.show()
    
    
    






######
#save a ton of stuff
######
defaultRes = 0.04
defaultSurfPoints = 300


#alternateSurfPoints1 = 300
#alternateSurfPoints2 = 500
#alternateRes = 0.04

alternateSurfPoints1 = 500
alternateSurfPoints2 = 700
alternateRes = 0.02

#alternateSurfPoints1 = 50
#alternateSurfPoints2 = 50
#alternateRes = 0.04

saveFolder = r'/Users/kollar2/Documents/HouckLab/EuclideanLatticePlanning/DOSs/'


 

def save_DOS(saveCell, numSurfPoints = 5, res = 0.04, folder = '', split_layout = False):
    DOS, FBs, binfreqs, res = compute_DOS(saveCell, numSurfPoints = numSurfPoints, modeType = 'FW', freq_res = res, detectFlatBands = True, split_layout = split_layout)
    saveDict = {}
    saveDict['DOS'] = DOS
    saveDict['FBs'] = FBs
    saveDict['binfreqs'] = binfreqs
    saveDict['res'] = res
    saveDict['numSurfPoints'] = numSurfPoints
    if split_layout:
        saveDict['name'] = saveCell.type + '_splitLayout_' + str(numSurfPoints) + 'x' + str(numSurfPoints) + '_res_' + str(res)
    else:
        saveDict['name'] = saveCell.type + '_' + str(numSurfPoints) + 'x' + str(numSurfPoints) + '_res_' + str(res)
    fileName = saveDict['name'] + '.pkl'
    filePath = os.path.join(folder, fileName)
    pickle.dump(saveDict, open(filePath, 'wb'))
    return

if saveLots:
    print '   '
    print 'num points 1 = ' + str(alternateSurfPoints1)
    print 'num points 2 = ' + str(alternateSurfPoints2)
    print 'res = ' + str(alternateRes)
    print '  '
    
    
    
    
    #make all the cells that I want
    kagomeCell = UnitCell('kagome') #graphene layout, kagome effective
    
    kapheneCell = kagomeCell.split_cell(name = 'kaphene', splitIn = 2) #split graphene layout, kaphene effective
    
    lineKapheneCell = kapheneCell.line_graph_cell(name = 'lineKaphene', resonatorsOnly = False) #kaphene layout, line kaffene effective
    
    hpkCell = UnitCell('Huse') #HPG layout, HPK effective
    
    oskCell = UnitCell('84Huse')
    
    hpKapheneCell = hpkCell.split_cell(name = 'HPKaphene', splitIn = 2) #split hpg layout, heptagon-pentagon kaphene effective
    
    lineHPKapheneCell = hpKapheneCell.line_graph_cell(name = 'lineHPKaphene') #HPKaphene layout, line HP kapheen effective
    
    CharonCell = lineKapheneCell.split_cell(name = 'Charon', splitIn = 2) #split kaphene layout, Charon effective
    
    AliciumCell = CharonCell.line_graph_cell(name = 'Alicium', resonatorsOnly = False) #charonLayout, Alicium effective
    
    lineKagomeCell = kagomeCell.line_graph_cell(name = 'lineKagome', resonatorsOnly = False) #kagome layout, soemthing effective

    t0 = time.time()
    
#    print 'kagome'
#    save_DOS(kagomeCell, numSurfPoints = alternateSurfPoints1 , res = alternateRes, folder = saveFolder)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#    
#    print 'kaphene'
#    save_DOS(kapheneCell, numSurfPoints = alternateSurfPoints1 , res = alternateRes, folder = saveFolder)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#    
#    print 'linekaphene'
#    save_DOS(lineKapheneCell, numSurfPoints = alternateSurfPoints1 , res = alternateRes, folder = saveFolder)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#    
#    
#    
#    
#    print 'hpk'
#    save_DOS(hpkCell, numSurfPoints = alternateSurfPoints2, res = alternateRes, folder = saveFolder)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#    
#    print 'osk'
#    save_DOS(oskCell, numSurfPoints = alternateSurfPoints2, res = alternateRes, folder = saveFolder)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#    
#    print 'hpkaphene'
#    save_DOS(hpKapheneCell, numSurfPoints = alternateSurfPoints2, res = alternateRes, folder = saveFolder)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#    
#    
#    
#    
#    
#    print 'Charon'
#    save_DOS(CharonCell, numSurfPoints = alternateSurfPoints2, res = alternateRes, folder = saveFolder)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#
#    print 'Alicum'
#    save_DOS(AliciumCell, numSurfPoints = alternateSurfPoints2, res = alternateRes, folder = saveFolder)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#    
#    
#    print 'splitGraphene'
#    save_DOS(kagomeCell, numSurfPoints = alternateSurfPoints1, res = alternateRes, folder = saveFolder, split_layout = True)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#
#    print 'splitHPG'
#    save_DOS(hpkCell, numSurfPoints = alternateSurfPoints2, res = alternateRes, folder = saveFolder, split_layout = True)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1
#    
#    print 'splitKaphene'
#    save_DOS(lineKapheneCell, numSurfPoints = alternateSurfPoints1, res = alternateRes, folder = saveFolder, split_layout = True)
#    t1 = time.time()
#    print 'elapsed time = ' + str(t1-t0) + ' s'
#    t0 = t1



#DOS, FBs, binfreqs, res = compute_DOS(kagomeCell, numSurfPoints = 50, modeType = 'FW', freq_res = res, detectFlatBands = True, split_layout = True)


file_path = os.path.join(saveFolder, 'kagome_splitLayout_500x500_res_0.02.pkl')
#file_path = os.path.join(saveFolder, 'Charon_500x500_res_0.04.pkl')
#file_path = os.path.join(saveFolder, 'Alicium_500x500_res_0.04.pkl')
pickledict = pickle.load(open(file_path, "rb" ) )

DOS2 = pickledict['DOS']
FBs2 = pickledict['FBs']
freqs = pickledict['binfreqs']
res2 = pickledict['res']
name = pickledict['name']

pylab.figure(5)
pylab.clf()
ax = pylab.subplot(1,1,1) 
#pylab.plot(binfreqs, 1.*DOS2, color =  'deepskyblue', label = str('fullDOS'), alpha = 1, zorder = 1)
#ax.fill_between(binfreqs, 0, DOS2, color = 'deepskyblue', zorder = 2)
#unmaskinds = numpy.where(DOS > 0)[0]
#mask = numpy.ones(len(DOS))*numpy.max(DOS)
#mask[unmaskinds] = 0
#ax.fill_between(binfreqs, 0, mask, color = 'w', zorder = 3)

#pylab.bar(binfreqs, 1.*DOS2, width = res2, color =  'mediumslateblue', label = str('fullDOS'), alpha = 1, align = 'center')
#pylab.bar(binfreqs, 1.*DOS2, width = res2, color =  'cornflowerblue', label = str('fullDOS'), alpha = 1, align = 'center')
pylab.bar(binfreqs, 1.*DOS2, width = res2, color =  'dodgerblue', label = str('fullDOS'), alpha = 1, align = 'center')
#pylab.bar(binfreqs, 1.*DOS2, width = res2, color =  'goldenrod', label = str('fullDOS'), alpha = 1, align = 'center')
pylab.bar(binfreqs, 1*FBs2, width = res2, color =  'mediumblue', label = str('fullDOS'), alpha = 1, align = 'center', zorder = 4)


ax.set_ylim([0,0.03])
pylab.title(name)
fig2.set_size_inches(9.2,2.75)
pylab.xlabel('Energy (|t|)')
pylab.ylabel('DOS')
pylab.tight_layout()



