#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 19:06:59 2018

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


from scipy.sparse import coo_matrix



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





def makeMinPickle(gon,vertex,modeType, maxItter):
    '''
    Old function used for making minimal pickles of the energies for effective graphs.
    '''
    test = PlanarLayout(gon = gon, vertex = vertex, side =1, radius_method = 'lin', modeType = modeType)
    test.populate(maxItter = maxItter)
    
    pickleDict = {}
    pickleDict['gon'] = test.gon
    pickleDict['vertex'] = test.vertex
    pickleDict['itter'] = test.itter
    pickleDict['modeType'] = test.modeType
    pickleDict['Es'] = test.Es
    pickleDict['Eorder'] = test.Eorder

    
    pickleName = str(test.gon) + 'gon_' + str(test.vertex) + 'vertex_' + str(test.itter) + '_' + str(test.modeType) + '_DOSminimal'
    pickleDict['name'] = pickleName
    
    pickle.dump(pickleDict, open(pickleName + '.pkl', 'wb'))


def generate_vertex_dict_min(layout, roundDepth = 3):
    '''
    custom function so I can get vertex dict without having to run the full populate of general layout
    and thereby having to also diagonalize the effective Hamiltonian.
    
    Returns a dictionary of which vertices are connected to which.
    '''
    resonators = layout.get_all_resonators()
    plusEnds = numpy.round(resonators[:,0:2],roundDepth)
    minusEnds = numpy.round(resonators[:,2:4],roundDepth)
    
    vertexDict = {}
    
    #loop over the vertices.
    for vind in range(0, layout.coords.shape[0]):
        vertex = numpy.round(layout.coords[vind, :],roundDepth)
        
        startMatch = numpy.where((plusEnds == (vertex[0], vertex[1])).all(axis=1))[0]
        endMatch = numpy.where((minusEnds == (vertex[0], vertex[1])).all(axis=1))[0]
        
        matchList = []
        for rind in startMatch:
            matchList.append(int(rind))
        for rind in endMatch:
            matchList.append(int(rind))
         
        #store the results
        vertexDict[vind] = numpy.asarray(matchList)
    
    return vertexDict

#@profile
def generate_layout_Hamiltonian(layout, roundDepth = 3, t = 1, verbose = True, sparse = True, flags = 5):
    '''
    custom function so I can get vertex dict without having to run the full populate of general layout
    and thereby having to also diagonalize the effective Hamiltonian.
    '''
    resonators = layout.get_all_resonators()
    resonators = numpy.round(resonators, roundDepth)
    
    numVerts = layout.coords.shape[0]
    if sparse:
        rowVec = numpy.zeros(numVerts*4+flags)
        colVec = numpy.zeros(numVerts*4+flags)
        Hvec = numpy.zeros(numVerts*4+flags)
    else:
        Hmat = numpy.zeros((numVerts, numVerts))
    
    coords_complex = numpy.round(layout.coords[:,0] + 1j*layout.coords[:,1], roundDepth)
    
    currInd = 0
    for rind in range(0, resonators.shape[0]):
        resPos = resonators[rind,:]
        startPos = numpy.round(resPos[0],roundDepth)+ 1j*numpy.round(resPos[1],roundDepth)
        stopPos = numpy.round(resPos[2],roundDepth)+ 1j*numpy.round(resPos[3],roundDepth)
        
        startInd = numpy.where(startPos == coords_complex)[0][0]
        stopInd = numpy.where(stopPos == coords_complex)[0][0]

        if sparse:
            rowVec[currInd] = startInd
            colVec[currInd] = stopInd
            Hvec[currInd] = t #will end up adding t towhatever this entry was before.
            currInd = currInd +1
            
            rowVec[currInd] = stopInd
            colVec[currInd] = startInd
            Hvec[currInd] = t #will end up adding t towhatever this entry was before.
            currInd = currInd +1
            
        else:
            Hmat[startInd, stopInd] = Hmat[startInd, stopInd] + t
            Hmat[stopInd, startInd] = Hmat[stopInd, startInd] + t
    
    #finish making the sparse matrix if we are in sparse matrix mode.
    if sparse:
        #pad the end of the matrix with values so that I can see if one of those is the missing one
        for ind in range(0, flags):
#            rowVec[currInd] = currInd
#            colVec[currInd] = currInd
            rowVec[currInd] = numVerts+ind
            colVec[currInd] = numVerts+ind
            Hvec[currInd] = -7.5 #will end up adding t towhatever this entry was before.
            currInd = currInd +1

#        Hmat = coo_matrix((Hvec,(rowVec,colVec)), shape = (currInd,currInd), dtype = 'd')
        Hmat = coo_matrix((Hvec,(rowVec,colVec)), shape = (numVerts+flags,numVerts+flags), dtype = 'd')
        Hmat.eliminate_zeros() #removed the unused spots since this is not a regular graph
        
    if verbose:
        temp = numpy.sum(Hmat)/numVerts
        print 'average degree = ' + str(temp)
        
    return Hmat

#@profile
def makeMinLayoutPickle_hyperbolic(gon,vertex, maxItter,modeType = 'FW', verbose = True, saveFolder = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/layoutDOSs', sparse = True, flags = 5):
    '''
    Old function used for making minimal pickles of the energies for effective graphs.
    '''
    test = PlanarLayout(gon = gon, vertex = vertex, side =1, radius_method = 'lin', modeType = modeType)
    test.populate(maxItter = maxItter, resonatorsOnly = True)
    
    if sparse:
        Hmat = generate_layout_Hamiltonian(test, verbose = verbose, sparse = True, flags = flags)
#        Es,Psis = scipy.sparse.linalg.eigsh(Hmat, k = test.coords.shape[0]-1 + flags) #will get all but one eigenvalue
        Es= scipy.sparse.linalg.eigsh(Hmat, k = test.coords.shape[0]-1 + flags, return_eigenvectors=False) #will get all but one eigenvalue
        Eorder = numpy.argsort(Es)
        
        if Es[Eorder[flags-1]]< 7.4:
            #all the flag values are present, so an actual eigenvalue is missing
            print 'actual eigenvalue missing'
            Es = Es[flags:]
        else:
            print 'all real eigenvalues found'
            Es = Es[flags-1:]
    else:
        Hmat = generate_layout_Hamiltonian(test, verbose = verbose, sparse = False)
        Es,Psis = scipy.linalg.eigh(Hmat)
        Eorder = numpy.argsort(Es)
    
    if verbose:
        print 'max eigenvalue = ' + str(numpy.max(Es))
    
    pickleDict = {}
    pickleDict['gon'] = test.gon
    pickleDict['vertex'] = test.vertex
    pickleDict['itter'] = test.itter
    pickleDict['modeType'] = test.modeType
    pickleDict['Es'] = Es
    pickleDict['Eorder'] = Eorder
    if sparse:
        pickleDict['sparse'] = True
    else:
        pickleDict['sparse'] = False

    
#    pickleName = str(test.gon) + 'gon_' + str(test.vertex) + 'vertex_' + str(test.itter) + '_' + str(test.modeType) + '_layoutDOSminimal'
    pickleName1 = str(test.gon) + 'gon_' + str(test.vertex) + 'vertex_' + str(test.itter) + '_' + str(test.modeType) 
    if sparse:
        pickleName = pickleName1 +  '_splayoutDOSminimal'
    else:
        pickleName = pickleName1 +  '_layoutDOSminimal'
    pickleDict['name'] = pickleName
    
#    pickle.dump(pickleDict, open(pickleName + '.pkl', 'wb'))
    filePath = os.path.join(saveFolder, pickleName + '.pkl')
    pickle.dump(pickleDict, open(filePath, 'wb'))
    return Es
    

    
#    pickleName = str(test.gon) + 'gon_' + str(test.vertex) + 'vertex_' + str(test.itter) + '_' + str(test.modeType) + '_layoutDOSminimal'
    pickleName1 = str(test.gon) + 'gon_' + str(test.vertex) + 'vertex_' + str(test.itter) + '_' + str(test.modeType) 
    if sparse:
        pickleName = pickleName1 +  '_splayoutDOSminimal'
    else:
        pickleName = pickleName1 +  '_layoutDOSminimal'
    pickleDict['name'] = pickleName
    
#    pickle.dump(pickleDict, open(pickleName + '.pkl', 'wb'))
    filePath = os.path.join(saveFolder, pickleName + '.pkl')
    pickle.dump(pickleDict, open(filePath, 'wb'))
    return Es
    
def makeMinLayoutPickle_gen(genLayout, verbose = True, saveFolder = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/layoutDOSs', sparse = True, flags = 5):
    '''
    Old function used for making minimal pickles of the energies for effective graphs.
    
    '''
    if sparse:
        Hmat = generate_layout_Hamiltonian(genLayout, verbose = verbose,sparse = True)
        Es,Psis = scipy.sparse.linalg.eigsh(Hmat, k = genLayout.coords.shape[0]-1+flags)
        Eorder = numpy.argsort(Es)
        
        if Es[Eorder[flags-1]]< -7.4:
            #all the flag values are present, so an actual eigenvalue is missing
            print 'actual eigenvalue missing'
            Es = Es[flags:]
        else:
            print 'all real eigenvalues found'
            Es = Es[flags-1:]
    else:
        Hmat = generate_layout_Hamiltonian(genLayout, verbose = verbose, sparse = False)
        Es,Psis = scipy.linalg.eigh(Hmat)
        Eorder = numpy.argsort(Es)
    
    if verbose:
        print 'max eigenvalue = ' + str(numpy.max(Es))
    
    pickleDict = {}
    pickleDict['name'] = genLayout.name
    pickleDict['gon'] = 'NA'
    pickleDict['vertex'] = 'NA'
    pickleDict['itter'] = 'NA'
    pickleDict['modeType'] = 'NA'
    pickleDict['Es'] = Es
    pickleDict['Eorder'] = Eorder
    if sparse:
        pickleDict['sparse'] = True
    else:
        pickleDict['sparse'] = False

    if sparse:
        pickleName = genLayout.name +  '_splayoutDOSminimal'
    else:
        pickleName = genLayout.name +  '_layoutDOSminimal'
    pickleDict['name'] = pickleName
    
    filePath = os.path.join(saveFolder, pickleName + '.pkl')
    pickle.dump(pickleDict, open(filePath, 'wb'))
    return Es

def loadMinLayoutPickle(name, verbose = True, saveFolder = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/layoutDOSs'):
    '''
    Old function used for making minimal pickles of the energies for effective graphs.
    
    '''
    filePath = os.path.join(saveFolder, name)
    pickleDict = pickle.load(open(filePath, "rb" ) )
    
    Es = pickleDict['Es']
    name = pickleDict['name']
    
    if verbose:
        if 'sparse' in pickleDict.keys():
            print 'sparse = ' + str(pickleDict['sparse'])
        else:
            print 'sparse = False'
        print 'max eigenvalue = ' + str(numpy.max(Es))
        print 'min eigenvalue = ' + str(numpy.min(Es)) 


    return Es , name, pickleDict

def plot_DOS(ax, Es, freq_res = 0.04, freq_max = 3, color = 'mediumblue', alpha = 1):
    pylab.sca(ax)
    
    #set up frequency sweep
    freq_range = freq_max  + freq_res/2
    
#    freqs = scipy.arange(-freq_range, freq_range, freq_res) + freq_res/2.
    freq_bins = scipy.arange(-freq_range, freq_range+freq_res, freq_res)
    
    [DOS, bins_out] = numpy.histogram(Es, freq_bins)
    
    bins_centers = (bins_out[0:-1] + bins_out[1:])/2
    binWidth = bins_out[1] - bins_out[0]
    
    pylab.bar(bins_centers, 1.*DOS/len(Es), width = binWidth, color = color, label = '', alpha = alpha)
    
    pylab.xlabel('Eigenvalue')
    pylab.ylabel('DOS (AU)')



def make_layoutGraph_vid(GenLayout, startInd =0, stopInd = -1, figNum = 8, xsize = 14, ysize = 4):
    
    
    if stopInd <0:
        #set a break point that will not be reached
        breakPoint = GenLayout.coords.shape[0]+1
    else:
        breakPoint = stopInd
        
    if startInd > GenLayout.coords.shape[0]:
        raise ValueError , 'dont have this many eigenvectors'
        
    
    folder = GenLayout.name + '_layoutGraph'

    
    currDir = os.getcwd()
    saveDir = os.path.join(currDir,folder)
    if os.path.isdir(saveDir):
        pass
    else:
        os.mkdir(saveDir)
    
    
    print folder + '\n'
    
    ####
    #compute eigenvectors
    ###
    Hmat = generate_layout_Hamiltonian(GenLayout, verbose = False, sparse = False)
    Es,Psis = scipy.linalg.eigh(Hmat)
    Eorder = numpy.argsort(Es)
    
        
    fig = pylab.figure(figNum)
    #######################################################
    for eigNum in range(startInd, GenLayout.coords.shape[0]):
        print eigNum
        
        eigAmps = Psis[:,Eorder[eigNum]]
        
        
        
        fig.clf()
#        fig.figsize = (100,100)
        ax1 = pylab.subplot(1,2,1)
        xs = scipy.arange(0,GenLayout.coords.shape[0],1)
        pylab.plot(xs, Es[Eorder], 'b.')
        
        pylab.plot(xs[eigNum],Es[Eorder[eigNum]], color = 'firebrick' , marker = '.', markersize = '10' )
        pylab.title(' full tight-binding Spectrum')
        ax1.set_ylim([-3.1, 3.1])
        
#        ax1 = pylab.subplot(1,2,2)
#        titleStr = 'eigenvector weight : ' + str(eigNum) + ' ; E = ' + str(Es[Eorder[eigNum]])
###        GenLayout.plot_layout_state(eigAmps, ax1, title = titleStr, colorbar = True, plot_links = True, cmap = 'Wistia')
#        GenLayout.plot_layout_state(eigAmps, ax1, title = titleStr, colorbar = True, plot_links = False, cmap = 'Wistia')
#        GenLayout.draw_SDlinks(ax1, color = 'firebrick', linewidth = 0.5, minus_links = True, minus_color = 'dodgerblue', NaNs = True)
#        
#        
#        ax1.set_aspect('auto')
        
        
        
        ax1 = pylab.subplot(1,2,2)
        titleStr = 'layout eigenvector weight : ' + str(eigNum) + ' ; E = ' + str(Es[Eorder[eigNum]])
        #####
        #plot the eigenvector
        #####
#        GenLayout.plot_layout_state(eigAmps, ax1, title = titleStr, colorbar = True, plot_links = False, cmap = 'Wistia')
        Amps = eigAmps
        Probs = numpy.abs(Amps)**2
        mSizes = Probs * len(Probs)*30
        mColors = numpy.angle(Amps)/numpy.pi
        cm = pylab.cm.get_cmap('Wistia')
        pylab.scatter(GenLayout.coords[:,0], GenLayout.coords[:,1],c =  mColors, s = mSizes, marker = 'o', edgecolors = 'k', cmap = cm, vmin = -0.5, vmax = 1.5, zorder = 1)
#        if colorbar:
#            cbar = pylab.colorbar(fraction=0.046, pad=0.04)
#            cbar.set_label('phase (pi radians)', rotation=270)
        pylab.title(titleStr, fontsize=8)


#        GenLayout.draw_SDlinks(ax1, color = 'firebrick', linewidth = 0.5, minus_links = True, minus_color = 'dodgerblue', NaNs = True)
        GenLayout.draw_resonator_lattice(ax1, color = 'firebrick', linewidth = 0.5)
        ax1.set_aspect('auto')
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
#        ax1.set_aspect('equal')
        

        fig.set_size_inches(xsize, ysize)
        
        pylab.tight_layout()
    
        fig_name = GenLayout.name + '_' + GenLayout.modeType +  '_eig_' + str(eigNum) + '.png' 
        fig_path = os.path.join(saveDir, fig_name)
        pylab.savefig(fig_path)
        
        if eigNum == breakPoint:
            break
        
    pylab.show()
    
    return









#saveLots = True
saveLots = False

#showCompiled = True
showCompiled = False


if saveLots:
#    #hexagonal one
#    makeMinLayoutPickle_hyperbolic(gon = 6,vertex = 3,modeType = 'FW', maxItter = 12)
#    makeMinLayoutPickle_hyperbolic(gon = 6,vertex = 3,modeType = 'FW', maxItter = 24)
#    
#    #hyperbolics
#    makeMinLayoutPickle_hyperbolic(gon = 7,vertex = 3,modeType = 'FW', maxItter = 6)
#    
#    makeMinLayoutPickle_hyperbolic(gon = 8,vertex = 3,modeType = 'FW', maxItter = 3)
#    makeMinLayoutPickle_hyperbolic(gon = 8,vertex = 3,modeType = 'FW', maxItter = 4)
#    makeMinLayoutPickle_hyperbolic(gon = 8,vertex = 3,modeType = 'FW', maxItter = 5)
#    
#    makeMinLayoutPickle_hyperbolic(gon = 9,vertex = 3,modeType = 'FW', maxItter = 3)
#    makeMinLayoutPickle_hyperbolic(gon = 9,vertex = 3,modeType = 'FW', maxItter = 4)
#    
#    makeMinLayoutPickle_hyperbolic(gon = 12,vertex = 3,modeType = 'FW', maxItter = 3)
#    makeMinLayoutPickle_hyperbolic(gon = 14,vertex = 3,modeType = 'FW', maxItter = 3)
    
#    print 'starting the simulations'
#    t0 = time.time()
#    makeMinLayoutPickle_hyperbolic(gon = 8,vertex = 7,modeType = 'FW', maxItter = 3)
#    t1 = time.time()
#    print 'octagon_7 elapsed time = ' + str(t1-t0)
#    t0 = t1
#    makeMinLayoutPickle_hyperbolic(gon = 14,vertex = 4,modeType = 'FW', maxItter = 3)
#    t1 = time.time()
#    print '14-gon_4 elapsed time = ' + str(t1-t0)
#    t0 = t1
#    makeMinLayoutPickle_hyperbolic(gon = 8,vertex = 6,modeType = 'FW', maxItter = 3)
#    t1 = time.time()
#    print 'octagon_6 elapsed time = ' + str(t1-t0)
#    t0 = t1
    
#    res = TreeResonators(degree = 3, iterations = 3).resonators
#    layout = GeneralLayout(resonators = res, resonatorsOnly = True)
    
    pass
    






########################################
#compile stuff
########



#gon = 7
#itt = 3
#sparse = False
#temp = PlanarLayout(gon = gon, vertex = 3)
#temp.populate(itt, resonatorsOnly = True)
#res = temp.get_all_resonators()
#nameStr = str(gon) + '3_' + str(itt)
#layout = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr)
#Es = makeMinLayoutPickle_gen(layout, verbose = True, sparse = sparse)



#########make a tree
#deg = 3
#itt = 4
#res = TreeResonators(degree = deg, iterations = itt).resonators
#nameStr = str(deg) + 'tree_' + str(itt)
#layout = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr)
#Es = makeMinLayoutPickle_gen(layout, verbose = False)

#
#deg = 5
#itt = 4
#sparse = True
#res = TreeResonators(degree = deg, iterations = itt).resonators
#nameStr = str(deg) + 'tree_' + str(itt)
#layout = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr)
#Es = makeMinLayoutPickle_gen(layout, verbose = False, sparse = sparse)

#deg = 6
#itt = 5
#sparse = False
#res = TreeResonators(degree = deg, iterations = itt).resonators
#nameStr = str(deg) + 'tree_' + str(itt)
#layout = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr)
#Es = makeMinLayoutPickle_gen(layout, verbose = False, sparse = sparse)




##########make a hyperbolic
#gon = 7
#itt = 3
#sparse = True
#nameStr = str(gon) + '3_' + str(itt)
#Es = makeMinLayoutPickle_hyperbolic(gon = gon,vertex = 3,modeType = 'FW', maxItter = itt, sparse = sparse)
##Es = makeMinLayoutPickle_hyperbolic(gon = gon,vertex = 3,modeType = 'FW', maxItter = itt, sparse = sparse, saveFolder = '')


#####temp, i need this
#makeMinLayoutPickle_hyperbolic(gon = gon,vertex = 3,modeType = 'FW', maxItter = 1, sparse = sparse)

##########Mclauhlin
#deg = 3
#itt = 8
#sparse = False
#treeRes = TreeResonators(degree = deg, iterations = itt).resonators
#sRes = split_resonators(treeRes)
#res = generate_line_graph(sRes)
#nameStr = 'McLaughlin_' + str(itt)
#layout = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr, roundDepth = 6)
#Es = makeMinLayoutPickle_gen(layout, verbose = False, sparse = sparse)


########## trimmed Mclauhlin
#deg = 3
#itt = 3
#sparse = False
#treeRes = TreeResonators(degree = deg, iterations = itt, roundDepth = 6).resonators
#sRes = split_resonators(treeRes)
#
#xs = sRes[:,2]
#ys = sRes[:,3]
#mags1 = numpy.sqrt(xs**2 + ys**2)
#xs = sRes[:,0]
#ys = sRes[:,1]
#mags2 = numpy.sqrt(xs**2 + ys**2)
#outerRadius1 = numpy.max(mags1)
#outerRadius2 = numpy.max(mags2)
#outerRadius = numpy.max([outerRadius1, outerRadius2])
#cutInds1 = numpy.where(mags1 > 0.99*outerRadius)[0]
##cutInds2 = numpy.where(mags2 > 0.99*outerRadius)[0]
##cutInds = numpy.concatenate((cutInds1, cutInds2))
#cutInds = cutInds1
#sRes2 = numpy.copy(sRes)
#sRes2= numpy.delete(sRes2, cutInds, axis = 0)
#res = generate_line_graph(sRes2)
#nameStr = 'McLaughlin_trimmed_' + str(itt)
#res= numpy.round(res,6)
#layout = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr, roundDepth = 6)
#Es = makeMinLayoutPickle_gen(layout, verbose = False, sparse = sparse)
#####make_layoutGraph_vid(layout, startInd = 0, stopInd = -1, figNum = 1, xsize = 10.8, ysize = 5.7)



#######big ones to run later
#deg = 3
#itt = 12
#res = TreeResonators(degree = deg, iterations = itt).resonators
#nameStr = str(deg) + 'tree_' + str(itt)
#layout = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr)
#Es = makeMinLayoutPickle_gen(layout, verbose = False)
#
#deg = 4
#itt = 8
#res = TreeResonators(degree = deg, iterations = itt).resonators
#nameStr = str(deg) + 'tree_' + str(itt)
#layout = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr)
#Es = makeMinLayoutPickle_gen(layout, verbose = False)



#####timing test
#t0 = time.time()
#gon = 7
#itt = 4
#Es = makeMinLayoutPickle_hyperbolic(gon = gon,vertex = 3,modeType = 'FW', maxItter = itt, sparse= False)
#t1 =time.time()
#print 'non-sparse elapsed time = ' + str(t1-t0)
#t0 = t1
#gon = 7
#itt = 4
#Es = makeMinLayoutPickle_hyperbolic(gon = gon,vertex = 3,modeType = 'FW', maxItter = itt, sparse= True)
#t1 =time.time()
#print 'sparse elapsed time = ' + str(t1-t0)





########################################
#load stuff
######
if sparse:
    Es , name, pickleDict = loadMinLayoutPickle(nameStr + '_splayoutDOSminimal.pkl')
else:
    Es , name, pickleDict = loadMinLayoutPickle(nameStr + '_layoutDOSminimal.pkl')

##euclidean
#Es , name, pickleDict = loadMinLayoutPickle('6gon_3vertex_12_FW_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('6gon_3vertex_24_FW_layoutDOSminimal.pkl')

##hyperbolics
#Es , name, pickleDict = loadMinLayoutPickle('7gon_3vertex_3_FW_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('7gon_3vertex_6_FW_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('8gon_3vertex_5_FW_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('9gon_3vertex_4_FW_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('12gon_3vertex_3_FW_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('14gon_3vertex_3_FW_layoutDOSminimal.pkl')

#trees
#Es , name, pickleDict = loadMinLayoutPickle('3tree_9_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('3tree_10_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('3tree_11_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('4tree_7_layoutDOSminimal.pkl')
#Es , name, pickleDict = loadMinLayoutPickle('4tree_8_layoutDOSminimal.pkl')






print name + ': number of sites = ' + str(len(Es))

pylab.figure(1)
pylab.clf()
ax = pylab.subplot(1,1,1)
plot_DOS(ax, Es, freq_res = 0.04, freq_max = 4, color = 'mediumblue')
pylab.title(name)
pylab.show()


#pylab.figure(2)
#pylab.clf()
#ax = pylab.subplot(1,1,1)
##layout.draw_resonator_lattice(ax)
##layout.draw_resonator_end_points(ax)
#layout.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 2.5)
#layout.draw_resonator_end_points(ax, color = layoutCapColor, edgecolor = 'k',  marker = 'o' , size = smallCdefault, zorder = 5)
#ax.set_aspect('equal')
#ax.axis('off')
#pylab.show()







###########
##try to find lambda_max another way.
##########
#
#roundDepth1 = 6
#
##gon = 7
##itt = 8
##test = PlanarLayout(gon = gon, vertex = 3, side =1, radius_method = 'exp', modeType = 'FW')
##test.populate(maxItter = itt, resonatorsOnly = True)
###test = GeneralLayout(resonators = test.get_all_resonators(), resonatorsOnly = True, modeType = 'FW', name = nameStr, roundDepth = roundDepth1)
###Hfull = generate_layout_Hamiltonian(test, sparse = False, flags = 0)
###Es, Piss = scipy.linalg.eigh(Hfull)
##Es , name, pickleDict = loadMinLayoutPickle('7gon_3vertex_6_FW_layoutDOSminimal.pkl')
##lambda_max = numpy.max(Es)
##lambda_min = numpy.min(Es)
###tries  = len(test.SDx)
##tries = 10000
#
#gon = 8
#itt = 6
#test = PlanarLayout(gon = gon, vertex = 3, side =1, radius_method = 'exp', modeType = 'FW')
#test.populate(maxItter = itt, resonatorsOnly = True)
##Hfull = generate_layout_Hamiltonian(test, sparse = False, flags = 0)
##Es, Piss = scipy.linalg.eigh(Hfull)
#Es , name, pickleDict = loadMinLayoutPickle('8gon_3vertex_5_FW_layoutDOSminimal.pkl')
#lambda_max = numpy.max(Es)
#lambda_min = numpy.min(Es)
##tries  = len(test.SDx)
#tries = 10000
#
#
##deg = 3
##itt = 11
##res = TreeResonators(degree = deg, iterations = itt, roundDepth = roundDepth1).resonators
##nameStr = str(deg) + 'tree_' + str(itt)
##test = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr, roundDepth = roundDepth1)
##Es , name, pickleDict = loadMinLayoutPickle('3tree_11_layoutDOSminimal.pkl')
##lambda_max = numpy.max(Es)
##lambda_min = numpy.min(Es)
##tries = test.coords.shape[0]*2
#
#
#
#
#Hsparse = generate_layout_Hamiltonian(test, sparse = True, flags = 0, roundDepth = roundDepth1)
#Hsparse.sum_duplicates()
#v0 = numpy.ones(test.coords.shape[0])
#v0 = v0/numpy.linalg.norm(v0)
#
#for tind in range(0, tries):
#    v1 = Hsparse.dot(v0)
#    v0 = v1/numpy.linalg.norm(v1)
#vfinal = v1/numpy.linalg.norm(v1)
#
####test the result
#vprime = Hsparse.dot(vfinal)
#mag = numpy.linalg.norm(vprime) #roughly lambda
#overlap = numpy.dot(vprime, vfinal)/numpy.linalg.norm(vprime) #overlap of lambda *eigvector dot eigvector
#print '    '
#print '    '
#print 'output magnitude = ' + str(mag)
#print 'output eigenness = ' + str(overlap)
#print 'old lambda max = ' + str(lambda_max)
#
#
#print Hsparse.max(axis = None)







#pylab.figure(2)
#pylab.clf()
#ax = pylab.subplot(1,1,1)
##layout.draw_resonator_lattice(ax)
##layout.draw_resonator_end_points(ax)
#test.draw_resonator_lattice(ax, color = layoutLineColor, alpha = 1 , linewidth = 0.25)
#ax.set_aspect('equal')
#ax.axis('off')
#pylab.show()





############
##make layout state vids
#######
#
#deg = 3
#itt = 4
#res = TreeResonators(degree = deg, iterations = itt).resonators
#nameStr = str(deg) + 'tree_' + str(itt)
#layout = GeneralLayout(resonators = res, resonatorsOnly = True, modeType = 'FW', name = nameStr)
#xsize = 10.8
#ysize = 5.7
#
#
#
#gon = 8
#itt = 3
#test = PlanarLayout(gon = gon, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
#test.populate(maxItter = itt, resonatorsOnly = True)
#nameStr = str(test.gon) + 'gon_' + str(test.vertex) + 'vertex_' + str(test.itter)
#layout = GeneralLayout(resonators = test.get_all_resonators(), resonatorsOnly = True, modeType = 'FW', name = nameStr, roundDepth = 3)
#xsize = 10.8
#ysize = 5.7
#
#make_layoutGraph_vid(layout, startInd = 0, stopInd = 3, figNum = 1, xsize = xsize, ysize = ysize)
#make_layoutGraph_vid(layout, startInd = test.coords.shape[0]-3, stopInd = -11, figNum = 1, xsize = xsize, ysize = ysize)
#
##make_layoutGraph_vid(layout, startInd = 0, stopInd = 3, figNum = 1, xsize = xsize, ysize = ysize)
##make_layoutGraph_vid(layout, startInd = 0, stopInd = -1, figNum = 1, xsize = xsize, ysize = ysize)



xsize = 10.7
ysize = 7.4
kagomeCell = UnitCell('kagome') #graphene layout, kagome effective
kapheneCell = kagomeCell.split_cell(name = 'kaphene', splitIn = 2) #split graphene layout, kaphene effective
lineKapheneCell = kapheneCell.line_graph_cell(name = 'lineKaphene', resonatorsOnly = False) #kaphene layout, line kaffene effective
splitKapheneLayoutCell =  lineKapheneCell.split_cell(name = 'splitKapheneLay_CharonEff') #split kaphene layout, Charon effective

EuclidLayout = EuclideanLayout(4,4,lattice_type = splitKapheneLayoutCell.type, modeType = 'FW', initialCell = splitKapheneLayoutCell )
layout = GeneralLayout(resonators = EuclidLayout.get_all_resonators(), resonatorsOnly = True, modeType = 'FW', name = EuclidLayout.lattice_type, roundDepth = 3)
#make_layoutGraph_vid(layout, startInd = 0, stopInd = 3, figNum = 1, xsize = xsize, ysize = ysize)
#make_layoutGraph_vid(layout, startInd = 0, stopInd = -1, figNum = 1, xsize = xsize, ysize = ysize)






#showGon = 14
#showInds = [3]
##showGon = 12
##showInds = [3]
##showGon = 9
##showInds = [4,3]
##showGon = 8
##showInds = [5,4,3]
##showGon = 7
##showInds = [6,5,4,3]
##showGon = 6
##showInds = [24, 12, 8] #showInds = [24, 12, 8,7,6]
#
#largestOnly = True
##largestOnly = False
#
#
#if showCompiled:
#    colors = {1:'r', 2:'b', 3:'deepskyblue', 4:'goldenrod', 5:'indigo', 6: 'maroon', 7:'r', 8:'b', 9:'deepskyblue', 10:'goldenrod'}
#    colors = {1:'mediumblue', 2:'firebrick', 3:'goldenrod', 4:'deepskyblue', 5:'indigo', 6:'firebrick', 8:'b', 9:'deepskyblue', 10:'goldenrod'}
#    
#    pylab.figure(4)
#    pylab.clf()
#    ax = pylab.subplot(1,1,1)
#    ind = 1
#    for itt in showInds:
#        name = str(showGon) + 'gon_3vertex_' + str(itt) +'_FW_layoutDOSminimal.pkl'
#
#        testDict = pickle.load( open(name, "rb"))
#        Evals = testDict['Es']
#        Eorder = testDict['Eorder']
#        xaxis = numpy.linspace(0,1, len(Eorder))
#        
#        pylab.plot(xaxis, Evals[Eorder], color = colors[ind], marker = '.', linestyle = '', label = str(itt), alpha = 1)
#        ind = ind +1
#        if largestOnly:
#            break
#    pylab.xlabel('Eigenvalue number')
#    pylab.ylabel('Energy (|t|)')
#
#    ax.legend(loc = 'upper left')
#    titleStr = str(showGon) + ':layout energy spectra versus system size'
#    pylab.title(titleStr)
#    pylab.show()
#    
#    
#    
#    
#    #set up frequency sweep
#    freq_res = 0.04
##    freq_res = 0.01
#    freq_range = 3.0 + freq_res/2
#    
#    freqs = scipy.arange(-freq_range, freq_range, freq_res) + freq_res/2.
#    freq_bins = scipy.arange(-freq_range, freq_range+freq_res, freq_res)
#    
#    
#    pylab.figure(5)
#    pylab.clf()
#    ax = pylab.subplot(1,1,1)
#    ind = 1
#    for itt in showInds:
#        name = str(showGon) + 'gon_3vertex_' + str(itt) +'_FW_layoutDOSminimal.pkl'
#
#        testDict = pickle.load( open(name, "rb"))
#        Evals = testDict['Es']
#        Eorder = testDict['Eorder']
#        xaxis = numpy.linspace(0,1, len(Eorder))
#        
#        [DOS, bins_out] = numpy.histogram(Evals, freq_bins)
#        
#        bins_centers = (bins_out[0:-1] + bins_out[1:])/2
#        binWidth = bins_out[1] - bins_out[0]
#        
#        pylab.bar(bins_centers, 1.*DOS/len(Evals), width = binWidth, color =  colors[ind], label = str(itt), alpha = 0.8)
##        pylab.bar(bins_centers, 1.*DOS/1, width = binWidth, color =  colors[ind], label = str(itt), alpha = 0.8)
#        ind = ind+1
#        if largestOnly:
#            break
#    pylab.xlabel('Eigenvalue')
#    pylab.ylabel('DOS (AU)')
#
#    ax.legend(loc = 'upper left')
#    titleStr = str(showGon) + ':layout DOS versus system size'
#    pylab.title(titleStr)
#    pylab.show()
#    
#    
##    pylab.figure(6)
##    pylab.clf()
##    ax = pylab.subplot(1,1,1)
##
##    #from the layout
##    name = '7gon_3vertex_6_FW_layoutDOSminimal.pkl'
##
##    testDict = pickle.load( open(name, "rb"))
##    Evals = testDict['Es']
##    Eorder = testDict['Eorder']
##    xaxis = numpy.linspace(0,1, len(Eorder))
##    
##    [DOS, bins_out] = numpy.histogram(Evals, freq_bins)
##    bins_centers = (bins_out[0:-1] + bins_out[1:])/2
##    binWidth = bins_out[1] - bins_out[0]
##    
##    pylab.bar(bins_centers+1, 1.*DOS/len(Evals), width = binWidth, color =  'mediumblue', label = 'layout+1', alpha = 0.8)
##    
##    #from the effective
##    testDict = pickle.load( open('7gon_3vertex_6_FW_DOSminimal.pkl', "rb"))
##    Evals = testDict['Es']
##    Eorder = testDict['Eorder']
##    xaxis = numpy.linspace(0,1, len(Eorder))
##    
##    freq_bins2 = scipy.arange(-freq_range+1, freq_range+1+freq_res, freq_res)
##    [DOS, bins_out] = numpy.histogram(Evals, freq_bins2)
##    bins_centers = (bins_out[0:-1] + bins_out[1:])/2
##    binWidth = bins_out[1] - bins_out[0]
##    
##    pylab.bar(bins_centers, 1.*DOS/len(Evals), width = binWidth, color =  'deepskyblue', label = 'effective', alpha = 0.8)
##
##    pylab.xlabel('Eigenvalue')
##    pylab.ylabel('DOS (AU)')
##    
##    ax.set_ylim([0,0.05])
##
##    ax.legend(loc = 'upper right')
##    titleStr = '7: layout v effective'
##    pylab.title(titleStr)
##    pylab.show()
##
##    pylab.figure(7)
##    pylab.clf()
##    ax = pylab.subplot(1,1,1)
##
##    #from the layout
##    name = '8gon_3vertex_5_FW_layoutDOSminimal.pkl'
##
##    testDict = pickle.load( open(name, "rb"))
##    Evals = testDict['Es']
##    Eorder = testDict['Eorder']
##    xaxis = numpy.linspace(0,1, len(Eorder))
##    
##    [DOS, bins_out] = numpy.histogram(Evals, freq_bins)
##    bins_centers = (bins_out[0:-1] + bins_out[1:])/2
##    binWidth = bins_out[1] - bins_out[0]
##    
##    pylab.bar(bins_centers+1, 1.*DOS/len(Evals), width = binWidth, color =  'mediumblue', label = 'layout+1', alpha = 0.8)
##    
##    #from the effective
##    testDict = pickle.load( open('8gon_3vertex_5_FW_DOSminimal.pkl', "rb"))
##    Evals = testDict['Es']
##    Eorder = testDict['Eorder']
##    xaxis = numpy.linspace(0,1, len(Eorder))
##    
##    freq_bins2 = scipy.arange(-freq_range+1, freq_range+1+freq_res, freq_res)
##    [DOS, bins_out] = numpy.histogram(Evals, freq_bins2)
##    bins_centers = (bins_out[0:-1] + bins_out[1:])/2
##    binWidth = bins_out[1] - bins_out[0]
##    
##    pylab.bar(bins_centers, 1.*DOS/len(Evals), width = binWidth, color =  'deepskyblue', label = 'effective', alpha = 0.8)
##
##    pylab.xlabel('Eigenvalue')
##    pylab.ylabel('DOS (AU)')
##    
##    ax.set_ylim([0,0.05])
##
##    ax.legend(loc = 'upper right')
##    titleStr = '8: layout v effective'
##    pylab.title(titleStr)
##    pylab.show()





























