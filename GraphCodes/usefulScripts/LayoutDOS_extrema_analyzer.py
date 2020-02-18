#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:57:25 2018

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


#def generate_layout_Hamiltonian(layout, roundDepth = 3, t = 1, verbose = True, sparse = True, flags = 0):
#    '''
#    custom function so I can get vertex dict without having to run the full populate of general layout
#    and thereby having to also diagonalize the effective Hamiltonian.
#    '''
#    resonators = layout.get_all_resonators()
#    resonators = numpy.round(resonators, roundDepth)
#    
#    numVerts = layout.coords.shape[0]
#    if sparse:
#        rowVec = numpy.zeros(numVerts*4+flags)
#        colVec = numpy.zeros(numVerts*4+flags)
#        Hvec = numpy.zeros(numVerts*4+flags)
#    else:
#        Hmat = numpy.zeros((numVerts, numVerts))
#    
#    coords_complex = numpy.round(layout.coords[:,0] + 1j*layout.coords[:,1], roundDepth)
#    
#    currInd = 0
#    for rind in range(0, resonators.shape[0]):
#        resPos = resonators[rind,:]
#        startPos = numpy.round(resPos[0],roundDepth)+ 1j*numpy.round(resPos[1],roundDepth)
#        stopPos = numpy.round(resPos[2],roundDepth)+ 1j*numpy.round(resPos[3],roundDepth)
#        
#        startInd = numpy.where(startPos == coords_complex)[0][0]
#        stopInd = numpy.where(stopPos == coords_complex)[0][0]
#
#        if sparse:
#            rowVec[currInd] = startInd
#            colVec[currInd] = stopInd
#            Hvec[currInd] = t #will end up adding t towhatever this entry was before.
#            currInd = currInd +1
#            
#            rowVec[currInd] = stopInd
#            colVec[currInd] = startInd
#            Hvec[currInd] = t #will end up adding t towhatever this entry was before.
#            currInd = currInd +1
#            
#        else:
#            Hmat[startInd, stopInd] = Hmat[startInd, stopInd] + t
#            Hmat[stopInd, startInd] = Hmat[stopInd, startInd] + t
#    
#    #finish making the sparse matrix if we are in sparse matrix mode.
#    if sparse:
#        #pad the end of the matrix with values so that I can see if one of those is the missing one
#        for ind in range(0, flags):
##            rowVec[currInd] = currInd
##            colVec[currInd] = currInd
#            rowVec[currInd] = numVerts+ind
#            colVec[currInd] = numVerts+ind
#            Hvec[currInd] = -7.5 #will end up adding t towhatever this entry was before.
#            currInd = currInd +1
#
##        Hmat = coo_matrix((Hvec,(rowVec,colVec)), shape = (currInd,currInd), dtype = 'd')
#        Hmat = coo_matrix((Hvec,(rowVec,colVec)), shape = (numVerts+flags,numVerts+flags), dtype = 'd')
#        Hmat.eliminate_zeros() #removed the unused spots since this is not a regular graph
#        
#    if verbose:
#        temp = numpy.sum(Hmat)/numVerts
#        print 'average degree = ' + str(temp)
#        
#    return Hmat



#def compile_lambda_extrema(maxItt, gon = 7, tree = False, degree = 3, roundDepth = 6, save = True, saveFolder = ''):
#    lambdaMax_vals = numpy.zeros(maxItt)
#    lambdaMin_vals = numpy.zeros(maxItt)
#    
#    sizes = numpy.zeros(maxItt)
#    
#    
#    
#    if tree:
#        stepSize = 2
#        itts = scipy.arange(1,maxItt+1, 1)*stepSize #go up ins steps of 3 generations at a time
#        maxItt = maxItt*stepSize
#        
#        for ind in range(0, len(itts)):
#            itt = itts[ind]
#            base = TreeResonators(degree = degree, iterations = itt)
#            baseName = '3tree'
#            print itt
#            
#            nameStr = baseName  + '_' + str(itt)
#            test = GeneralLayout(resonators = base.get_all_resonators(), resonatorsOnly = True, modeType = 'FW', name = nameStr, roundDepth = roundDepth)
#            
#            Hsparse = generate_layout_Hamiltonian(test, roundDepth = roundDepth, t = 1, verbose = False, sparse = True, flags = 0)
#            
#            almostMaxs = scipy.sparse.linalg.eigsh(Hsparse, sigma = 3, return_eigenvectors = False)
#            almostMins = scipy.sparse.linalg.eigsh(Hsparse, sigma = -3, return_eigenvectors = False)
#            
#            lambdaMax_vals[ind] = numpy.max(almostMaxs)
#            lambdaMin_vals[ind] = numpy.min(almostMins)
#            sizes[ind] = test.coords.shape[0]
#        
#    else:
#        itts = scipy.arange(1,maxItt+1, 1)
#        
#        base = PlanarLayout(gon = gon, vertex = degree, side =1, radius_method = 'exp', modeType = 'FW')
#        base.populate(maxItter = 1, resonatorsOnly = True)   
#        baseName = str(gon) + 'gon_' + str(degree) + 'vertex'
#        for itt in itts:
#            if itt > 1:
#                base.itter_generate()
#                base.generate_semiduals()
#                print itt
#            
#            nameStr = baseName  + '_' + str(itt)
#            test = GeneralLayout(resonators = base.get_all_resonators(), resonatorsOnly = True, modeType = 'FW', name = nameStr, roundDepth = roundDepth)
#            
#            Hsparse = generate_layout_Hamiltonian(test, roundDepth = roundDepth, t = 1, verbose = False, sparse = True, flags = 0)
#            
#            almostMaxs = scipy.sparse.linalg.eigsh(Hsparse, sigma = 3, return_eigenvectors = False)
#            almostMins = scipy.sparse.linalg.eigsh(Hsparse, sigma = -3, return_eigenvectors = False)
#            
#            lambdaMax_vals[itt-1] = numpy.max(almostMaxs)
#            lambdaMin_vals[itt-1] = numpy.min(almostMins)
#            sizes[itt-1] = test.coords.shape[0]
#        
#       
#    pickleDict = {}
#    pickleDict['itts'] = itts
#    pickleDict['lambdaMax'] = lambdaMax_vals
#    pickleDict['lambdaMin'] = lambdaMin_vals
#    pickleDict['name'] = baseName
#    pickleDict['systemSizes'] = sizes
#    
#    if save:
#        saveName = baseName + '_layoutLambdaExtrema_' + str(maxItt) + '.pkl'
#        save_path = os.path.join(saveFolder, saveName)
#        pickle.dump(pickleDict, open(save_path, 'wb'))
#    
#    
#    fig1 = pylab.figure(1)
#    pylab.clf()
#    ax = pylab.subplot(1,2,1)
#    pylab.plot(itts, lambdaMax_vals, 'b.-', label = 'sparse')
##    pylab.plot(itts_full[2:], lambdaMax_vals_full[2:], color = 'deepskyblue', linestyle = '-', marker = '.', label = 'full')
#    ax.legend(loc = 'lower right')
#    pylab.title('$\lambda_{max}$')
#    pylab.xlabel('itterations')
#    
#    ax = pylab.subplot(1,2,2)
#    pylab.plot(itts, lambdaMin_vals, 'r.-', label = 'sparse')
##    pylab.plot(itts_full[2:], lambdaMin_vals_full[2:], color = 'deepskyblue', linestyle = '-', marker = '.', label = 'full')
#    ax.legend(loc = 'upper right')
#    pylab.title('$\lambda_{min}$')
#    pylab.xlabel('itterations')
#    
#    pylab.suptitle(baseName)
#    
#    fig1.set_size_inches([9.7,5.1])
#    pylab.show()
#    
#    if save:
#        fig_saveName = baseName + '_layoutLambdaExtremaFig_' + str(maxItt) + '.png'
#        fig_save_path = os.path.join(saveFolder, fig_saveName)
#        fig1.savefig(fig_save_path,transparent= False, dpi = 200)
#        
#    return lambdaMax_vals, lambdaMin_vals, itts






    










plotVSysSize = True






saveFolder = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/layoutDOSs'
#######################################################################  
from scipy.optimize import curve_fit

def lorish_ff(x, offset, width, scaling, power):
    vals = offset - 1./(width +  (1.*x/scaling)**power)
    return vals

def lor_ff(x, offset, width, scaling):
    vals = offset - 1./(width +  (1.*x/scaling)**2)
    return vals

def lorish_bottom_ff(x, offset, width, scaling, power):
    vals = offset + 1./(width +  (1.*x/scaling)**power)
    return vals


def lorish_ff_V(V, offset, width, scaling, power):
    vals = offset - 1./(width +  (1.*numpy.log(V)/scaling)**power)
    return vals 

def lorish_bottom_ff_V(V, offset, width, scaling, power):
    vals = offset + 1./(width +  (1.*numpy.log(V)/scaling)**power)
    return vals 
    

#load up and old one
pickleName = '3tree_layoutLambdaExtrema_18.pkl'
filePath = os.path.join(saveFolder, pickleName)

pickleDict = pickle.load(open(filePath, "rb" ) )
lmax = pickleDict['lambdaMax']
lmin = pickleDict['lambdaMin']
itts = pickleDict['itts']
baseName = pickleDict['name']
sizes = pickleDict['systemSizes']

##        guess1 = 2.8 - numpy.exp(-1.*itts/3)
#        guess1 = 2.78 - numpy.exp(-1.*itts/2.7)
#        
#        
#        guess2 = 2.8 - 1./itts**1.2
#        guess3 = 2.8 - numpy.exp(-(1.*itts/4)**2)
##        guess4 = 2.8 - numpy.exp(-1.*itts*numpy.log(itts)/5)
#        guess4 = 2.8 - numpy.exp(-1.*itts/3)*numpy.log(itts+2)




#############################################################
fig1 = pylab.figure(1)
pylab.clf()


#####
#fit lambda_max
#####
init_guess = [2*numpy.sqrt(2), 2.5, 1.0, 1.9]
fit_out, pcov = curve_fit(lorish_ff, itts, lmax, p0 = init_guess)
guessedLambdas = lorish_ff(itts, init_guess[0], init_guess[1], init_guess[2], init_guess[3])
fittedLambdas = lorish_ff(itts, fit_out[0], fit_out[1], fit_out[2], fit_out[3])

print '   '
print baseName + ' upper:' 
print 'theoretical = ' + str(2*numpy.sqrt(2))
print 'fitted assymptotic = ' + str(fit_out[0])
print 'fitted width = ' + str(fit_out[1])
print 'fitted scaling = ' + str(fit_out[2])
print 'fitted power = ' + str(fit_out[3])


###alternate fit function
#init_guess = [2*numpy.sqrt(2), 2.5, 1.0]
#fit_out, pcov = curve_fit(lor_ff, itts, lmax, p0 = init_guess)
#guessedLambdas = lor_ff(itts, init_guess[0], init_guess[1], init_guess[2])
#fittedLambdas = lor_ff(itts, fit_out[0], fit_out[1], fit_out[2])
#
#print '   '
#print baseName + ' upper:' 
#print 'theoretical = ' + str(2*numpy.sqrt(2))
#print 'fitted assymptotic = ' + str(fit_out[0])
#print 'fitted width = ' + str(fit_out[1])
#print 'fitted scaling = ' + str(fit_out[2])
##print 'fitted power = ' + str(fit_out[3])



ax = pylab.subplot(1,2,1)
pylab.plot(itts, lmax, 'b.-', label = 'sparse')
pylab.plot(itts, guessedLambdas, color = 'firebrick', linestyle = '-', marker = '', label = 'init_guess', linewidth = 0.5)
pylab.plot(itts, fittedLambdas, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'fit_out')
ax.legend(loc = 'lower right')
pylab.title('$\lambda_{max}$')
pylab.xlabel('itterations')

#        pylab.plot(itts, guess1, color = 'deepskyblue', linestyle = '-', marker = '.', label = 'exp')
#        pylab.plot(itts, guess2, color = 'darkgoldenrod', linestyle = '-', marker = '.', label = 'power law')
#        pylab.plot(itts, guess3, color = 'firebrick', linestyle = '-', marker = '.', label = 'gaussian')
#        pylab.plot(itts, guess4, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'gaussian-log')
#        pylab.plot(itts_full[minItt:], lambdaMax_vals_full[minItt:], color = 'deepskyblue', linestyle = '-', marker = '.', label = 'full')


#####
#fit lambda_min
#####
init_guess = [-2*numpy.sqrt(2), 2.5, 1.0, 1.9]
fit_out, pcov = curve_fit(lorish_bottom_ff, itts, lmin, p0 = init_guess)
guessedLambdas = lorish_bottom_ff(itts, init_guess[0], init_guess[1], init_guess[2], init_guess[3])
fittedLambdas = lorish_bottom_ff(itts, fit_out[0], fit_out[1], fit_out[2], fit_out[3])

print '   '
print baseName + ' lower:' 
print 'theoretical = ' + str(2*numpy.sqrt(2))
print 'fitted assymptotic = ' + str(fit_out[0])
print 'fitted width = ' + str(fit_out[1])
print 'fitted scaling = ' + str(fit_out[2])
print 'fitted power = ' + str(fit_out[3])

ax = pylab.subplot(1,2,2)
pylab.plot(itts, lmin, 'r.-', label = 'sparse')
pylab.plot(itts, guessedLambdas, color = 'firebrick', linestyle = '-', marker = '', label = 'init_guess', linewidth = 0.5)
pylab.plot(itts, fittedLambdas, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'fit_out')
ax.legend(loc = 'upper right')
pylab.title('$\lambda_{min}$')
pylab.xlabel('itterations')
#        pylab.plot(itts_full[minItt:], lambdaMin_vals_full[minItt:], color = 'deepskyblue', linestyle = '-', marker = '.', label = 'full')


pylab.suptitle(baseName)

fig1.set_size_inches([9.7,5.1])
pylab.show()
#############################################################












#        pickleName = '7gon_3vertex_layoutLambdaExtrema_9.pkl'
pickleName = '7gon_3vertex_layoutLambdaExtrema_11.pkl'
filePath = os.path.join(saveFolder, pickleName)

pickleDict = pickle.load(open(filePath, "rb" ) )
lmax = pickleDict['lambdaMax']
lmin = pickleDict['lambdaMin']
itts = pickleDict['itts']
baseName = pickleDict['name']
sizes = pickleDict['systemSizes']

##        guess1 = 2.8 - numpy.exp(-1.*itts/3)
##        guess1 = 2.93 - numpy.exp(-1.*itts/5)
#        guess1 = 2.925 - numpy.exp(-1.*itts/1.5)
#        
#        
#        
#        guess2 =2.92- 1./(3*itts)**1.2
##        guess2 =2.92- 0.25/(itts)**2
#        guess2 =2.94- 0.37/(itts*numpy.log(itts))**1.0
#        
##        guess3 = 2.94 - 1/(2. +  (itts/1.0)**2)
#        guess3 = 2.94 - 1/(3. +  (itts/1.5)**2)
#        guess3 = 2.9375 - 1/(2.5 +  (itts/1.0)**1.9)
#        
#        
##        guess4 = 2.9 - numpy.exp(-1.*itts*numpy.log(itts)/5)
##        guess4 = 2.9 - numpy.exp(-1.*itts)*numpy.log(itts+1.2)
##        guess4 = 2.9 - numpy.exp(-1.*itts)*itts**0.5


#############################################################
fig3 = pylab.figure(3)
pylab.clf()

#####
#fit lambda_max
#####
if pickleName[0] == '7':
    init_guess = [2.9375, 2.5, 1.0, 1.9]
if pickleName[0] == '8':
    init_guess = [2.9375, 2.5, 1.0, 1.9]
fit_out, pcov = curve_fit(lorish_ff, itts, lmax, p0 = init_guess)
guessedLambdas = lorish_ff(itts, init_guess[0], init_guess[1], init_guess[2], init_guess[3])
fittedLambdas = lorish_ff(itts, fit_out[0], fit_out[1], fit_out[2], fit_out[3])

print '   '
print baseName + ' upper:' 
print 'fitted assymptotic = ' + str(fit_out[0])
print 'fitted width = ' + str(fit_out[1])
print 'fitted scaling = ' + str(fit_out[2])
print 'fitted power = ' + str(fit_out[3])

ax = pylab.subplot(1,2,1)
pylab.plot(itts, lmax, 'b.-', label = 'sparse')
pylab.plot(itts, guessedLambdas, color = 'firebrick', linestyle = '-', marker = '', label = 'init_guess', linewidth = 0.5)
pylab.plot(itts, fittedLambdas, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'fit_out')
ax.legend(loc = 'lower right')
pylab.title('$\lambda_{max}$')
pylab.xlabel('itterations')

#        pylab.plot(itts, guess1, color = 'deepskyblue', linestyle = '-', marker = '.', label = 'exp')
#        pylab.plot(itts, guess2, color = 'darkgoldenrod', linestyle = '-', marker = '.', label = '1/(r*log(r))')
#        pylab.plot(itts, guess3, color = 'firebrick', linestyle = '-', marker = '.', label = 'other')
#        pylab.plot(itts, guess4, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'exp-log')


#####
#fit lambda_min
#####
if pickleName[0] == '7':
    init_guess = [-2.57, 2.5, 1.0, 1.9]
if pickleName[0] == '8':
    init_guess = [-2.9375, 2.5, 1.0, 1.9]
fit_out, pcov = curve_fit(lorish_bottom_ff, itts, lmin, p0 = init_guess)
guessedLambdas = lorish_bottom_ff(itts, init_guess[0], init_guess[1], init_guess[2], init_guess[3])
fittedLambdas = lorish_bottom_ff(itts, fit_out[0], fit_out[1], fit_out[2], fit_out[3])

print '   '
print baseName + ' lower:' 
print 'fitted assymptotic = ' + str(fit_out[0])
print 'fitted width = ' + str(fit_out[1])
print 'fitted scaling = ' + str(fit_out[2])
print 'fitted power = ' + str(fit_out[3])

ax = pylab.subplot(1,2,2)
pylab.plot(itts, lmin, 'r.-', label = 'sparse')
pylab.plot(itts, guessedLambdas, color = 'firebrick', linestyle = '-', marker = '', label = 'init_guess', linewidth = 0.5)
pylab.plot(itts, fittedLambdas, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'fit_out')
ax.legend(loc = 'upper right')
pylab.title('$\lambda_{min}$')
pylab.xlabel('itterations')

pylab.suptitle(baseName)

fig3.set_size_inches([9.7,5.1])
pylab.show()
#############################################################








pickleName = '8gon_3vertex_layoutLambdaExtrema_8.pkl'
filePath = os.path.join(saveFolder, pickleName)

pickleDict = pickle.load(open(filePath, "rb" ) )
lmax = pickleDict['lambdaMax']
lmin = pickleDict['lambdaMin']
itts = pickleDict['itts']
baseName = pickleDict['name']
sizes = pickleDict['systemSizes']



#############################################################
fig5 = pylab.figure(5)
pylab.clf()

#####
#fit lambda_max
#####
if pickleName[0] == '7':
    init_guess = [2.9375, 2.5, 1.0, 1.9]
if pickleName[0] == '8':
    init_guess = [2.9375, 2.5, 1.0, 1.9]
fit_out, pcov = curve_fit(lorish_ff, itts, lmax, p0 = init_guess)
guessedLambdas = lorish_ff(itts, init_guess[0], init_guess[1], init_guess[2], init_guess[3])
fittedLambdas = lorish_ff(itts, fit_out[0], fit_out[1], fit_out[2], fit_out[3])

print '   '
print baseName + ' upper:' 
print 'fitted assymptotic = ' + str(fit_out[0])
print 'fitted width = ' + str(fit_out[1])
print 'fitted scaling = ' + str(fit_out[2])
print 'fitted power = ' + str(fit_out[3])

ax = pylab.subplot(1,2,1)
pylab.plot(itts, lmax, 'b.-', label = 'sparse')
pylab.plot(itts, guessedLambdas, color = 'firebrick', linestyle = '-', marker = '', label = 'init_guess', linewidth = 0.5)
pylab.plot(itts, fittedLambdas, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'fit_out')
ax.legend(loc = 'lower right')
pylab.title('$\lambda_{max}$')
pylab.xlabel('itterations')


#####
#fit lambda_min
#####
if pickleName[0] == '7':
    init_guess = [-2.57, 2.5, 1.0, 1.9]
if pickleName[0] == '8':
    init_guess = [-2.9375, 2.5, 1.0, 1.9]
fit_out, pcov = curve_fit(lorish_bottom_ff, itts, lmin, p0 = init_guess)
guessedLambdas = lorish_bottom_ff(itts, init_guess[0], init_guess[1], init_guess[2], init_guess[3])
fittedLambdas = lorish_bottom_ff(itts, fit_out[0], fit_out[1], fit_out[2], fit_out[3])

print '   '
print baseName + ' lower:' 
print 'fitted assymptotic = ' + str(fit_out[0])
print 'fitted width = ' + str(fit_out[1])
print 'fitted scaling = ' + str(fit_out[2])
print 'fitted power = ' + str(fit_out[3])

ax = pylab.subplot(1,2,2)
pylab.plot(itts, lmin, 'r.-', label = 'sparse')
pylab.plot(itts, guessedLambdas, color = 'firebrick', linestyle = '-', marker = '', label = 'init_guess', linewidth = 0.5)
pylab.plot(itts, fittedLambdas, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'fit_out')
ax.legend(loc = 'upper right')
pylab.title('$\lambda_{min}$')
pylab.xlabel('itterations')

pylab.suptitle(baseName)

fig5.set_size_inches([9.7,5.1])
pylab.show()
#############################################################














if plotVSysSize:
    #############################################################
#    fig2 = pylab.figure(2)
#    pylab.clf()
#    ax = pylab.subplot(1,2,1)
#    pylab.plot(sizes, lmax, 'b.-', label = 'sparse')
#    ax.legend(loc = 'lower right')
#    pylab.title('$\lambda_{max}$')
#    pylab.xlabel('|V|')
#    
#    ax = pylab.subplot(1,2,2)
#    pylab.plot(sizes, lmin, 'r.-', label = 'sparse')
#    ax.legend(loc = 'upper right')
#    pylab.title('$\lambda_{min}$')
#    pylab.xlabel('|V|')
#    
#    pylab.suptitle(baseName)
#    
#    fig2.set_size_inches([9.7,5.1])
#    pylab.show()
#    
    #############################################################
    
    #        pickleName = '7gon_3vertex_layoutLambdaExtrema_9.pkl'
    pickleName = '7gon_3vertex_layoutLambdaExtrema_11.pkl'
    filePath = os.path.join(saveFolder, pickleName)
    
    pickleDict = pickle.load(open(filePath, "rb" ) )
    lmax = pickleDict['lambdaMax']
    lmin = pickleDict['lambdaMin']
    itts = pickleDict['itts']
    baseName = pickleDict['name']
    sizes = pickleDict['systemSizes']
    
    fig2 = pylab.figure(2)
    pylab.clf()
    
    
    #####
    #fit lambda_max
    #####
    init_guess = [2.9525, 2.5, 2., 1.9]
    fit_out, pcov = curve_fit(lorish_ff_V, sizes, lmax, p0 = init_guess)
    guessedLambdas = lorish_ff_V(sizes, init_guess[0], init_guess[1], init_guess[2], init_guess[3])
    fittedLambdas = lorish_ff_V(sizes, fit_out[0], fit_out[1], fit_out[2], fit_out[3])
    
    print '   '
    print baseName + ' v size upper:' 
    print 'fitted assymptotic = ' + str(fit_out[0])
    print 'fitted width = ' + str(fit_out[1])
    print 'fitted scaling = ' + str(fit_out[2])
    print 'fitted power = ' + str(fit_out[3])
    

    
    
    ax = pylab.subplot(1,2,1)
    pylab.plot(sizes, lmax, 'b.-', label = 'sparse')
    pylab.plot(sizes, guessedLambdas, color = 'firebrick', linestyle = '-', marker = '', label = 'init_guess', linewidth = 0.5)
    pylab.plot(sizes, fittedLambdas, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'fit_out')
    ax.legend(loc = 'lower right')
    pylab.title('$\lambda_{max}$')
    pylab.xlabel('|V|')
    
    
    #####
    #fit lambda_min
    #####
#    init_guess = [-2.575, 2.5, 1.0, 1.9]
    init_guess = [-2.6, 2.5, 2., 1.9]
    fit_out, pcov = curve_fit(lorish_bottom_ff_V, sizes, lmin, p0 = init_guess)
    guessedLambdas = lorish_bottom_ff_V(sizes, init_guess[0], init_guess[1], init_guess[2], init_guess[3])
    fittedLambdas = lorish_bottom_ff_V(sizes, fit_out[0], fit_out[1], fit_out[2], fit_out[3])
    
    print '   '
    print baseName + ' v size lower:' 
    print 'fitted assymptotic = ' + str(fit_out[0])
    print 'fitted width = ' + str(fit_out[1])
    print 'fitted scaling = ' + str(fit_out[2])
    print 'fitted power = ' + str(fit_out[3])
    
    ax = pylab.subplot(1,2,2)
    pylab.plot(sizes, lmin, 'r.-', label = 'sparse')
    pylab.plot(sizes, guessedLambdas, color = 'firebrick', linestyle = '-', marker = '', label = 'init_guess', linewidth = 0.5)
    pylab.plot(sizes, fittedLambdas, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'fit_out')
    ax.legend(loc = 'upper right')
    pylab.title('$\lambda_{min}$')
    pylab.xlabel('|V|')
    
    pylab.suptitle(baseName)
    
    fig1.set_size_inches([9.7,5.1])
    pylab.show()
    #############################################################
    
    
    
    #############################################################
    
    
    
    
    #        guess1 = 2.91 - numpy.exp(-1.*sizes/1000.)
    #        guess2 = 2.9 - 1./(3*sizes)**1.2
    #        guess3 = 2.9 - numpy.exp(-(1.*sizes)**2)
    ##        guess4 = 2.9 - numpy.exp(-1.*itts*numpy.log(itts)/5)
    #        guess4 = 2.91 - numpy.exp(-1.*itts)*numpy.log(itts+1.
    
    #        guess1 = 2.91 - numpy.exp(-(1.*(sizes+1000)/10000)**3.)
    
    
    
    #############################################################
    fig4 = pylab.figure(4)
    pylab.clf()
    ax = pylab.subplot(1,2,1)
    pylab.plot(sizes, lmax, 'b.-', label = 'sparse')
    #        pylab.plot(sizes, guess1, color = 'deepskyblue', linestyle = '-', marker = '.', label = 'testing')
    #        pylab.plot(itts, guess2, color = 'darkgoldenrod', linestyle = '-', marker = '.', label = 'power law')
    #        pylab.plot(itts, guess3, color = 'firebrick', linestyle = '-', marker = '.', label = 'gaussian')
    #        pylab.plot(itts, guess4, color = 'dodgerblue', linestyle = '-', marker = '.', label = 'gaussian')
    ax.legend(loc = 'lower right')
    pylab.title('$\lambda_{max}$')
    pylab.xlabel('|V|')
    
    ax = pylab.subplot(1,2,2)
    pylab.plot(sizes, lmin, 'r.-', label = 'sparse')
    #        pylab.plot(numpy.log(sizes), lmin, 'r.-', label = 'sparse')
    ax.legend(loc = 'upper right')
    pylab.title('$\lambda_{min}$')
    pylab.xlabel('|V|')
    
    pylab.suptitle(baseName)
    
    fig4.set_size_inches([9.7,5.1])
    pylab.show()
    #############################################################





        
    











