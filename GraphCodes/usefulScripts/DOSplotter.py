#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:31:10 2018

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



hyperbolicFolderPath = r'/Users/kollar2/Documents/HouckLab/HyperbolicPlanning/'
euclideanFolderPath = r'/Users/kollar2/Documents/HouckLab/EuclideanLatticePlanning/'
generalFolderPath = r'/Users/kollar2/Documents/HouckLab/GeneralLayoutCode/'
if not hyperbolicFolderPath in sys.path:
    sys.path.append(hyperbolicFolderPath)
if not euclideanFolderPath in sys.path:
    sys.path.append(euclideanFolderPath)
if not generalFolderPath in sys.path:
    sys.path.append(generalFolderPath)

from EuclideanLayoutGenerator2 import UnitCell
from EuclideanLayoutGenerator2 import EuclideanLayout

from LayoutGenerator5 import PlanarLayout


from GeneralLayoutGenerator import GeneralLayout
from GeneralLayoutGenerator import TreeResonators

from GeneralLayoutGenerator import split_resonators
from GeneralLayoutGenerator import generate_line_graph
from GeneralLayoutGenerator import decorate_layout



def makeMinPickle(gon,vertex,modeType, maxItter):
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






#showRegularized = True
#showHardWall = True
showRegularized = True
showHardWall = True



##makeMinPickle(gon = 8, vertex = 3, modeType = 'FW', maxItter = 4)
#makeMinPickle(gon = 8, vertex = 3, modeType = 'FW', maxItter = 5)
#makeMinPickle(gon = 8, vertex = 3, modeType = 'HW', maxItter = 5)
#
#makeMinPickle(gon = 7, vertex = 3, modeType = 'FW', maxItter = 5)
#makeMinPickle(gon = 7, vertex = 3, modeType = 'FW', maxItter = 6)
#makeMinPickle(gon = 7, vertex = 3, modeType = 'HW', maxItter = 6)




#makeMinPickle = True
##makeMinPickle = False
#
#if makeMinPickle:
#    test = PlanarLayout(gon = 8, vertex = 3, side =1, radius_method = 'lin', modeType = 'FW')
#    test.populate(maxItter = 4)
#    
#    pickleDict = {}
#    pickleDict['gon'] = test.gon
#    pickleDict['vertex'] = test.vertex
#    pickleDict['itter'] = test.itter
#    pickleDict['modeType'] = test.modeType
#    pickleDict['Es'] = test.Es
#    pickleDict['Eorder'] = test.Eorder
#
#    
#    pickleName = str(test.gon) + 'gon_' + str(test.vertex) + 'vertex_' + str(test.itter) + '_' + str(test.modeType) + '_DOSminimal'
#    pickleDict['name'] = pickleName
#    
#    pickle.dump(pickleDict, open(pickleName + '.pkl', 'wb'))
    




minimal = True
#minimal = False
if showHardWall:
    if minimal:
#        testDict = pickle.load( open('8gon_3vertex_4_FW_DOSminimal.pkl', "rb"))
#        testDict = pickle.load( open('8gon_3vertex_5_FW_DOSminimal.pkl', "rb"))
        
#        testDict = pickle.load( open('7gon_3vertex_6_FW_DOSminimal.pkl', "rb"))
        testDict = pickle.load( open('7gon_3vertex_6_HW_DOSminimal.pkl', "rb"))
    else:
    #    test = PlanarLayout(file_path = '7gon_3vertex_ 3.pkl')
    #    test = PlanarLayout(file_path = '7gon_3vertex_ 4.pkl')
    #    test = PlanarLayout(file_path = '7gon_3vertex_ 5.pkl')
        #test = PlanarLayout(file_path = '7gon_3vertex_ 6.pkl')
        
        
        test = PlanarLayout(file_path = '8gon_3vertex_ 3.pkl')
    #    test = PlanarLayout(file_path = '8gon_3vertex_ 4.pkl')
    
    
    ######
    #trying DOS comparison versus sytem size instead
    ######
    
    #set up frequency sweep
    freq_range = 4.01
    freq_res = 0.04
#    freq_res = 0.06
    #freq_res = 0.12
    
    freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
    freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
    
    
    
    pylab.figure(7)
    pylab.clf()
    
    ax1 = pylab.subplot(1,1,1)

    if minimal:
        Evals = testDict['Es']
        Eorder = testDict['Eorder']
        itt = testDict['itter']
        
        [fullDOS, bins_out] = numpy.histogram(Evals, freq_bins)
    else:
        [Evals, Eorder] =  [test.Es, test.Eorder]
        itt = test.itter
        
        [fullDOS, bins_out] = numpy.histogram(test.Es, freq_bins)

    
    [DOS, bins_out] = numpy.histogram(Evals, freq_bins)

    bins_centers = (bins_out[0:-1] + bins_out[1:])/2
    binWidth = bins_out[1] - bins_out[0]

    
    pylab.bar(bins_centers, 1.*DOS/len(Evals), width = binWidth, color =  'deepskyblue', label = str(itt), alpha = 1)
#    pylab.bar(bins_centers, 1.*DOS/len(Evals), width = binWidth, color =  'deepskyblue', label = '6 shells', alpha = 1)
    
        
    pylab.xlabel('Energy (|t|)')
    pylab.ylabel('DOS')
    ax1.legend(loc = 'upper right')
    pylab.title('Finite Hepatgon-Kagome DOS')
    ax1.set_xlim([-2.5,4])
    ax1.set_ylim([0,40])
    ax1.set_ylim([0,0.05])
    
    pylab.show()






if showRegularized:
    ###sample figure of FW v HW regu;arized
#    loaded = pickle.load(open('7gon_3vertex_4_FW_reg.pkl', "rb" ) )
#    loaded2 = pickle.load(open('7gon_3vertex_4_HW_reg.pkl', "rb" ) )
    loaded = pickle.load(open('7gon_3vertex_5_FW_reg.pkl', "rb" ) )
    loaded2 = pickle.load(open('7gon_3vertex_5_HW_reg.pkl', "rb" ) )
    pylab.figure(44)
    pylab.clf()
    ax = pylab.subplot(1,1,1)
    
    xs = numpy.linspace(0,1,len(loaded['regEs']))
    xs2 = numpy.linspace(0,1,len(loaded2['regEs']))
    pylab.plot(xs, loaded['regEs'], color =  'mediumblue', marker = '.', linestyle = '', label = loaded['name'])
    pylab.plot(xs2, loaded2['regEs'], color =  'deepskyblue', marker = '.', linestyle = '', label = loaded2['name'])
    
    ax.set_ylim([-2.5, 4.2])
    pylab.ylabel('Energy (|t|)')
    pylab.xlabel('normalized eigenvector index')
    pylab.title('regularized layout')
    
    ax.legend(loc = 'upper left')
    pylab.show()
        
        
    ####   sample DOS plot
    #set up frequency sweep
    freq_range = 4.01
    freq_res = 0.04
    
    freqs = scipy.arange(-freq_range/2, freq_range, freq_res) + freq_res/2.
    freq_bins = scipy.arange(-freq_range/2, freq_range+freq_res, freq_res)
    
    #[fullDOS, bins_out] = numpy.histogram(loaded['regEs'], freq_bins)
    #[fullDOS2, bins_out] = numpy.histogram(loaded2['regEs'], freq_bins)
    
    pylab.figure(45)
    pylab.clf()
    ax1 = pylab.subplot(1,1,1)
    
    
    [DOS_unnorm, bins_out] = numpy.histogram(loaded['regEs'], freq_bins)
    DOS= DOS_unnorm*1.0/len(loaded['regEs'])
    
    [DOS2_unnorm, bins_out] = numpy.histogram(loaded2['regEs'], freq_bins)
    DOS2= DOS2_unnorm*1.0/len(loaded2['regEs'])
    
    bins_centers = (bins_out[0:-1] + bins_out[1:])/2
    binWidth = bins_out[1] - bins_out[0]
    
    pylab.bar(bins_centers, 1.*DOS, width = binWidth, color = 'mediumblue', label = loaded['name'], alpha = 0.6)
    pylab.bar(bins_centers, 1.*DOS2, width = binWidth, color = 'deepskyblue', label = loaded2['name'], alpha = 0.6)   
#    pylab.bar(bins_centers, 1.*DOS, width = binWidth, color = 'mediumblue', label = 'full-wave', alpha = 0.6)
#    pylab.bar(bins_centers, 1.*DOS2, width = binWidth, color = 'deepskyblue', label = 'half-wave', alpha = 0.6)
        
    pylab.xlabel('Energy (|t|)')
    pylab.ylabel('DOS')
    pylab.title('Regularized Heptagon-Kagome DOS')
    ax1.set_xlim([-2.5,4.1])
    #ax1.set_xlim([-2.5,4.0])
    ax1.set_ylim([0,0.035])
    ax1.legend(loc= 'upper right')
    
    pylab.show()
        
