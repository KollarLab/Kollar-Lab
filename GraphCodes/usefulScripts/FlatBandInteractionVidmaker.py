#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 18:27:52 2017

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

#from LayoutGenerator2 import PlanarLayout
#from LayoutGenerator3 import PlanarLayout
from LayoutGenerator4 import PlanarLayout


def interaction_make_vid(source, startInd =0, stopInd = -1, figNum = 8, colorbar = True, plot_links = True, autoscale = False):
    
    test = PlanarLayout(file_path = source)
    startInd = numpy.mod(startInd, len(test.SDx))

    if stopInd <0:
        #set a break point that will not be reached
        breakPoint = len(test.Es)+1
    else:
        breakPoint = stopInd
        
    if startInd > len(test.Es):
        raise ValueError , 'dont have this many eigenvectors'
        
    

    folder = 'Heptagon_' + source[-5] + '_FB_interactions'

    
    currDir = os.getcwd()
    saveDir = os.path.join(currDir,folder)
    if os.path.isdir(saveDir):
        pass
    else:
        os.mkdir(saveDir)
    
    
    print folder + '\n'
    
    #find the flat band
    flat_band = numpy.where(test.Es < -1.95*test.t)[0]
    
    
        
    fig = pylab.figure(figNum)

    #######################################################
    for source_ind in range(startInd, len(test.SDx)):
        print source_ind
        
        map_vect = test.V_int_map(source_ind, flat_band)
        
        fig.clf()
        ax1 = pylab.subplot(1,1,1)

        test.plot_map_state(map_vect, ax1, title = 'flat-band interactions', colorbar = colorbar, plot_links = plot_links, cmap = 'winter', autoscale = autoscale)
        pylab.scatter([test.SDx[source_ind]], [test.SDy[source_ind]], c =  'gold', s = 150, edgecolors = 'k')
        
        fig.set_size_inches(8.0, 7.4)
    
        fig_name = str(test.gon) + '_' + str(test.vertex) + '_' + str(test.itter) + '_site_' + str(source_ind) + '.png' 
        fig_path = os.path.join(saveDir, fig_name)
        pylab.savefig(fig_path)
        
        if source_ind == breakPoint:
            break
        
    pylab.show()
    
    return 








source = '7gon_3vertex_ 2.pkl'
#source = '7gon_3vertex_ 3.pkl'










#interaction_make_vid('7gon_3vertex_ 2.pkl', startInd = -3, stopInd = -1)
#interaction_make_vid('7gon_3vertex_ 2.pkl', startInd = 0, stopInd = 2)
#interaction_make_vid('7gon_3vertex_ 2.pkl', startInd = 0, stopInd = -1)

#make_vid('7gon_3vertex_ 3.pkl', startInd = 397, stopInd = -1, phaseScale = True)


#interaction_make_vid('7gon_3vertex_ 3.pkl', startInd = -3, stopInd = -1)
#interaction_make_vid('7gon_3vertex_ 3.pkl', startInd = 100, stopInd = 103)
interaction_make_vid('7gon_3vertex_ 3.pkl', startInd = 0, stopInd = 2)
#interaction_make_vid('7gon_3vertex_ 3.pkl', startInd = 7, stopInd = 9)
#interaction_make_vid('7gon_3vertex_ 3.pkl', startInd = 35, stopInd = 37)
#interaction_make_vid('7gon_3vertex_ 3.pkl', startInd = 125, stopInd = 127)
#interaction_make_vid('7gon_3vertex_ 3.pkl', startInd = 250, stopInd = 252)







