# -*- coding: utf-8 -*-
"""
Created on Wed May 13 11:59:51 2020

@author: Kollarlab
"""

import pylab
import pickle
import os
#from mpldatacursor import datacursor
#import mplcursors.cursor as datacursor
from mplcursors import cursor as datacursor


#saveFolder = r'Z:\Data\MeasurementSetupTests\random'
#saveName = 'testfig.png'
#pathStr= os.path.join(saveFolder, saveName) 
#fig.savefig(pathStr )
#    
#saveFolder = r'Z:\Data\MeasurementSetupTests\random'
#saveName = 'testfig.pkl'
#pathStr= os.path.join(saveFolder, saveName) 
##pickle.dump(fig, open('FigureObject.fig.pickle', 'wb'))
#pickle.dump(fig, open(pathStr, 'wb'))
    


saveFolder = r'Z:\Data\MeasurementSetupTests\random'
saveName = 'testfig.pkl'
pathStr= os.path.join(saveFolder, saveName) 
figx = pickle.load(open(pathStr, 'rb'))    


#datacursor()
datatips = datacursor(multiple = True)
figx.show()


#datatips.remove() ##clears all cursors.