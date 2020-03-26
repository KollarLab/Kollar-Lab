# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 18:45:42 2020

@author: Kollarlab
"""

import numpy
import scipy
import pylab


dataSize = 512000

array = scipy.arange(0, dataSize,1)
array2 = scipy.arange(0, dataSize*10000,10000)
array3 = scipy.arange(0, dataSize,1.)*11000

pylab.figure(11)
pylab.clf()

pylab.plot(array*1000)
pylab.plot(array*10000)
pylab.plot(array2)
pylab.plot(array3)

pylab.show()