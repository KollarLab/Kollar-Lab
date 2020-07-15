# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:36:46 2020

@author: Kollarlab
"""

import Calibration.sgsIQcal as SGSiq
import Calibration.mixerIQcal as mixeriq
import userfuncs as uf
import pylab
import numpy
from datetime import datetime

instruments              = {}
instruments['AWG']       = hdawg
instruments['Digitizer'] = card
instruments['LO']        = logen
instruments['RFgen']     = rfgen
instruments['Trigger']   = triggergen

settings = SGSiq.GetDefaultSettings()

def runIQcal(instruments, settings, angle, imp):
    settings['NumPoints'] = 100
    settings['AngleImp']  = angle
    settings['AmpImp']    = imp
    
    SGSa, SGSc, SGSp, SGSe = SGSiq.calibrate_SGS_IQ(instruments, settings)
    
    SGSxx, SGSyy = uf.make_elipse(SGSa, SGSc, SGSp,150)
    
    maxSGSx = abs(max(SGSxx))
    maxSGSy = abs(max(SGSyy))
    maxdim  = max([maxSGSx,maxSGSy])
    
    figtitle = 'IQ ellipse relative to perfect circle.\n Alpha: {}, Phi: {}'.format(settings['AmpImp'], settings['AngleImp'])
    myfig = pylab.figure(1)
    pylab.clf()
    pylab.plot(SGSxx, SGSyy, color='r')
    
    circle2 = pylab.Circle((SGSc[0], SGSc[1]), maxdim, color='b', fill=False)
    ax = pylab.gca()
    ax.add_artist(circle2)
    ax.set_xlim([-1.1*maxdim, 1.1*maxdim])
    ax.set_ylim([-1.1*maxdim, 1.1*maxdim])
    ax.set_title(figtitle)
    ax.set_aspect('equal')
    
    date  = datetime.now()
    stamp = date.strftime('%Y%m%d_%H%M%S')
    
    figname = 'IQcal_{}'.format(stamp)
    figpath = r'C:\Users\Kollarlab\Desktop\IQcal'
    uf.savefig(myfig, figname, figpath, png=True)

def IQcal(instruments, settings, angle, imp):
    settings['NumPoints'] = 100
    settings['AngleImp']  = angle
    settings['AmpImp']    = imp
    settings['verbose']   = False
    
    SGSa, SGSc, SGSp, SGSe = SGSiq.calibrate_SGS_IQ(instruments, settings)
    return SGSe

#angles  = numpy.linspace(-0.10,.025,num = 20)
#impairs = numpy.linspace(0.87,.93,num = 20)
#
#ecc = numpy.zeros((len(angles), len(impairs)))
#for i in range(len(angles)):
#    print('{} out of {} angle measurements, {} amp measures per setting'.format(i, len(angles), len(impairs)))
#    for j in range(len(impairs)):
#        if j%5 == 0:
#            print('{} measurements completed'.format(j))
#        ecc[i,j] = IQcal(instruments, settings, angles[i], impairs[j])
#
#saveMe = pylab.figure(1)
#pylab.imshow(ecc, origin = 'lower', extent = [impairs[0], impairs[-1], angles[0], angles[-1]])
#pylab.colorbar()
#ax = pylab.gca()
#ax.set_ylabel('Angles')
#ax.set_xlabel('Amplitudes')
#ax.set_title('Eccentricity vs angle impairement and amplitude impairement')
#
#date  = datetime.now()
#stamp = date.strftime('%Y%m%d_%H%M%S')
#
#filename = 'IQcal_{}points_{}'.format(len(angles)*len(impairs),stamp)
#path = r'C:\Users\Kollarlab\Desktop\IQcal'
#uf.SaveFull(path, filename, 'ecc', locals(), settings, instruments, saveMe)
#uf.savefig(saveMe, filename, path, png=True)
    
settings['NumPoints'] = 20
settings['MinAngle'] = -0.2
settings['MaxAngle'] = 0.2
settings['AngleNum'] = 1
settings['Amp_points'] = 10

ecc = SGSiq.calibrate_SGS_IQ_fast(instruments, settings)
impairs = [0.8,1.2]
angles  = [settings['MinAngle'], settings['MaxAngle']]
pylab.figure()
pylab.imshow(ecc, origin = 'lower', extent = [impairs[0], impairs[-1], angles[0], angles[-1]])
pylab.colorbar()
ax = pylab.gca()
ax.set_ylabel('Angles')
ax.set_xlabel('Amplitudes')
ax.set_title('Eccentricity vs angle impairement and amplitude impairement')
