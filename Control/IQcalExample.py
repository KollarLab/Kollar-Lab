# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 11:28:59 2020

@author: Kollarlab
"""
import numpy as np
import Calibration.mixerIQcal as mixer
import Calibration.sgsIQcal as SGSiq
import pylab
import userfuncs as uf

identity = np.array([[1,0],[0,1]])

def correction_factor(instruments, settings, mixermat, minecc = 0.15):
    initmat = identity
    ecc=1
    print('Current ecc: {}'.format(ecc))
    while(ecc>minecc):
        Is, Qs = SGSiq.calibrate_SGS_IQ_basic(instruments, settings, initmat.T)
        
        Is = Is-np.mean(Is)
        Qs = Qs-np.mean(Qs)
        
        Is, Qs = remove_mixer(Is, Qs, mixermat)
        
        a,c,p,ecc = uf.fit_ell_martin(Is, Qs)
        
        absmax = max([max(Is),max(Qs)])
        finalPos = np.array([(idat, qdat) for idat, qdat in zip(Is, Qs)])
        initPos = np.array([(np.cos(th)*absmax, np.sin(th)*absmax) for th in np.linspace(0,2*np.pi,len(Is))])
        
        A, res, rank, s = np.linalg.lstsq(finalPos, initPos, rcond=None)
        initmat = initmat@A
        
        print('Current ecc: {}'.format(ecc))
    
    return initmat

def correction_factor_mixer(instruments, settings):
    Is, Qs, ecc = mixer.calibrate_mixer_IQ(instruments, settings)
    
    Is = Is-np.mean(Is)
    Qs = Qs-np.mean(Qs)
    
    absmax = max([max(Is),max(Qs)])
    finalPos = np.array([(idat, qdat) for idat, qdat in zip(Is, Qs)])
    initPos  = np.array([(np.cos(th)*absmax, np.sin(th)*absmax) for th in np.linspace(0,2*np.pi,len(Is))])
    
    correction, res, rank, s = np.linalg.lstsq(finalPos, initPos, rcond=None)
    
    return correction, Is, Qs

def remove_mixer(Is, Qs, mixermat):
    finalPos = np.array([(idat, qdat) for idat, qdat in zip(Is, Qs)])
    corrected = [(mixermat.T)@p for p in finalPos]
    fixed = np.array(list(zip(*corrected)))
    return fixed[0], fixed[1]

mixerSet = mixer.GetDefaultSettings()
mixerSet ['rfpower'] = -4

mixermat, Id, Qd = correction_factor_mixer(instruments, mixerSet)
Ic, Qc = remove_mixer(Id, Qd, mixermat)

axes1, center1, phi1, ecc1 = uf.fit_ell_martin(Id, Qd)
axes2, center2, phi2, ecc2 = uf.fit_ell_martin(Ic, Qc)

pylab.figure()
pylab.plot(Id, Qd)
pylab.plot(Ic, Qc)
ax = pylab.gca()
ax.set_aspect('equal')
ax.set_title('Removing mixer ellipse from data \n Initial ecc: {}, final ecc: {}'.format(numpy.round(ecc1,3), numpy.round(ecc2,3)))
for i in range(3):
    ax.annotate(i,(Id[i],Qd[i]))
    ax.annotate(i,(Ic[i],Qc[i]))

sgsSet  = SGSiq.GetDefaultSettings()
sgsSet['rfpower'] = -4
sgsSet['NumPoints'] = 20

calibmat = correction_factor(instruments, sgsSet, mixermat)

Is1, Qs1 = SGSiq.calibrate_SGS_IQ_basic(instruments, sgsSet, identity.T)
Is2, Qs2 = SGSiq.calibrate_SGS_IQ_basic(instruments, sgsSet, calibmat.T)
        
Is1 = Is1-np.mean(Is1)
Qs1 = Qs1-np.mean(Qs1)
Is2 = Is2-np.mean(Is2)
Qs2 = Qs2-np.mean(Qs2)

Is1, Qs1 = remove_mixer(Is1, Qs1, mixermat)
Is2, Qs2 = remove_mixer(Is2, Qs2, mixermat)

a,c,p,ecc1 = uf.fit_ell_martin(Is1, Qs1)
a,c,p,ecc2 = uf.fit_ell_martin(Is2, Qs2)

pylab.figure()
pylab.plot(Is1,Qs1, 'x', label='No fix')
pylab.plot(Is2,Qs2, 'o', label='Removing aberrations')
pylab.title('Removing HDAWG/ SGS imperfections \n Init ecc: {}, final ecc: {}'.format(ecc1,ecc2))
ax = pylab.gca()
ax.set_aspect('equal')
ax.legend(loc='upper right')
for i in range(3):
    ax.annotate(i,(Is1[i],Qs1[i]))
    ax.annotate(i,(Is2[i],Qs2[i]))





