# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:36:46 2020

@author: Kollarlab
"""

import Calibration.sgsIQcal as SGSiq
import Calibration.mixerIQcal as mixeriq
import userfuncs as uf
import pylab

instruments = {}
instruments['AWG'] = hdawg
instruments['Digitizer'] = card
instruments['LO'] = logen
instruments['RFgen'] = rfgen

settings = mixeriq.GetDefaultSettings()
settings['NumPoints'] = 50

SGSa, SGSc, SGSp = SGSiq.calibrate_SGS_IQ(instruments, settings)
mixera, mixerc, mixerp = mixeriq.calibrate_mixer_IQ(instruments, settings)

SGSxx, SGSyy = uf.make_elipse(SGSa, SGSc, SGSp,150)
mixerxx, mixeryy = uf.make_elipse(mixera, mixerc, mixerp, 150)

maxSGSx = abs(max(SGSxx))
maxSGSy = abs(max(SGSyy))
maxMixerx = abs(max(mixerxx))
maxMixery = abs(max(mixeryy))

normSGSx = SGSxx/maxSGSx
normSGSy = SGSyy/maxSGSy
normMixerx = mixerxx/maxMixerx
normMixery = mixeryy/maxMixery

pylab.figure(1)
pylab.plot(SGSxx, SGSyy)
pylab.plot(mixerxx, mixeryy)
pylab.title('Raw comparison between IQ mod on SGS and phase control on SGS')
pylab.legend(['SGS','Mixer'])
pylab.show()

pylab.figure(2)
pylab.plot(normSGSx, normSGSy)
pylab.plot(normMixerx, normMixery)
pylab.title('Normalized traces')
pylab.legend(['SGS','Mixer'])
pylab.show()