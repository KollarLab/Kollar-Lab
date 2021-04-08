# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 10:45:36 2020

@author: Kollarlab
"""

import pylab
import userfuncs

saveDir = r'Z:\Data\BF1\20201117'

stamp = userfuncs.timestamp()

name = 'circulator_at3from1_throughdevice'

filename = name + '_' + stamp

freqs = vna.get_channel_axis(1)
mag = vna.get_trace(1,'Trc1')

plt.figure()
ax = pylab.subplot(1,1,1)
pylab.plot(freqs/1e9, mag)
pylab.xlabel('frequency (GHz)')
pylab.ylabel('logmag')

pylab.suptitle(filename)

pylab.show()

varsToSave = ['mag', 'freqs','filename']

userfuncs.SaveFull( saveDir, filename, varsToSave, locals())