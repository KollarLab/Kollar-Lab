# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 14:36:59 2020

@author: Kollarlab

example code for reading traces from the VNA.

-assumes the hardware is basically set up for S21

-assumes hardware is preinitialized.


"""

init_dir = dir()

import pylab
import userfuncs

finished_imports = dir()

import_list = [imp for imp in finished_imports if imp not in set(init_dir+['init_dir'])]

vars_to_save = vars_to_save + import_list
userfuncs.reset_local_vars(locals(), globals(), vars_to_save+import_list)

##########
#save location
###########

saveDir = r'Z:\Data\BF1\InitialCalibrations_10-20'

stamp = userfuncs.timestamp()

name = 'Output3'

filename = name + '_' + stamp


settings = vna.TransDefaultSettings()


settings['start'] = 1e9
settings['stop'] = 15e9
settings['sweep_points'] = 2001
settings['RFpower'] = -20
settings['averages'] = 10
settings['measurement'] = 'S21'

mag, phase, freqs = vna.TransMeas(settings)

vna.output = 'OFF'

fig = pylab.figure(1)
pylab.clf()

ax = pylab.subplot(1,2,1)
pylab.plot(freqs/1e9, mag)
pylab.xlabel('frequency (GHz)')
pylab.ylabel('logmag')


ax = pylab.subplot(1,2,2)
pylab.plot(freqs/1e9, phase)
pylab.xlabel('frequency(GHz)')
pylab.ylabel('unwrapped phase')

pylab.suptitle(filename)

pylab.show()

#####
#try to save
#####

varsToSave = ['mag', 'phase', 'freqs','filename','fig']
#figsToSave = [fig]

#userfuncs.SaveFull( saveDir, filename, varsToSave, locals(), expsettings = settings)

#######
##read the file
#######
#import pickle
#def tempLoadFull(path):
#    fullData    = pickle.load(open(path,'rb'))
#
##    figures     = fullData['Figures']
#    expsettings = fullData['ExpSettings']
#    hwsettings  = fullData['HWSettings']
#    data        = fullData['Data']
#
##    return [data, expsettings, hwsettings, figures]
#    return [data, expsettings, hwsettings]





##figures = dataFile[3]
#
#
#
#dataFile = userfuncs.LoadFull(r'Z:\Data\BF1\InitialCalibrations_10-20\lineCal_LinesToFridge_20201015_111108.pkl')
#data = dataFile[0]
#
#readDataMag = data['mag']
#readDataPhase = data['phase']
#
#fig2 = pylab.figure(2)
#pylab.clf()
#
#ax = pylab.subplot(1,2,1)
#pylab.plot(readDataMag)
#pylab.title('reread data')
#
#ax = pylab.subplot(1,2,2)
#pylab.plot(readDataPhase)
#pylab.title('reread data')
#
#
#pylab.tight_layout()
#pylab.show()





#    #Save png of figure for easy viewing
#    uf.savefig(fig1,filename,savepath, png = True)
#    
#    print("std(angles) = " , numpy.std(Angles))
#
#    rfgen.power_Off()
#    logen.power_Off()    
#
#    dataTosave = ['Is','Qs','Amps','Angles', 'actualTimes','xx','yy','Idata', 'Qdata', 'Iav', 'Qav']
#    figsTosave = [fig1]
#
#    if settings['save']:    
#        uf.SaveFull(savepath, filename, dataTosave, locals(), settings, instruments, figsTosave)
