# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 14:36:59 2020

@author: Kollarlab

example code for reading traces from the VNA.

-assumes the hardware is basically set up for S21

-assumes hardware is preinitialized.


"""

#init_dir = dir()

import pylab
import userfuncs
import os

#finished_imports = dir()
#
#import_list = [imp for imp in finished_imports if imp not in set(init_dir+['init_dir'])]
#
#vars_to_save = vars_to_save + import_list
#userfuncs.reset_local_vars(locals(), globals(), vars_to_save+import_list)

##########
#save location
###########

saveDir = r'Z:\Data\Fluxonium_Raman\CRF01_A3\Trans\20210305'
if not os.path.isdir(saveDir):
    #directory does not exist. We must make it.
    os.mkdir(saveDir)

stamp = userfuncs.timestamp()

name = 'cavity_search'

filename = name + '_' + stamp


settings = vna.trans_default_settings()


settings['start_freq'] = 7e9
settings['stop_freq'] = 8e9
settings['sweep_points'] = 501
settings['RFpower'] = -
settings['avg_time'] = 10
settings['measurement'] = 'S21'

#mag, phase, freqs = vna.trans_meas(settings)
data_dict = vna.trans_meas(settings)
mag = data_dict['mag']
phase = data_dict['phase']
freqs = data_dict['xaxis']

vna.output = 'OFF'

fig = pylab.figure(1)
pylab.clf()

#ax = pylab.subplot(1,2,1)
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
figsToSave = [fig]

userfuncs.SaveFull( saveDir, filename, varsToSave, locals(), expsettings = settings)
userfuncs.savefig(fig, filename, path = saveDir,  png = True)

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
