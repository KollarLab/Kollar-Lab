# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 16:20:18 2020

@author: Kollarlab
"""
import time
import os
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
#import VNAplottingTools as VNAplots
import  utility.plotting_tools as VNAplots


def general_colormap_subplot(ax,xaxis, yaxis, data, cmap = 'jet_r', vmax = np.NaN, vmin = np.NaN):
    '''given a subplot object, it should make a decent imshow plot of the data.
    
    You have to handle subtracting out any attenuation and do the labeling after.
    
    Trying to make this the most generic function possible.'''
    
    plt.sca(ax)
    
    #this version will handle the color plot if there is only one row.
    if (type(xaxis) == list) or (type(xaxis) == np.ndarray):
        if (len(xaxis) == 1): 
            xlimits = [xaxis[0]-1e-3, xaxis[0]+1e-3]
        else:
            xlimits = [xaxis[0], xaxis[-1]]
    else:
        xlimits =  [xaxis-1e-3, xaxis+1e-3]
        
    if (type(yaxis) == list) or (type(yaxis) == np.ndarray):
        if (len(yaxis)) == 1:
            ylimits = [yaxis[0]-1e-3, yaxis[0]+1e-3]
        else:
            ylimits = [yaxis[0], yaxis[-1]]  
    else:
        ylimits = [yaxis-1e-3, yaxis+1e-3]
        
    limits = np.concatenate((xlimits, ylimits))
    
    if (not np.isnan(vmax)) and (not np.isnan(vmin)):
#        print('version with color scale limits')
        plt.imshow(data, extent = limits, origin='lower', aspect='auto', cmap=cmap, vmin = vmin, vmax = vmax)
    else:
        plt.imshow(data, extent = limits, origin='lower', aspect='auto', cmap=cmap)
    plt.colorbar()
        
    return


#userfuncs.SaveFull(saveDir, filename, ['mags', 'phases', 'freqs', 
#                                            'CAV_Attenuation', 
#                                            'trans_freqs', 'trans_mags', 'trans_phases'], 
#                                            locals(), expsettings=fullsettings)

#pull up a data file, and reprocess the data into the same structure that the acquisition code
#dat = userfuncs.LoadFull(r'Z:\Data\Fluxonium_Raman\CRF01_A3\spec_flux_scan\20210117\diagnosticOvernight_belowCav_20210117_091112.pkl')
#dat = userfuncs.LoadFull(r'Z:\Data\Fluxonium_Raman\CRF01_A3\spec_flux_scan\20210116\diagnosticOvernight_aboveCav_20210116_205523.pkl')
#dat = userfuncs.LoadFull(r'Z:\Data\Fluxonium_Raman\CRF01_A3\spec_flux_scan\20210118\plasmonCrosingFollowUp_20210118_115030.pkl')
#dat = userfuncs.LoadFull(r'Z:\Data\Fluxonium_Raman\CRF01_A3\spec_flux_scan\20210118\plasmonCrosingFollowUp2_20210118_132316.pkl')
#dat = userfuncs.LoadFull(r'Z:\Data\HouckDualHangerFluxonium\spec_flux_scan\20210119\20kHzfilter_finescan_fluxon_tracking_20210119_234214.pkl')
saveDir = r'Z:\Data\Fluxonium_Raman\WTF01_B3\spec_flux_scan\20210607'
filename = 'cavity1_spec_scan_0_spec_20210607_170153.pkl'
dat = userfuncs.LoadFull(os.path.join(saveDir, filename))

data = dat[0]
expsettings = dat[1]
#voltages = expsettings['voltages']
#yaxis = data['autler_freqs']
yaxis = data['voltages']
#yaxis = data['powers']
scanname = expsettings['scanname']

#testData = {}
#testData['xaxis'] = specdata['xaxis']
#testData['mags'] = specdata['mags']
#testData['phases'] = specdata['phases']
#testData = data['full_data']
testData = data['specdata']
#testData['xaxis'] = data['freqs']
#testData['mags'] = data['mags']
#testData['phases'] = data['phases']

specdata = data['specdata']
##specdata['xaxis'] = data['freqs']
##specdata['mags'] = data['mags']
##specdata['phases'] = data['phases']
#
#
singledata = data['singledata']
##singledata['xaxis'] = data['freqs']
##singledata['mag'] = data['mags'][-1,:]
##singledata['phase'] = data['phases'][-1,:]
#
#
transdata = data['transdata']
##transdata['xaxis'] = data['trans_freqs']
##transdata['mags'] = data['trans_mags']
##transdata['phases'] = data['trans_phases']

trans_labels = ['Freq (GHz)','Voltage (V)']
spec_labels = ['Freq (GHz)','Voltage (V)']
#labels = ['Freq (GHz)', 'Freq (GHz)']
labels = ['Freq (GHz)','Voltage (V)']


mat = np.copy(testData['mags'])
for ind in range(0, mat.shape[0]):
    mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
    
testData['mags'] = mat


mat = np.copy(testData['phases'])
for ind in range(0, mat.shape[0]):
    mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
    
testData['phases'] = mat


#VNAplots.spec_fluxscanplot(transdata, testData, singledata, voltages, scanname, trans_labels, spec_labels, identifier, fig_num = 4)
#spec_fluxscanplot(transdata, testData, singledata, voltages, scanname, trans_labels, spec_labels, identifier, fig_num = 4)



#defaultcmap = 'jet'
#defaultcmap = 'jet_r'
#defaultcmap = 'bwr'
#defaultcmap = 'RdYlBu'
#defaultcmap = 'hot_r'
defaultcmap = 'hot'
#defaultcmap = 'afmhot_r'
#defaultcmap = 'CMRmap'
#defaultcmap = 'bone'
#defaultcmap = 'cool'
#defaultcmap = 'cividis'
#defaultcmap = 'viridis'
#defaultcmap = 'Blues'


#yaxis = voltages
#
fig = plt.figure(4, figsize=(13,8))
plt.clf()

ax = plt.subplot(3,2,1)
general_colormap_subplot(ax,transdata['xaxis'], yaxis, transdata['mags'], cmap = defaultcmap)
plt.xlabel(trans_labels[0])
plt.ylabel(trans_labels[1])
plt.title('Trans mag')

ax = plt.subplot(3,2,2)
general_colormap_subplot(ax,transdata['xaxis'], yaxis, transdata['phases'], cmap = defaultcmap)
plt.xlabel(trans_labels[0])
plt.ylabel(trans_labels[1])
plt.title('Trans phase')



ax = plt.subplot(3,2,3)
ampMin = -2
ampMax = 1

ampMin = -5
ampMax = 2

#ampMin = -1
#ampMax = 0
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, testData['mags'], cmap = 'viridis')
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, testData['mags'], vmin = -0.5, vmax = 0, cmap = 'viridis')
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, testData['mags'], vmin = -2, vmax = 0., cmap = 'jet')
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, testData['mags'], vmin = -0.5, vmax = 0., cmap = 'hot')
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, testData['mags'], vmin = -2, vmax = 0., cmap = 'hot_r')
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, testData['mags'], vmin = ampMin, vmax = ampMax, cmap = 'jet_r')


general_colormap_subplot(ax,testData['xaxis'], yaxis, testData['mags'], vmin = ampMin, vmax = ampMax, cmap = defaultcmap)
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, testData['mags'],  cmap = defaultcmap)
plt.xlabel(labels[0])
plt.ylabel(labels[1])
plt.title('mag')

ax = plt.subplot(3,2,4)
phaseMin = -30
phaseMax = 30
general_colormap_subplot(ax,testData['xaxis'], yaxis, testData['phases'], vmin = phaseMin, vmax = phaseMax, cmap = defaultcmap)
plt.xlabel(labels[0])
plt.ylabel(labels[1])
plt.title(' phase')

ax = plt.subplot(3,2,5)
plt.plot(singledata['xaxis'], singledata['mag'])
plt.xlabel(spec_labels[0])
plt.title('Single shot spec mag')

ax = plt.subplot(3,2,6)
plt.plot(singledata['xaxis'], singledata['phase'])
plt.xlabel(spec_labels[0])
plt.title('Single shot spec mag')


#plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
plt.suptitle('Filename: {}, {}'.format(scanname, 'replot'))

fig.canvas.draw()
fig.canvas.flush_events()


from scipy.ndimage import gaussian_filter

plt.figure(7)
plt.clf()
ax = plt.subplot(1,1,1)

mat3  = np.copy(testData['mags'])
mat4 = gaussian_filter(mat3, sigma = 0.5)
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, mat4, vmin = ampMin, vmax = ampMax, cmap = 'jet_r')
general_colormap_subplot(ax,specdata['xaxis'], yaxis, mat4, vmin = ampMin, vmax = ampMax, cmap = defaultcmap)
plt.xlabel(spec_labels[0])
plt.ylabel(spec_labels[1])
plt.title('Spec mag')
#plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
plt.suptitle('Filename: {}, {}'.format(scanname, 'replot'))
plt.show()



plt.figure(8)
plt.clf()
ax = plt.subplot(1,1,1)

mat3  = np.copy(testData['phases'])
mat4 = gaussian_filter(mat3, sigma =0.5)
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, mat4, vmin = phaseMin, vmax =phaseMax, cmap = 'jet_r')
general_colormap_subplot(ax,specdata['xaxis'], yaxis, mat4, vmin = phaseMin, vmax =phaseMax, cmap = defaultcmap)
plt.xlabel(spec_labels[0])
plt.ylabel(spec_labels[1])
plt.title('Spec phase')
#plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
plt.suptitle('Filename: {}, {}'.format(scanname, 'replot'))
plt.show()




#mat3  = np.copy(testData['mags'])
#general_colormap_subplot(ax,specdata['xaxis'], yaxis, mat3, vmin = -1.5, vmax = -0.1, cmap = 'jet_r')
#plt.xlabel(spec_labels[0])
#plt.ylabel(spec_labels[1])
#plt.title('Spec mag')
#plt.suptitle('Filename: {}, {}'.format(scanname, identifier))
#plt.show()






