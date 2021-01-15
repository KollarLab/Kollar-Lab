# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 12:02:58 2020

@author: Kollarlab
"""

import time
import os
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
import VNAplottingTools as VNAplots

project_dir = r'Z:\Data\HouckDualHangerFluxonium'
meas_type = 'Spec'
save_Dir = userfuncs.saveDir(project_dir, meas_type)


#def SRS_voltage_ramp(newV, step_size = 0.001, step_time = 0.02):
#    deltaV = newV - SRS.Volt
#    numSteps = int(np.abs(np.ceil(deltaV/step_size)))
#    vsteps = np.linspace(SRS.Volt, newV, numSteps)
#    for vstep in vsteps:
#        SRS.Volt = np.round(vstep,6)
#        time.sleep(step_time)
#    return

CAV_Attenuation = 30

#start_power = -40
#stop_power = 0
#power_points = 21
#powers = np.linspace(start_power, stop_power, power_points)

settings = vna.spec_default_settings()



settings['scanname'] = 'qubit2_overnight'

settings['start_voltage'] = -2
settings['stop_voltage'] =  2
settings['voltage_points'] = 251

settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1

settings['channel'] = 1
settings['avg_time'] = 200
settings['measurement'] = 'S21'
settings['start'] = 6.5e9
settings['stop'] = 9.5e9
settings['sweep_points'] = 3001
settings['CAVpower'] = -55 + CAV_Attenuation
settings['RFpower'] = -15
settings['ifBW'] = .2e3

autoscan_settings = vna.trans_default_settings()
cavity_center =  7.577e9
cavity_span = 50e6
autoscan_settings['sweep_points'] = 501
autoscan_settings['ifBW'] = settings['ifBW']
autoscan_settings['avg_time'] = 15
autoscan_settings['start'] = cavity_center-cavity_span/2
autoscan_settings['stop'] = cavity_center+cavity_span/2
autoscan_settings['RFpower'] = settings['CAVpower']





#settings['scanname'] = 'qubit2_speccheck'
#
#settings['start_voltage'] = 0.525
#settings['stop_voltage'] =  0.575
#settings['voltage_points'] = 10
#
#settings['RFport'] = 3
#settings['Mport'] = 2
#settings['CAVport'] = 1
#
#settings['channel'] = 1
#settings['avg_time'] = 20
#settings['measurement'] = 'S21'
#settings['start'] = 6.4e9
#settings['stop'] = 7.1e9
#settings['sweep_points'] = 501
#settings['CAVpower'] = -55 + CAV_Attenuation
#settings['RFpower'] = -10
#settings['ifBW'] = .2e3
#
#autoscan_settings = vna.trans_default_settings()
#cavity_center =  7.577e9
#cavity_span = 50e6
#autoscan_settings['sweep_points'] = 501
#autoscan_settings['ifBW'] = 0.5e3
#autoscan_settings['avg_time'] = 7
#autoscan_settings['start'] = cavity_center-cavity_span/2
#autoscan_settings['stop'] = cavity_center+cavity_span/2
#autoscan_settings['RFpower'] = -55 + CAV_Attenuation


#settings['scanname'] = 'qubit2_fluxonCheck'
#
#settings['start_voltage'] = 0.464
#settings['stop_voltage'] =  0.466847
#settings['voltage_points'] = 2
#
#settings['RFport'] = 3
#settings['Mport'] = 2
#settings['CAVport'] = 1
#
#settings['channel'] = 1
#settings['avg_time'] = 40
#settings['measurement'] = 'S21'
#settings['start'] = 6.062e9 - 15e6
#settings['stop'] = 6.062e9 + 15e6
#settings['sweep_points'] = 251
#settings['CAVpower'] = -55 + CAV_Attenuation
#settings['RFpower'] = -10
#settings['ifBW'] = 0.5e3
#
#autoscan_settings = vna.trans_default_settings()
#cavity_center =  7.576e9
#cavity_span = 50e6
#autoscan_settings['sweep_points'] = 501
#autoscan_settings['ifBW'] = 0.5e3
#autoscan_settings['avg_time'] = 7
#autoscan_settings['start'] = cavity_center-cavity_span/2
#autoscan_settings['stop'] = cavity_center+cavity_span/2












#set voltage sweep
voltages = np.round(np.linspace(settings['start_voltage'], settings['stop_voltage'], settings['voltage_points']),6)
max_voltage = 3
if np.max(voltages) > 3:
    raise ValueError('max voltage too! large')
else:
    settings['voltages'] = voltages

#set file name
stamp = userfuncs.timestamp()
scanname = 'specFluxScan_' + settings['scanname'] + '_' + stamp




trans_mags = np.zeros((settings['voltage_points'], autoscan_settings['sweep_points']))
trans_phases = np.zeros((settings['voltage_points'], autoscan_settings['sweep_points']))

mags = np.zeros((settings['voltage_points'], settings['sweep_points']))
phases = np.zeros((settings['voltage_points'], settings['sweep_points']))

SRS.Volt = 0
SRS.Output = 'On'

t0 = time.time()

identifier = 'Cav Power : ' + str(settings['CAVpower'] - CAV_Attenuation) + ' dB'
#identifier = 'RF Power : ' + str(settings['RFpower']) + ' dB'

for vind in range(len(voltages)):
    voltage = voltages[vind]
    print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
    
    SRS.voltage_ramp(voltage)
    time.sleep(0.1)
    
    vna.reset()
    vna.output = 'on'
    
    print('trans')
    trans_data = vna.trans_meas(autoscan_settings)
    trans_freqs = trans_data['xaxis']
    trans_mags[vind] = trans_data['mag']
    trans_phases[vind] = trans_data['phase']
    
#    settings['CAVfreq'] = trans_freqs[np.argmin(trans_data['mag'])]
    settings['CAVfreq'] = trans_freqs[np.argmax(trans_data['mag'])]
    print('spec')
    data = vna.spec_meas(settings)
    
    vna.autoscale()
    
    ##took out this normalization
    #it was here from transmission power scans (and used to be devision.)
    #I've switched it to subtraction, but basically I do the real offset subtraction that I want lower down.
    
#    normalization = np.mean(data['mag'][-10:]) #I think that this is only the -10th data point. And for log data we might want to add or subtract it, not devide by
#    mags[vind] = data['mag'] - normalization
#    
    
    mags[vind] = data['mag']
    phases[vind] = data['phase']
    freqs = data['xaxis']  
    
    if vind==0:
        t1=time.time()
        tdiff = t1-t0
        ttotal = tdiff*len(voltages)
        
        estimatedTime = tdiff*len(voltages)
        
        print('    ')
        print('Single run time: {}, estimated total time: {}'.format(tdiff, ttotal))
        print('    ')
        print('estimated time for this scan : ' + str(np.round(estimatedTime/60, 1)) + ' minutes')
        print('estimated time for this scan : ' + str(np.round(estimatedTime/60/60, 2)) + ' hours')
        print('    ')
        
#    VNAplots.general_VNAplot(freqs, mags[0:vind+1,:], phases[0:vind+1,:], voltages[0:vind+1], scanname, 
#                         xlabel = 'Frequency (GHz)', ylabel = 'Voltage (V)', identifier = identifier,
#                         fig_num = 1)
    transdata = {}
    transdata['xaxis'] = trans_freqs
    transdata['mags'] = trans_mags[0:vind+1,:]
    transdata['phases'] = trans_phases[0:vind+1,:]
    
    specdata = {}
    specdata['xaxis'] = freqs
    specdata['mags'] = mags[0:vind+1,:]
    specdata['phases'] = phases[0:vind+1,:]
    
    singledata = {}
    singledata['xaxis'] = freqs
    singledata['mag'] = data['mag'] 
    singledata['phase'] = data['phase']
    
    trans_labels = ['Freq (GHz)','Voltage (V)']
    spec_labels = ['Freq (GHz)','Voltage (V)']
    
    
    #modify the spec data to subtract the offset in amp and phase
    #and then plot the modified version
    specplotdata = {}
    specplotdata['xaxis'] = specdata['xaxis']
    specplotdata['mags'] = specdata['mags']
    specplotdata['phases'] = specdata['phases']
    
    mat = np.copy(specplotdata['mags'])
    for ind in range(0, mat.shape[0]):
        mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
    specplotdata['mags'] = mat
    
    mat = np.copy(specplotdata['phases'])
    for ind in range(0, mat.shape[0]):
        mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
    specplotdata['phases'] = mat
    
    VNAplots.spec_fluxscanplot(transdata, specplotdata, singledata, voltages[0:vind+1], scanname, trans_labels, spec_labels, identifier, fig_num = 1)
        
    if np.mod(vind,10) ==4:
        userfuncs.SaveFull(save_Dir, scanname, ['mags', 'phases', 'freqs', 'powers', 
                                        'CAV_Attenuation', 
                                        'trans_freqs', 'trans_mags', 'trans_phases', 
                                        'autoscan_settings'], locals(), expsettings=settings)
        

t2 = time.time()
print('Elapsed time: {}'.format(t2-t0))

#return to zero voltage
SRS_voltage_ramp(0)
SRS.Output = 'Off'
vna.output = 'off'

#identifier = 'Cav Power : ' + str(settings['RFpower']) + ' dB'
#VNAplots.general_VNAplot(freqs, mags, phases, voltages, scanname, 
#                         xlabel = 'Frequency (GHz)', ylabel = 'Voltage (V)', identifier = identifier,
#                         fig_num = 1)

transdata = {}
transdata['xaxis'] = trans_freqs
transdata['mags'] = trans_mags
transdata['phases'] = trans_phases

specdata = {}
specdata['xaxis'] = freqs
specdata['mags'] = mags
specdata['phases'] = phases

singledata = {}
singledata['xaxis'] = freqs
singledata['mag'] = data['mag']
singledata['phase'] = data['phase']

trans_labels = ['Freq (GHz)','Voltage (V)']
spec_labels = ['Freq (GHz)','Voltage (V)']

#identifier = 'Cav Power : ' + str(settings['RFpower']) + ' dB'

VNAplots.spec_fluxscanplot(transdata, specplotdata, singledata, voltages, scanname, trans_labels, spec_labels, identifier, fig_num = 1)


userfuncs.SaveFull(save_Dir, scanname, ['mags', 'phases', 'freqs', 'powers', 
                                        'CAV_Attenuation', 
                                        'trans_freqs', 'trans_mags', 'trans_phases', 
                                        'autoscan_settings'], locals(), expsettings=settings)
plt.savefig(os.path.join(save_Dir, scanname+'.png'), dpi = 150)


