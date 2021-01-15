# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 14:07:57 2020

@author: Kollarlab
"""
import time
import os
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
import VNAplottingTools as VNAplots

project_dir = r'Z:\Data\HouckDualHangerFluxonium'
meas_type = 'Trans'
save_Dir = userfuncs.saveDir(project_dir, meas_type)


def SRS_voltage_ramp(newV, step_size = 0.005, step_time = 0.001):
    deltaV = newV - SRS.Volt
    numSteps = int(np.abs(np.ceil(deltaV/step_size)))
    vsteps = np.linspace(SRS.Volt, newV, numSteps)
    for vstep in vsteps:
        SRS.Volt = np.round(vstep,6)
        time.sleep(step_time)
    return




CAV_Attenuation = 30
SRS.Range = '10 V'



#
#settings = vna.trans_default_settings()
#
#settings['scanname'] = 'qubit2_overnight'
#
#settings['start_voltage'] = -1.5
#settings['stop_voltage'] = 3.4
#settings['voltage_points'] = 275
#
#settings['powers'] = np.array([-70,-60,-50])  + CAV_Attenuation
#settings['avg_times'] = np.array([60, 30, 15])
#
#settings['start'] = 7.576e9 -40e6
#settings['stop'] = 7.576e9 + 40e6
#
#settings['channel'] = 1
#settings['measurement'] = 'S21'
#settings['sweep_points'] = 701
#settings['ifBW'] = 0.5e3






settings = vna.trans_default_settings()

settings['scanname'] = 'qubit2_avoidedCrossingMap'

settings['start_voltage'] = 0.4
settings['stop_voltage'] = 0.9
settings['voltage_points'] = 15

settings['powers'] = np.array([-55])  + CAV_Attenuation
settings['avg_times'] = np.array([8])

settings['start'] = 7.576e9 -40e6
settings['stop'] = 7.576e9 + 40e6

settings['channel'] = 1
settings['measurement'] = 'S21'
settings['sweep_points'] = 501
settings['ifBW'] = 4e3


















#settings['avg_time'] = 15 #25
#settings['RFpower'] = -60 + CAV_Attenuation

if len(settings['avg_times']) != len(settings['powers']):
    raise ValueError('incorrect number of averaging times specified')

#set voltage sweep
voltages = np.round(np.linspace(settings['start_voltage'], settings['stop_voltage'], settings['voltage_points']),6)
max_voltage = 3.5
if np.max(voltages) > max_voltage:
    raise ValueError('max voltage too! large')
else:
    settings['voltages'] = voltages

##set file name
#stamp = userfuncs.timestamp()
#scanname = 'transFluxScan_' + settings['scanname'] + '_' + stamp




mags = np.zeros((settings['voltage_points'], settings['sweep_points']))
phases = np.zeros((settings['voltage_points'], settings['sweep_points']))


SRS.Volt = 0
SRS.Output = 'On'



for pind in range(0, len(settings['powers'])):
    power = settings['powers'][pind]
    settings['RFpower'] = power
    settings['avg_time'] = settings['avg_times'][pind]
    identifier = 'Cav Power : ' + str(settings['RFpower'] - CAV_Attenuation) + ' dB'
    
    stamp = userfuncs.timestamp()
    scanname = 'transFluxScan_' + settings['scanname'] + '_Power'+str(settings['RFpower'] - CAV_Attenuation) + '_' + stamp
    
    print('Power: {}'.format(settings['RFpower'] - CAV_Attenuation))
    
    for vind in range(len(voltages)):
        voltage = voltages[vind]
#        print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
        
        SRS_voltage_ramp(voltage)
    #    SRS.Volt = voltage
        time.sleep(0.1)
        t0 = time.time()
        
        data = vna.trans_meas(settings)
        
        vna.autoscale()
        
        mags[vind] = data['mag']
        phases[vind] = data['phase']
        freqs = data['xaxis']  
        
        if vind==0 and pind == 0:
            t1=time.time()
            tdiff = t1-t0
#            ttotal = tdiff*len(voltages)*len(settings['powers'])
            
            timePer = tdiff/settings['avg_time'] #time needed pwer second of avergaing time
            totalAveraging = np.sum(settings['avg_times'])
            estimatedTime = timePer*len(voltages)*totalAveraging
            
            print('    ')
#            print('Single run time: {}, estimated total time: {}'.format(tdiff, ttotal))
            print('Single run time: {}'.format(np.round(tdiff,1)))
            print('    ')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60, 1)) + ' minutes')
            print('estimated time for this scan : ' + str(np.round(estimatedTime/60/60, 2)) + ' hours')
            print('    ')
            
        
        #plot the data as it comes in
        VNAplots.general_VNAplot(freqs, mags[0:vind+1,:], phases[0:vind+1,:], voltages[0:vind+1], scanname, 
                                 xlabel = 'Frequency (GHz)', ylabel = 'Voltage (V)', identifier = identifier,
                                 fig_num = 1)
        
    #end loop over voltages
    userfuncs.SaveFull(save_Dir, scanname, ['mags', 'phases', 'freqs', 'powers', 'CAV_Attenuation'], locals(), expsettings=settings)
    plt.savefig(os.path.join(save_Dir, scanname+'.png'), dpi = 150)

t2 = time.time()
print('Elapsed time: {}'.format(t2-t0))

#return to zero voltage
SRS_voltage_ramp(0)
#SRS.Volt = 0
SRS.Output = 'Off'


 
####Martin's original version
#VNAplots.power_plot(freqs, mags, phases, voltages, scanname, -CAV_Attenuation)

##Alicia's next try
VNAplots.general_VNAplot(freqs, mags, phases, voltages, scanname, 
                         xlabel = 'Frequency (GHz)', ylabel = 'Voltage (V)', identifier = identifier,
                         fig_num = 1)


####manual plotting example
#fig = plt.figure(1,figsize=(13,8))
#fig.clf()
#
#ax = plt.subplot(1,2,1)
#VNAplots.general_colormap_subplot(ax,freqs, voltages, mags)
#plt.xlabel('Frequency (GHz)')
#plt.ylabel('Voltage (V)')
#plt.title('S21 mag')
#
#ax = plt.subplot(1,2,2)
#VNAplots.general_colormap_subplot(ax,freqs, voltages, phases)
#plt.xlabel('Frequency (GHz)')
#plt.ylabel('Voltage (V)')
#plt.title('S21 phase')
#
#plt.suptitle('Filename: {}'.format(scanname))
#plt.show()


#userfuncs.SaveFull(save_Dir, scanname, ['mags', 'phases', 'freqs', 'powers', 'CAV_Attenuation'], locals(), expsettings=settings)
#plt.savefig(os.path.join(save_Dir, scanname+'.png'), dpi = 150)




