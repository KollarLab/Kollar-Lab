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


def SRS_voltage_ramp(newV, step_size = 0.001, step_time = 0.02):
    deltaV = newV - SRS.Volt
    numSteps = int(np.abs(np.ceil(deltaV/step_size)))
    vsteps = np.linspace(SRS.Volt, newV, numSteps)
    for vstep in vsteps:
        SRS.Volt = np.round(vstep,6)
        time.sleep(step_time)
    return




CAV_Attenuation = 30

#start_power = -70 + CAV_Attenuation
#stop_power = -18 + CAV_Attenuation
#power_points = 21


settings = vna.trans_default_settings()

settings['scanname'] = 'qubit2_basicScan'

settings['start_voltage'] = -0.150
settings['stop_voltage'] = 0.350
settings['voltage_points'] = 35

settings['channel'] = 1
settings['avg_time'] = 15 #25
settings['measurement'] = 'S21'
settings['start'] = 7.576e9 -40e6
settings['stop'] = 7.576e9 + 40e6
settings['sweep_points'] = 701
settings['RFpower'] = -60 + CAV_Attenuation
#settings['RFpower'] = np.array([-60,-59])  + CAV_Attenuation
settings['ifBW'] = 0.5e3


#set voltage sweep
voltages = np.round(np.linspace(settings['start_voltage'], settings['stop_voltage'], settings['voltage_points']),6)
max_voltage = 3
if np.max(voltages) > 3:
    raise ValueError('max voltage too! large')
else:
    settings['voltages'] = voltages

#set file name
stamp = userfuncs.timestamp()
scanname = 'transFluxScan_' + settings['scanname'] + '_' + stamp




mags = np.zeros((settings['voltage_points'], settings['sweep_points']))
phases = np.zeros((settings['voltage_points'], settings['sweep_points']))


SRS.Volt = 0
SRS.Output = 'On'

t0 = time.time()

identifier = 'Cav Power : ' + str(settings['RFpower'] - CAV_Attenuation) + ' dB'

for vind in range(len(voltages)):
    voltage = voltages[vind]
    print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
    
    SRS_voltage_ramp(voltage)
#    SRS.Volt = voltage
    time.sleep(0.1)
    
    data = vna.trans_meas(settings)
    
    vna.autoscale()
    
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
        
    
    #plot the data as it comes in
#    pulsed_debug(fig, freqs, powers[0:powerind+2], powerdat[0:powerind+1], phasedat[0:powerind+1], powerslice, phaseslice, filename, power)
    VNAplots.general_VNAplot(freqs, mags[0:vind+1,:], phases[0:vind+1,:], voltages[0:vind+1], scanname, 
                         xlabel = 'Frequency (GHz)', ylabel = 'Voltage (V)', identifier = identifier,
                         fig_num = 1)
        
        

t2 = time.time()
print('Elapsed time: {}'.format(t2-t0))

#return to zero voltage
SRS_voltage_ramp(0)
#SRS.Volt = 0
SRS.Output = 'Off'


 

#VNAplots.power_plot(freqs, mags, phases, voltages, scanname, -CAV_Attenuation)

#identifier = 'Cav Power : ' + str(settings['RFpower']) + ' dB'
VNAplots.general_VNAplot(freqs, mags, phases, voltages, scanname, 
                         xlabel = 'Frequency (GHz)', ylabel = 'Voltage (V)', identifier = identifier,
                         fig_num = 1)


userfuncs.SaveFull(save_Dir, scanname, ['mags', 'phases', 'freqs', 'powers', 'CAV_Attenuation'], locals(), expsettings=settings)
plt.savefig(os.path.join(save_Dir, scanname+'.png'), dpi = 150)




