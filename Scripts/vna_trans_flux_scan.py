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
import VNAplottingTools as plots

def get_default_settings():

    settings = {}
    
    #Save location
    settings['scanname']    = 'flux_Scan'
    settings['meas_type']   = 'trans_flux_scan'
    settings['project_dir'] = r'Z:\Data'
    
    #Sweep parameter
    settings['CAV_Attenuation'] = 30

    settings['start_volt']  = 0.4
    settings['stop_volt']   = 0.9
    settings['volt_points'] = 15

    settings['start_power'] = -50
    settings['stop_power'] = -48
    settings['power_points'] = 1
    
    settings['avg_times'] = np.array([0.8])
    
    #VNA settings
    settings['channel']  = 1
    settings['avg_time'] = 1
    settings['measurement'] = 'S21'
    settings['start_freq']  = 7.576e9-40e6 
    settings['stop_freq']   = 7.576e9+40e6 
    settings['freq_points'] = 501
    settings['ifBW'] = 4e3

    return settings

def vna_trans_flux_scan(instruments, settings):
    ##Instruments used
    vna = instruments['VNA']
    SRS = instruments['DCsupply']

    ##Data saving and naming
    saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp = userfuncs.timestamp()
    filename = settings['scanname'] + '_' + stamp

    CAV_Attenuation = settings['CAV_attenuation']
    scanname = settings['scanname']
    
    SRS.Range = '10 V'
    
    #set voltage sweep
    start_voltage = settings['start_voltage']
    stop_voltage  = settings['stop_voltage']
    voltage_points = settings['voltage_points']
    voltages = np.round(np.linspace(settings['start_voltage'], settings['stop_voltage'], settings['voltage_points']),6)
    max_voltage = 3.5
    if np.max(voltages) > max_voltage:
        raise ValueError('max voltage too! large')
    else:
        settings['voltages'] = voltages
    
    #optional power sweep
    start_power = settings['start_power'] + CAV_Attenuation
    stop_power = settings['stop_power'] + CAV_Attenuation
    power_points = settings['power_points']
    powers = np.linspace(start_power, stop_power, power_points)
    
    avg_times = settings['avg_times']
    
    if len(settings['avg_times']) != len(settings['powers']):
        raise ValueError('incorrect number of averaging times specified')
    
    mags = np.zeros((settings['voltage_points'], settings['freq_points']))
    phases = np.zeros((settings['voltage_points'], settings['freq_points']))
    
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
            print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
            
            SRS.voltage_ramp(voltage)
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
                print('Single run time: {}'.format(np.round(tdiff,1)))
                print('    ')
                print('estimated time for this scan : ' + str(np.round(estimatedTime/60, 1)) + ' minutes')
                print('estimated time for this scan : ' + str(np.round(estimatedTime/60/60, 2)) + ' hours')
                print('    ')
                
            freqs = data['xaxis']   

            full_data = {}
            full_data['xaxis'] = freqs
            full_data['mags'] = mags[0:vind+1]
            full_data['phases'] = phases[0:vind+1]
    
            single_data = data
    
            labels = ['Freq (GHz)', 'Voltage (V)']
            yaxis = voltages[0:vind+1]-CAV_Attenuation
            plots.simplescan_plot(full_data, single_data, yaxis, scanname, labels, identifier=identifier, fig_num=2)
    
            userfuncs.SaveFull(saveDir, filename, ['mags', 'phases', 'freqs', 'powers'], locals(), expsettings=settings)
        
#            #plot the data as it comes in
#            plots.general_VNAplot(freqs, mags[0:vind+1,:], phases[0:vind+1,:], voltages[0:vind+1], scanname, 
#                                     xlabel = 'Frequency (GHz)', ylabel = 'Voltage (V)', identifier = identifier,
#                                     fig_num = 1)
            
        #end loop over voltages
        userfuncs.SaveFull(saveDir, scanname, ['mags', 'phases', 'freqs', 'powers', 'CAV_Attenuation'], locals(), expsettings=settings)
        plt.savefig(os.path.join(saveDir, scanname+'.png'), dpi = 150)
    
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-t0))
    
    #return to zero voltage
    SRS.voltage_ramp(0)
    SRS.Output = 'Off'

#    plots.general_VNAplot(freqs, mags, phases, voltages, scanname, 
                             xlabel = 'Frequency (GHz)', ylabel = 'Voltage (V)', identifier = identifier,
                             fig_num = 1)