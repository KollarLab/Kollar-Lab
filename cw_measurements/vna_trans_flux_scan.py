# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 14:07:57 2020

@author: Kollarlab
"""
import time
import os
from utility.measurement_helpers import estimate_time
import numpy as np
import matplotlib.pyplot as plt

import userfuncs
import utility.plotting_tools as plots

def get_default_settings():
    '''
    get_default_settings _summary_
    '''
    settings = {}
    
    #Save location
    settings['scanname']    = 'flux_Scan'
    settings['meas_type']   = 'trans_flux_scan'
    #settings['project_dir'] = r'Z:\Data'
    settings['DC_mode'] = 'Voltage'
    #Sweep parameter
#    settings['CAV_Attenuation'] = 30

    settings['start_voltage']  = 0.4
    settings['stop_voltage']   = 0.9
    settings['voltage_points'] = 15
    
    settings['start_current'] = -0.02#In amps
    settings['stop_current'] =  0.02
    settings['current_points'] = 251

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
    settings['unwrap_phase'] = False

    return settings

def vna_trans_flux_scan(instruments, settings):
    '''
    vna_trans_flux_scan _summary_

    :param instruments: _description_
    :type instruments: _type_
    :param settings: _description_
    :type settings: _type_
    '''   
    ##Instruments used
    vna = instruments['VNA']
    DCS = instruments['DCsupply']

#    vna.reset() #####!!! warning put me back
    print('WARNING" Auto reset has been overridden')

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    ##Data saving and naming
    #saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename_template = exp_settings['scanname'] + '_power{}dBm_' + stamp

    CAV_Attenuation = exp_globals['CAV_Attenuation']
    
    #set voltage sweep
    # start_voltage  = exp_settings['start_voltage']
    # stop_voltage   = exp_settings['stop_voltage']
    # voltage_points = exp_settings['voltage_points']
    # voltages = np.round(np.linspace(start_voltage, stop_voltage, voltage_points),6)
    # max_voltage = 9.1
    # if np.max(voltages) > max_voltage:
    #     raise ValueError('max voltage too! large')
    # else:
    #     pass
    #set DC sweep
    sweep_steps = 0
    if exp_settings['DC_mode'] == 'Voltage' :
        start_voltage = exp_settings['start_voltage']
        stop_voltage  = exp_settings['stop_voltage']
        voltage_points = exp_settings['voltage_points']
        sweep_steps = voltage_points
        voltages = np.round(np.linspace(start_voltage, stop_voltage, voltage_points),6)
        max_voltage = 10#3.5
        if np.max(voltages) > max_voltage:
            raise ValueError('max voltage too large!')
        else:
            settings['voltages'] = voltages
            DCS.Output = 'On'
    elif exp_settings['DC_mode'] == 'Current' :
        start_current = exp_settings['start_current']
        stop_current  = exp_settings['stop_current']
        current_points = exp_settings['current_points']
        sweep_steps = current_points
        currents = np.round(np.linspace(start_current, stop_current, current_points),7)
        max_current = 0.1#3.5
        if np.max(currents) > max_current:
            raise ValueError('max current too large!')
        else:
            settings['currents'] = currents
            DCS.output = 1
    else:
        print("DC_mode must be either Voltage or Current")
    #optional power sweep
    start_power  = exp_settings['start_power'] + CAV_Attenuation
    stop_power   = exp_settings['stop_power'] + CAV_Attenuation
    power_points = exp_settings['power_points']
    powers = np.linspace(start_power, stop_power, power_points)
    
    if len(exp_settings['avg_times']) != len(powers):
        raise ValueError('incorrect number of averaging times specified')
    
    mags   = np.zeros((sweep_steps, exp_settings['freq_points']))
    phases = np.zeros((sweep_steps, exp_settings['freq_points']))
    
    
    tstart = time.time()
    for pind in range(len(powers)):
        power = powers[pind]
        exp_settings['RFpower']  = power
        exp_settings['avg_time'] = exp_settings['avg_times'][pind]
        identifier = 'Cav Power : ' + str(exp_settings['RFpower'] - CAV_Attenuation) + ' dB'
        
        stamp    = userfuncs.timestamp()
        filename = filename_template.format(power - CAV_Attenuation)
        #filename = 'transFluxScan_' + settings['scanname'] + '_Power'+str(exp_settings['RFpower'] - CAV_Attenuation) + '_' + stamp
        
        print('Power: {}'.format(exp_settings['RFpower'] - CAV_Attenuation))
        
        for vind in range(sweep_steps):
            if exp_settings['DC_mode'] == 'Voltage' :
                voltage = voltages[vind]
                print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
                
                DCS.voltage_ramp(voltage)
                time.sleep(0.1)
            elif exp_settings['DC_mode'] == 'Current' :
                current = currents[vind]
                print('Current: {}, final current: {}'.format(current, currents[-1]))
                
                DCS.current_ramp(current)
                time.sleep(0.1)
            
            vna.reset()  
            vna.output = 'on'
            
            t0 = time.time()
            
            data = vna.trans_meas(exp_settings)
            
            vna.autoscale()
            
            mags[vind]   = data['mag']
            phases[vind] = data['phase']

            freqs = data['xaxis']  
            
            if vind==0 and pind == 0:
                tstop=time.time()
                mean_avg_time = np.mean(exp_settings['avg_times'])/exp_settings['avg_times'][0]
                estimate_time(tstart, tstop, len(powers)*(sweep_steps)*mean_avg_time)
                
            freqs = data['xaxis']   

            full_data = {}
            full_data['xaxis']  = freqs/1e9
            full_data['mags']   = mags[0:vind+1]
            full_data['phases'] = phases[0:vind+1]
    
            single_data = data
            single_data['xaxis'] = freqs/1e9
            
            if exp_settings['DC_mode'] == 'Voltage' :
                labels = ['Freq (GHz)','Voltage (V)']
                yaxis  = voltages[0:vind+1]
            elif exp_settings['DC_mode'] == 'Current' :
                labels = ['Freq (GHz)','Current (A)']
                yaxis  = currents[0:vind+1]
            
            
            plots.simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=2)
            
            if exp_settings['DC_mode'] == 'Voltage' :
                userfuncs.SaveFull(saveDir, filename, ['full_data', 'single_data', 'powers', 'voltages', 'filename', 'labels'], 
                                    locals(), expsettings=settings, instruments=instruments)
            elif exp_settings['DC_mode'] == 'Current' :
                userfuncs.SaveFull(saveDir, filename, ['full_data', 'single_data', 'powers', 'currents', 'filename', 'labels'], 
                                    locals(), expsettings=settings, instruments=instruments)
            plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))
    
    # #return to zero voltage
    # SRS.voltage_ramp(0)
    # SRS.Output = 'Off'
    if exp_settings['DC_mode'] == 'Voltage' :
        data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data, 'powers': powers, 'voltages': voltages}
    elif exp_settings['DC_mode'] == 'Current' :
        data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data, 'powers': powers, 'currents': currents}
    return data
