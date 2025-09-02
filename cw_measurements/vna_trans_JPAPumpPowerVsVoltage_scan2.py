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
    
    #Sweep parameter
#    settings['CAV_Attenuation'] = 30

    settings['start_voltage']  = 0.4
    settings['stop_voltage']   = 0.9
    settings['voltage_points'] = 15

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

def vna_trans_JPAPumpPowerVsVoltage_scan2(instruments, settings):
    '''
    vna_trans_flux_scan _summary_

    :param instruments: _description_
    :type instruments: _type_
    :param settings: _description_
    :type settings: _type_
    '''   
    ##Instruments used
    vna = instruments['VNA']
    SRS = instruments['DCsupply']
    pump_gen = instruments['pump_gen'] 

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    vna.reset() 

    if exp_settings['Electrical Length'][0]: 
        electrical_length = exp_settings['Electrical Length'][1]
        el_str = f'SENS1:CORR:EDEL1:ELEN {electrical_length}'
        vna.inst.write(el_str)
   
    ##Data saving and naming
    #saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename_template = exp_settings['scanname'] + '_Pump{}GHz_' + stamp

    CAV_Attenuation = exp_globals['CAV_Attenuation']

    # set pump power settings 
    pump_amp = exp_settings['pump_amp']
    line_attenuation = exp_settings['line_attenuation']
    
    start_pump_power  = exp_settings['start_pump_power']
    stop_pump_power   = exp_settings['stop_pump_power']
    pump_power_points = exp_settings['pump_power_points']
    pump_powers = np.round(np.linspace(start_pump_power - pump_amp - line_attenuation, stop_pump_power - pump_amp - line_attenuation, pump_power_points),6)
    max_power = 13
    min_power = -25
    if np.max(pump_powers) > max_power:
        raise ValueError('max power too large!')
    elif np.min(pump_powers) < min_power:
        raise ValueError('min power too low!')
    else:
        pass

    # set pump frequency
    pump_freq = 2 * exp_settings['CAV_freq'] 
    pump_gen.set_frequency(pump_freq)

    # pump voltage sweep 
    start_voltage  = exp_settings['start_voltage']
    stop_voltage   = exp_settings['stop_voltage']
    voltage_points = exp_settings['voltage_points']
    voltages = np.round(np.linspace(start_voltage, stop_voltage, voltage_points),6)
    if np.max(voltages) > 1:
        raise ValueError('SRS voltage too high')
    else:
        pass


    # RF settings
    rf_power = exp_settings['RFpower'] 

    
    mags   = np.zeros((voltage_points, exp_settings['freq_points']))
    phases = np.zeros((voltage_points, exp_settings['freq_points']))
    
    tstart = time.time()

    for vind in range(len(voltages)):
        voltage = voltages[vind]
        SRS.voltage_ramp(voltage)
        SRS.Output = 'On'
        exp_settings['RFpower'] = rf_power + CAV_Attenuation
        exp_settings['avg_time'] = exp_settings['avg_times'][0]
        identifier = 'Cav Power : ' + str(exp_settings['RFpower'] - CAV_Attenuation) + ' dB'

        stamp    = userfuncs.timestamp()
        filename = filename_template.format(pump_freq/1e9)
        #filename = 'transFluxScan_' + settings['scanname'] + '_Power'+str(exp_settings['RFpower'] - CAV_Attenuation) + '_' + stamp
        
        print(f'voltage: {voltage}, final voltage: {voltages[-1]}')
        # print('Power: {}'.format(exp_settings['RFpower'] - CAV_Attenuation))
        
        for pind in range(len(pump_powers)):
            pump_power = pump_powers[pind]
            pump_gen.set_level(pump_power)
            pump_gen.set_output(0)
            
            t0 = time.time()
            data = vna.trans_meas(exp_settings)   
            vna.autoscale()

            # get data for pump off so that gain is calculated
            if exp_settings['subtract_PumpOff']:
                pump_gen.set_output(0)
                data_off = vna.trans_meas(exp_settings)
                vna.autoscale()
                mags[pind]   = data['mag'] - data_off['mag']
            else:
                mags[pind]   = data['mag']
            phases[pind] = data['phase']
            
            if vind==0 and pind == 0:
                tstop=time.time() 
                mean_avg_time = np.mean(exp_settings['avg_times'])/exp_settings['avg_times'][0]
                estimate_time(tstart, tstop, len(pump_powers)*len(voltages)*mean_avg_time)
                
            # freqs = data['xaxis']   

            full_data = {}
            full_data['xaxis']  = data['xaxis']# pump_powers#voltages #pump_freqs/1e9
            full_data['mags']   = mags[0:vind+1]
            full_data['phases'] = phases[0:vind+1]
    
            single_data = data
            single_data['xaxis'] = data['xaxis']/1e9
    
            labels = ['Pump Power (dBm)', 'Voltage (v)']
            yaxis  = voltages[0:vind+1]
            plots.simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=2)
            
            userfuncs.SaveFull(saveDir, filename, ['full_data', 'single_data', 'pump_powers', 'voltages', 'filename', 'labels'], 
                                locals(), expsettings=settings) #, instruments=instruments)
            plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))
    
    #return to zero voltage and turn pump gen off 
    SRS.voltage_ramp(0)
    SRS.Output = 'Off'
    pump_gen.set_output(0)


    data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data, 'pump_powers': pump_powers, 'voltages': voltages}

    return data
