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
    settings['meas_type']   = 'multi_trans_flux_scan'
    
    settings['start_voltage']  = 0.4
    settings['stop_voltage']   = 0.9
    settings['voltage_points'] = 15

    settings['start_pump_power'] = -50
    settings['stop_pump_power'] = -48
    settings['power_points'] = 1

    settings['start_pump_freq'] = 11e9
    settings['start_pump_freq'] = 13e9
    settings['pump_freq_points'] = 10

    settings['subtract_Pump'] = 'On'

    
    #VNA settings
    settings['channel']  = 1
    settings['avg_time'] = 1
    settings['measurement'] = 'S21'
    settings['start_freq']  = 7.576e9-40e6 
    settings['stop_freq']   = 7.576e9+40e6 
    settings['CAVpower'] = -50
    settings['freq_points'] = 501
    settings['ifBW'] = 4e3
    settings['unwrap_phase'] = False
    settings['eletrical_length'] = 30

    return settings






def vna_trans_JPAPumpPower_scan(instruments, settings):
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

    # offset embed
    electrical_length = exp_settings['electrical_length']
    vna.inst.write('SENS1:CORR:EDEL1:ELEN {}'.format(electrical_length))
   
    ##Data saving and naming
    #saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename_template = exp_settings['scanname'] + '_Vbias_{}V_PumpFreq_{}GHz' + stamp

    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = exp_settings['CAVpower']



    # pump voltage sweep 
    start_voltage  = exp_settings['start_voltage']
    stop_voltage   = exp_settings['stop_voltage']
    voltage_points = exp_settings['voltage_points']
    voltages = np.round(np.linspace(start_voltage, stop_voltage, voltage_points),6)
    if np.max(voltages) > 1:
        raise ValueError('SRS voltage too high')
    else:
        pass


    # set pump power settings 
    pump_amp = exp_globals['pump_amp']
    line_attenuation = exp_globals['line_attenuation']
    
    start_pump_power  = exp_settings['start_pump_power']
    stop_pump_power   = exp_settings['stop_pump_power']
    pump_power_points = exp_settings['power_points']
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
    # start_pump_freq = exp_settings['start_pump_freq']
    # stop_pump_freq = exp_settings['stop_pump_freq']
    # pump_freq_points = exp_settings['pump_freq_points']
    # pump_freqs = np.linspace(start_pump_freq, stop_pump_freq, pump_freq_points)
    pump_freq = exp_settings['pump_freq']   
    pump_gen.set_frequency(pump_freq)

    # start settings
    SRS.voltage_ramp(0)
    SRS.output = 'On'
    pump_gen.set_output(0)
    
    tstart = time.time()

    for vind in range(len(voltages)):
        mags = np.zeros((pump_power_points, exp_settings['freq_points']))
        phases = np.zeros((pump_power_points, exp_settings['freq_points']))

        voltage = voltages[vind]
        SRS.voltage_ramp(voltage)
        time.sleep(0.1)
        exp_settings['RFpower'] = CAV_power + CAV_Attenuation
        exp_settings['avg_time'] = exp_settings['avg_time']
        identifier = 'Cav Power : ' + str(exp_settings['RFpower'] - CAV_Attenuation) + ' dB'

        stamp    = userfuncs.timestamp()
        filename = filename_template.format(voltage, pump_freq/1e9)
        #filename = 'transFluxScan_' + settings['scanname'] + '_Power'+str(exp_settings['RFpower'] - CAV_Attenuation) + '_' + stamp
        
        print('voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
        # print('Power: {}'.format(exp_settings['RFpower'] - CAV_Attenuation))

        # get the pump off data to substract
        pump_gen.set_output(0)
        data_off = vna.trans_meas(exp_settings)
        vna.autoscale()
        pump_gen.set_output(1)

        # this code is assuming 1 pump power point, later erase the pump power for loop.  
        for pind in range(len(pump_powers)):
            pump_power = pump_powers[pind]
            pump_gen.set_level(pump_power)
            print('pump power: {}, final pump power: {}'.format(pump_power, pump_powers[-1]))

            t0 = time.time()
            data = vna.trans_meas(exp_settings)   
            vna.autoscale()

            # get data for pump off so that gain is calculated
            if exp_settings['subtract_Pump'] == 'On':
                #print(data['mag'].shape)
                mags[pind]   = data['mag'] - data_off['mag']
                phases[pind] = data['phase'] - data_off['phase']
            else:
                mags[pind]   = data['mag']
                phases[pind] = data['phase']
            
            # if find==0 and pind == 0:
            #     tstop=time.time() 
            #     mean_avg_time = exp_settings['avg_time']
            #     estimate_time(tstart, tstop, power_points*pump_freq_points*mean_avg_time)
                
            freqs = data['xaxis']   

            full_data = {}
            full_data['xaxis']  = freqs/1e9
            full_data['mags']   = mags[0:pind+1]
            full_data['phases'] = phases[0:pind+1]
    
            single_data = data
            single_data['xaxis'] = data['xaxis']/1e9
    
            labels = ['Freq (GHz)', 'Pump power (dBm)']
            yaxis  = pump_powers[0:pind+1]
            plots.simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=21)
            
            userfuncs.SaveFull(saveDir, filename, ['full_data', 'single_data', 'pump_powers', 'pump_freq', 'filename', 'labels'], 
                                locals(), expsettings=settings) #, instruments=instruments)
            plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))
    
    #return to zero voltage and turn pump gen off 
    SRS.voltage_ramp(0)
    SRS.Output = 'Off'
    pump_gen.set_output(0)


    data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data, 'pump_powers': pump_powers, 'voltages': voltages, 'pump_freq': pump_freq}

    return data
