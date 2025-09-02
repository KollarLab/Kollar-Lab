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

def vna_trans_JPAPump_scan(instruments, settings):
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
    filename_template = exp_settings['scanname'] + '_VoltBias{}v_' + stamp

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
    
    # pump frequency sweep
    start_pump_freq = exp_settings['start_pump_freq']
    stop_pump_freq = exp_settings['stop_pump_freq']
    pump_freq_points = exp_settings['pump_freq_points']
    pump_freqs  = np.linspace(start_pump_freq, stop_pump_freq, pump_freq_points)

    # RF settings
    rf_power = exp_settings['RFpower'] 
    
    # if len(exp_settings['avg_times']) != len(powers):
    #     raise ValueError('incorrect number of averaging times specified')
    
    mags   = np.zeros((pump_freq_points, pump_power_points))
    phases = np.zeros((pump_freq_points, pump_power_points))
    
    # set bias voltage
    bias_voltage = exp_settings['bias_voltage'] 
    if np.abs(bias_voltage) > 1:
        raise ValueError('SRS voltage too high')
    else:
        pass
    SRS.voltage_ramp(bias_voltage)
    SRS.Output = 'On'
    
    tstart = time.time()

    # get zero pump reference 
    exp_settings['RFpower'] = rf_power + CAV_Attenuation
    exp_settings['avg_time'] = exp_settings['avg_times'][0]
    pump_gen.set_output(0)
    time.sleep(0.2)
    zero_pump_data = vna.trans_meas(exp_settings)
    vna.autoscale()
    zero_pump_mag = zero_pump_data['mag']

    for find in range(len(pump_freqs)):
        
        pump_freq = pump_freqs[find]
        pump_gen.set_frequency(pump_freq) 
        print(f'pump freq: {np.round(pump_freq/1e9,3)} GHz, final pump freq: {np.round(pump_freqs[-1]/1e9,3)} GHz')
        identifier = 'Cav Power : ' + str(exp_settings['RFpower'] - CAV_Attenuation) + ' dB'
        
        stamp    = userfuncs.timestamp()
        filename = filename_template.format(bias_voltage)
        #filename = 'transFluxScan_' + settings['scanname'] + '_Power'+str(exp_settings['RFpower'] - CAV_Attenuation) + '_' + stamp
        
        
        
        for pind in range(len(pump_powers)):
            pump_power = pump_powers[pind]
            
            # set signal core pump power
            pump_gen.set_level(pump_power)
            pump_gen.set_output(1)
            
            print(f'pump power: {pump_power}, final pump power: {pump_powers[-1]}')
            
            t0 = time.time()
            data = vna.trans_meas(exp_settings)
            
            vna.autoscale()
            
            if exp_settings['subtract_pump']:
                mag = data['mag'] - zero_pump_mag
                phase = data['phase']
            else:
                mag = data['mag']
                phase = data['phase']
                
            # Now get the mag and phase data for the resonant point (center)
            center_index = exp_settings['freq_points']//2
            mags[find, pind]   = mag[center_index]
            phases[find, pind] = phase[center_index]
            
            if find==0 and pind == 0:
                tstop=time.time() 
                mean_avg_time = np.mean(exp_settings['avg_times'])/exp_settings['avg_times'][0]
                estimate_time(tstart, tstop, len(pump_powers)*len(pump_freqs)*mean_avg_time)
                
            # freqs = data['xaxis']   

    full_data = {}
    full_data['xaxis']  = pump_freqs/1e9
    full_data['mags']   = mags
    full_data['phases'] = phases

    single_data = data
    single_data['xaxis'] = data['xaxis']/1e9

    labels = ['Freq (GHz)', 'Pump Power (dBm)']
    yaxis  = pump_powers
    plots.simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=2)
    
    userfuncs.SaveFull(saveDir, filename, ['full_data', 'single_data', 'pump_powers', 'pump_freqs', 'filename', 'labels'], 
                        locals(), expsettings=settings) #, instruments=instruments)
    plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))
    
    #return to zero voltage and turn pump gen off 
    SRS.voltage_ramp(0)
    SRS.Output = 'Off'
    pump_gen.set_output(0)


    data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data, 'pump_powers': pump_powers, 'pump_freqs': pump_freqs}

    return data
