# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 15:53:30 2024

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

    settings = {}
    
    #Save location
    settings['scanname']    = 'differential_phase_scan'
    settings['meas_type']   = 'trans_differential_phase_scan'
    #settings['project_dir'] = r'Z:\Data'
    
    #Sweep parameter
#    settings['CAV_Attenuation'] = 30


    settings['start_power'] = -50
    settings['stop_power'] = -48
    settings['power_points'] = 1
    
    settings['avg_times'] = np.array([0.8])
    
    #VNA settings
    settings['channel']  = 1
    settings['avg_time'] = 1
    settings['measurement'] = 'S21'
    settings['measurement_types'] = ['S41', 'S23']
    settings['start_freq']  = 7.576e9-40e6 
    settings['stop_freq']   = 7.576e9+40e6 
    settings['freq_points'] = 501
    settings['ifBW'] = 4e3
    settings['unwrap_phase'] = True

    #Dual_gen settings
    settings['modulation_frequency'] = 20e6
    settings['modulation_amplitude_1'] = 0.001
    settings['modulation_amplitude_2'] = 0.001
    settings['dc_offset_voltage_1'] = 0
    settings['dc_offset_voltage_2'] = 0
    settings['start_phi'] = 0
    settings['stop_phi'] = 180
    settings['phi_points'] = 181

    return settings

def vna_differential_phase_scan(instruments, settings, avg_span, detection_span):
    ##Instruments used
    vna = instruments['VNA']
    Dual_gen = instruments['DCsupply']

    vna.reset()
    Dual_gen.reset()
    time.sleep(1)

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    ##Data saving and naming
    #saveDir = userfuncs.saveDir(settings['project_dir'], settings['meas_type'])
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename_template = exp_settings['scanname'] + '_power{}dBm_' + stamp

    CAV_Attenuation = exp_globals['CAV_Attenuation']
    
    #set differential phase sweep
    start_phi  = exp_settings['start_phi']
    stop_phi   = exp_settings['stop_phi']
    phi_points = exp_settings['phi_points']
    phis = np.round(np.linspace(start_phi, stop_phi, phi_points), 3)

    
    #optional power sweep
    start_power  = exp_settings['start_power'] + CAV_Attenuation
    stop_power   = exp_settings['stop_power'] + CAV_Attenuation
    power_points = exp_settings['power_points']
    powers = np.linspace(start_power, stop_power, power_points)
    
    if len(exp_settings['avg_times']) != len(powers):
        raise ValueError('incorrect number of averaging times specified')
    
    # I know the elegant way of doing this is to use dictionary and for loop,
    # but here I just want have somehting that can quickly work

    if 'S41' in exp_settings['measurement_types']:
        S41_mags = np.zeros((phi_points, exp_settings['freq_points']))
        S41_phases = np.zeros((phi_points, exp_settings['freq_points']))
        S41_mags_unmodulated = np.zeros((power_points, exp_settings['freq_points']))
        S41_phases_unmodulated = np.zeros((power_points, exp_settings['freq_points']))

    if 'S23' in exp_settings['measurement_types']:
        S23_mags = np.zeros((phi_points, exp_settings['freq_points']))
        S23_phases = np.zeros((phi_points, exp_settings['freq_points']))
        S23_mags_unmodulated = np.zeros((power_points, exp_settings['freq_points']))
        S23_phases_unmodulated = np.zeros((power_points, exp_settings['freq_points']))

    if 'S21' in exp_settings['measurement_types']:
        S21_mags = np.zeros((phi_points, exp_settings['freq_points']))
        S21_phases = np.zeros((phi_points, exp_settings['freq_points']))
        S21_mags_unmodulated = np.zeros((power_points, exp_settings['freq_points']))
        S21_phases_unmodulated = np.zeros((power_points, exp_settings['freq_points']))

    if 'S43' in exp_settings['measurement_types']:
        S43_mags = np.zeros((phi_points, exp_settings['freq_points']))
        S43_phases = np.zeros((phi_points, exp_settings['freq_points']))
        S43_mags_unmodulated = np.zeros((power_points, exp_settings['freq_points']))
        S43_phases_unmodulated = np.zeros((power_points, exp_settings['freq_points']))

    transmission_mag_offset = np.zeros(power_points)
    center_isolation_contrast_raw = np.zeros((power_points, phi_points))
    center_isolation_contrast_calibrated = np.zeros((power_points, phi_points))

    Dual_gen.phase_unit = 'DEG'
    time.sleep(0.1)
    
    Dual_gen.Ch1_dc_voltage_ramp(0)
    Dual_gen.Ch1_output = 1
    Dual_gen.Ch1_dc_voltage_ramp(exp_settings['dc_offset_voltage_1'])

    Dual_gen.Ch2_dc_voltage_ramp(0)
    Dual_gen.Ch2_output = 1
    Dual_gen.Ch2_dc_voltage_ramp(exp_settings['dc_offset_voltage_2'])
    

    t0 = time.time()
    for pind in range(len(powers)):
        power = powers[pind]
        exp_settings['RFpower']  = power
        exp_settings['avg_time'] = exp_settings['avg_times'][pind]
        identifier = 'Cav Power : ' + str(exp_settings['RFpower'] - CAV_Attenuation) + ' dB'
        
        stamp    = userfuncs.timestamp()
        filename = filename_template.format(power - CAV_Attenuation)
        
        print('Power: {}'.format(exp_settings['RFpower'] - CAV_Attenuation))

        if 'S41' in exp_settings['measurement_types']:

            exp_settings['measurement'] = 'S41'
            data_S41_unmodulated = vna.trans_meas(exp_settings)
            vna.autoscale()

            S41_mags_unmodulated[pind] = data_S41_unmodulated['mag']
            S41_phases_unmodulated[pind] = data_S41_unmodulated['phase']

        if 'S23' in exp_settings['measurement_types']:

            exp_settings['measurement'] = 'S23'
            data_S23_unmodulated = vna.trans_meas(exp_settings)
            vna.autoscale()

            S23_mags_unmodulated[pind] = data_S23_unmodulated['mag']
            S23_phases_unmodulated[pind] = data_S23_unmodulated['phase']

        if 'S21' in exp_settings['measurement_types']:

            exp_settings['measurement'] = 'S21'
            data_S21_unmodulated = vna.trans_meas(exp_settings)
            vna.autoscale()

            S21_mags_unmodulated[pind] = data_S21_unmodulated['mag']
            S21_phases_unmodulated[pind] = data_S21_unmodulated['phase']

        if 'S43' in exp_settings['measurement_types']:

            exp_settings['measurement'] = 'S43'
            data_S43_unmodulated = vna.trans_meas(exp_settings)
            vna.autoscale()

            S43_mags_unmodulated[pind] = data_S43_unmodulated['mag']
            S43_phases_unmodulated[pind] = data_S43_unmodulated['phase']

        index_center_freq = int((len(data_S41_unmodulated['xaxis'])-1)/2)
        ydata_1_sliced = S41_mags_unmodulated[pind][index_center_freq-int(avg_span/2):index_center_freq+int(avg_span/2)] 
        ydata_2_sliced = S23_mags_unmodulated[pind][index_center_freq-int(avg_span/2):index_center_freq+int(avg_span/2)]
        
        transmission_mag_offset[pind] = np.mean(ydata_1_sliced - ydata_2_sliced)
        
        t1 = time.time() #to account for the time taken by the calibration traces

        for phi_ind in range(len(phis)):
            phi = phis[phi_ind]
            print('phi: {}, final phi: {}'.format(phi, phis[-1]))

            Dual_gen.Ch1_sin_gen(exp_settings['modulation_amplitude_1'], exp_settings['modulation_frequency'], phase=0, offset=exp_settings['dc_offset_voltage_1'])
            Dual_gen.Ch2_sin_gen(exp_settings['modulation_amplitude_2'], exp_settings['modulation_frequency'], phase=phi, offset=exp_settings['dc_offset_voltage_2'])
            Dual_gen.phase_sync()
            time.sleep(0.1)

            if 'S41' in exp_settings['measurement_types']:
                exp_settings['measurement'] = 'S41'
                S41_data = vna.trans_meas(exp_settings)
                vna.autoscale()
            
                S41_mags[phi_ind]   = S41_data['mag']
                S41_phases[phi_ind] = S41_data['phase']

                freqs = S41_data['xaxis']  
            
            if 'S23' in exp_settings['measurement_types']:
                exp_settings['measurement'] = 'S23'
                S23_data = vna.trans_meas(exp_settings)
                vna.autoscale()
            
                S23_mags[phi_ind]   = S23_data['mag']
                S23_phases[phi_ind] = S23_data['phase']

                freqs = S23_data['xaxis']

            if 'S21' in exp_settings['measurement_types']:
                exp_settings['measurement'] = 'S21'
                S21_data = vna.trans_meas(exp_settings)
                vna.autoscale()
            
                S21_mags[phi_ind]   = S21_data['mag']
                S21_phases[phi_ind] = S21_data['phase']

                freqs = S21_data['xaxis']

            if 'S43' in exp_settings['measurement_types']:
                exp_settings['measurement'] = 'S43'
                S43_data = vna.trans_meas(exp_settings)
                vna.autoscale()
            
                S43_mags[phi_ind]   = S43_data['mag']
                S43_phases[phi_ind] = S43_data['phase']

                freqs = S43_data['xaxis']
            
            if phi_ind==0 and pind == 0:
                t2=time.time()

                avg_interval_1 = (t1 - t0)*np.mean(exp_settings['avg_times'])/exp_settings['avg_times'][0]
                avg_interval_2 = (t2 - t1)*np.mean(exp_settings['avg_times'])/exp_settings['avg_times'][0]
                avg_time_per_power = avg_interval_2*len(phis) + avg_interval_1
                avg_total_time = np.round((len(powers)*avg_time_per_power)/60, 1) # in mins

                print('Estimated time for the scan is ' + str(avg_total_time) + ' mins')
                
            S41_mag_sliced = S41_mags[phi_ind][index_center_freq-int(detection_span/2):index_center_freq+int(detection_span/2)]
            S23_mag_sliced = S23_mags[phi_ind][index_center_freq-int(detection_span/2):index_center_freq+int(detection_span/2)]
            freqs_sliced = freqs[index_center_freq-int(detection_span/2):index_center_freq+int(detection_span/2)]

            if S41_mag_sliced[int((len(S41_mag_sliced)-1)/2)] > S23_mag_sliced[int((len(S23_mag_sliced)-1)/2)]:
                index_S23_min = np.argmin(S23_mag_sliced) # detecting where the S23 has the minimum transmission

                center_isolation_contrast_raw[pind, phi_ind] = S41_mag_sliced[index_S23_min] - S23_mag_sliced[index_S23_min]
                center_isolation_contrast_calibrated[pind, phi_ind] = center_isolation_contrast_raw[pind, phi_ind] - transmission_mag_offset[pind]

             

            else:
                index_S41_min = np.argmin(S41_mag_sliced) # detecting where the S41 has the minimum transmission

                center_isolation_contrast_raw[pind, phi_ind] = S41_mag_sliced[index_S41_min] - S23_mag_sliced[index_S41_min]
                center_isolation_contrast_calibrated[pind, phi_ind] = center_isolation_contrast_raw[pind, phi_ind] - transmission_mag_offset[pind]

                  

            full_data = {}
            full_data['xaxis']  = freqs/1e9
            
            if 'S41' in exp_settings['measurement_types']:
                full_data['S41_mags']   = S41_mags[0:phi_ind+1]
                full_data['S41_phases'] = S41_phases[0:phi_ind+1]
                full_data['S41_mags_unmodulated'] = S41_mags_unmodulated[0:pind+1]
                full_data['S41_phases_unmodulated'] = S41_phases_unmodulated[0:pind+1]


            if 'S23' in exp_settings['measurement_types']:
                full_data['S23_mags']   = S23_mags[0:phi_ind+1]
                full_data['S23_phases'] = S23_phases[0:phi_ind+1]
                full_data['S23_mags_unmodulated'] = S23_mags_unmodulated[0:pind+1]
                full_data['S23_phases_unmodulated'] = S23_phases_unmodulated[0:pind+1]

            if 'S21' in exp_settings['measurement_types']:
                full_data['S21_mags']   = S21_mags[0:phi_ind+1]
                full_data['S21_phases'] = S21_phases[0:phi_ind+1]
                full_data['S21_mags_unmodulated'] = S21_mags_unmodulated[0:pind+1]
                full_data['S21_phases_unmodulated'] = S21_phases_unmodulated[0:pind+1]

            if 'S43' in exp_settings['measurement_types']:
                full_data['S43_mags']   = S43_mags[0:phi_ind+1]
                full_data['S43_phases'] = S43_phases[0:phi_ind+1]
                full_data['S43_mags_unmodulated'] = S43_mags_unmodulated[0:pind+1]
                full_data['S43_phases_unmodulated'] = S43_phases_unmodulated[0:pind+1]

            
            full_data['center_isolation_contrast_raw'] = center_isolation_contrast_raw
            full_data['center_isolation_contrast_calibrated'] = center_isolation_contrast_calibrated
            full_data['transmission_mag_offset'] = transmission_mag_offset



            fig = plt.figure(25, figsize=(12,9))
            plt.clf()

            # Adjust subplot spacing
            plt.subplots_adjust(hspace=0.5, wspace=0.3)

            phi_limits = phis[0:phi_ind+1]

            ax = plt.subplot(2,2,1)
            plt.imshow(full_data['S41_mags'], extent=[full_data['xaxis'][0], full_data['xaxis'][-1], phi_limits[0], phi_limits[-1]], origin='lower', aspect='auto', cmap='viridis')
            plt.colorbar(label='S41 mag (dB)')
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('phi (degree)')
            plt.title('S41_phi_scan')

            ax = plt.subplot(2,2,2)
            plt.imshow(full_data['S23_mags'], extent=[full_data['xaxis'][0], full_data['xaxis'][-1], phi_limits[0], phi_limits[-1]], origin='lower', aspect='auto', cmap='viridis')
            plt.colorbar(label='S23 mag (dB)')
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('phi (degree)')
            plt.title('S23_phi_scan')

            ax = plt.subplot(2,2,3)
            plt.plot(full_data['xaxis'], S41_data['mag'], label='S41')
            plt.plot(full_data['xaxis'], S23_data['mag'], label='S23')
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('Mag (dB)')
            plt.title('Live traces of S41 and S23')
            plt.legend(loc='lower right')

            ax = plt.subplot(2,2,4)
            plt.plot(phi_limits, center_isolation_contrast_calibrated[pind][0:phi_ind+1], label='S41-S23')
            plt.xlabel('phi (degree)')
            plt.ylabel('S41-S23 (dB)')
            plt.title('Center isolation contrast calibrated')

            
            # Set supertitle with appropriate spacing
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            plt.suptitle('Filename: {}'.format(filename))
    
            fig.canvas.draw()
            fig.canvas.flush_events()

            
            
            
            userfuncs.SaveFull(saveDir, filename, ['full_data', 'powers', 'phis', 'filename'], 
                                locals(), expsettings=settings, instruments=instruments)
            plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-t0))
    
    #return to zero voltage
    Dual_gen.Ch1_sin_gen(0.001, exp_settings['modulation_frequency'], phase=0, offset=0)
    Dual_gen.Ch2_sin_gen(0.001, exp_settings['modulation_frequency'], phase=0, offset=0)
    Dual_gen.Ch1_output = 0
    Dual_gen.Ch2_output = 0
    vna.reset()

    data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data, 'powers': powers, 'phis': phis}

    return data
