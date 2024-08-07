# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 15:53:30 2024

@author: Kollarlab
"""
import time
import os
import numpy as np
import matplotlib.pyplot as plt

import userfuncs

def get_default_settings():

    settings = {}
    
    #Save location
    settings['scanname']    = 'fm_and_phi_scan'
    settings['meas_type']   = 'trans_fm_and_phi_scan'
    #settings['project_dir'] = r'Z:\Data'
    
    #Sweep parameter
#    settings['CAV_Attenuation'] = 30


    settings['power'] = -50

    settings['avg_times'] = 10
    
    
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
    settings['modulation_amplitude_1'] = 0.001
    settings['modulation_amplitude_2'] = 0.001
    settings['dc_offset_voltage_1'] = 0
    settings['dc_offset_voltage_2'] = 0
    settings['start_phi'] = 0
    settings['stop_phi'] = 180
    settings['phi_points'] = 181
    settings['start_fm'] = 19e6
    settings['stop_fm'] = 26e6
    settings['fm_points'] = 40

    return settings

def vna_fm_and_phi_scan(instruments, settings, avg_span, detection_span):
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
    filename_template = exp_settings['scanname'] + '{}MHz_fm_power{}dBm_' + stamp

    CAV_Attenuation = exp_globals['CAV_Attenuation']
    
    #set differential phase sweep
    start_phi  = exp_settings['start_phi']
    stop_phi   = exp_settings['stop_phi']
    phi_points = exp_settings['phi_points']
    phis = np.round(np.linspace(start_phi, stop_phi, phi_points), 3)

    
    #set fm sweep
    start_fm  = exp_settings['start_fm']
    stop_fm   = exp_settings['stop_fm']
    fm_points = exp_settings['fm_points']
    fms = np.linspace(start_fm, stop_fm, fm_points)

    
    # I know the elegant way of doing this is to use dictionary and for loop,
    # but here I just want have somehting that can quickly work

    if 'S41' in exp_settings['measurement_types']:
        S41_mags = np.zeros((phi_points, exp_settings['freq_points']))
        S41_phases = np.zeros((phi_points, exp_settings['freq_points']))
        S41_mags_unmodulated = np.zeros((1, exp_settings['freq_points']))
        S41_phases_unmodulated = np.zeros((1, exp_settings['freq_points']))

    if 'S23' in exp_settings['measurement_types']:
        S23_mags = np.zeros((phi_points, exp_settings['freq_points']))
        S23_phases = np.zeros((phi_points, exp_settings['freq_points']))
        S23_mags_unmodulated = np.zeros((1, exp_settings['freq_points']))
        S23_phases_unmodulated = np.zeros((1, exp_settings['freq_points']))

    if 'S21' in exp_settings['measurement_types']:
        S21_mags = np.zeros((phi_points, exp_settings['freq_points']))
        S21_phases = np.zeros((phi_points, exp_settings['freq_points']))
        S21_mags_unmodulated = np.zeros((1, exp_settings['freq_points']))
        S21_phases_unmodulated = np.zeros((1, exp_settings['freq_points']))

    if 'S43' in exp_settings['measurement_types']:
        S43_mags = np.zeros((phi_points, exp_settings['freq_points']))
        S43_phases = np.zeros((phi_points, exp_settings['freq_points']))
        S43_mags_unmodulated = np.zeros((1, exp_settings['freq_points']))
        S43_phases_unmodulated = np.zeros((1, exp_settings['freq_points']))

    transmission_mag_offset = np.zeros(1)
    center_isolation_contrast_raw = np.zeros((fm_points, phi_points))
    center_isolation_contrast_calibrated = np.zeros((fm_points, phi_points))

    Dual_gen.phase_unit = 'DEG'
    time.sleep(0.1)
    
    Dual_gen.Ch1_dc_voltage_ramp(0)
    Dual_gen.Ch1_output = 1
    Dual_gen.Ch1_dc_voltage_ramp(exp_settings['dc_offset_voltage_1'])

    Dual_gen.Ch2_dc_voltage_ramp(0)
    Dual_gen.Ch2_output = 1
    Dual_gen.Ch2_dc_voltage_ramp(exp_settings['dc_offset_voltage_2'])
    

    # start measurement
    t0 = time.time()

    power = exp_settings['power'] + CAV_Attenuation
    exp_settings['RFpower']  = power
    exp_settings['avg_time'] = exp_settings['avg_times'][0]
    
    if 'S41' in exp_settings['measurement_types']:

        exp_settings['measurement'] = 'S41'
        data_S41_unmodulated = vna.trans_meas(exp_settings)
        vna.autoscale()

        S41_mags_unmodulated[0] = data_S41_unmodulated['mag']
        S41_phases_unmodulated[0] = data_S41_unmodulated['phase']

    if 'S23' in exp_settings['measurement_types']:

        exp_settings['measurement'] = 'S23'
        data_S23_unmodulated = vna.trans_meas(exp_settings)
        vna.autoscale()

        S23_mags_unmodulated[0] = data_S23_unmodulated['mag']
        S23_phases_unmodulated[0] = data_S23_unmodulated['phase']

    if 'S21' in exp_settings['measurement_types']:

        exp_settings['measurement'] = 'S21'
        data_S21_unmodulated = vna.trans_meas(exp_settings)
        vna.autoscale()

        S21_mags_unmodulated[0] = data_S21_unmodulated['mag']
        S21_phases_unmodulated[0] = data_S21_unmodulated['phase']

    if 'S43' in exp_settings['measurement_types']:

        exp_settings['measurement'] = 'S43'
        data_S43_unmodulated = vna.trans_meas(exp_settings)
        vna.autoscale()

        S43_mags_unmodulated[0] = data_S43_unmodulated['mag']
        S43_phases_unmodulated[0] = data_S43_unmodulated['phase']

    index_center_freq = int((len(data_S41_unmodulated['xaxis'])-1)/2)
    ydata_1_sliced = S41_mags_unmodulated[0][index_center_freq-int(avg_span/2):index_center_freq+int(avg_span/2)] 
    ydata_2_sliced = S23_mags_unmodulated[0][index_center_freq-int(avg_span/2):index_center_freq+int(avg_span/2)]
    
    transmission_mag_offset[0] = np.mean(ydata_1_sliced - ydata_2_sliced)

    t1 = time.time() #to account for the time taken by the calibration traces

    for fm_ind, fm in enumerate(fms):

        stamp    = userfuncs.timestamp()
        filename = filename_template.format(fm/1e6, power - CAV_Attenuation)
        
        print('fm: {} MHz'.format(fm/1e6))


        for phi_ind in range(len(phis)):
            phi = phis[phi_ind]
            print('phi: {}, final phi: {}'.format(phi, phis[-1]))

            Dual_gen.Ch1_sin_gen(exp_settings['modulation_amplitude_1'], fm, phase=0, offset=exp_settings['dc_offset_voltage_1'])
            Dual_gen.Ch2_sin_gen(exp_settings['modulation_amplitude_2'], fm, phase=phi, offset=exp_settings['dc_offset_voltage_2'])
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
            
            if phi_ind==0 and fm_ind == 0:
                t2=time.time()

                avg_interval_1 = t1 - t0
                avg_interval_2 = t2 - t1
                avg_time_per_fm = avg_interval_2*len(phis)
                avg_total_time = np.round((len(fms)*avg_time_per_fm + avg_interval_1)/3600, 1) # in hrs

                print('Estimated time for the scan is ' + str(avg_total_time) + ' hrs')
                
            S41_mag_sliced = S41_mags[phi_ind][index_center_freq-int(detection_span/2):index_center_freq+int(detection_span/2)]
            S23_mag_sliced = S23_mags[phi_ind][index_center_freq-int(detection_span/2):index_center_freq+int(detection_span/2)]
            freqs_sliced = freqs[index_center_freq-int(detection_span/2):index_center_freq+int(detection_span/2)]

            if S41_mag_sliced[int((len(S41_mag_sliced)-1)/2)] > S23_mag_sliced[int((len(S23_mag_sliced)-1)/2)]:
                index_S23_min = np.argmin(S23_mag_sliced) # detecting where the S23 has the minimum transmission

                center_isolation_contrast_raw[fm_ind, phi_ind] = S41_mag_sliced[index_S23_min] - S23_mag_sliced[index_S23_min]
                center_isolation_contrast_calibrated[fm_ind, phi_ind] = center_isolation_contrast_raw[fm_ind, phi_ind] - transmission_mag_offset[0]

             

            else:
                index_S41_min = np.argmin(S41_mag_sliced) # detecting where the S41 has the minimum transmission

                center_isolation_contrast_raw[fm_ind, phi_ind] = S41_mag_sliced[index_S41_min] - S23_mag_sliced[index_S41_min]
                center_isolation_contrast_calibrated[fm_ind, phi_ind] = center_isolation_contrast_raw[fm_ind, phi_ind] - transmission_mag_offset[0]

                  

            full_data = {}
            full_data['xaxis']  = freqs/1e9
            
            if 'S41' in exp_settings['measurement_types']:
                full_data['S41_mags']   = S41_mags[0:phi_ind+1]
                full_data['S41_phases'] = S41_phases[0:phi_ind+1]
                full_data['S41_mags_unmodulated'] = S41_mags_unmodulated[0:1]
                full_data['S41_phases_unmodulated'] = S41_phases_unmodulated[0:1]


            if 'S23' in exp_settings['measurement_types']:
                full_data['S23_mags']   = S23_mags[0:phi_ind+1]
                full_data['S23_phases'] = S23_phases[0:phi_ind+1]
                full_data['S23_mags_unmodulated'] = S23_mags_unmodulated[0:1]
                full_data['S23_phases_unmodulated'] = S23_phases_unmodulated[0:1]

            if 'S21' in exp_settings['measurement_types']:
                full_data['S21_mags']   = S21_mags[0:phi_ind+1]
                full_data['S21_phases'] = S21_phases[0:phi_ind+1]
                full_data['S21_mags_unmodulated'] = S21_mags_unmodulated[0:1]
                full_data['S21_phases_unmodulated'] = S21_phases_unmodulated[0:1]

            if 'S43' in exp_settings['measurement_types']:
                full_data['S43_mags']   = S43_mags[0:phi_ind+1]
                full_data['S43_phases'] = S43_phases[0:phi_ind+1]
                full_data['S43_mags_unmodulated'] = S43_mags_unmodulated[0:1]
                full_data['S43_phases_unmodulated'] = S43_phases_unmodulated[0:1]

            
            full_data['center_isolation_contrast_raw'] = center_isolation_contrast_raw[fm_ind]
            full_data['center_isolation_contrast_calibrated'] = center_isolation_contrast_calibrated[fm_ind]
            full_data['transmission_mag_offset'] = transmission_mag_offset[0]



            fig = plt.figure(55, figsize=(12,9))
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
            plt.imshow(center_isolation_contrast_calibrated, extent=[phis[0], phis[-1], fms[0]/1e6, fms[-1]/1e6], origin='lower', aspect='auto', cmap='RdBu_r')
            plt.colorbar(label='Center isolation contrast (dB)')
            plt.xlabel('phi (degree)')
            plt.ylabel('fm (MHz)')
            plt.title('fm and phi sweep')

            
            # Set supertitle with appropriate spacing
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            plt.suptitle('Filename: {}'.format(filename))
    
            fig.canvas.draw()
            fig.canvas.flush_events()

            
            
            
            userfuncs.SaveFull(saveDir, filename, ['full_data', 'fms', 'phis', 'filename'], 
                                locals(), expsettings=settings, instruments=instruments)
            plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
    t2 = time.time()
    print('Elapsed time: {} hrs'.format(np.round((t2-t0)/3600, 2)))
    
    #return to zero voltage
    Dual_gen.Ch1_sin_gen(0.001, 20e6, phase=0, offset=0)
    Dual_gen.Ch2_sin_gen(0.001, 20e6, phase=0, offset=0)
    Dual_gen.Ch1_output = 0
    Dual_gen.Ch2_output = 0
    vna.reset()

    data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data, 'fms': fms, 'phis': phis}

    return data
