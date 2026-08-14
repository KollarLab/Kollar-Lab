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
    settings['scanname']    = 'flux_Scan'
    settings['meas_type']   = 'JPA_calibration'
    
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



def vna_JPA_calibration(instruments, settings):

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
    filename_template = exp_settings['scanname'] + '_{}V_{}dBm_' + stamp

    # Cavity settings
    CAV_Attenuation = exp_globals['CAV_Attenuation']
    CAV_power = exp_settings['CAVpower']
    exp_settings['RFpower'] = CAV_power + CAV_Attenuation
    exp_settings['avg_time'] = exp_settings['avg_time']
    identifier = 'Cav Power : ' + str(exp_settings['RFpower'] - CAV_Attenuation) + ' dB'

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
    pump_attenuation = exp_globals['pump_attenuation']
    
    start_pump_power  = exp_settings['start_pump_power']
    stop_pump_power   = exp_settings['stop_pump_power']
    power_points = exp_settings['power_points']
    pump_powers = np.round(np.linspace(start_pump_power - pump_amp + pump_attenuation, stop_pump_power - pump_amp + pump_attenuation, power_points),6)
    pump_powers_panel = np.round(np.linspace(start_pump_power, stop_pump_power, power_points),6)
    max_power = 13
    # min_power = -30
    if np.max(pump_powers) > max_power:
        raise ValueError('max power too large!')
    else:
        pass

    # set pump frequency
    start_pump_freq = exp_settings['start_pump_freq']
    stop_pump_freq = exp_settings['stop_pump_freq']
    pump_freq_points = exp_settings['pump_freq_points']
    pump_freqs = np.linspace(start_pump_freq, stop_pump_freq, pump_freq_points)

    # start settings
    SRS.voltage_ramp(0)
    SRS.output = 'On'
    
    tstart = time.time()

    for vind in range(len(voltages)):
        mags = np.zeros((pump_freq_points, exp_settings['freq_points']))
        phases = np.zeros((pump_freq_points, exp_settings['freq_points']))

        SRS.voltage_ramp(voltages[vind])
        time.sleep(0.1)
        
        print('voltage: {}, final voltage: {}'.format(voltages[vind], voltages[-1]))

        # get the pump off data to substract
        pump_gen.Output = 'Off'
        time.sleep(2)
        data_off = vna.trans_meas(exp_settings)
        vna.autoscale()
        

        for pind in range(len(pump_powers)):
            pump_gen.Power = pump_powers[pind]
            print('pump power: {}, final pump power: {}'.format(pump_powers[pind], pump_powers[-1]))
            filename = filename_template.format(voltages[vind], pump_powers_panel[pind])

            for find in range(len(pump_freqs)):
                pump_gen.Freq = pump_freqs[find]
                pump_gen.Output = 'On'
                time.sleep(5)
                
                
                t0 = time.time()
                data = vna.trans_meas(exp_settings)   
                vna.autoscale()

                # subtract pump to get gain
                mags[find]   = data['mag'] - data_off['mag']
                phases[find] = data['phase'] - data_off['phase']
                freqs = data['xaxis']   

                full_data = {}
                full_data['xaxis']  = freqs/1e9
                full_data['mags']   = mags[0:find+1]
                full_data['phases'] = phases[0:find+1]
        
                single_data = {}
                single_data['mag'] = data['mag'] - data_off['mag']
                single_data['phase'] = data['phase'] - data_off['phase']
                single_data['xaxis'] = data['xaxis']/1e9
        
                labels = ['Freq (GHz)', 'Pump freq (GHz)']
                yaxis  = pump_freqs[0:find+1]
                plots.simplescan_plot(full_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=2)
                
                userfuncs.SaveFull(saveDir, filename, ['full_data', 'single_data', 'pump_powers', 'pump_freqs', 'filename', 'labels'], 
                                    locals(), expsettings=settings) #, instruments=instruments)
                plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
            
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))
    
    #return to zero voltage and turn pump gen off 
    SRS.voltage_ramp(0)
    SRS.Output = 'Off'
    pump_gen.Output = 'Off'


    data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data, 'pump_powers': pump_powers, 'voltages': voltages, 'pump_freqs': pump_freqs}

    return data