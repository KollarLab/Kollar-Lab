# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 18:24:45 2021

@author: Kollarlab
"""
import matplotlib.pyplot as plt
import userfuncs
import os

saveDir = r'Z:\Data\Fluxonium_Raman\CRF01_A3\Spec\20210406'
scanname = 'VNA_optimization_test_{}_points_{}_mode_{}_ifBW_{}_averaging'

average_modes = ['AUTO', 'FLAT', 'RED', 'MOV']
ifBWs = [1e2, 2e2, 5e2, 1e3]
average_times = [10,30,60]
points = [51,101,201,501]

settings = vna.spec_default_settings()

settings['CAV_attenuation'] = 30
settings['Qbit_attenuation'] = 20

settings['CAVpower'] = -60
settings['RFpower'] = -20

center = 3.473e9
span = 50e6
#VNA settings
settings['channel']  = 1
settings['measurement'] = 'S21'
settings['start_freq']  = center - span/2
settings['stop_freq']   = center + span/2

settings['CAVfreq'] = 6.558872e9
settings['CAVpower'] = settings['CAV_attenuation'] + settings['CAVpower']
settings['RFpower'] = settings['Qbit_attenuation'] + settings['RFpower']

plt.figure()
for point in points:
    print('Points: {}'.format(point))
    for a_time in average_times:
        print('Averaging time: {}'.format(a_time))
        for mode in average_modes:
            print('Averaging mode: {}'.format(mode))
            for ifBW in ifBWs:
                print('ifBW:{}'.format(ifBW))
                filename = scanname.format(point, mode, ifBW, a_time)
                settings['ifBW'] = ifBW
                settings['mode'] = mode
                settings['avg_time'] = a_time
                settings['freq_points'] = point
                
                data = vna.spec_meas(settings)
                userfuncs.SaveFull(saveDir, filename, ['data', 'ifBW', 'mode', 'a_time', 'point'], locals(), expsettings=settings)
                plt.clf()
                plt.plot(data['xaxis'], data['mag'])
                plt.xlabel('Freq (GHz)')
                plt.ylabel('Mag (dB)')
                plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)
                

from vna_spec import get_default_settings, vna_spec

settings = get_default_settings()
#Save location
settings['scanname']    = 'fluxon_0mV_line_drift_1'
settings['meas_type']   = 'Spec'
settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'

#Sweep parameters
settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 20

settings['start_power']  = -20.1
settings['stop_power']   = -20
settings['power_points'] = 121

#VNA settings
center = 3.473e9
span = 50e6
settings['channel'] = 1
settings['avg_time'] = 30
settings['measurement'] = 'S21'
settings['start_freq'] = center - span/2
settings['stop_freq'] = center + span/2
settings['freq_points'] = 201
settings['RFpower'] = -45
settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1
settings['CAVpower'] = -65
settings['CAVfreq'] = 6.558872e9
settings['ifBW'] = 1e3
settings['mode'] = 'MOV'

vna_spec(instruments, settings)

vna.output = 'Off'
vna.power = -10

settings = get_default_settings()
#Save location
settings['scanname']    = 'fluxon_0mV_line_drift_2'
settings['meas_type']   = 'Spec'
settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'

#Sweep parameters
settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 20

settings['start_power']  = -20.1
settings['stop_power']   = -20
settings['power_points'] = 121

#VNA settings
center = 3.473e9
span = 50e6
settings['channel'] = 1
settings['avg_time'] = 30
settings['measurement'] = 'S21'
settings['start_freq'] = center - span/2
settings['stop_freq'] = center + span/2
settings['freq_points'] = 201
settings['RFpower'] = -45
settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1
settings['CAVpower'] = -65
settings['CAVfreq'] = 6.558872e9
settings['ifBW'] = 1e3
settings['mode'] = 'MOV'

vna_spec(instruments, settings)

vna.output = 'Off'
vna.power = -10