# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:19:18 2024

@author: Kollarlab
"""

exp_globals = {}

#Project naming
exp_globals['root_folder'] = r'K:\Data\Oxford'
exp_globals['project_name'] = 'PDH'
exp_globals['device_name'] = 'PDH_v5'

#HW config
exp_globals['CAV_Attenuation'] = 20 #0
exp_globals['Qbit_Attenuation'] = 20 #0
exp_globals['Gain'] = '70dB'

#SW config
exp_globals['trigger_rate'] = 40e3
exp_globals['hanger'] = True
exp_globals['IF'] = 1e6

#for safety try to set the trigger period
try:
    triggergen.Freq = exp_globals['trigger_rate']
except:
    print('WARNING: Trigger period not set. Generator may not exist yet.')


#Pulse config
measurement_pulse = {}
measurement_pulse['meas_pos'] = 150e-6
measurement_pulse['init_buffer'] = 1e-6
measurement_pulse['emp_delay'] = 1.9e-6
measurement_pulse['meas_window'] = 4.5e-6
measurement_pulse['post_buffer'] = 3e-6

qubit_pulse = {}
qubit_pulse['sigma'] = 4e-9 #4e-9
qubit_pulse['num_sigma'] = 4
qubit_pulse['delay'] = 5e-8
qubit_pulse['piAmp'] = 1
qubit_pulse['hold_time'] = 150e-9

exp_globals['measurement_pulse'] = measurement_pulse
exp_globals['qubit_pulse'] = qubit_pulse

#HW corrections and mixer corrections
mixer_config = {}
mixer_config['axes'] = [1,1]
mixer_config['center'] = [0,0]
mixer_config['phi'] = 0
exp_globals['mixer_config'] = mixer_config

cavitygen_config = {}
#cavitygen_config['Ileak'] = -0.16   #holzworth is cavitygen. There are no IQ settings.
#cavitygen_config['Qleak'] = -0.1
exp_globals['cavitygen_config'] = cavitygen_config

qubitgen_config = {}
qubitgen_config['Ileak'] = 0
qubitgen_config['Qleak'] = 0
exp_globals['qubitgen_config'] = qubitgen_config

try:
    qubitgen.IQ.Imp = 'On'
    qubitgen.IQ.Ileak = exp_globals['qubitgen_config']['Ileak']
    qubitgen.IQ.Qleak = exp_globals['qubitgen_config']['Qleak']
except:
    print('WARNING. Qubitgen IQ settings not applied.')

card_config = {}
card_config['activeChannels'] = [1,2]
card_config['sampleRate'] = 2e9/8
card_config['channelRange'] = 2.5
card_config['timeout'] = 30
card_config['triggerSlope'] = 'Falling'
exp_globals['card_config'] = card_config

hdawg_config = {}
hdawg_config['trigger_slope'] = 'falling'
hdawg_config['samplerate'] = '2.4GHz'
hdawg_config['channelgrouping'] = '1x4'
hdawg_config['amplitude'] = 0.5 #NOTE: THE SGS PORT ONLY TAKES IN +-0.5V MAX      
exp_globals['hdawg_config'] = hdawg_config

#DDC configuration
ddc = {}
ddc['method'] = 'low_pass'
ddc['order'] = 5
ddc['stop_atten'] = 40
ddc['cutoff'] = exp_globals['IF']
ddc['sample_rate'] = card_config['sampleRate']
exp_globals['ddc_config'] = ddc