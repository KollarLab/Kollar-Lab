'''
This script should contains all the big settings in a measurement that do not change
much from run to run (think HW configurations, cavity attenuation etc.). It also contains
the project name and directory so that all the subsequent scripts save in the same
folder. This is an early version so feel free to add more settings to the list. You
will have to go through the other scripts though to implement those changes as well
'''

"""Hello World """

exp_globals = {}

#Project naming
exp_globals['root_folder'] = r'Z:\Data'
exp_globals['project_name'] = 'Qubit_Development'
exp_globals['device_name'] = 'FPS01_01_B2'

#HW config
exp_globals['CAV_Attenuation'] = 20
exp_globals['Qbit_Attenuation'] = 0
exp_globals['Input_filters'] = ''     #'LPF 6.4GHz'
exp_globals['Output_filters'] =  'LPF 8GHz, HPF 4GHz'
exp_globals['Gain'] = '50dB'

#SW config
exp_globals['trigger_rate'] = 1e3
exp_globals['hanger'] = False
exp_globals['IF'] = 1e6

#Pulse config
measurement_pulse = {}
measurement_pulse['meas_pos'] = 50e-6
measurement_pulse['init_buffer'] = 1e-6
measurement_pulse['emp_delay'] = 0#1.3e-6
measurement_pulse['meas_window'] = 10e-6
measurement_pulse['post_buffer'] = 1e-6

qubit_pulse = {}
qubit_pulse['sigma'] = 200e-9
qubit_pulse['num_sigma'] = 4
qubit_pulse['delay'] = 1e-7

exp_globals['measurement_pulse'] = measurement_pulse
exp_globals['qubit_pulse'] = qubit_pulse

#HW corrections and mixer corrections
mixer_config = {}
mixer_config['axes'] = [1,1]
mixer_config['center'] = [0,0]
mixer_config['phi'] = 0
exp_globals['mixer_config'] = mixer_config

cavitygen_config = {}
cavitygen_config['Ileak'] = 0
cavitygen_config['Qleak'] = 0
exp_globals['cavitygen_config'] = cavitygen_config

qubitgen_config = {}
qubitgen_config['Ileak'] = 0
qubitgen_config['Qleak'] = 0
exp_globals['qubitgen_config'] = qubitgen_config

card_config = {}
card_config['activeChannels'] = [1,2]
card_config['sampleRate'] = 2e9/8
card_config['channelRange'] = 2.5
card_config['timeout'] = 40
card_config['triggerSlope'] = 'Falling'
exp_globals['card_config'] = card_config

hdawg_config = {}
hdawg_config['trigger_slope'] = 'falling'
hdawg_config['samplerate'] = '2.4GHz'
hdawg_config['channelgrouping'] = '1x4'
hdawg_config['amplitude'] = 0.5 #NOTE: THE SGS PORT ONLY TAKES IN +-0.5V MAX      
exp_globals['hdawg_config'] = hdawg_config