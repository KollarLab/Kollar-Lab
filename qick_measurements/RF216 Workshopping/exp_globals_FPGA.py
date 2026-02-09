exp_globals = {}

#LO settings
exp_globals['LO']       = True #Set to False for debugging purposes
exp_globals['LO_freq']  = 6e9
exp_globals['LO_power'] = 22

#Channel config
exp_globals['cav_channel'] = {'ID' : 6, 'Atten_1' : 30, 'Atten_2': 30}
exp_globals['qub_channel'] = {'ID' : 5, 'Atten_1' : 15, 'Atten_2': 30}
exp_globals['ro_channel']   = {'ID' : 0, 'Atten' : 30} #Might not work
# exp_globals['ro_channels']   = [0]

#Project naming
# exp_globals['root_folder'] =  r'Z:\Users\Kellen'  
# exp_globals['project_name'] = 'FPGA_Control\TestsWithNDAvg_2023_05'
# exp_globals['device_name'] = 'Loopback'

exp_globals['root_folder'] = r'Z:\Data'
exp_globals['project_name'] = 'Hot_Qubits'
exp_globals['device_name'] = 'TXQ_v1'

#HW config (this currently does NOTHING)
#exp_globals['CAV_Attenuation']  = 0
#exp_globals['Qbit_Attenuation'] = 0
#exp_globals['Input_filters']    = '' 
#exp_globals['Output_filters']   =  ''
#exp_globals['Gain']             = ''

#Clock config
exp_globals['relax_delay'] = 100 

#Pulse config 
measurement_pulse = {} #ALL IN MICROSECONDS
measurement_pulse['meas_pos']    = 80
measurement_pulse['init_buffer'] = 0.2
measurement_pulse['emp_delay']   = 0.44 #0.8
measurement_pulse['meas_window'] = 2
measurement_pulse['post_buffer'] = 0.5
measurement_pulse['cav_phase']   = 30 #0
#ADC always triggers at zero. If the ADC and DAC were both triggered at t=0 with reference to the master clock,
#the DAC takes more time to set up its pulse. This additional time, along with travel time for the pulse, is factored
#into the emp_delay. One difference from old code is that this is in us not s
#Configuration notes:
#meas_pulse: t = meas_pos + init_buffer, convert from us to clock cycles
#meas_start = init_buffer, convert from us
#adc_trig_offset = emp_delay
#readout_length = init_buffer + meas_window + post_buffer



qubit_pulse = {} #ALL IN MICROSECONDS
qubit_pulse['sigma'] = 0.04 #0.02 
qubit_pulse['num_sigma'] = 4
qubit_pulse['delay'] = 0.1
qubit_pulse['piAmp'] = 1 #?
qubit_pulse['hold_time'] = 10#20 #THIS ONLY AFFECTS CERTAIN SCRIPTS. Pulsed spec, for example, is always just a Gaussian
qubit_pulse['qub_phase'] = 0
#qubit_pulse['mode']      = 'None'


exp_globals['measurement_pulse'] = measurement_pulse
exp_globals['qubit_pulse'] = qubit_pulse