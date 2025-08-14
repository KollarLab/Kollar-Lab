#import numpy as np

exp_globals = {}

#LO settings
exp_globals['LO']       = True #Set to False for debugging purposes
exp_globals['LO_freq']  = 0
exp_globals['LO_power'] = 22

#Channel config
exp_globals['cav_channel']   = 4 #4 is correct
exp_globals['qub_channel'] = 6 #6 #6 is correct
exp_globals['ro_channels']   = [0]

#Project naming
# exp_globals['root_folder'] =  r'K:\Data'  
# exp_globals['project_name'] = 'FPGA_Loopback'
# exp_globals['device_name'] = 'CW_Testing'

exp_globals['root_folder'] = r'K:\Data'
exp_globals['project_name'] = 'FPGA_Loopback'
exp_globals['device_name'] = 'CW_Flux_Testing'

exp_globals['hanger'] = True #This now does something for spec_flux_scan

#HW config (this currently does NOTHING)
#exp_globals['CAV_Attenuation']  = 0
#exp_globals['Qbit_Attenuation'] = 0
#exp_globals['Input_filters']    = '' 
#exp_globals['Output_filters']   =  ''
#exp_globals['Gain']             = ''

#Clock config
exp_globals['relax_delay'] = 60 

#Pulse config 
measurement_pulse = {} #ALL IN MICROSECONDS 
measurement_pulse['meas_pos']    = 45
measurement_pulse['init_buffer'] = 0.2
measurement_pulse['emp_delay']   = 0.48 + 0.046 #0.48+0.046
measurement_pulse['meas_window'] = 6 ##900 #20 #60#20#20#60 #6
measurement_pulse['post_buffer'] = 0.4
measurement_pulse['cav_phase']   = 0
measurement_pulse['side_buff']   = 0


#ADC always triggers at zero. If the ADC and DAC were both triggered at t=0 with reference to the master clock,
#the DAC takes more time to set up its pulse. This additional time, along with travel time for the pulse, is factored
#into the emp_delay. One difference from old code is that this is in us not s
#Configuration notes:
#meas_pulse: t = meas_pos + init_buffer, convert from us to clock cycles
#meas_start = init_buffer, convert from us
#adc_trig_offset = emp_delay
#readout_length = init_buffer + meas_window + post_buffer



qubit_pulse = {} #ALL IN MICROSECONDS
qubit_pulse['sigma'] = 0.01
qubit_pulse['num_sigma'] = 4
qubit_pulse['delay'] = 0.3
qubit_pulse['piAmp'] = 1 #?
qubit_pulse['Hold_time'] = 0 #THIS ONLY AFFECTS CERTAIN SCRIPTS. Pulsed spec, for example, is always just a Gaussian
qubit_pulse['qub_phase'] = 0
#qubit_pulse['mode']      = 'None'





exp_globals['measurement_pulse'] = measurement_pulse
exp_globals['qubit_pulse'] = qubit_pulse
