from vna_spec import get_default_settings, vna_spec

instruments = {}
instruments['VNA'] = vna

settings = get_default_settings()

#SRS.output = 'On'    
#SRS.voltage_ramp(0.0)

##Save location
#settings['scanname']    = 'plasmon_recheck'
#settings['meas_type']   = 'Spec'
#settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'
#
##Sweep parameters
#settings['CAV_Attenuation'] = 30
#settings['Qbit_Attenuation'] = 20
#
#settings['start_power']  = -50
#settings['stop_power']   = -40
#settings['power_points'] = 3
#
##VNA settings
#center = 6.925e9
#span = 100e6
#settings['channel'] = 1
#settings['avg_time'] = 20
#settings['measurement'] = 'S21'
#settings['start_freq'] = center - span/2
#settings['stop_freq'] = center + span/2
#settings['freq_points'] = 251
#settings['RFpower'] = -45
#settings['RFport'] = 3
#settings['Mport'] = 2
#settings['CAVport'] = 1
#settings['CAVpower'] = -65
#settings['CAVfreq'] = 6.56063e9
#settings['ifBW'] = 1e3
#
#vna_spec(instruments, settings)
#
#vna.output = 'Off'
#vna.power = -10

settings = get_default_settings()
#Save location
settings['scanname']    = 'fluxon_0mV_line_drift'
settings['meas_type']   = 'Spec'
settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'

#Sweep parameters
settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 20

settings['start_power']  = -20.1
settings['stop_power']   = -20
settings['power_points'] = 121

#VNA settings
center = 3.4715e9
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