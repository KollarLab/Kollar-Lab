from vna_trans import get_default_settings, vna_trans

instruments = {}
instruments['VNA'] = vna

settings = get_default_settings()

SRS.output = 'On'    
SRS.voltage_ramp(0.0)

#Save location
settings['scanname']    = 'cavity1_0mV_check'
settings['meas_type']   = 'Trans'
settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'
    
#Sweep parameter
settings['CAV_attenuation'] = 30

settings['start_power']  = -70
settings['stop_power']   = -60
settings['power_points'] = 3

center = 6.56e9
span = 40e6
#VNA settings
settings['channel']  = 1
settings['avg_time'] = 20
settings['measurement'] = 'S21'
settings['start_freq']  = center - span/2
settings['stop_freq']   = center + span/2
settings['freq_points'] = 801
settings['ifBW'] = 1e3
settings['mode'] = 'FLAT'
vna_trans(instruments, settings)

#SRS.voltage_ramp(0)
#SRS.output = 'Off'    