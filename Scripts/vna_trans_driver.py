from vna_trans import get_default_settings, vna_trans

instruments = {}
instruments['VNA'] = vna

settings = get_default_settings()

SRS.output = 'On'    
SRS.voltage_ramp(-0.50)

#Save location
settings['scanname']    = 'cavity1_trans_sweep'
settings['meas_type']   = 'Trans'
settings['project_dir'] = r'Z:\Data\HouckDualHangerFluxonium'
    
#Sweep parameter
settings['CAV_attenuation'] = 30

settings['start_power']  = -60
settings['stop_power']   = -50
settings['power_points'] = 5

center = 8.5e9
span = 600e6
#VNA settings
settings['channel']  = 1
settings['avg_time'] = 30
settings['measurement'] = 'S21'
settings['start_freq']  = 8.65e9
settings['stop_freq']   = 8.7e9
settings['freq_points'] = 201
settings['ifBW'] = 1e3

vna_trans(instruments, settings)

SRS.voltage_ramp(0)
SRS.output = 'Off'    