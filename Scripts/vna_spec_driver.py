from vna_spec import get_default_settings, vna_spec

instruments = {}
instruments['VNA'] = vna_spec

settings = get_default_settings()
#Save location
settings['scanname']    = 'initial_power_scan_q4'
settings['meas_type']   = 'Spec'
settings['project_dir'] = r'Z:\Data\HouckQuadTransmon'

#Sweep parameters
settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10

settings['start_power']  = -20
settings['stop_power']   = 10
settings['power_points'] = 31

#VNA settings
settings['channel'] = 1
settings['avg_time'] = 30
settings['measurement'] = 'S21'
settings['start_freq'] = 3.5e9
settings['stop_freq'] = 4.5e9
settings['freq_points'] = 501
settings['RFpower'] = -25
settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1
settings['CAVpower'] = -55
settings['CAVfreq'] = 8.12555e9
settings['ifBW'] = 2e2

vna_spec(instruments, settings)