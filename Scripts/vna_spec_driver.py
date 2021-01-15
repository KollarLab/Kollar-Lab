from vna_spec import get_default_settings, vna_spec

instruments = {}
instruments['VNA'] = vna

settings = get_default_settings()

#Save location
settings['scanname']    = 'testingcode'
settings['meas_type']   = 'Spec'
settings['project_dir'] = r'Z:\Data\deleteme'

#Sweep parameters
settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10

settings['start_power']  = -50
settings['stop_power']   = -20
settings['power_points'] = 5

#VNA settings
settings['channel'] = 1
settings['avg_time'] = 10
settings['measurement'] = 'S21'
settings['start_freq'] = 3.5e9
settings['stop_freq'] = 4.5e9
settings['freq_points'] = 501
settings['RFpower'] = -25
settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1
settings['CAVpower'] = -48
settings['CAVfreq'] = 8.126e9
settings['ifBW'] = 2e2

vna_spec(instruments, settings)