from vna_spec import get_default_settings, vna_spec

instruments = {}
instruments['VNA'] = vna

settings = get_default_settings()

SRS.output = 'On'
SRS.voltage_ramp(0.112)

#Save location
settings['scanname']    = 'fluxon_powersweep_fine'
settings['meas_type']   = 'Spec'
settings['project_dir'] = r'Z:\Data\HouckDualHangerFluxonium'

#Sweep parameters
settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10

settings['start_power']  = -20
settings['stop_power']   = 0
settings['power_points'] = 21

#VNA settings
settings['channel'] = 1
settings['avg_time'] = 20
settings['measurement'] = 'S21'
settings['start_freq'] = 4.585e9
settings['stop_freq'] = 4.635e9
settings['freq_points'] = 201
settings['RFpower'] = -25
settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1
settings['CAVpower'] = -55
settings['CAVfreq'] = 7.5765e9
settings['ifBW'] = 1e3

vna_spec(instruments, settings)

vna.output = 'Off'
vna.power = -10
SRS.voltage_ramp(0.)
SRS.output = 'Off'
