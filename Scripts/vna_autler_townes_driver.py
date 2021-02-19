from vna_autler_townes import get_default_settings, vna_autler_townes

instruments = {}
instruments['VNA'] = vna
instruments['RFsource'] = qubitgen
instruments['SRS'] = SRS

settings = get_default_settings()

#Save location
settings['scanname']    = 'fluxon_meas_AT_low_cav'
settings['meas_type']   = 'Autler_Townes'
settings['project_dir'] = r'Z:\Data\HouckDualHangerFluxonium'

#Sweep parameters
settings['CAV_Attenuation'] = 30
settings['Qbit_Attenuation'] = 10
settings['Autler_Attenuation'] = 10

settings['ext_flux'] = 0.112
settings['autler_power'] = 10
settings['start_autler_freq']  = 4.19e9
settings['stop_autler_freq']   = 4.23e9
settings['autler_points'] = 21

#VNA settings
settings['channel'] = 1
settings['avg_time'] = 20
settings['measurement'] = 'S21'
settings['start_freq'] = 4.57e9
settings['stop_freq'] = 4.65e9
settings['freq_points'] = 101
settings['RFpower'] = -5
settings['RFport'] = 3
settings['Mport'] = 2
settings['CAVport'] = 1
settings['CAVpower'] = -55
settings['CAVfreq'] = 7.5765e9
settings['ifBW'] = 1e3

vna_autler_townes(instruments, settings)

vna.output = 'Off'
vna.power = -10
#SRS.voltage_ramp(0.)
#SRS.output = 'Off'
