from vna_trans import get_default_settings, vna_trans

instruments = {}
instruments['VNA'] = vna

settings = get_default_settings()
    
#Save location
settings['scanname']    = 'testingcode'
settings['meas_type']   = 'Trans'
settings['project_dir'] = r'Z:\Data\deleteme'
    
#Sweep parameter
settings['CAV_attenuation'] = 30

settings['start_power']  = -50
settings['stop_power']   = -20
settings['power_points'] = 7

#VNA settings
settings['channel']  = 1
settings['avg_time'] = 6
settings['measurement'] = 'S21'
settings['start_freq']  = 7.7*1e9  
settings['stop_freq']   = 7.8*1e9 
settings['freq_points'] = 501
settings['ifBW'] = 1e3

vna_trans(instruments, settings)