from vna_autler_townes import get_default_settings, vna_autler_townes

instruments = {}
instruments['VNA'] = vna
instruments['RFsource'] = qubitgen
instruments['SRS'] = SRS

settings = get_default_settings()

fluxon = [5.32483e9]
plasmon = [6.9715e9]
fluxon_power = -25
plasmon_power = -40

#for ffreq in fluxon:
#    for pfreq in plasmon:
#        #Save location
#        settings['scanname']    = 'AT_plasmon_search'
#        settings['meas_type']   = 'Autler_Townes'
#        settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'
#        
#        #Sweep parameters
#        settings['CAV_Attenuation'] = 30
#        settings['Qbit_Attenuation'] = 20
#        settings['Autler_Attenuation'] = 10
#        
#        plasmon_freq = pfreq
#        fluxon_freq = ffreq
#        spec_freq = plasmon_freq
#        spec_power = plasmon_power
#        
#        autler_center = plasmon_freq - fluxon_freq
#        autler_span = 20e6
#        settings['ext_flux'] = 0.087
#        settings['autler_power'] = -10
#        settings['start_autler_freq']  = autler_center - autler_span/2
#        settings['stop_autler_freq']   = autler_center + autler_span/2
#        settings['autler_points'] = 61
#        
#        #VNA settings
#        
#        span = 20e6
#        settings['channel'] = 1
#        settings['avg_time'] = 60
#        settings['measurement'] = 'S21'
#        settings['start_freq'] = spec_freq - span/2
#        settings['stop_freq'] = spec_freq + span/2
#        settings['freq_points'] = 101
#        settings['RFpower'] = spec_power
#        settings['RFport'] = 3
#        settings['Mport'] = 2
#        settings['CAVport'] = 1
#        settings['CAVpower'] = -65
#        settings['CAVfreq'] = 6.56063e9
#        settings['ifBW'] = 1e3
#        
#        settings['reverse'] = False
#        settings['random'] = False
#        vna_autler_townes(instruments, settings)
        
for ffreq in fluxon:
    for pfreq in plasmon:
        #Save location
        settings['scanname']    = 'AT_fluxon_debugging'
        settings['meas_type']   = 'Autler_Townes'
        settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'
        
        #Sweep parameters
        settings['CAV_Attenuation'] = 30
        settings['Qbit_Attenuation'] = 20
        settings['Autler_Attenuation'] = 10
        
        plasmon_freq = pfreq
        fluxon_freq = ffreq
        spec_freq = fluxon_freq
        spec_power = fluxon_power
        
        autler_center = plasmon_freq - fluxon_freq
        autler_span = 60e6
        settings['ext_flux'] = 0.087
        settings['autler_power'] = 0
        settings['start_autler_freq']  = autler_center - autler_span/2
        settings['stop_autler_freq']   = autler_center + autler_span/2
        settings['autler_points'] = 61
        
        #VNA settings
        
        span = 20e6
        settings['channel'] = 1
        settings['avg_time'] = 20
        settings['measurement'] = 'S21'
        settings['start_freq'] = spec_freq - span/2
        settings['stop_freq'] = spec_freq + span/2
        settings['freq_points'] = 101
        settings['RFpower'] = spec_power
        settings['RFport'] = 3
        settings['Mport'] = 2
        settings['CAVport'] = 1
        settings['CAVpower'] = -65
        settings['CAVfreq'] = 6.56063e9
        settings['ifBW'] = 1e3
        
        settings['reverse'] = False
        settings['random'] = False
        vna_autler_townes(instruments, settings)

for ffreq in fluxon:
    for pfreq in plasmon:
        #Save location
        settings['scanname']    = 'AT_fluxon_AToff_debugging'
        settings['meas_type']   = 'Autler_Townes'
        settings['project_dir'] = r'Z:\Data\Fluxonium_Raman\CRF01_A3'
        
        #Sweep parameters
        settings['CAV_Attenuation'] = 30
        settings['Qbit_Attenuation'] = 20
        settings['Autler_Attenuation'] = 10
        
        plasmon_freq = pfreq
        fluxon_freq = ffreq
        spec_freq = fluxon_freq
        spec_power = fluxon_power
        
        autler_center = plasmon_freq - fluxon_freq
        autler_span = 60e6
        settings['ext_flux'] = 0.087
        settings['autler_power'] = -130
        settings['start_autler_freq']  = autler_center - autler_span/2
        settings['stop_autler_freq']   = autler_center + autler_span/2
        settings['autler_points'] = 61
        
        #VNA settings
        
        span = 20e6
        settings['channel'] = 1
        settings['avg_time'] = 20
        settings['measurement'] = 'S21'
        settings['start_freq'] = spec_freq - span/2
        settings['stop_freq'] = spec_freq + span/2
        settings['freq_points'] = 101
        settings['RFpower'] = spec_power
        settings['RFport'] = 3
        settings['Mport'] = 2
        settings['CAVport'] = 1
        settings['CAVpower'] = -65
        settings['CAVfreq'] = 6.56063e9
        settings['ifBW'] = 1e3
        
        settings['reverse'] = False
        settings['random'] = False
        vna_autler_townes(instruments, settings)