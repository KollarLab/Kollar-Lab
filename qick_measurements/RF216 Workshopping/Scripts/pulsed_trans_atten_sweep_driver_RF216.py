from pulsed_trans_atten_sweep_RF216 import pulsed_trans, get_trans_settings


instruments = {}
#instruments['LO'] = logen


settings = get_trans_settings()
settings['scanname']  = 'Loopback_Check'
 
cav_center = 4e9 
span = 10e6 
freq_points = 51

settings['atten_start']     = 1
settings['atten_step']      = 0.1
settings['atten_points']    = 1

settings['meas_gain'] = 0.5

#ADC settings
settings['reps']      = 1
settings['rounds']  = 1

settings['freq_start']      = cav_center-span/2 
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 
settings['mixer_detuning']  = 200e6

settings['filter'] = False


#settings['nqz_c'] = 2

#settings['subtract_background'] = True 
#Currently no background subtraction

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = pulsed_trans(soc,soccfg,instruments,fullsettings)



