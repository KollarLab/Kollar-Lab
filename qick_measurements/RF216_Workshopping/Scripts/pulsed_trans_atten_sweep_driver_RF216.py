from pulsed_trans_atten_sweep_RF216 import pulsed_trans, get_trans_settings


instruments = {}
#instruments['LO'] = logen


settings = get_trans_settings()
settings['scanname']  = 'cavity_atten_sweep'
settings['meas_type']   = 'pulsed_trans'

#offset=-0.0035
SRS2.voltage_ramp(1.449)#0.093+offset)
SRS3.voltage_ramp(0.21)

cav_center = 7.177228e9#6.10208e9 
span = 4e6 
freq_points = 11

settings['atten_start']     = 40
settings['atten_step']      = -5
settings['atten_points']    = 1

settings['meas_gain'] = 0.5

#ADC settings
settings['reps']      = 1
settings['rounds']  = 500

settings['freq_start']      = cav_center-span/2 
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 
settings['mixer_detuning']  = 200e6

settings['filter'] = True


#settings['nqz_c'] = 2

#settings['subtract_background'] = True 
#Currently no background subtraction

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = pulsed_trans(soc,soccfg,instruments,fullsettings)


