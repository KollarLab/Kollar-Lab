from pulsed_trans import pulsed_trans, get_trans_settings


instruments = {}
instruments['LO'] = logen


settings = get_trans_settings()
settings['scanname']  = 'test_scan_take_3_fine'
 
cav_center = 7.31272e9 
span = 5e6 
freq_points = 51

settings['gain_start']     = 1000
settings['gain_step']      = 1000
settings['gain_points']    = 1


#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 500

settings['freq_start']      = cav_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 

#settings['subtract_background'] = True 
#Currently no background subtraction

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = pulsed_trans(soc,soccfg,instruments,fullsettings)


