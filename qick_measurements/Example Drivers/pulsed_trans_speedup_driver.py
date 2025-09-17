from qick_measurements.pulsed_trans_speedup import pulsed_trans, get_trans_settings


instruments = {}
#instruments['LO'] = logen


settings = get_trans_settings()
settings['scanname']  = 'Low_Pow_Search'
cav_center = 5e9 #7.2309e9 #7.08979e9 #7.1995e9 #7.19962e9 #q2 7.11736e9 #q3 7.1409e9
span = 10e6 
freq_points = 2

settings['gain_start']     = 10000#20000#8000 #7000
settings['gain_step']      = 2000
settings['gain_points']    = 5


#ADC settings
settings['reps']      = 101
settings['soft_avgs']  = 1
# settings['phases']    = -phases[0]
settings['dphi_df']   = 0 #5.96965198e7


settings['freq_start']      = cav_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = pulsed_trans(soc,soccfg,instruments,fullsettings)


