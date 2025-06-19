# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 16:28:05 2025

@author: jhyang
"""
# exp_globals are acquired by running exp_globals_FPGA in shared namespace
from CW_trans import get_CW_trans_settings, CW_trans

settings = get_CW_trans_settings()
instruments = {}
settings['scanname']  = 'Test'
 
cav_center = 2e9 #7.08979e9 #7.1995e9 #7.19962e9 #q2 7.11736e9 #q3 7.1409e9
span = 2e9 
freq_points = 201

# what are these used for? necessary in CW?
settings['gain_start']     = 4000#8000 #7000
settings['gain_step']      = 2000
settings['gain_points']    = 1


#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 80
# change to 0.2
settings['dphi_df']   = 5.96965198e7


settings['freq_start']      = cav_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = CW_trans(soc,soccfg,instruments,fullsettings)