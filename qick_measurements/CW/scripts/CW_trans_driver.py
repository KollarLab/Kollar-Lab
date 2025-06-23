# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 16:28:05 2025

@author: jhyang
"""
# exp_globals are acquired by running exp_globals_FPGA in shared namespace
from CW_trans import get_CW_trans_settings, CW_trans

settings = get_CW_trans_settings()
instruments = {}
settings['scanname']  = 'Applied_Phase_Correction_zoom'
 
cav_center = 4.389545e9 + 1e6 #7.08979e9 #7.1995e9 #7.19962e9 #q2 7.11736e9 #q3 7.1409e9
span = 2e6 
freq_points = 201


settings['gain_start']     = 5000#8000 #7000
settings['gain_step']      = 2000
settings['gain_points']    = 1


#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 1
settings['meas_window'] = 1000
settings['initial_phase']   = -0.001041 -6.860519651869944e-07 # 0.13 # rad


settings['freq_start']      = cav_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = CW_trans(soc,soccfg,instruments,fullsettings)