# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 16:28:05 2025

@author: jhyang
"""
# exp_globals are acquired by running exp_globals_FPGA in shared namespace
from qick_measurements_CW.CW_trans import get_CW_trans_settings, CW_trans

settings = get_CW_trans_settings()
instruments = {}
settings['scanname']  = 'E_Delay_Calibration_Check2'
 
cav_center = 8.069e9#7.230796e9#4.389545e9 + 1e6 #7.08979e9 #7.1995e9 #7.19962e9 #q2 7.11736e9 #q3 7.1409e9
span = 5e6
freq_points = 801


settings['gain_start']     = 15000#8000 #7000
settings['gain_step']      = 12000
settings['gain_points']    = 1


#ADC settings
settings['reps']      = 1
settings['soft_avgs']  = 1
settings['meas_window'] = 10000
settings['initial_phase']   = 0.019509930789210244 #-0.001041 -8.860519651869944e-07 + 1.2801928797568243e-06 # 0.13 # rad


settings['freq_start']      = cav_center-span/2
settings['freq_step']       = span/(freq_points-1)
settings['freq_points']     = freq_points 


fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

full_data = CW_trans(soc,soccfg,instruments,fullsettings)