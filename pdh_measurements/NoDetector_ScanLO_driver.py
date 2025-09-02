# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:05:51 2024

@author: Kollarlab
"""

import os
import time
from pdh_measurements.pdh_funcs import NoDetector_LOScan

# instruments used
settings = {}

exp_globals = {}
exp_globals['root_folder'] = r"K:\Data\Oxford"
exp_globals['project_name'] = 'PDH'
exp_globals['device_name'] = 'PDH_v5'
settings['exp_globals'] = exp_globals

exp_settings = {}
exp_settings['meas_type'] = "PDH" 
exp_settings['extra_comment'] = "NoDetectporPDH_NewCode"
settings['exp_settings'] = exp_settings

instruments = {}
instruments["Rf_gen"] = gen
instruments["AWG"] = hdawg
instruments["digitizer"] = card
instruments['lo'] = holz

###################### ONLY THINGS TO CHANGE ###########################
averages = 1000
card.averages = averages
card.SetParams()

norm_amp = [1,2,1]
set1_1 = [0.38,0.24,0.38]
set1_2 = [0.12,0.76,0.12]
set2_1 = [2,1,2]
set2_2 = [0.5, 4, 0.5]
amp_list = norm_amp #set4:[0.5, 4, 0.5]#set3:[2,1,2]#[1, 2, 1] #[1, 2, 1] set1=[0.38,0.24,0.38], set2=[0.12,0.76,0.12]
phase_list = [1, 0, 0]


settings['mod_freq'] = 10e6 
settings['ddc_freq'] = 10e6 
settings['lo_power'] = 12
settings['IF_freq'] = 30e6
Normalize = True 
settings['gen_power'] = -30
settings['cavity_freq'] = 6.09e9 #7.28e9 #6.808e9
settings['num_points']  = 501
span = 40e6
settings['start_freq']  =   settings['cavity_freq'] - span/2
settings['stop_freq'] = settings['cavity_freq'] + span/2
settings['plot_data'] = True
settings['narrow_span'] = False # if true, 60MHz span
settings['amp_list'] = amp_list 
settings['phase_list'] = phase_list
settings['normalize'] = Normalize

settings['buffer_start'] = 3000
settings['buffer_end'] = 6000

    
    

########################################################################

start_time = time.time()
# measurement
data_dict = NoDetector_LOScan(instruments, settings)

end_time = time.time() - start_time
print(f'time_elasped:{end_time}')


################################################################################

