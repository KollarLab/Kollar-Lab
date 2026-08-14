import numpy as np
from pdh_measurements.pulsed_measurements.pulsed_DPA_4wm_cal import get_default_settings, pulsed_DPA_4wm_cal

instruments = {}
instruments['cavitygen'] = cavitygen
instruments['pumpgen'] = qubitgen
instruments['LO'] = holz.ch1
instruments['card'] = card
instruments['AWG'] = hdawg
instruments['pumpbias'] = SRS1

settings = get_default_settings()

settings['scanname']  = 'highfreqtarget_powerscan'
settings['meas_type'] = 'pulsed_DPA_cal'

# vbias parameters
settings['start_voltage']     = 0
settings['stop_voltage']      = 0
settings['voltage_points']    = 1

# pump parameters
settings['total_pump_power']  = -23
settings['start_pump_power']  = -43
settings['stop_pump_power']   = -33
settings['pump_power_points'] = 5

settings['detuning']          = 40e6

settings['pump_window']       = 100e-6

settings['start_phase']       = [0,0,0]
settings['stop_phase']        = [0,1,0]
settings['phase_points']      = 41

# cavity parameters
settings['cav_power']         = -45
settings['start_cav_freq']    = 5.85e9
settings['stop_cav_freq']     = 5.86e9
settings['cav_freq_points']   = 11


# card settings
settings['segments']          = 1
settings['reads']             = 1
settings['averages']          = 5e2

fullsettings = {}
fullsettings['exp_globals'] = exp_globals
fullsettings['exp_settings'] = settings

data, data_off = pulsed_DPA_4wm_cal(instruments, fullsettings)