# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 15:33:02 2023

@author: Kollarlab

File to hold a lot of the qubit calibration stuff that we like to do:
    -cavity frequency scan
    -qubit frequency scan 
    -mixer calibration
    -hold time sweep to get the pi pulse amplitude 
For now, I'm just going to copy paste a lot of the drivers and adapt them to
the "new" format but longer term I think we'll want to change how we do the 
lower level stuff
"""
from utility.userfits_v2 import fit_model
import numpy as np
import matplotlib.pyplot as plt

def get_cavity_freq(instruments, calibration, exp_globals, freq_range=10e6, freq_points=51, ax=None):
    from pulsed_measurements.pulsed_trans import get_default_settings, pulsed_trans
    settings = get_default_settings()
    settings['scanname'] = 'cavity_scan_calibration'
    
    cav_center = calibration['CAV_Freq']
    span = freq_range
    settings['start_freq']      = cav_center-span/2
    settings['stop_freq']       = cav_center+span/2
    settings['freq_points']     = freq_points
    
    settings['start_power']     = -50
    settings['power_points']    = 1
    
    #Card settings
    settings['averages']         = 2e3
    
    settings['subtract_background'] = True
    
    fullsettings = {}
    fullsettings['exp_globals'] = exp_globals
    fullsettings['exp_settings'] = settings

    data = pulsed_trans(instruments, fullsettings)
    params = fit_model(data['xaxis']*1e9, data['mags'][0], 'lorenz', plot=True, ax=ax)
    return np.round(params['center'])

def calibrate_mixer(instruments, calibration, exp_globals, ax=None):
    from pulsed_measurements.mixer_cal import mixer_cal, get_default_settings
    settings = get_default_settings()
    settings['CAV_Freq'] = calibration['CAV_Freq']
    settings['CAV_Power'] = calibration['CAV_Power']
    settings['plot'] = True
    fullsettings = {}
    fullsettings['exp_globals'] = exp_globals
    fullsettings['exp_settings'] = settings
    
    mix_config = mixer_cal(instruments, fullsettings, ax)
    mix_config['center'] = [0,0]
    return mix_config

def get_qubit_freq(instruments, calibration, exp_globals, freq_range=50e6, freq_points=101, power=-20, ax=None):
    from pulsed_measurements.pulsed_spec import get_default_settings, pulsed_spec
    exp_settings = get_default_settings()
    exp_settings['scanname'] = 'qubit_scan_calibration'
    
    hold_time = exp_globals['qubit_pulse']['hold_time']
    exp_globals['qubit_pulse']['hold_time'] = 30e-6
    
    #Cavity parameters
    exp_settings['CAV_Power']       = calibration['CAV_Power']
    exp_settings['CAV_freq']        = calibration['CAV_Freq']
    
    #Qubit parameters
    q_freq = calibration['Q_Freq']
    span = freq_range
    exp_settings['start_freq']      = q_freq - span/2
    exp_settings['stop_freq']       = q_freq + span/2
    exp_settings['freq_points']     = freq_points
    
    exp_settings['start_power']     = power
    exp_settings['power_points']    = 1
    exp_settings['averages']         = 2e3
    
    exp_settings['subtract_background'] = False
    
    settings = {}
    settings['exp_globals'] = exp_globals
    settings['exp_settings'] = exp_settings
    
    spec_data = pulsed_spec(instruments, settings)
    mag = spec_data['mags'][0]
    try:
        params = fit_model(spec_data['xaxis']*1e9, abs(mag-max(mag)), 'lorenz', plot=True, ax=ax)
    except:
        params = {'center':0}
    exp_globals['qubit_pulse']['hold_time'] = hold_time
    return np.round(params['center']), abs(mag-max(mag)), spec_data['xaxis']
    
def get_pi_pulse(instruments, calibration, exp_globals, ax=None):
    from pulsed_measurements.hold_time_sweep import get_default_settings, hold_time_sweep
    exp_settings = get_default_settings()
    exp_settings['scanname'] = 'pi_pulse_calibration'
    
    exp_settings['Q_Freq']    = calibration['Q_Freq']
    exp_settings['Q_Power']   = calibration['Q_Power']
    exp_settings['CAV_Power'] = calibration['CAV_Power']
    exp_settings['CAV_Freq']  = calibration['CAV_Freq']
    
    #Card settings
    exp_settings['segments']      = 1
    exp_settings['reads']         = 1
    exp_settings['averages']      = 1e3
    
    exp_settings['Tau_min']    = 0e-9
    exp_settings['Tau_max']    = 150e-9
    exp_settings['Tau_points'] = 101
    exp_settings['spacing'] = 'lin'
    exp_settings['subtract_background'] = True
    exp_settings['DRAG'] = False
    
    fullsettings = {}
    fullsettings['exp_globals'] = exp_globals
    fullsettings['exp_settings'] = exp_settings
    t_ax, amps = hold_time_sweep(instruments, fullsettings)
    params = fit_model(t_ax, amps, 'cos',True, ax=ax)
    period = 1/params['freq']
    t_ax_fine = np.linspace(0, exp_settings['Tau_max'],201)
    p_len_points = int(period/t_ax_fine[1])
    osc = params['amp']*np.cos(params['freq']*t_ax_fine[:p_len_points]*2*np.pi+params['phi'])+params['offset']
    pi_time = int(t_ax_fine[np.argmax(osc)]*1e9)/1e9
    if not ax:
        plt.figure()
        ax = plt.subplot(111)
        ax.plot(t_ax*1e6, amps)
        ax.vlines(pi_time*1e6,min(amps), max(amps))
        ax.plot(t_ax_fine[:p_len_points]*1e6, osc)
    else:
        ax.vlines(pi_time*1e6,min(amps), max(amps))
        ax.set_title('Pi pulse cal')
    return pi_time
