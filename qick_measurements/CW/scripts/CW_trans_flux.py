# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: jhyang

"""

from qick.averager_program import AveragerProgram


import numpy as np
import time
import userfuncs
from utility.measurement_helpers import estimate_time
from scipy.stats import linregress
import logging
from utility.plotting_tools import simplescan_plot

#######################################
# Taken from CW_trans LoopBackProgram #
#######################################

# Plays a constant tone, swept under cw_spec_flux()
class CavitySweepFluxTrans(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=1) # nqz zone fixed

        #configure the readout lengths and downconversion frequencies
        readout = self.us2cycles(cfg["meas_window"],ro_ch=cfg["ro_channels"][0])
        self.declare_readout(ch=cfg["ro_channels"][0], length=readout,
                             freq=self.cfg["cav_freq"], gen_ch=cfg["cav_channel"])

        freq=self.freq2reg(cfg["cav_freq"], gen_ch=cfg["cav_channel"], 
                                 ro_ch=cfg["ro_channels"][0])
        self.set_pulse_registers(ch=cfg["cav_channel"], style="const", freq=freq,
                                 # converts phase degrees to QICK register val
                                 phase=self.deg2reg(cfg["cav_phase"]), 
                                 gain=cfg["cav_gain"], 
                                 length=self.us2cycles(cfg["cav_pulse_len"],gen_ch=self.cfg["cav_channel"]), mode = "periodic")
        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        self.measure(pulse_ch=self.cfg["cav_channel"], 
             adcs=[self.cfg["ro_channels"][0]],
             adc_trig_offset=self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"]),
             wait=True,
             syncdelay=self.us2cycles(self.cfg["relax_delay"],gen_ch=self.cfg["cav_channel"]))
        

def get_cw_trans_flux_settings():
    settings = {}
    
    settings['scanname'] = 'continuous_power_scan'
    settings['meas_type'] = 'CW_Trans_Flux_Testing'
    
    settings['start_voltage']  = 0
    settings['stop_voltage']   = 0.1
    settings['voltage_points'] = 5
    
    
    settings['cav_gain'] = 1000
    settings['meas_window'] = 900
    settings['cav_pulse_len'] = 10
    settings['initial_phase'] = 0
    
    #Sweep parameters
    settings['freq_start']   = 4e9  
    settings['freq_step']    = 0.5e9
    settings['freq_points']  = 6

    #Card settings
    settings['reps'] = 1
    settings['soft_avgs'] = 1#5e3
    
    return settings

def cw_spec_flux(soc,soccfg,instruments,settings):
    
    print("Initializing...")
    
    # suppresses sum buffer overflow warning
    logging.getLogger("qick").setLevel(logging.ERROR)

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    
    SRS = instruments['DCsupply']
    soc.reset_gens()
    

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel'],
        'ro_channels'     : exp_globals['ro_channels'],

        'nqz_c'           : 1,
        'cav_phase'       : m_pulse['cav_phase'],
        'cav_pulse_len'   : exp_settings['cav_pulse_len'],
        'meas_time'       : m_pulse['meas_pos'],
        #'meas_gain'       : exp_settings['meas_gain'],
        'cav_gain'        : exp_settings['cav_gain'],

        'freq_start'      : exp_settings['freq_start']/1e6,
        'freq_step'       : exp_settings['freq_step']/1e6,
        'freq_points'     : exp_settings['freq_points'],
        'qub_gain'        : exp_settings['qub_gain'],

        'meas_window'  : exp_settings['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'],  #+ m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs']
        }

    #############################
    # Taken from spec_flux_scan #
    #############################
    
    #set voltage sweep
    start_voltage = exp_settings['start_voltage']
    stop_voltage  = exp_settings['stop_voltage']
    voltage_points = exp_settings['voltage_points']
    voltages = np.round(np.linspace(start_voltage, stop_voltage, voltage_points),6)
    max_voltage = 3.5
    if np.max(voltages) > max_voltage:
        raise ValueError('max voltage too large!')
    else:
        settings['voltages'] = voltages

    SRS.Output = 'On'
    
    #Making array of cavity frequencies for transmission scan (to be looped through later)
    start_freq = exp_settings['freq_start']
    stop_freq = exp_settings['freq_stop']
    freq_points = exp_settings['freq_points']
    trans_fpts = np.linspace(start_freq,stop_freq,freq_points)

    #Dummy arrays for cavity scan
    trans_mags   = np.zeros((voltage_points, exp_settings['freq_points']))
    trans_phases = np.zeros((voltage_points, exp_settings['freq_points']))

    #Defining file variables
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp
    
    phase_slope = (exp_settings['initial_phase']) * (180/np.pi) # degrees
    def freq2phase(fpts):
        return (exp_settings["freq_start"] - fpts) * phase_slope # note the negative sign!
    # phase correcting array
    phase_fpts=freq2phase(trans_fpts)
    
    first_it = True
    t_start = time.time()
    phase_result = {}
    for vind in range(len(voltages)):
        
        if exp_settings['stability'] == True: # If running a stability scan, will ramp to starting voltage on the first iteration but then never ramp again.
            if first_it == True:
                print('Stability scan enabled.')
                voltage = voltages[vind]
                print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
                
                SRS.voltage_ramp(voltage)
                time.sleep(0.1)

                voltages = np.linspace(0,len(voltages)-1,len(voltages))
            voltage = voltages[vind]
        else: # If not, will ramp to appropriate voltage every loop in standard fashion.
            voltage = voltages[vind]
            print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
            
            SRS.voltage_ramp(voltage)
            time.sleep(0.1)
        
        print('Sweeping transmission.')
        
        for find in range(len(trans_fpts)):
            tstart = time.time()
            config['cav_freq'] = trans_fpts[find]/1e6 #Update the frequency, board wants it in MHz so converting now
            config['cav_phase'] = phase_fpts[find]

            trans_prog = CavitySweepFluxTrans(soccfg,config) #Make transmission pulse sequence object

            trans_I, trans_Q = trans_prog.acquire(soc,load_pulses=True,progress=False) #Transmission data acquisition occurs here
            soc.reset_gens()
            
            mag = np.sqrt(trans_I[0][0]**2 + trans_Q[0][0]**2)
            phase = np.arctan2(trans_Q[0][0], trans_I[0][0])*180/np.pi

            trans_mags[vind][find] = mag
            trans_phases[vind][find] = phase
            
            if find == 0:
                t_stop = time.time()
                estimate_time(t_start, t_stop, len(trans_fpts))
                
        if vind == 0:
            phase_result.append([trans_I, trans_Q])
        
    print("Transmission sweep complete. ")
    
    phase_result = np.transpose(phase_result)
    print("Slope: ")
    print(linregress(trans_fpts, np.unwrap
                    (np.arctan2((phase_result[0][0][0]),(phase_result[0][0][1]))))[0])
    
    full_data = {}

    full_data['xaxis']  = trans_fpts/1e9
    full_data['mags']   = trans_mags
    full_data['phases'] = trans_phases
    full_data['Is']     = trans_I
    full_data['Qs']     = trans_Q
    
    plot_data = {}
    plot_data['xaxis']  = trans_fpts/1e9
    plot_data['mags']   = trans_mags[0:vind+1]
    plot_data['phases'] = trans_phases[0:vind+1]

    single_data = {}
    single_data['xaxis'] = trans_fpts/1e9
    single_data['mag']   = trans_mags[vind]
    single_data['phase'] = trans_phases[vind]
    
    yaxis  = voltages[0:vind+1]
    labels = ['Freq (GHz)', 'Voltages (V)']
    
    simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=1, IQdata = False)
    
    userfuncs.SaveFull(saveDir, filename, ['voltages','trans_fpts','full_data','filename'],
    locals(), expsettings=settings, instruments={})
    
    return full_data