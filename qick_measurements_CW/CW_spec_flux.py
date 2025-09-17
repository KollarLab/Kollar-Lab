# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: kollarlab

Modified by: jhyang
"""

from qick.averager_program import AveragerProgram


import numpy as np
import matplotlib.pyplot as plt
import time
import userfuncs
import os
from utility.measurement_helpers import estimate_time
import logging
import utility.plotting_tools as plots
from utility.userfits_v2 import fit_model

class CW_trans(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=1) # nqz zone fixed

        #configure the readout lengths and downconversion frequencies
        readout = self.us2cycles(cfg["meas_window"],ro_ch=cfg["ro_channels"][0])
        self.declare_readout(ch=cfg["ro_channels"][0], length=readout,
                             freq=self.cfg["pulse_freq"], gen_ch=cfg["cav_channel"])

        freq=self.freq2reg(cfg["pulse_freq"], gen_ch=cfg["cav_channel"], 
                                 ro_ch=cfg["ro_channels"][0])
        self.set_pulse_registers(ch=cfg["cav_channel"], style="const", freq=freq,
                                 # converts phase degrees to QICK register val
                                 phase=self.deg2reg(cfg["cav_phase"]), 
                                 gain=cfg["pulse_gain"], 
                                 length=self.us2cycles(cfg["pulse_length"],gen_ch=self.cfg["cav_channel"]), mode = "periodic")
        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        self.measure(pulse_ch=self.cfg["cav_channel"], 
             adcs=[self.cfg["ro_channels"][0]],
             adc_trig_offset=self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"]),
             wait=True,
             syncdelay=self.us2cycles(self.cfg["relax_delay"],gen_ch=self.cfg["cav_channel"]))


class CW_spec(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   
        gen_ch = cfg["cav_channel"]
        qub_ch = cfg["qub_channel"]

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=cfg["nqz_c"])
        self.declare_gen(ch=cfg["qub_channel"], nqz=cfg["nqz_q"])
        
        # configure the readout lengths and downconversion frequencies (ensuring it is an available DAC frequency)
        readout = self.us2cycles(cfg["readout_length"],ro_ch=cfg["ro_channels"][0])
        for ch in cfg["ro_channels"]:
            self.declare_readout(ch=ch, length=readout,
                                 freq=self.cfg["cav_freq"], gen_ch=cfg["cav_channel"])

        # Configure cavity DAC
        freq_c  = self.freq2reg(cfg["cav_freq"],gen_ch=gen_ch, ro_ch=cfg["ro_channels"][0])
        phase_c = self.deg2reg(cfg["cav_phase"], gen_ch=gen_ch)
        gain_c  = cfg["cav_gain"]
        self.default_pulse_registers(ch=gen_ch, freq=freq_c, phase=phase_c, gain=gain_c)
        self.set_pulse_registers(ch=gen_ch, style="const", 
                                 length=self.us2cycles(self.cfg["cav_pulse_len"],
                                 gen_ch=gen_ch), mode='periodic')
        
        # Configure qubit DAC
        freq_q  = self.freq2reg(cfg["qub_freq"],gen_ch=qub_ch)
        phase_q = self.deg2reg(cfg["qub_phase"], gen_ch=qub_ch)
        gain_q  = cfg["qub_gain"]
        
        self.default_pulse_registers(ch=qub_ch, freq=freq_q, phase=phase_q, gain=gain_q)
        self.set_pulse_registers(ch=qub_ch, style="const",
                                 length=self.us2cycles(cfg['qub_len'],
                                 gen_ch=qub_ch), mode= 'periodic')
        
        self.synci(200)   
    
    def body(self):
        
        #############################
        # LEGACY CODE FROM QUASI_CW #
        #############################
        
        #self.reset_phase(gen_ch = self.cfg['cav_channel'], t=0)
        #self.reset_phase(gen_ch = self.cfg['qub_channel'], t=0)
        
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.
        #sigma = self.us2cycles(self.cfg["qub_sigma"])
        #num_sigma = self.cfg["num_sigma"]
        #pulse_len = self.us2cycles(self.cfg['qub_len'],gen_ch=self.cfg['qub_channel']) + int(num_sigma*sigma)
        
        offset = self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"]) + self.us2cycles(20,gen_ch = self.cfg["cav_channel"])
        
        #meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=self.cfg["cav_channel"])
        #ex_time = meas_time - self.us2cycles(self.cfg['qub_delay'], gen_ch=self.cfg["qub_channel"]) - pulse_len
        #Sets off the ADC
        
        self.trigger(adcs=self.ro_chs,
                    pins=[0],
                    adc_trig_offset=offset)
        
        # Both measure() and pulse/pulse/wait/sync work.
        # Send pulses and trigger measurement
        #self.measure(pulse_ch=[self.cfg["cav_channel"], self.cfg['qub_channel']], 
             # adcs=[self.cfg["ro_channels"][0]],
             # adc_trig_offset=self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"]),
             # wait=True,
             # syncdelay=self.us2cycles(self.cfg["relax_delay"],gen_ch=self.cfg["cav_channel"]))
        self.pulse(ch=self.cfg["cav_channel"],t=0)
        self.pulse(ch=self.cfg["qub_channel"],t=0)
        self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        self.sync_all(self.us2cycles(self.cfg["relax_delay"])) #Syncs to an offset time after the final pulse is sent


def get_cw_spec_flux_settings():
    fullsettings = {}
    settings = {}
    autoscan_settings = {}

    settings['scanname'] = 'continuous_power_scan'
    settings['meas_type'] = 'CW_Spec_Flux'
    
    settings['cav_gain'] = 1000
    settings['meas_window'] = 900

    settings['qub_gain'] = 0
    
    #Freq sweep parameters
    settings['freq_start']   = 4e9  
    settings['freq_stop']    = 4.5e9
    settings['freq_points']  = 6

    #Voltage sweep parameters
    settings['start_voltage']   = 0
    settings['stop_voltage']    = 1
    settings['voltage_points']  = 2

    #Card settings
    settings['reps'] = 1
    settings['soft_avgs'] = 1#5e3

    autoscan_settings['reps'] = 1
    autoscan_settings['soft_avgs'] = 1 

    autoscan_settings['freq_start'] = 5e6
    autoscan_settings['freq_stop']  = 5.5e6
    autoscan_settings['freq_points'] = 6

    fullsettings['spec'] = settings
    fullsettings['autoscan'] = autoscan_settings

    return fullsettings

def cw_spec_flux(soc,soccfg,instruments,settings):
    
    # suppresses sum buffer overflow warning
    logging.getLogger("qick").setLevel(logging.ERROR)

    SRS = instruments['DCsupply']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    spec_set = exp_settings['spec']
    autoscan_set = exp_settings['autoscan']

    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = spec_set['scanname'] + '_' + stamp
    config_trans = {
        'cav_channel'     : exp_globals['cav_channel'],
        'ro_channels'     : exp_globals['ro_channels'],
        'relax_delay'     : exp_globals['relax_delay'],

        'nqz_c'           : 1,
        'cav_phase'       : m_pulse['cav_phase'],
        'pulse_length'    : 1,
        'meas_window'     : spec_set['meas_window'],
        'pulse_gain'      : spec_set['cav_gain'],
        'pulse_freq'      : 0, #Placeholder
        'adc_trig_offset' : m_pulse['emp_delay'],
        'reps'            : autoscan_set['reps'],
        'soft_avgs'       : autoscan_set['soft_avgs']
    }

    config_qubit = {
        'cav_channel'     : exp_globals['cav_channel'],
        'qub_channel'     : exp_globals['qub_channel'],
        'ro_channels'     : exp_globals['ro_channels'],

        'nqz_c'           : 1,
        'cav_phase'       : m_pulse['cav_phase'],
        'cav_pulse_len'   : 1,

        'meas_time'       : m_pulse['meas_pos'],
        'cav_gain'        : spec_set['cav_gain'],
        'cav_freq'        : 0, #exp_settings['cav_freq']/1e6,
        
        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],
        'qub_gain'        : spec_set['qub_gain'], 

        #'qub_sigma'       : q_pulse['sigma'],
        #'qub_delay'       : q_pulse['delay'],
        #'num_sigma'       : q_pulse['num_sigma'],
        'qub_len'         : 1,

        'readout_length'  : spec_set['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'],  #+ m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : spec_set['reps'],
        'soft_avgs'       : spec_set['soft_avgs']
        }

    
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = spec_set['scanname'] + '_' + stamp


    #set voltage sweep
    start_voltage = spec_set['start_voltage']
    stop_voltage  = spec_set['stop_voltage']
    voltage_points = spec_set['voltage_points']
    voltages = np.round(np.linspace(start_voltage, stop_voltage, voltage_points),6)
    max_voltage = 10#3.5
    if np.max(voltages) > max_voltage:
        raise ValueError('max voltage too large!')
    else:
        settings['voltages'] = voltages

    SRS.Output = 'On'
    
    f_start_trans = autoscan_set['freq_start']
    f_stop_trans  = autoscan_set['freq_stop']
    expts_trans   = autoscan_set['freq_points']

    fpts_trans = np.linspace(f_start_trans,f_stop_trans,expts_trans)

    #fpts_trans = np.arange(0,expts_trans)*f_step_trans+f_start_trans

    trans_mags   = np.zeros((voltage_points, autoscan_set['freq_points']))
    trans_phases = np.zeros((voltage_points, autoscan_set['freq_points']))

    raw_trans = {}
    raw_trans['Is'] = np.zeros((voltage_points, autoscan_set['freq_points']))
    raw_trans['Qs'] = np.zeros((voltage_points, autoscan_set['freq_points']))   

###################

    f0_start = spec_set['freq_start']
    f0_stop = spec_set['freq_stop']
    expts = spec_set['freq_points']
    
    fpts = np.linspace(f0_start,f0_stop,expts) #np.arange(0,expts)*f0_step+f0_start


    mags = np.zeros((voltage_points, spec_set['freq_points']))
    phases = np.zeros((voltage_points, spec_set['freq_points']))

    raw_spec = {}
    raw_spec['Is'] = np.zeros((voltage_points, spec_set['freq_points']))
    raw_spec['Qs'] = np.zeros((voltage_points, spec_set['freq_points']))
    
    t_start = time.time()
    
    identifier = 'Cav Gain : ' + str(spec_set['cav_gain'])  + ' au'

    # qubit gain sweep
    for vind in range(len(voltages)):
        voltage = voltages[vind]
        print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
        
        SRS.voltage_ramp(voltage)
        time.sleep(0.1)

        print('trans')

        for tfind in range(len(fpts_trans)):
            config_trans['pulse_freq'] = fpts_trans[tfind]/1e6
            prog = CW_trans(soccfg,config_trans)

            trans_I, trans_Q = prog.acquire(soc,progress=False)
            soc.reset_gens()

            mag = np.sqrt(trans_I[0][0]**2 + trans_Q[0][0]**2)
            phase = np.arctan2(trans_Q[0][0], trans_I[0][0])*180/np.pi

            trans_mags[vind,tfind] = mag
            trans_phases[vind,tfind] = phase
            raw_trans['Is'][vind,tfind] = trans_I[0][0]
            raw_trans['Qs'][vind,tfind] = trans_Q[0][0]

        hanger = exp_globals['hanger'] #"Fitting" cav freq

        try:
            print("Fitting Lorenzian to Cavity")
            config_qubit['cav_freq'] = fit_model(fpts_trans, trans_mags[vind], 'lorenz')['center']/1e6
        except:
            print("Fitting Lorenzian Failed, taking extrema...")
            if hanger:
                config_qubit['cav_freq'] = fpts_trans[np.argmin(trans_mags[vind])]/1e6
            else:
                config_qubit['cav_freq'] = fpts_trans[np.argmax(trans_mags[vind])]/1e6


        print('spec, cav freq: {}'.format(config_qubit['cav_freq']))

        # qubit frequency sweep
        for find in range(0,len(fpts)):
            
            config_qubit["qub_freq"]=fpts[find]/1e6 # convert to MHz
            #print(config_qubit)
            prog = CW_spec(soccfg, config_qubit)
            I, Q = prog.acquire(soc,progress=False)
            soc.reset_gens()
            
            mag = np.sqrt(I[0][0]**2 + Q[0][0]**2)
            phase = np.arctan2(Q[0][0], I[0][0])*180/np.pi
            
            mags[vind,find] = mag
            phases[vind,find] = phase
            raw_spec['Is'][vind,find] = I[0][0]
            raw_spec['Qs'][vind,find] = Q[0][0]
            
        if vind == 0:
            t_stop = time.time()
            estimate_time(t_start, t_stop, len(voltages))
    
        transdata = {}
        transdata['xaxis'] = fpts_trans
        transdata['mags'] = trans_mags[0:vind+1,:]
        transdata['phases'] = trans_phases[0:vind+1,:]

        specdata = {}
        specdata['xaxis'] = fpts/1e9 # Convert to GHz
        specdata['mags'] = mags[0:vind+1,:]
        specdata['phases'] = phases[0:vind+1,:]

        singledata = {}
        singledata['xaxis'] = fpts/1e9
        singledata['mag']   = specdata['mags'][vind]
        singledata['phase'] = specdata['phases'][vind]

        trans_labels = ['Freq (GHz)','Voltage (V)']
        spec_labels  = ['Freq (GHz)','Voltage (V)']
        
        specplotdata = {}
        specplotdata['xaxis']  = specdata['xaxis']
        specplotdata['mags']   = specdata['mags']
        specplotdata['phases'] = specdata['phases']
        
        mat = np.copy(specplotdata['mags'])
        for ind in range(0, mat.shape[0]):
            mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
        specplotdata['mags'] = mat
        
        mat = np.copy(specplotdata['phases'])
        for ind in range(0, mat.shape[0]):
            mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
        specplotdata['phases'] = mat
        
        plots.autoscan_plot(transdata, specplotdata, singledata, voltages[0:vind+1], filename, trans_labels, spec_labels, identifier, fig_num = 1)
        userfuncs.SaveFull(saveDir, filename, ['transdata', 'raw_trans', 'specdata', 'raw_spec', 'singledata', 'voltages', 
                                       'filename', 'trans_labels', 'spec_labels'], 
                                       locals(), expsettings=settings, instruments=instruments)
        plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)    
    

    data = {'saveDir': saveDir, 'filename': filename, 'specdata': specdata}

    return data,prog
