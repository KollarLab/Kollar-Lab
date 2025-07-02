# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: kollarlab

Modified by: jhyang
"""

from qick.averager_program import QickSweep, AveragerProgram

import numpy as np
import os
import time
import matplotlib.pyplot as plt
import userfuncs
from utility.measurement_helpers import estimate_time
#from utility.plotting_tools import simplescan_plot
      
#Heavily considering getting rid of the initial and post buffers for the speedup classes...
#Don't see the use when we can't acquire_decimated` anyway.

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
        # does cav_phase matter here since it does not change?
        phase_c = self.deg2reg(cfg["cav_phase"], gen_ch=gen_ch)
        gain_c  = cfg["cav_gain"]
        self.default_pulse_registers(ch=gen_ch, freq=freq_c, phase=phase_c, gain=gain_c)
        self.set_pulse_registers(ch=gen_ch, style="const", 
                                 length=self.us2cycles(self.cfg["cav_pulse_len"],
                                 gen_ch=gen_ch), mode='periodic')
        
        # Configure qubit DAC
        freq_q  = self.freq2reg(cfg["qub_freq"],gen_ch=qub_ch)
        phase_q = self.deg2reg(cfg["qub_phase"], gen_ch=qub_ch)
        # possibly add a qubit gain
        gain_q  = cfg["qub_gain"]
        
        self.default_pulse_registers(ch=qub_ch, freq=freq_q, phase=phase_q, gain=gain_q)
        self.set_pulse_registers(ch=qub_ch, style="const",
                                 length=self.us2cycles(cfg['qub_len'],
                                 gen_ch=qub_ch), mode= 'periodic')
        
        #self.default_pulse_registers(ch=qub_ch, phase=phase_q, freq=freq_q, gain=gain_q)
        #sigma = self.us2cycles(cfg["qub_sigma"],gen_ch=qub_ch) 
        #num_sigma = cfg["num_sigma"]
        #self.add_gauss(ch=qub_ch, name="ex", sigma=sigma,length=int(sigma*num_sigma))        
        #self.set_pulse_registers(ch=qub_ch, style="flat_top",
        #                         waveform="ex",length=self.us2cycles(cfg['qub_len'],
        #                         gen_ch=qub_ch))

        ###Start sweep definition
        
        #self.qub_r_freq = self.get_gen_reg(qub_ch,"freq")
        
        
        #freq_start_reg = self.freq2reg(cfg["freq_start"],gen_ch=qub_ch)
        #freq_stop_reg = self.freq2reg(cfg["freq_stop"],gen_ch=qub_ch)
        
        #self.add_sweep(QickSweep(self, self.qub_r_freq, cfg['freq_start'], cfg['freq_stop'], cfg["freq_points"]))
        
        # give processor some time to configure pulses, I believe it sets the offset time 200 cycles into the future?
        # Try varying this and seeing if it moves. Also try putting a synci after the trigger in the body
        self.synci(200)   
    
    def body(self):
        
        #self.reset_phase(gen_ch = self.cfg['cav_channel'], t=0)
        #self.reset_phase(gen_ch = self.cfg['qub_channel'], t=0)
        
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.
        #sigma = self.us2cycles(self.cfg["qub_sigma"])
        #num_sigma = self.cfg["num_sigma"]
        #pulse_len = self.us2cycles(self.cfg['qub_len'],gen_ch=self.cfg['qub_channel']) + int(num_sigma*sigma)
        #offset = self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"])
        #meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=self.cfg["cav_channel"])
        #ex_time = meas_time - self.us2cycles(self.cfg['qub_delay'], gen_ch=self.cfg["qub_channel"]) - pulse_len
        #Sets off the ADC
        #self.trigger(adcs=self.ro_chs,
        #            pins=[0],
        #            adc_trig_offset=offset)
        
        #Sends measurement pulse
        self.measure(pulse_ch=[self.cfg["cav_channel"], self.cfg['qub_channel']], 
             adcs=[self.cfg["ro_channels"][0]],
             adc_trig_offset=self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"]),
             wait=True,
             syncdelay=self.us2cycles(self.cfg["relax_delay"],gen_ch=self.cfg["cav_channel"]))
        #self.pulse(ch=self.cfg["cav_channel"],t=meas_time)
        #self.pulse(ch=self.cfg["qub_channel"],t=ex_time)
        #self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        #self.sync_all(self.us2cycles(self.cfg["relax_delay"])) #Syncs to an offset time after the final pulse is sent


def get_cw_spec_settings():
    settings = {}
    
    settings['scanname'] = 'continuous_power_scan'
    settings['meas_type'] = 'CW_Qubit_Testing'
    
    settings['cav_freq'] = 1e9
    settings['cav_gain'] = 1000
    settings['cav_pulse_len'] = 10
    settings['meas_window'] = 900

    settings['qub_gain'] = 1000
    settings['qub_len'] = 20
    
    #Sweep parameters
    settings['freq_start']   = 4e9  
    settings['freq_step']    = 0.5e9
    settings['freq_points']  = 6

    #Card settings
    settings['reps'] = 1
    settings['soft_avgs'] = 1#5e3
    
    return settings

def cw_spec(soc,soccfg,instruments,settings):

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    
    #soc_reset_gens() necessary?

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel'],
        'qub_channel'     : exp_globals['qub_channel'],
        'ro_channels'     : exp_globals['ro_channels'],

        'nqz_c'           : 1,
        'cav_phase'       : m_pulse['cav_phase'],
        'cav_pulse_len'   : exp_settings['cav_pulse_len'],
        #'meas_window'     : exp_settings['meas_window'],
        # set meas_pos to 0?
        'meas_time'       : m_pulse['meas_pos'],
        #'meas_gain'       : exp_settings['meas_gain'],
        'cav_gain'        : exp_settings['cav_gain'],
        'cav_freq'        : exp_settings['cav_freq']/1e6,
        
        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],

        'freq_start'      : exp_settings['freq_start']/1e6,
        'freq_step'       : exp_settings['freq_step']/1e6,
        'freq_points'     : exp_settings['freq_points'],
        'qub_gain'        : exp_settings['qub_gain'],

        #'qub_sigma'       : q_pulse['sigma'],
        #'qub_delay'       : q_pulse['delay'],
        #'num_sigma'       : q_pulse['num_sigma'],
        'qub_len'         : exp_settings['qub_len'],

        'readout_length'  : exp_settings['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs']
        }


    #prog = CW_spec(soccfg,config)
    #rep_period = config['adc_trig_offset'] + config['readout_length'] + config['relax_delay']
    #projected_time = config['reps']*config['soft_avgs']*config['freq_points']*rep_period/1e6
    #print("Projected Time: " + str(projected_time))
    
    # it was only two in CW_trans, why three here?
    #avg_di, avg_dq = prog.acquire(soc, load_pulses=True, progress=False)
    #Is = avg_di[0][0]
    #Qs = avg_dq[0][0]
    
    f0_start = exp_settings['freq_start']
    f0_step = exp_settings['freq_step']
    expts = exp_settings['freq_points']
    
    fpts = np.arange(0,expts)*f0_step+f0_start
    end_freq = exp_settings["freq_start"] + exp_settings["freq_step"]*exp_settings["freq_points"]
    
    powerdat = np.zeros(len(fpts))
    phasedat = np.zeros(len(fpts))
    Is = np.zeros(len(fpts))
    Qs = np.zeros(len(fpts))
    
    t_start = time.time()
    
    for f in range(0,len(fpts)):
        config["qub_freq"]=fpts[f]/1e6 # convert to MHz
        prog = CW_spec(soccfg, config)
        trans_I, trans_Q = prog.acquire(soc,progress=False)
        mag = np.sqrt(trans_I[0][0]**2 + trans_Q[0][0]**2)
        phase = np.arctan2(trans_Q[0][0], trans_I[0][0])*180/np.pi
        powerdat[f] = mag
        phasedat[f] = phase
        Is[f] = trans_I[0][0]
        Qs[f] = trans_Q[0][0]
        
        if f == 0:
            t_stop = time.time()
            estimate_time(t_start, t_stop, len(fpts))

    full_data = {}

    full_data['xaxis']  = fpts/1e9
    full_data['mags']   = powerdat
    full_data['phases'] = phasedat
    full_data['Is']     = Is
    full_data['Qs']     = Qs

    fig1 = plt.figure(1)
    plt.clf()
    plt.plot(fpts/1e9, powerdat)
    plt.title('Mag {}'.format(filename))
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Amplitude')
    fig1.canvas.draw()
    fig1.canvas.flush_events()
    plt.savefig(os.path.join(saveDir, filename+'_mag.png'), dpi=150)
    
    fig2 = plt.figure(2)
    plt.clf()
    plt.plot(fpts/1e9, phasedat)
    plt.title('Phase {}'.format(filename))
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Phase (Degrees)')
    fig1.canvas.draw()
    fig1.canvas.flush_events()
    plt.savefig(os.path.join(saveDir, filename+'_phase.png'), dpi=150)
    
    userfuncs.SaveFull(saveDir, filename, ['fpts','full_data','filename'],
    locals(), expsettings=settings, instruments={})

    #t_final    = time.time()
    #t_single = t_final - t_start
    
    #print("Elapsed Time: " + str(t_single))

    data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data}

    return data,prog
