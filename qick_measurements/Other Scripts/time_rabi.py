# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: kollarlab
"""

from qick.averager_program import AveragerProgram

import time
import numpy as np
import matplotlib.pyplot as plt
import userfuncs
import os
# from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import estimate_time
from utility.plotting_tools import general_colormap_subplot
      
class RabiSequence(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   
        gen_ch = cfg["cav_channel"]
        qub_ch = cfg["qub_channel"]

        # set the nyquist zone
        # Implement an if statement here to catch? 
        self.declare_gen(ch=cfg["cav_channel"], nqz=cfg["nqz_c"])
        self.declare_gen(ch=cfg["qub_channel"], nqz=cfg["nqz_q"])
        
        # configure the readout lengths and downconversion frequencies (ensuring it is an available DAC frequency)
        
        readout = self.us2cycles(cfg["readout_length"],ro_ch=cfg["ro_channels"][0])
        for ch in cfg["ro_channels"]:
            self.declare_readout(ch=ch, length=readout,
                                 freq=self.cfg["cav_freq"], gen_ch=cfg["cav_channel"])

        # convert frequency to DAC freqency (ensuring it is an available ADC frequency)
        freq_c  = self.freq2reg(cfg["cav_freq"],gen_ch=gen_ch, ro_ch=cfg["ro_channels"][0])
        phase_c = self.deg2reg(cfg["cav_phase"], gen_ch=gen_ch)
        gain_c  = cfg["meas_gain"]
        
        self.default_pulse_registers(ch=gen_ch, freq=freq_c, phase=phase_c, gain=gain_c, mode = "oneshot")
        self.set_pulse_registers(ch=gen_ch, style="const", length=self.us2cycles(self.cfg["meas_window"],gen_ch=gen_ch))
        
        freq_q  = self.freq2reg(cfg["qub_freq"],gen_ch=qub_ch)
        phase_q = self.deg2reg(cfg["qub_phase"], gen_ch=gen_ch)
        gain_q  = cfg["qub_gain"]
        
        self.default_pulse_registers(ch=qub_ch, phase=phase_q, freq=freq_q, gain=gain_q)#, mode = "oneshot")
        
        sigma = self.us2cycles(cfg["qub_sigma"],gen_ch=qub_ch)
        num_sigma = cfg["num_sigma"]
        hold = self.us2cycles(self.cfg['length'],gen_ch = self.cfg["qub_channel"])
        
        self.add_gauss(ch=qub_ch, name="ex", sigma=sigma,length=int(sigma*num_sigma))        
        self.set_pulse_registers(ch=qub_ch, style="flat_top", waveform="ex",length=hold)

        
        # give processor some time to configure pulses, I believe it sets the offset time 200 cycles into the future?
        # Try varying this and seeing if it moves. Also try putting a synci after the trigger in the body
        self.synci(200)   
    
    def body(self):
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.
        sigma = self.us2cycles(self.cfg["qub_sigma"])
        num_sigma = self.cfg["num_sigma"]
        hold = self.us2cycles(self.cfg['length'],gen_ch = self.cfg["qub_channel"])
        
        offset = self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"])
        meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=self.cfg["cav_channel"])
        ex_time = meas_time - self.us2cycles(self.cfg['qub_delay'],gen_ch=self.cfg["qub_channel"]) - int(num_sigma*sigma) - hold
        #Sets off the ADC
        self.trigger(adcs=self.ro_chs,
                    pins=[0],
                    adc_trig_offset=offset)
        
        #Sends measurement pulse
        self.pulse(ch=self.cfg["qub_channel"],t=ex_time)
        self.pulse(ch=self.cfg["cav_channel"],t=meas_time)
        self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        self.sync_all(self.us2cycles(self.cfg["relax_delay"])) #Syncs to an offset time after the final pulse is sent


def get_time_Rabi_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'RabiRate'
    
    settings['cav_freq'] = 1e9
    settings['cav_gain'] = 1000
    
    #Sweep parameters
    # settings['freq_start']   = 4e9  
    # settings['freq_step']    = 100e6
    # settings['freq_points']  = 6
    
    settings['qub_freq']       = 5.332e9
    settings['qub_gain']       = 3000
    
    settings['hold_start']     = 0.1e-6
    settings['hold_step']      = 0.5e-6
    settings['hold_points']    = 7

    # settings['gain_start']  = 500
    # settings['gain_step']   = 100
    # settings['gain_points'] = 11
    #Card settings
    settings['reps'] = 1
    settings['soft_avgs'] = 5e3
    
    return settings

def time_Rabi_sweep(soc,soccfg,instruments,settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    
    soc.reset_gens()
    
    if exp_globals['LO']:
        logen = instruments['LO']
        
        logen.freq   = exp_globals['LO_freq']
        logen.power  = exp_globals['LO_power']
        logen.output = 1

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel'],
        'qub_channel'     : exp_globals['qub_channel'],
        'ro_channels'     : exp_globals['ro_channels'],

        'nqz_c'           : 1,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'meas_gain'       : exp_settings['cav_gain'],
        'cav_freq'        : (exp_settings['cav_freq']-exp_globals['LO_freq'])/1e6,
        
        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],
        'qub_freq'        : exp_settings['qub_freq']/1e6,
        'qub_gain'        : exp_settings['qub_gain'], 
        'qub_sigma'       : q_pulse['sigma'],
        'qub_delay'       : q_pulse['delay'],
        'num_sigma'       : q_pulse['num_sigma'],
        'length'          : 0.1, #Placeholder
        
        'readout_length'  : m_pulse['init_buffer'] + m_pulse['meas_window'] + m_pulse['post_buffer'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'] - m_pulse['init_buffer'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs']
        }

    hpts = (exp_settings["hold_start"]+exp_settings["hold_step"]*np.arange(exp_settings["hold_points"]))*1e6

    prog = RabiSequence(soccfg, config)
    meas_start = prog.us2cycles(m_pulse["init_buffer"],ro_ch=0)
    meas_end = meas_start+prog.us2cycles(m_pulse["meas_window"],ro_ch=0)
    total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)
   

    # powerdat = np.zeros((len(gpts), len(fpts)))
    # phasedat = np.zeros((len(gpts), len(fpts)))
    
    amp_int = np.zeros(len(hpts))
    ang_int = np.zeros(len(hpts))
    amps    = np.zeros((len(hpts),total_samples))
    angles  = np.zeros((len(hpts),total_samples))

    #drive_powers_lin = 10**(powers/10) ?
    #drive_amps_lin = np.sqrt(drive_powers_lin)
    
    tstart = time.time()
    first_it = True

    for h in range(0,len(hpts)):
        print("Current Hold Time: " + str(hpts[h]) + ", Max Gain: " + str(hpts[-1]))
        config["length"] = hpts[h]

        total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)
        
        prog = RabiSequence(soccfg,config)
        
        #Need to assign Iwindow, Qwindow, Ifull, Qfull, xaxis (which should just be timeus)
        holder = prog.acquire_decimated(soc, load_pulses=True, progress=False, debug=False)
        I_full = holder[0][0]
        Q_full = holder[0][1]
        I_window = I_full[meas_start:meas_end]
        Q_window = Q_full[meas_start:meas_end]
       
        I_final = np.mean(I_window)
        Q_final = np.mean(Q_window)

        amps[h] = np.sqrt(I_full**2+Q_full**2)
        angles[h] = np.arctan2(Q_full,I_full)*180/np.pi
        amp_int[h] = np.sqrt(I_final**2+Q_final**2)
        ang_int[h] = np.arctan2(Q_final, I_final)*180/np.pi
        
        if first_it:
            xaxis = np.linspace(0,len(I_full)-1,len(I_full))
            
            tstop = time.time()
            estimate_time(tstart, tstop, len(hpts))

            for x in range(0,len(xaxis)):
                xaxis[x] = prog.cycles2us(xaxis[x],ro_ch=0)
                
            first_it = False
            
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(hpts, amp_int)
        plt.xlabel('Hold Time (us)')
        plt.ylabel('Amplitude')  
        plt.subplot(122)
        plt.plot(hpts, ang_int)
        plt.xlabel('Hold Time (us)')
        plt.ylabel('Phase')  
        plt.title('Live Rabi data \n'+filename)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'_live.png'), dpi = 150)
    
        fig2 = plt.figure(2,figsize=(13,8))
        plt.clf()
    
        ax = plt.subplot(1,1,1)
        general_colormap_subplot(ax, xaxis, hpts, amps, ['Time (us)', 'Hold Time (us)'], 'Raw data\n'+filename)
    
        plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
        userfuncs.SaveFull(saveDir, filename, ['hpts','xaxis', 'amps', 'amp_int'], locals(), expsettings=settings, instruments=instruments)

    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))
        




       
    
    if exp_globals['LO']:
        logen.output = 0
        
    return hpts, amp_int