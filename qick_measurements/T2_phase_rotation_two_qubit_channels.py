# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:16:54 2023

@author: kollarlab
"""

from qick.averager_program import AveragerProgram

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
import time

#from utility.userfits import fit_T2
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import estimate_time

      
class T2_sequence(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   
        cav_ch = cfg["cav_channel"]
        qub_ch_D1 = cfg["qub_channel_D1"]
        qub_ch_D2 = cfg["qub_channel_D2"]

        # set the nyquist zone
        # Implement an if statement here to catch? 
        self.declare_gen(ch=cav_ch, nqz=cfg["nqz_c"])
        self.declare_gen(ch=qub_ch_D1, nqz=cfg["nqz_q"])
        self.declare_gen(ch=qub_ch_D2, nqz=cfg["nqz_q"])
        
        # configure the readout lengths and downconversion frequencies (ensuring it is an available DAC frequency)
        
        readout = self.us2cycles(cfg["readout_length"],ro_ch=cfg["ro_channels"][0])
        for ch in cfg["ro_channels"]:
            self.declare_readout(ch=ch, length=readout,
                                 freq=self.cfg["cav_freq"], gen_ch=cfg["cav_channel"])

        # convert frequency to DAC freqency (ensuring it is an available ADC frequency)
        freq_c  = self.freq2reg(cfg["cav_freq"],gen_ch=cav_ch, ro_ch=cfg["ro_channels"][0])
        phase_c = self.deg2reg(cfg["cav_phase"], gen_ch=cav_ch)
        gain_c  = cfg["meas_gain"]
        
        self.default_pulse_registers(ch=cav_ch, freq=freq_c, phase=phase_c, gain=gain_c)
        self.set_pulse_registers(ch=cav_ch, style="const", length=self.us2cycles(self.cfg["meas_window"],gen_ch=cav_ch))
        
        freq_q_D1  = self.freq2reg(cfg["qub_freq_D1"],gen_ch=qub_ch_D1)
        #phase_q_D1 = self.deg2reg(cfg["qub_phase_D1"], gen_ch=qub_ch_D1)
        gain_q_D1  = cfg["qub_gain_D1"]
        
        freq_q_D2  = self.freq2reg(cfg["qub_freq_D2"],gen_ch=qub_ch_D2)
        #phase_q_D2 = self.deg2reg(cfg["qub_phase_D2"], gen_ch=qub_ch_D2)
        gain_q_D2  = cfg["qub_gain_D2"]

        self.default_pulse_registers(ch=qub_ch_D1, freq=freq_q_D1, gain=gain_q_D1) #removed the phase settings here to avoid the conflict with later phase settings
        self.default_pulse_registers(ch=qub_ch_D2, freq=freq_q_D2, gain=gain_q_D2)

        sigma_D1 = self.us2cycles(cfg["qub_sigma_D1"],gen_ch=qub_ch_D1)
        num_sigma_D1 = cfg["num_sigma_D1"]

        sigma_D2 = self.us2cycles(cfg["qub_sigma_D2"],gen_ch=qub_ch_D2)
        num_sigma_D2 = cfg["num_sigma_D2"]
        
        self.add_gauss(ch=qub_ch_D1, name="ex_D1", sigma=sigma_D1,length=int(sigma_D1*num_sigma_D1),maxv=int(self.soccfg['gens'][qub_ch_D1]['maxv'])) # Here I didn't use pi pulse voltage divided by 2, instead I will directly use the gain for pi/2 pulse 
        self.add_gauss(ch=qub_ch_D2, name="ex_D2", sigma=sigma_D2,length=int(sigma_D2*num_sigma_D2),maxv=int(self.soccfg['gens'][qub_ch_D2]['maxv']))        
        #self.set_pulse_registers(ch=qub_ch, style="arb", waveform="ex")

        base_phase_reg_D1 = self.deg2reg(self.cfg["qub_phase_D1"], gen_ch=qub_ch_D1)
        base_phase_reg_D2 = self.deg2reg(self.cfg["qub_phase_D2"], gen_ch=qub_ch_D2)

        self.set_pulse_registers(ch=qub_ch_D1, style="arb", waveform="ex_D1", phase=base_phase_reg_D1)
        self.set_pulse_registers(ch=qub_ch_D2, style="arb", waveform="ex_D2", phase=base_phase_reg_D2)

        
        # give processor some time to configure pulses, I believe it sets the offset time 200 cycles into the future?
        # Try varying this and seeing if it moves. Also try putting a synci after the trigger in the body
        self.synci(200)   
    
    def body(self):
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.

        cfg    = self.cfg
        cav_ch = cfg["cav_channel"]
        qub_ch_D1  = cfg["qub_channel_D1"]
        qub_ch_D2  = cfg["qub_channel_D2"]


        if cfg["phase_reset"]:
            self.reset_phase(gen_ch = [qub_ch_D1, qub_ch_D2], t=0)
        else:
            pass

        
        # Choose which channel is pulse 1 and pulse 2
        order = cfg["pulse_order"]  # tuple/list like ("D1","D2"), ("D1","D1"), ("D2","D1")
        ch_map = {"D1": qub_ch_D1, "D2": qub_ch_D2}
        wf_map = {"D1": "ex_D1", "D2": "ex_D2"} # waveform map: links named waveforms to data in FPGA memory
        base_phase_deg = {"D1": cfg["qub_phase_D1"], "D2": cfg["qub_phase_D2"]}
        sigma_cyc = {
            "D1": self.us2cycles(cfg["qub_sigma_D1"], gen_ch=qub_ch_D1),
            "D2": self.us2cycles(cfg["qub_sigma_D2"], gen_ch=qub_ch_D2),
        }
        len_cyc = {
            "D1": int(sigma_cyc["D1"]*cfg["num_sigma_D1"]),
            "D2": int(sigma_cyc["D2"]*cfg["num_sigma_D2"]),
        }

        ch1 = ch_map[order[0]]
        ch2 = ch_map[order[1]]
        wf1 = wf_map[order[0]]
        wf2 = wf_map[order[1]]

        # convert base phases to regs for the *appropriate* channels
        base1_reg = self.deg2reg(base_phase_deg[order[0]], gen_ch=ch1)
        base2_reg = self.deg2reg(base_phase_deg[order[1]], gen_ch=ch2)
        
        # Measurement Timing
        offset    = self.us2cycles(cfg["adc_trig_offset"], gen_ch=cav_ch)
        meas_time = self.us2cycles(cfg["meas_time"],       gen_ch=cav_ch)
        

        # Place the 2nd pulse at a fixed offset before measurement; 1st is placed len1 (+ gap) earlier.
        gap   = cfg.get("gap", 0.0)  # default 0 us
        gap_cyc  = self.us2cycles(gap, gen_ch=ch1)

        # schedule by start times:
        t2_start = meas_time - self.us2cycles(cfg["qub_delay_fixed"], gen_ch=ch2) - len_cyc[order[1]]
        
        # Two sweep modes:
        mode = cfg.get("sweep_mode")
        
        if mode == 'phase':
            # 1st pulse immediately before 2nd plus optional gap
            #t1_start = t2_start - gap_cyc - len_cyc[order[0]]  # immediate: gap_cyc=0 -> adjacent pulses
            t1_start = t2_start - gap_cyc

            # Phase rotation for 2nd pulse: φ 
            phase2_deg  = float(cfg.get("phase2_deg", 0.0))  # degrees, set per sweep point
            phase2_reg  = self.deg2reg(phase2_deg, gen_ch=ch2)
            
        elif mode == 'interval':
            # Sweep the interval between the STARTS of pulse 1 and pulse 2.
            # tau is in microseconds; we translate to cycles and shift pulse 1 back.
            tau_us   = float(cfg.get("tau_us", 0.0))               ### NEW/CHANGED
            tau_cyc  = self.us2cycles(tau_us, gen_ch=ch1)
            t1_start = t2_start - tau_cyc
            
            # Allow a fixed user-set extra phase on pulse #2 (deg -> reg)
            phase2_deg = float(cfg.get("phase2_deg", 0.0))
            phase2_reg = self.deg2reg(phase2_deg, gen_ch=ch2)
            
        else:
            raise ValueError(f"Unsupported sweep_mode: {mode}")


        #Sets off the ADC
        self.trigger(adcs=self.ro_chs,
                    pins=[0],
                    adc_trig_offset=offset)
        
# =============================================================================
#         # ---- reset phases just before the pulses ----
#         # 3 cycles min; use a small safety margin (e.g. 8 cycles)
#         fabric_margin = 8  # cycles
#         t1_reset = max(0, t1_start - fabric_margin)
#         t2_reset = max(0, t2_start - fabric_margin)
#         
#         if t1_start == t2_start:
#             self.reset_phase(gen_ch=[ch1, ch2], t=t1_reset)
#         else:
#             self.reset_phase(gen_ch=[ch1,ch2], t=t1_reset)
#             self.reset_phase(gen_ch=[ch1,ch2], t=t2_reset)
# =============================================================================
        
        # π/2 #1 on chosen channel
        self.set_pulse_registers(ch=ch1, style="arb", waveform=wf1, phase=base1_reg)
        self.pulse(ch=ch1, t=t1_start)

        # π/2 #2 on chosen channel, with phase advance
        self.set_pulse_registers(ch=ch2, style="arb", waveform=wf2, phase=base2_reg + phase2_reg)
        self.pulse(ch=ch2, t=t2_start)
        
        # measurement pulse
        self.pulse(ch=cav_ch, t=meas_time)
        self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        self.sync_all(self.us2cycles(cfg["relax_delay"])) #Syncs to an offset time after the final pulse is sent


class CavitySweep(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   
        gen_ch = cfg["cav_channel"]

        # set the nyquist zone
        # Implement an if statement here to catch? 
        self.declare_gen(ch=cfg["cav_channel"], nqz=cfg["nqz_c"])
        
        # configure the readout lengths and downconversion frequencies (ensuring it is an available DAC frequency)
        
        readout = self.us2cycles(cfg["readout_length"],ro_ch=cfg["ro_channels"][0])
        for ch in cfg["ro_channels"]:
            self.declare_readout(ch=ch, length=readout,
                                 freq=self.cfg["cav_freq"], gen_ch=cfg["cav_channel"])

        # convert frequency to DAC freqency (ensuring it is an available ADC frequency)
        freq = self.freq2reg(cfg["cav_freq"],gen_ch=gen_ch, ro_ch=cfg["ro_channels"][0])
        phase = self.deg2reg(cfg["cav_phase"], gen_ch=gen_ch)
        gain = cfg["meas_gain"]
        self.default_pulse_registers(ch=gen_ch, freq=freq, phase=phase, gain=gain)

        
        self.set_pulse_registers(ch=gen_ch, style="const", length=self.us2cycles(self.cfg["meas_window"],gen_ch=gen_ch))

        
        # give processor some time to configure pulses, I believe it sets the offset time 200 cycles into the future?
        # Try varying this and seeing if it moves. Also try putting a synci after the trigger in the body
        self.synci(200)   
    
    def body(self):
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.
        offset = self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"])
        meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=self.cfg["cav_channel"])
        #Sets off the ADC
        self.trigger(adcs=self.ro_chs,
                    pins=[0],
                    adc_trig_offset=offset)
        
        #Sends measurement pulse
        self.pulse(ch=self.cfg["cav_channel"],t=meas_time)
        self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        self.sync_all(self.us2cycles(self.cfg["relax_delay"]))

def get_phase_rotation_settings():
    settings = {}
    
    settings['scanname'] = 'T2_phase_rotation_two_qubit_channels'
    settings['meas_type'] = 'Drive_crosstalk'
    settings['phase_reset'] = False
    
    settings['cav_freq'] = 6e9
    settings['meas_gain'] = 1000
    
    # qubit (D1)
    settings['qub_freq_D1']  = 4e9
    settings['qub_gain_D1']  = 1000
    settings['qub_phase_D1'] = 0.0  # deg
    settings['qub_sigma_D1'] = 0.02e-6
    settings['num_sigma_D1'] = 4

    # qubit (D2)
    settings['qub_freq_D2']  = 4e9
    settings['qub_gain_D2']  = 1000
    settings['qub_phase_D2'] = 0.0  # deg
    settings['qub_sigma_D2'] = 0.02e-6
    settings['num_sigma_D2'] = 4
    
    # --- Sweep controls ---
    settings['sweep_mode'] = 'phase'   ###: 'phase' or 'interval'
    
    # for 'phase' mode    
    settings['spacing']    = 'Linear'
    settings['phase2_min_deg']  = -180.0
    settings['phase2_max_deg']  =  +180.0
    settings['phase2_points']   =  73  
    
    # for 'interval' mode (units: microseconds)
    settings['phase2_deg']      = 0.0   # fixed extra phase for the 2nd pulse (deg)
    settings['tau_min_us']      = 0.0
    settings['tau_max_us']      = 2.0
    settings['tau_points']      = 51 


    settings['gap']              = 0.0            # immediate by default
    settings['pulse_order']         = ('D1','D2')    # ('D1','D1'), ('D2','D1'), etc.

    settings['qub_delay_fixed'] = 0.0  # example: 0.2 us
    settings['subtract_background'] = False

    #Card settings
    settings['reps'] = 1
    settings['soft_avgs'] = 5e3
    
    return settings

def meas_phase_rotation(soc,soccfg,instruments,settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse_D1   = exp_globals['qubit_pulse_D1']
    q_pulse_D2   = exp_globals['qubit_pulse_D2']
    
   # if exp_globals['LO']:
    if False:
        logen = instruments['LO']
        
        logen.freq   = exp_globals['LO_freq']
        logen.power  = exp_globals['LO_power']
        logen.output = 1
    else:
        print("Warning: No LO has been initialized.")

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel'],
        'qub_channel_D1'  : exp_globals['qub_channel_D1'],
        'qub_channel_D2'  : exp_globals['qub_channel_D2'],
        'ro_channels'     : exp_globals['ro_channels'],

        'nqz_c'           : 1,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'meas_gain'       : exp_settings['meas_gain'],
        'cav_freq'        : (exp_settings['cav_freq']-exp_globals['LO_freq']*exp_globals['LO'])/1e6,
        
        # qubits
        'nqz_q'           : 2,
        # D1
        'qub_phase_D1'    : q_pulse_D1['qub_phase'],
        'qub_freq_D1'     : (exp_settings['qub_freq_D1'])/1e6,  # MHz
        'qub_gain_D1'     : exp_settings['qub_gain_D1'],
        'qub_sigma_D1'    : q_pulse_D1['sigma'],
        'num_sigma_D1'    : q_pulse_D1['num_sigma'],
        # D2
        'qub_phase_D2'    : q_pulse_D2['qub_phase'],
        'qub_freq_D2'     : (exp_settings['qub_freq_D2'])/1e6,  # MHz
        'qub_gain_D2'     : exp_settings['qub_gain_D2'],
        'qub_sigma_D2'    : q_pulse_D2['sigma'],
        'num_sigma_D2'    : q_pulse_D2['num_sigma'],
        
        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],
        
        'sweep_mode'      : exp_settings['sweep_mode'],

        'phase2_deg'      : exp_settings.get('phase2_deg', 0.0), # can be assigned in the driver for 'interval' mode, otherwise default to zero; updated in the 'phase' mode
        'tau_us'          : 0.0, # us, updated in loop
        
        'qub_delay_fixed' : exp_globals['qub_delay_fixed'], # minimum delay between 2nd pi/2 pulse and measurement pulse, us
        'pulse_order'     : exp_settings['pulse_order'], # ('D1','D2'), ('D1','D1'), ('D2','D1'), etc.
        'gap'             : exp_settings['gap'], # us, physical time gap between the two pi/2 pulses


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs'],

        'phase_reset'     : exp_settings['phase_reset']
        }


    prog = T2_sequence(soccfg,config)
    meas_start = prog.us2cycles(m_pulse["init_buffer"],ro_ch=0)
    meas_end = meas_start+prog.us2cycles(m_pulse["meas_window"],ro_ch=0)
    total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)
    
    # Choose sweep based on mode
    sweep_mode = exp_settings['sweep_mode']
    
    if sweep_mode == 'phase':
        if exp_settings['spacing'] == 'Log':
            raise ValueError("Log spacing not meaningful for phase; use Linear.")
        xvals = np.linspace(
            exp_settings['phase2_min_deg'],
            exp_settings['phase2_max_deg'],
            exp_settings['phase2_points']
        )
        xvals = np.round(xvals, 3)
        xlabel = 'Phase2 (deg)'
    
    elif sweep_mode == 'interval':
        # tau in microseconds
        if exp_settings['spacing'] == 'Log':
            raise ValueError("Log spacing not meaningful for tau; use Linear.")
        tau_list = np.linspace(
            exp_settings['tau_min_us'],
            exp_settings['tau_max_us'],
            exp_settings['tau_points']
        )
        xvals  = np.round(tau_list, 3)
        xlabel = 'τ (us)'
    
    else:
        raise ValueError(f"Unsupported sweep_mode: {sweep_mode}")

    
    # Allocate arrays
    N = len(xvals)
    amp_int  = np.zeros(N); amp_orig = np.zeros(N)
    ang_int  = np.zeros(N); ang_orig = np.zeros(N)
    
    
    
    
    first_it = True
    
    if exp_settings['subtract_background']:
        #Acquire background trace
#            qubitgen.freq=3.8e9
#            time.sleep(0.1)
        print('Starting Background Trace')
        bprog = CavitySweep(soccfg,config)
        holder = bprog.acquire(soc, load_pulses=True, progress=False)
        print('Background Trace Complete')
        I_back = holder[0][0][0]
        Q_back = holder[1][0][0]
    else:
        I_back, Q_back = 0,0
    
    tstart = time.time()
    
    for index, x in enumerate(xvals):
        if sweep_mode == 'phase':
            config['phase2_deg'] = float(x)
        else:  # 'interval'
            config['tau_us'] = float(x)
    
        prog = T2_sequence(soccfg, config)
        holder = prog.acquire(soc, load_pulses=True, progress=False)
        I_sig = holder[0][0][0]
        Q_sig = holder[1][0][0]
    
        I_final = I_sig - I_back
        Q_final = Q_sig - Q_back
    
        amp_int[index]  = np.sqrt(I_final**2 + Q_final**2)
        ang_int[index]  = np.degrees(np.arctan2(Q_final, I_final))
        amp_orig[index] = np.sqrt(I_sig**2 + Q_sig**2)
        ang_orig[index] = np.degrees(np.arctan2(Q_sig, I_sig))
    
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, N)
            first_it = False
    
        # live plot
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121); plt.plot(xvals, amp_int)
        plt.suptitle('Live sweep (no fit)')
        plt.xlabel(xlabel); plt.ylabel('Amplitude')
        plt.subplot(122); plt.plot(xvals, ang_int)
        plt.xlabel(xlabel); plt.ylabel('Phase (deg)')
        fig.canvas.draw(); fig.canvas.flush_events()
    
    
        #last save at the end    
        plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)
        userfuncs.SaveFull(saveDir, filename, ['xvals','amp_int','ang_int','amp_orig','ang_orig','sweep_mode','xlabel'], locals(), expsettings=settings, instruments=instruments)
    
    
    


    t2 = time.time()
    
    print('Elapsed time: {}'.format(t2-tstart))
    

    userfuncs.SaveFull(saveDir, filename, ['xvals','amp_int','ang_int','amp_orig','ang_orig','sweep_mode','xlabel'],
                         locals(), expsettings=settings, instruments=instruments)
    
    
    if exp_globals['LO']:
        pass
        #logen.output = 0
        
    return prog

    


   
    

        