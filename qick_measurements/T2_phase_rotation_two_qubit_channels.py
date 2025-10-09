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
            self.reset_phase(gen_ch = [cav_ch, qub_ch_D1, qub_ch_D2], t=0)
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
        
        # Timing – 2nd pulse starts immediately after 1st ends (plus optional tiny gap)
        offset    = self.us2cycles(cfg["adc_trig_offset"], gen_ch=cav_ch)
        meas_time = self.us2cycles(cfg["meas_time"],       gen_ch=cav_ch)

        # Place the 2nd pulse at a fixed offset before measurement; 1st is placed len1 (+ gap) earlier.
        gap   = cfg.get("gap", 0.0)  # default 0 us
        gap_cyc  = self.us2cycles(gap, gen_ch=ch2)

        # end time of 2nd pulse is at (meas_time - qub_delay_fixed), like before
        # we schedule by start times:
        t2_start = meas_time - self.us2cycles(cfg["qub_delay_fixed"], gen_ch=ch2) - len_cyc[order[1]]
        t1_start = t2_start - gap_cyc - len_cyc[order[0]]  # immediate: gap_cyc=0 -> adjacent pulses

        # Phase rotation for 2nd pulse: φ 
        phase2_deg  = float(cfg.get("phase2_deg", 0.0))  # degrees, set per sweep point
        phase2_reg  = self.deg2reg(phase2_deg, gen_ch=ch2)


        #Sets off the ADC
        self.trigger(adcs=self.ro_chs,
                    pins=[0],
                    adc_trig_offset=offset)
        
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
    
    #Sweep parameters    
    settings['spacing']    = 'Linear'
    settings['phase2_min_deg']  = -180.0
    settings['phase2_max_deg']  =  +180.0
    settings['phase2_points']   =  73  


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

        'phase2_deg'      : 0.0, # deg, updated in loop
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
    
    # Build phase sweep (degrees)
    if exp_settings['spacing'] == 'Log':
        raise ValueError("Log spacing not meaningful for phase; use Linear.")
    else:
        phase_list = np.linspace(
            exp_settings['phase2_min_deg'],
            exp_settings['phase2_max_deg'],
            exp_settings['phase2_points']
        )
    phases = np.round(phase_list, 6)
    
    amp_int = np.zeros(len(phases))
    amp_orig = np.zeros(len(phases))

    
    ang_int = np.zeros(len(phases))
    ang_orig = np.zeros(len(phases))

    
    
    
    
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
    
    for pind, phi in enumerate(phases):
        
        print(f'Phase2 (deg): {phi}')
        config['phase2_deg'] = float(phi)
        
        prog = T2_sequence(soccfg,config)
        holder = prog.acquire(soc, load_pulses=True, progress=False)
        I_sig = holder[0][0][0]
        Q_sig = holder[1][0][0]
        
              
#        I_final, Q_final   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back
        
#        amps[tind] = np.sqrt(I_full**2+Q_full**2)
#        angles[tind] = np.arctan2(Q_full,I_full)*180/np.pi
#        amp_int[tind] = np.sqrt(I_final**2+Q_final**2)
#        ang_int[tind] = np.arctan2(Q_final, I_final)*180/np.pi
        
        

        amp_int[pind] = np.sqrt(I_final**2 + Q_final**2)
        ang_int[pind] = np.arctan2(Q_final, I_final)*180/np.pi

        amp_orig[pind] = np.sqrt(I_sig**2 + Q_sig**2)
        ang_orig[pind] = np.arctan2(Q_sig, I_sig)*180/np.pi
        
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(phases))
            
                
            first_it = False  
            
        #_______________________________________#
    
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(phases, amp_int)
        plt.suptitle('Live phase rotation data (no fit)')#, {} pi pulses'.format(exp_settings['pulse_count']))
        plt.xlabel('Phase2 (deg)')
        plt.ylabel('Amplitude')
        plt.subplot(122)
        plt.plot(phases, ang_int)
        plt.xlabel('Phase2 (deg)')
        plt.ylabel('Phase (deg)')
        fig.canvas.draw()
        fig.canvas.flush_events()

        # plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)
        # userfuncs.SaveFull(saveDir, filename, ['taus', 'amp_int', 'ang_int'], locals(), expsettings=settings, instruments=instruments)

    #last save at the end    
    plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)
    userfuncs.SaveFull(saveDir, filename, ['phases', 'amp_int', 'ang_int'], locals(), expsettings=settings, instruments=instruments)
    
    
    # if exp_settings['debug']:

    #     config['readout_length'] = 3
    #     config['adc_trig_offset'] = m_pulse['emp_delay'] + m_pulse['meas_pos'] + exp_settings['debug_time']
   
    #     total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)

    #     amp_intd = np.zeros(len(phases))
    #     ang_intd = np.zeros(len(phases))
    #     ampsd    = np.zeros((len(phases),total_samples))
    #     anglesd  = np.zeros((len(phases),total_samples))      
   
    #     if exp_settings['subtract_background']:
    #         #Acquire background trace
    # #            qubitgen.freq=3.8e9
    # #            time.sleep(0.1)
    #         print('Starting Background Trace')
    #         bprog = CavitySweep(soccfg,config)
    #         holder = bprog.acquire_decimated(soc, load_pulses=True, progress=False)
    #         print('Background Trace Complete')
    #         I_full_bd = holder[0][0]
    #         Q_full_bd = holder[0][1]
    #         I_window_bd = I_full_bd[meas_start:meas_end]
    #         Q_window_bd = Q_full_bd[meas_start:meas_end]     
    #     else:
    #         I_window_bd, Q_window_bd, I_full_bd, Q_full_bd = 0,0,0,0
            
    #     I_backd, Q_backd = [np.mean(I_window_bd), np.mean(Q_window_bd)]
        
    #     first_it = True
   
    #     for tind in indices:
            
    #         tau = taus[tind]
    #         print('Tau: {}'.format(tau))
            
    #         config['qub_delay_t2'] = tau*1e6
    #         prog = T2_sequence(soccfg,config)
    #         holder = prog.acquire_decimated(soc, load_pulses=True, progress=False)
    #         I_fulld = holder[0][0]
    #         Q_fulld = holder[0][1]
    #         I_windowd = I_fulld[meas_start:meas_end]
    #         Q_windowd = Q_fulld[meas_start:meas_end]
            
    #         ##No background subtraction here!!!!
    #         I_sigd, Q_sigd   = [np.mean(I_windowd), np.mean(Q_windowd)] #<I>, <Q> for signal trace
            
    #         I_finald = I_sigd-I_backd #compute <I_net> in the data window
    #         Q_finald = Q_sigd-Q_backd
   
    #         ampsd[tind] = np.sqrt((I_fulld-I_full_bd)**2+(Q_fulld-Q_full_bd)**2)
    #         anglesd[tind] = np.arctan2((Q_fulld-Q_full_bd), (I_fulld-I_full_bd))*180/np.pi
    #         amp_intd[tind] = np.sqrt(I_finald**2+Q_finald**2)
    #         ang_intd[tind] = np.arctan2(Q_finald, I_finald)*180/np.pi

    #         if first_it:
    #             tstop = time.time()
    #             estimate_time(tstart, tstop, len(taus))
                
    #             xaxis = np.linspace(0,len(I_fulld)-1,len(I_fulld))

    #             for x in range(0,len(xaxis)):
    #                 xaxis[x] = prog.cycles2us(xaxis[x],ro_ch=0)
                    
    #             first_it = False

    #         fig = plt.figure(3, figsize=(13,8))
    #         plt.clf()
    #         plt.subplot(121)
    #         plt.plot(taus*1e6, amp_intd)
    #         plt.xlabel('Tau (us)')
    #         plt.ylabel('Amplitude')  
    #         plt.subplot(122)
    #         plt.plot(taus*1e6, ang_intd)
    #         plt.xlabel('Tau (us)')
    #         plt.ylabel('Phase')  
    #         plt.suptitle('Live T2 data debug (no fit)\n'+filename)
    #         fig.canvas.draw()
    #         fig.canvas.flush_events()
    #         plt.savefig(os.path.join(saveDir, filename+'_debug_no_fit.png'), dpi = 150)
        
    #         fig2 = plt.figure(4,figsize=(13,8))
    #         plt.clf()
        
    #         ax = plt.subplot(121)
    #         general_colormap_subplot(ax, xaxis*1e6, taus*1e6, ampsd, ['Time (us)', 'Tau (us)'], 'Raw data\n'+filename)
    #         ax = plt.subplot(122)
    #         general_colormap_subplot(ax, xaxis*1e6, taus*1e6, anglesd, ['Time (us)', 'Tau (us)'], 'Raw data\n'+filename)
    #         plt.savefig(os.path.join(saveDir, filename+'_debug_fulldata.png'), dpi = 150)


    t2 = time.time()
    
    print('Elapsed time: {}'.format(t2-tstart))
    
    # T2_guess     = exp_settings['T2_guess']
    # amp_guess    = max(amp_int)-min(amp_int)
    # offset_guess = np.mean(amp_int[-10:])
    # freq_guess   = exp_settings['detuning'] # No detuning mode, always does this
    # phi_guess    = 0
    
    # fit_guess = [T2_guess, amp_guess, offset_guess, freq_guess, phi_guess]
    # T2, amp, offset, freq, phi, fit_xvals, fit_yvals = fit_T2(taus, amp_int, fit_guess)
    # fig3 = plt.figure(3)
    # plt.clf()
    # plt.plot(taus*1e6, amp_int)
    # plt.plot(fit_xvals*1e6, fit_yvals)
    # plt.title('T2:{}us freq:{}MHz. \n {}'.format(np.round(T2*1e6,3), np.round(freq/1e6, 3),  filename)) # {} pi pulses, exp_settings['pulse_count'],
    # plt.xlabel('Time (us)')
    # plt.ylabel('Amplitude')
    # fig3.canvas.draw()
    # fig3.canvas.flush_events()
    # plt.savefig(os.path.join(saveDir, filename+'_fit.png'), dpi=150)

    userfuncs.SaveFull(saveDir, filename, ['phases','ang_int', 'amp_int', 'amp_orig','ang_orig'],
                         locals(), expsettings=settings, instruments=instruments)
    
    
    if exp_globals['LO']:
        pass
        #logen.output = 0

    


   
    

        