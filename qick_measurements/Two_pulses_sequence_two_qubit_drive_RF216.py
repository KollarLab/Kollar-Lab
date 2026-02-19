# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:16:54 2023

@author: kollarlab
"""

from qick.asm_v2 import AveragerProgramV2

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
import time

from utility.userfits import fit_T2
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import estimate_time

class T2_sequence(AveragerProgramV2):
    def _initialize(self,cfg): 
        ro_ch  = cfg['ro_channel']
        gen_ch = cfg["cav_channel"]
        qub_ch_D1 = cfg["qub_channel_D1"]
        qub_ch_D2 = cfg["qub_channel_D2"]

        # set the nyquist zone
        # Implement an if statement here to catch? 
        self.declare_gen(ch=cfg["cav_channel"], nqz=cfg["nqz_c"], mixer_freq=cfg['cav_mixer_freq'],ro_ch=ro_ch)
        self.declare_gen(ch=qub_ch_D1, nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq'])
        self.declare_gen(ch=qub_ch_D2, nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq'])
        
        self.declare_readout(ch=ro_ch, length=cfg['readout_length'])
        
        
        self.add_readoutconfig(ch=ro_ch, name="myro",
                           freq=cfg['cav_freq'],
                           gen_ch=gen_ch,
                           outsel='product')
        self.add_pulse(ch=gen_ch, name="cav_pulse", ro_ch=ro_ch,
               style="const",
               freq=cfg['cav_freq'],
               length= cfg["meas_window"],
               phase=cfg['cav_phase'],
               gain=cfg['cav_gain'],
              )
        
        # Envelopes
        sigma1 = cfg["qub_sigma_D1"]; ns1 = cfg["num_sigma_D1"]
        sigma2 = cfg["qub_sigma_D2"]; ns2 = cfg["num_sigma_D2"]
        
        self.add_gauss(ch=qub_ch_D1, name="env_D1", sigma=sigma1, length=sigma1*ns1)
        self.add_gauss(ch=qub_ch_D2, name="env_D2", sigma=sigma2, length=sigma2*ns2)
        
        # pulse #1
        self.add_pulse(ch=qub_ch_D1, name="ex1_D1", ro_ch=ro_ch,
                       style="flat_top", envelope="env_D1",
                       freq=cfg["qub_freq"],
                       length=cfg.get("hold_length_D1", 0.0),
                       phase=cfg["qub_phase_D1"],
                       gain=cfg["qub_gain_D1"])
        
        self.add_pulse(ch=qub_ch_D2, name="ex1_D2", ro_ch=ro_ch,
                       style="flat_top", envelope="env_D2",
                       freq=cfg["qub_freq"],
                       length=cfg.get("hold_length_D2", 0.0),
                       phase=cfg["qub_phase_D2"],
                       gain=cfg["qub_gain_D2"])
        
        # pulse #2
        self.add_pulse(ch=qub_ch_D1, name="ex2_D1", ro_ch=ro_ch,
                       style="flat_top", envelope="env_D1",
                       freq=cfg["qub_freq"],
                       length=cfg.get("hold_length_D1", 0.0),
                       phase=cfg.get("qub_phase2_D1", cfg["qub_phase_D1"]), #find the definition of qub_phase2_Dx later
                       gain=cfg["qub_gain_D1"])
        
        self.add_pulse(ch=qub_ch_D2, name="ex2_D2", ro_ch=ro_ch,
                       style="flat_top", envelope="env_D2",
                       freq=cfg["qub_freq"],
                       length=cfg.get("hold_length_D2", 0.0),
                       phase=cfg.get("qub_phase2_D2", cfg["qub_phase_D2"]),
                       gain=cfg["qub_gain_D2"])


        
        # Phase reset pulses (phrst=1 resets phase accumulator)
        self.add_pulse(ch=qub_ch_D1, name="phrst_D1", ro_ch=ro_ch,
                       style="const",
                       freq=cfg["qub_freq"], phase=0, gain=0,
                       length=0.015, phrst=1)
        
        self.add_pulse(ch=qub_ch_D2, name="phrst_D2", ro_ch=ro_ch,
                       style="const",
                       freq=cfg["qub_freq"], phase=0, gain=0,
                       length=0.015, phrst=1)
        
        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)
        
        
    
    def _body(self, cfg):
        # ---- Select which drive line is pulse1 and pulse2 ----
        order = tuple(cfg["pulse_order"])   # e.g. ("D1","D2"), ("D1","D1"), ("D2","D1")
        ch_map = {"D1": cfg["qub_channel_D1"], "D2": cfg["qub_channel_D2"]}
        p1_map = {"D1": "ex1_D1", "D2": "ex1_D2"}
        p2_map = {"D1": "ex2_D1", "D2": "ex2_D2"}

        
        ch1 = ch_map[order[0]]
        ch2 = ch_map[order[1]]
        p1 = p1_map[order[0]]
        p2 = p2_map[order[1]]
        
        meas_time = cfg["meas_time"]           # (us)
        gap      = float(cfg.get("gap", 0.0))  # (us)
        
        # Estimate pulse lengths (us) for scheduling
        def pulse_len(which):
            if which == "D1":
                return cfg.get("hold_length_D1", 0.0) + cfg["qub_sigma_D1"]*cfg["num_sigma_D1"]
            else:
                return cfg.get("hold_length_D2", 0.0) + cfg["qub_sigma_D2"]*cfg["num_sigma_D2"]
        
        len2 = pulse_len(order[1])
        
        # Pulse #2 ends qub_delay_fixed before measurement pulse
        t2_start = meas_time - float(cfg["qub_delay_fixed"]) - len2
        
        mode = cfg.get("sweep_mode", "phase")

        len1 = pulse_len(order[0])  # <-- add this
        
        if mode == "phase":
            # start-to-start gap (fixed)
            t1_start = t2_start - gap
        
        elif mode == "interval":
            # start-to-start interval (sweep)
            tau_us = float(cfg.get("tau_us", 0.0))
            t1_start = t2_start - tau_us
                       
        else:
            raise ValueError(f"Unsupported sweep_mode: {mode}")

        # ---------------------------------------------------------
        # SAFEGUARD: ensure BOTH pulses end before cav_pulse starts
        # (keep at least qub_delay_fixed margin before measurement)
        # ---------------------------------------------------------
        end1 = t1_start + len1
        end2 = t2_start + len2

        latest_end = max(end1, end2)

        # Require all qubit activity ends before (meas_time - qub_delay_fixed)
        latest_allowed_end = float(meas_time) - float(cfg.get("qub_delay_fixed", 0.0))

        if latest_end > latest_allowed_end:
            shift_earlier = latest_end - latest_allowed_end
            t1_start -= shift_earlier
            t2_start -= shift_earlier

        # Optional: fail early if we shifted before t=0
        earliest_start = min(t1_start, t2_start)
        if earliest_start < 0:
            raise ValueError(
                f"Cannot safely schedule pulses: earliest_start={earliest_start:.6f} us < 0 after safeguard.\n"
                f"Increase meas_time/meas_pos, or reduce pulse lengths, or reduce required qub_delay_fixed."
            )


        if order[0] == order[1]:
            # same channel: must be sequential (no overlap)
            if (t1_start + len1) > t2_start:
                raise ValueError(
                    f"Same-channel pulses overlap: order={order}, "
                    f"t1_start={t1_start:.6f} us, len1={len1:.6f} us, t2_start={t2_start:.6f} us"
                )

        
        # optional phase reset
        if cfg.get("phase_reset", False):
            self.pulse(ch=cfg["qub_channel_D1"], name="phrst_D1", t=0)
            self.pulse(ch=cfg["qub_channel_D2"], name="phrst_D2", t=0)
        
        # trigger ADC
        self.trigger(ros=[cfg["ro_channel"]], pins=[0], t=cfg["adc_trig_offset"])

        
        # ---- Play pulses ----
        self.pulse(ch=ch1, name=p1, t=t1_start)
        self.pulse(ch=ch2, name=p2, t=t2_start)
        
        # ---- Measurement ----
        self.pulse(ch=cfg["cav_channel"], name="cav_pulse", t=meas_time)
        self.wait_auto()
        self.delay(cfg["relax_delay"])

      


def get_T2_settings():
    settings = {}
    
    settings['scanname'] = 'T2_meas'
    settings['meas_type'] = 'T2_meas'
    # settings['phase_reset'] = True
    
    # settings['cav_freq'] = 1e9
    # settings['meas_gain'] = 1000
    
    # settings['qub_freq'] = 4e9
    # settings['qub_gain'] = 1000
    
    # #Sweep parameters    
    # settings['Tau_min']    = 200e-9
    # settings['Tau_max']    = 30e-6
    # settings['Tau_points'] = 5
    settings['spacing']    = 'Linear'
    
    # settings['T2_guess'] = 10e-6
    # settings['detuning'] = 2e6
    # #Card settings
    # settings['reps'] = 1
    # settings['soft_avgs'] = 5e3
    
    return settings

def meas_T2(soc,soccfg,instruments,settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    
    
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
        'cav_channel'     : exp_globals['cav_channel']['ID'],
        'qub_channel_D1'  : exp_globals['qub_channel_1']['ID'],
        'qub_channel_D2'  : exp_globals['qub_channel_2']['ID'],
        'ro_channel'      : exp_globals['ro_channel']['ID'],

        'nqz_c'           : 2,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'cav_gain'        : exp_settings['cav_gain'],
        'cav_freq'        : (exp_settings['cav_freq']-exp_globals['LO_freq']*exp_globals['LO'])/1e6,
        'cav_mixer_freq'  : (exp_settings['cav_freq'] + exp_settings['cav_mixer_detuning'])/1e6,
        
        'nqz_q'           : 2,
        'qub_phase_D1'    : exp_globals['qubit_pulse_D1']['qub_phase'],
        'qub_phase_D2'    : exp_globals['qubit_pulse_D2']['qub_phase'],
        
        'qub_phase2_D1'   : exp_settings.get('qub_phase_D1', exp_globals['qubit_pulse_D1']['qub_phase']),
        'qub_phase2_D2'   : exp_settings.get('qub_phase_D2', exp_globals['qubit_pulse_D2']['qub_phase']), #just a placeholder here, we overwrite it in every sweep

        
        
        'qub_freq'        : exp_settings['qub_freq']/1e6,
        'qub_mixer_freq'  : (exp_settings['qub_freq']+exp_settings['qub_mixer_detuning'])/1e6,
        'qub_gain_D1'     : exp_settings['qub_gain_D1'],
        'qub_gain_D2'     : exp_settings['qub_gain_D2'],
        
        'qub_sigma_D1'    : exp_globals['qubit_pulse_D1']['sigma'],
        'num_sigma_D1'    : exp_globals['qubit_pulse_D1']['num_sigma'],
        'qub_sigma_D2'    : exp_globals['qubit_pulse_D2']['sigma'],
        'num_sigma_D2'    : exp_globals['qubit_pulse_D2']['num_sigma'],
        'qub_delay_fixed' : exp_globals['qub_delay_fixed'],
        'hold_length_D1'  : exp_settings['hold_time_D1_us'],
        'hold_length_D2'  : exp_settings['hold_time_D2_us'],
        
        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],

        'phase_reset'     : exp_settings['phase_reset'],
        
        'sweep_mode'      : exp_settings['sweep_mode'],     # 'phase' or 'interval'
        'pulse_order'     : exp_settings['pulse_order'],    # ('D1','D2'), etc.
        'gap'             : exp_settings.get('gap', 0.0),    # us
        
        'phase2_deg'      : exp_settings.get('phase2_deg', 0.0), # used in interval mode (fixed) or overwritten in loop
        'tau_us'          : 0.0,                                 # overwritten in interval loop

        }
    
    cav_ch = exp_globals['cav_channel']
    #qub_ch = exp_globals['qub_channel']
    qub_ch1 = exp_globals['qub_channel_1']
    qub_ch2 = exp_globals['qub_channel_2']
    

    ro_ch  = exp_globals['ro_channel']
    # Set attenuator on DAC.
    soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
    soc.rfb_set_gen_rf(qub_ch1['ID'], qub_ch1['Atten_1'], qub_ch1['Atten_2'])
    soc.rfb_set_gen_rf(qub_ch2['ID'], qub_ch2['Atten_1'], qub_ch2['Atten_2'])
    # Set attenuator on ADC.
    soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])
    
    if exp_settings['filter'] == 'all_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
        
        soc.rfb_set_gen_filter(config['qub_channel_D1'], fc=config['qub_freq']/1000, ftype='bandpass', bw=qub_ch1['BW'])
        soc.rfb_set_gen_filter(config['qub_channel_D2'], fc=config['qub_freq']/1000, ftype='bandpass', bw=qub_ch2['BW'])

        
    elif exp_settings['filter'] == 'no_qubit_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
    
        soc.rfb_set_gen_filter(config['qub_channel_D1'], fc=config['qub_freq']/1000, ftype='bypass')
        soc.rfb_set_gen_filter(config['qub_channel_D2'], fc=config['qub_freq']/1000, ftype='bypass')

        
    elif exp_settings['filter'] == 'no_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bypass')
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
        
        soc.rfb_set_gen_filter(config['qub_channel_D1'], fc=config['qub_freq']/1000, ftype='bypass')
        soc.rfb_set_gen_filter(config['qub_channel_D2'], fc=config['qub_freq']/1000, ftype='bypass')
        
    else:
        print('Please select one option from:')
        print('\'all_filter\', \'no_qubit_filter\', and \'no_filter\'')
        return


    # ----------------------------
    # Build sweep axis (script #1 style: ONE axis)
    # ----------------------------
    sweep_mode = exp_settings["sweep_mode"]  # 'phase' or 'interval'
    
    if sweep_mode == "phase":
        # phase sweep (deg)
        if exp_settings.get("spacing", "Linear") == "Log":
            raise ValueError("Log spacing not meaningful for phase; use Linear.")
        xvals = np.linspace(exp_settings["phase2_min_deg"],
                            exp_settings["phase2_max_deg"],
                            exp_settings["phase2_points"])
        xvals = np.round(xvals, 3)
        xlabel = "Phase2 (deg)"
    
    elif sweep_mode == "interval":
        # interval between STARTS (us)
        if exp_settings.get("spacing", "Linear") == "Log":
            raise ValueError("Log spacing not meaningful for tau; use Linear.")
        xvals = np.linspace(exp_settings["tau_min_us"],
                            exp_settings["tau_max_us"],
                            exp_settings["tau_points"])
        xvals = np.round(xvals, 6)
        xlabel = "τ (us)"
           
    else:
        raise ValueError(f"Unsupported sweep_mode: {sweep_mode}")
    
    # Allocate arrays (script #1 style)
    N = len(xvals)
    amp_int  = np.zeros(N)
    ang_int  = np.zeros(N)
    amp_orig = np.zeros(N)
    ang_orig = np.zeros(N)
    
    # Save base phases so we can restore and only rotate pulse #2
    base_phase_D1 = float(config["qub_phase_D1"])
    base_phase_D2 = float(config["qub_phase_D2"])
    
   
    first_it = True
# =============================================================================
#     if exp_settings["subtract_background"]:
#         print("Starting Background Trace")
#         bprog = CavitySweep(soccfg, config)
#         holder = bprog.acquire(soc, load_pulses=True, progress=False)
#         print("Background Trace Complete")
#         I_back = holder[0][0][0]
#         Q_back = holder[1][0][0]
#     else:
# =============================================================================
    I_back, Q_back = 0.0, 0.0
    
    # ----------------------------
    # Main sweep loop (script #1 style)
    # ----------------------------
    tstart = time.time()
    
    for idx, x in enumerate(xvals):
    
        if sweep_mode == "phase":
            dphi = float(x)  # this is "qub_phase2" (deg), the relative phase between pulse2 and pulse1
            order = tuple(config["pulse_order"])  # ("D1","D2"), etc.
        
            # phase of pulse #1 depends on which channel pulse #1 is on
            phi1 = base_phase_D1 if order[0] == "D1" else base_phase_D2
            phi2 = phi1 + dphi
        
            # pulse #1 phases stay at base (so ex1_D1/ex1_D2 are correct)
            config["qub_phase_D1"] = base_phase_D1
            config["qub_phase_D2"] = base_phase_D2
        
            # pulse #2 phase should be phi2 on whichever channel pulse #2 uses
        
            if order[1] == "D1":
                config["qub_phase2_D1"] = phi2
            else:  # "D2"
                config["qub_phase2_D2"] = phi2

    
        elif sweep_mode == "interval":
            config["tau_us"] = float(x)
            order = tuple(config["pulse_order"])
        
            phi1 = base_phase_D1 if order[0] == "D1" else base_phase_D2
            dphi_fixed = float(exp_settings.get("phase2_deg", 0.0))
            phi2 = phi1 + dphi_fixed
        
            config["qub_phase_D1"] = base_phase_D1
            config["qub_phase_D2"] = base_phase_D2
            config["qub_phase2_D1"] = base_phase_D1
            config["qub_phase2_D2"] = base_phase_D2
        
            if order[1] == "D1":
                config["qub_phase2_D1"] = phi2
            else:
                config["qub_phase2_D2"] = phi2
                     
    
        # Build + acquire (rebuild per point, like script #1)
        prog = T2_sequence(soccfg, reps=exp_settings["reps"], final_delay=None, final_wait=0, cfg=config)
        holder = prog.acquire(soc, rounds=exp_settings["rounds"], load_pulses=True, progress=False)
    
        # Read I/Q 
        I_sig = holder[0][0][0]
        Q_sig = holder[0][0][1]
    
        I_final = I_sig - I_back
        Q_final = Q_sig - Q_back
    
        amp_int[idx]  = np.sqrt(I_final**2 + Q_final**2)
        ang_int[idx]  = np.degrees(np.arctan2(Q_final, I_final))
    
        amp_orig[idx] = np.sqrt(I_sig**2 + Q_sig**2)
        ang_orig[idx] = np.degrees(np.arctan2(Q_sig, I_sig))
    
        # time estimate after first point
        if first_it:
            estimate_time(tstart, time.time(), N)
            first_it = False
    
        # live plot (script #1 style)
        fig = plt.figure(1, figsize=(13, 8))
        plt.clf()
        plt.subplot(121)
        plt.plot(xvals, amp_int)
        plt.xlabel(xlabel)
        plt.ylabel("Amplitude")
        plt.suptitle(f"{exp_settings['scanname']}_{stamp}")
        plt.subplot(122)
        plt.plot(xvals, ang_int)
        plt.xlabel(xlabel)
        plt.ylabel("Phase (deg)")
        fig.canvas.draw()
        fig.canvas.flush_events()
    
    # ----------------------------
    # Save at end (script #1 style)
    # ----------------------------
    plt.savefig(os.path.join(saveDir, filename + "_no_fit.png"), dpi=150)
    userfuncs.SaveFull(
        saveDir, filename,
        ["xvals", "amp_int", "ang_int", "amp_orig", "ang_orig", "sweep_mode", "xlabel"],
        locals(), expsettings=settings, instruments=instruments
    )
    
    print(f"Elapsed time: {time.time() - tstart:.2f} s")



   
    

        