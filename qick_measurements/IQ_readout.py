# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 16:53:27 2025

@author: Kollarlab
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

class IQReadoutProgram(AveragerProgram):
    """
    Versatile readout program:
      - pulse_scheme: 'none' (ground), 'one', 'two'
      - drive channels: 'D1' / 'D2' (same mapping as your code)
      - two-pulse controls: tau_us (interval), phase2_deg (extra phase for pulse #2)
      - always resets phase at beginning (both channels)
    """
    def initialize(self):
        cfg = self.cfg
        cav_ch    = cfg["cav_channel"]
        qub_ch_D1 = cfg["qub_channel_D1"]
        qub_ch_D2 = cfg["qub_channel_D2"]

        # Nyquist zones
        self.declare_gen(ch=cav_ch,    nqz=cfg["nqz_c"])
        self.declare_gen(ch=qub_ch_D1, nqz=cfg["nqz_q"])
        self.declare_gen(ch=qub_ch_D2, nqz=cfg["nqz_q"])

        # Readout declaration
        readout = self.us2cycles(cfg["readout_length"], ro_ch=cfg["ro_channels"][0])
        for ch in cfg["ro_channels"]:
            self.declare_readout(ch=ch, length=readout,
                                 freq=cfg["cav_freq"], gen_ch=cav_ch)

        # Cav registers
        freq_c  = self.freq2reg(cfg["cav_freq"], gen_ch=cav_ch, ro_ch=cfg["ro_channels"][0])
        phase_c = self.deg2reg(cfg["cav_phase"], gen_ch=cav_ch)
        self.default_pulse_registers(ch=cav_ch, freq=freq_c, phase=phase_c, gain=cfg["meas_gain"])
        self.set_pulse_registers(ch=cav_ch, style="const",
                                 length=self.us2cycles(cfg["meas_window"], gen_ch=cav_ch))

        # Qubit waveforms (we make both available; selection happens in body)
        # D1
        fq1   = self.freq2reg(cfg["qub_freq_D1"], gen_ch=qub_ch_D1)
        self.default_pulse_registers(ch=qub_ch_D1, freq=fq1)
        sigma1 = self.us2cycles(cfg["qub_sigma_D1"], gen_ch=qub_ch_D1)
        self.add_gauss(ch=qub_ch_D1, name="ex_D1",
                       sigma=sigma1,
                       length=int(sigma1*cfg["num_sigma_D1"]),
                       maxv=int(self.soccfg['gens'][qub_ch_D1]['maxv']))

        # D2
        fq2   = self.freq2reg(cfg["qub_freq_D2"], gen_ch=qub_ch_D2)
        self.default_pulse_registers(ch=qub_ch_D2, freq=fq2)
        sigma2 = self.us2cycles(cfg["qub_sigma_D2"], gen_ch=qub_ch_D2)
        self.add_gauss(ch=qub_ch_D2, name="ex_D2",
                       sigma=sigma2,
                       length=int(sigma2*cfg["num_sigma_D2"]),
                       maxv=int(self.soccfg['gens'][qub_ch_D2]['maxv']))

        # Cache flat-top lengths for convenience
        self.len_flat = {
            "D1": int(sigma1*cfg["num_sigma_D1"] + self.us2cycles(cfg["hold_time_D1_us"], gen_ch=qub_ch_D1)),
            "D2": int(sigma2*cfg["num_sigma_D2"] + self.us2cycles(cfg["hold_time_D2_us"], gen_ch=qub_ch_D2)),
        }
        self.synci(200)

    def body(self):
        cfg = self.cfg
        cav_ch    = cfg["cav_channel"]
        qub_ch_D1 = cfg["qub_channel_D1"]
        qub_ch_D2 = cfg["qub_channel_D2"]
        

        # Maps
        ch_map = {"D1": qub_ch_D1, "D2": qub_ch_D2}
        wf_map = {"D1": "ex_D1",   "D2": "ex_D2"}
        base_phase_deg = {"D1": cfg["qub_phase_D1"], "D2": cfg["qub_phase_D2"]}
        base_gain_map  = {"D1": cfg["qub_gain_D1"],  "D2": cfg["qub_gain_D2"]}

        # Timing
        offset    = self.us2cycles(cfg["adc_trig_offset"], gen_ch=cav_ch)
        meas_time = self.us2cycles(cfg["meas_time"],       gen_ch=cav_ch)

        # Unpack scheme
        scheme = cfg.get("pulse_scheme", "none")  # 'none'|'one'|'two'

        # ----- Build qubit pulses if needed -----
        if scheme == "none":
            # Ground measurement: no qubit pulses
            pass

        elif scheme == "one":
            drv = cfg["single_drive"]           # 'D1' or 'D2'
            ch  = ch_map[drv]
            wf  = wf_map[drv]
            ht  = self.us2cycles(cfg[f"hold_time_{drv}_us"], gen_ch=ch)
            base_reg = self.deg2reg(base_phase_deg[drv], gen_ch=ch)
            gain = base_gain_map[drv]

            # Start the single pulse a fixed time before meas_time
            t_start = meas_time - self.us2cycles(cfg["qub_delay_fixed"], gen_ch=ch) - self.len_flat[drv]

            self.set_pulse_registers(ch=ch, style="flat_top", waveform=wf,
                                     length=ht, phase=base_reg, gain=gain)
            
            # --- always reset phase at the start ---
            self.reset_phase(gen_ch=[ch], t=0)
            
            self.pulse(ch=ch, t=t_start)

        elif scheme == "two":
            # Use pulse_order = ('D1','D2') etc.
            order = cfg["pulse_order"]
            ch1, ch2 = ch_map[order[0]], ch_map[order[1]]
            wf1, wf2 = wf_map[order[0]], wf_map[order[1]]
            ht1 = self.us2cycles(cfg[f"hold_time_{order[0]}_us"], gen_ch=ch1)
            ht2 = self.us2cycles(cfg[f"hold_time_{order[1]}_us"], gen_ch=ch2)
            g1, g2 = base_gain_map[order[0]], base_gain_map[order[1]]

            base1 = self.deg2reg(base_phase_deg[order[0]], gen_ch=ch1)
            base2 = self.deg2reg(base_phase_deg[order[1]], gen_ch=ch2)

            # Schedule: anchor pulse #2 relative to meas_time, pulse #1 by tau_us
            tau_cyc   = self.us2cycles(cfg.get("tau_us", 0.0), gen_ch=ch1)
            phase2reg = self.deg2reg(cfg.get("phase2_deg", 0.0), gen_ch=ch2)

            L1 = self.len_flat[order[0]]
            L2 = self.len_flat[order[1]]

            # Default: place pulse #2 close to measurement
            t2_end   = meas_time - self.us2cycles(cfg["qub_delay_fixed"], gen_ch=ch2)
            t2_start = t2_end - L2
            t1_start = t2_start - tau_cyc  

            # Program pulses
            self.set_pulse_registers(ch=ch1, style="flat_top", waveform=wf1,
                                     length=ht1, phase=base1, gain=g1)
            

            self.set_pulse_registers(ch=ch2, style="flat_top", waveform=wf2,
                                     length=ht2, phase=base2 + phase2reg, gain=g2)
            
            # --- always reset phase at the start (both channels) ---
            self.reset_phase(gen_ch=[ch1, ch2], t=0)
            
            self.pulse(ch=ch1, t=t1_start)
            self.pulse(ch=ch2, t=t2_start)

        else:
            raise ValueError(f"Unknown pulse_scheme: {scheme}")

        # ----- Measurement -----
        self.trigger(adcs=self.ro_chs, pins=[0], adc_trig_offset=offset)
        self.pulse(ch=cav_ch, t=meas_time)
        self.wait_all()
        self.sync_all(self.us2cycles(cfg["relax_delay"]))
        

def run_single_shot_iq(soc, soccfg, base_config, nshots, state="ground"):
    """
    Collect nshots single-shot IQ points with minimal per-shot overhead.
      - state='ground': pulse_scheme='none' (no qubit pulses)
      - state='excited': use base_config['pulse_scheme'] ('one' or 'two') to prepare |e>
    Speed tips:
      * load_pulses once (first acquire), then set load_pulses=False
      * no plotting in the loop
    Returns: I (N,), Q (N,)
    """
    cfg = dict(base_config)  # shallow copy
    # Single-shot settings
    cfg["reps"]       = 1      # no internal averaging
    cfg["soft_avgs"]  = 1

    # Choose scheme based on state
    if state == "ground":
        cfg["pulse_scheme"] = "none"
    elif state == "excited":
        # Use whatever was requested in base_config ('one' or 'two')
        cfg["pulse_scheme"] = base_config.get("pulse_scheme", "one")
    else:
        raise ValueError("state must be 'ground' or 'excited'")

    I = np.zeros(nshots); Q = np.zeros(nshots)

    t_start = time.time()
    prog = IQReadoutProgram(soccfg, cfg)

    # 1st shot (loads pulses)
    holder = prog.acquire(soc, load_pulses=True, progress=False)
    I[0] = holder[0][0][0]; Q[0] = holder[1][0][0]

    # After a few shots, report a time estimate
    warmup = 100  # use ~100 shots for a decent per-shot estimate
    for k in range(1, nshots):
        holder = prog.acquire(soc, load_pulses=False, progress=False)
        I[k] = holder[0][0][0]; Q[k] = holder[1][0][0]

        if k == warmup:
            t_now = time.time()
            print('For ' + state + ':')
            estimate_time(t_start, t_now, int(nshots/k))
    return I, Q


def run_averaged_iq(soc, soccfg, base_config, navg, state="ground"):
    """
    Do one averaged readout for the given state by setting reps=navg.
    Returns: (I_avg, Q_avg)
    """
    cfg = dict(base_config)
    cfg["reps"]      = int(navg)  # average in-FPGA across reps
    cfg["soft_avgs"] = 1          # keep this 1 to avoid extra Python averaging

    if state == "ground":
        cfg["pulse_scheme"] = "none"
    elif state == "excited":
        cfg["pulse_scheme"] = base_config.get("pulse_scheme", "one")
    else:
        raise ValueError("state must be 'ground' or 'excited'")

    prog = IQReadoutProgram(soccfg, cfg)
    holder = prog.acquire(soc, load_pulses=True, progress=False)
    I = holder[0][0][0]; Q = holder[1][0][0]
    return float(I), float(Q)

def get_iq_readout_settings():
    s = {}
    # ---- legacy/global pieces you already have ----
    s['scanname']       = 'IQ_readout'
    s['meas_type']      = 'IQ_readout'
    s['phase_reset']    = True
    s['cav_freq']       = 6e9
    s['meas_gain']      = 1000

    # Qubit D1
    s['qub_freq_D1']    = 4e9
    s['qub_gain_D1']    = 1000
    s['qub_phase_D1']   = 0.0
    s['qub_sigma_D1']   = 0.02e-6
    s['num_sigma_D1']   = 4
    s['hold_time_D1_us']= 1.0

    # Qubit D2
    s['qub_freq_D2']    = 4e9
    s['qub_gain_D2']    = 1000
    s['qub_phase_D2']   = 0.0
    s['qub_sigma_D2']   = 0.02e-6
    s['num_sigma_D2']   = 4
    s['hold_time_D2_us']= 1.0

    # Measurement timing
    s['readout_length'] = 4.0    # us (same as meas_window typically)
    s['meas_window']    = 4.0    # us
    s['meas_time']      = 2.0    # us
    s['adc_trig_offset']= 0.2    # us after trigger
    s['relax_delay']    = 200.0  # us between shots
    s['qub_delay_fixed']= 0.2    # us: delay from qubit pulse #2 (or single pulse) to readout

    # -------- new high-level switches --------
    s['meas_mode']      = 'single'      # 'single' or 'averaged'
    s['n_points']       = 2000          # for 'single' mode only
    s['navg']           = 5000         # for 'averaged' mode only

    s['pulse_scheme']   = 'one'         # 'one' or 'two' for EXCITED state
    s['single_drive']   = 'D1'          # if pulse_scheme == 'one': 'D1' or 'D2'

    # Two-pulse controls (used when pulse_scheme == 'two')
    s['pulse_order']    = ('D1','D2')
    s['tau_us']         = 0.0
    s['phase2_deg']     = 0.0

    # Card/QICK
    s['nqz_c']          = 1
    s['nqz_q']          = 2
    s['ro_channels']    = [0]         # adapt to your hw
    s['cav_channel']    = 0
    s['qub_channel_D1'] = 2
    s['qub_channel_D2'] = 3
    s['cav_phase']      = 0.0

    return s



def measure_iq_states(soc, soccfg, instruments, settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    # Build QICK config from your globals/settings (same approach as before)
    cfg = {
        'cav_channel'     : exp_globals['cav_channel'],
        'qub_channel_D1'  : exp_globals['qub_channel_D1'],
        'qub_channel_D2'  : exp_globals['qub_channel_D2'],
        'ro_channels'     : exp_globals['ro_channels'],
        'nqz_c'           : 1,
        'nqz_q'           : 2,

        'cav_phase'       : exp_globals['measurement_pulse']['cav_phase'],
        'meas_window'     : exp_globals['measurement_pulse']['meas_window'],
        'meas_time'       : exp_globals['measurement_pulse']['meas_pos'],
        'meas_gain'       : exp_settings['meas_gain'],
        'cav_freq'        : (exp_settings['cav_freq']-exp_globals['LO_freq']*exp_globals['LO'])/1e6,

        # D1
        'qub_phase_D1'    : exp_globals['qubit_pulse_D1']['qub_phase'],
        'qub_freq_D1'     : (exp_settings['qub_freq_D1'])/1e6,
        'qub_gain_D1'     : exp_settings['qub_gain_D1'],
        'qub_sigma_D1'    : exp_globals['qubit_pulse_D1']['sigma'],
        'num_sigma_D1'    : exp_globals['qubit_pulse_D1']['num_sigma'],
        'hold_time_D1_us' : exp_globals['qubit_pulse_D1']['hold_time'],

        # D2
        'qub_phase_D2'    : exp_globals['qubit_pulse_D2']['qub_phase'],
        'qub_freq_D2'     : (exp_settings['qub_freq_D2'])/1e6,
        'qub_gain_D2'     : exp_settings['qub_gain_D2'],
        'qub_sigma_D2'    : exp_globals['qubit_pulse_D2']['sigma'],
        'num_sigma_D2'    : exp_globals['qubit_pulse_D2']['num_sigma'],
        'hold_time_D2_us' : exp_globals['qubit_pulse_D2']['hold_time'],

        'readout_length'  : exp_globals['measurement_pulse']['meas_window'],
        'adc_trig_offset' : exp_globals['measurement_pulse']['emp_delay'] + exp_globals['measurement_pulse']['meas_pos'],

        # high-level choices
        'pulse_scheme'    : exp_settings['pulse_scheme'],    # 'one' or 'two' (for EXCITED only)
        'single_drive'    : exp_settings.get('single_drive','D1'),
        'pulse_order'     : exp_settings.get('pulse_order',('D1','D2')),
        'tau_us'          : exp_settings.get('tau_us', 0.0),
        'phase2_deg'      : exp_settings.get('phase2_deg', 0.0),

        'qub_delay_fixed' : exp_globals['qub_delay_fixed'],
        'relax_delay'     : exp_globals['relax_delay'],
    }

    mode    = exp_settings['meas_mode']  # 'single' or 'averaged'
    nshots  = int(exp_settings.get('n_points', 1000))
    navg    = int(exp_settings.get('navg', 5000))

    t0 = time.time()
    if mode == "single":
        # Ground first (no pulses), then excited (with your chosen scheme)
        Ig, Qg = run_single_shot_iq(soc, soccfg, cfg, nshots, state="ground")
        Ie, Qe = run_single_shot_iq(soc, soccfg, cfg, nshots, state="excited")

        # Quick scatter (single figure; replaceable by fig number if you prefer)
        fig = plt.figure(1, figsize=(5,5))
        plt.clf()
        plt.scatter(Ig, Qg, s=6, label='|g>')
        plt.scatter(Ie, Qe, s=6, label='|e>')
        plt.xlabel('I'); plt.ylabel('Q'); plt.legend(); plt.axis('equal')
        plt.suptitle('Single-shot IQ\n'+filename)
        plt.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_singleshot_scatter.png'), dpi=150)

        userfuncs.SaveFull(saveDir, filename,
            ['Ig','Qg','Ie','Qe','mode','nshots','cfg'],
            locals(), expsettings=settings, instruments=instruments)

    elif mode == "averaged":
        Ig, Qg = run_averaged_iq(soc, soccfg, cfg, navg, state="ground")
        Ie, Qe = run_averaged_iq(soc, soccfg, cfg, navg, state="excited")
        
        g_mag = np.sqrt(Ig**2 + Qg**2)
        e_mag = np.sqrt(Ie**2 + Qe**2)

        fig = plt.figure(2, figsize=(5,5))
        plt.clf()
        plt.scatter([Ig],[Qg], s=60, label=f'|g>, mag={g_mag:.2f}')
        plt.scatter([Ie],[Qe], s=60, label=f'|e>, mag={e_mag:.2f}')
        plt.xlabel('I'); plt.ylabel('Q'); plt.legend(); plt.axis('equal')
        plt.suptitle(f'Averaged IQ (reps={navg})\n'+filename)
        plt.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_averaged_points.png'), dpi=150)

        userfuncs.SaveFull(saveDir, filename,
            ['Ig','Qg','Ie','Qe','mode','navg','cfg'],
            locals(), expsettings=settings, instruments=instruments)
    else:
        raise ValueError("meas_mode must be 'single' or 'averaged'")

    print(f"Elapsed time: {time.time()-t0:.2f} s")




    


   
    

        