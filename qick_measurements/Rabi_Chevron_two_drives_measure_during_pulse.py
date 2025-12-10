# -*- coding: utf-8 -*-
"""
Rabi Chevron (two Gaussian-square pulses on D1 & D2, shared hold sweep + common detuning)
- D1/D2 pulses: Gaussian-square via add_gauss + flat_top; sigma/num_sigma/max_amp (gain) from exp_globals
- Configurable start-to-start interval (tau_us) and extra phase on pulse #2 (phase2_deg)
- 2D sweep: detuning (MHz) x hold_time (us)
- Live colormap refresh after each hold row
- Total-time estimate printed after finishing the first hold row
- Timestamp included in saved figure file name and title
"""

from qick.averager_program import AveragerProgram

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
import time

from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import estimate_time


class RabiChevronSequence(AveragerProgram):
    def initialize(self):
        cfg=self.cfg
        cav_ch     = cfg["cav_channel"]
        qub_ch_D1  = cfg["qub_channel_D1"]
        qub_ch_D2  = cfg["qub_channel_D2"]

        # Nyquist zones
        self.declare_gen(ch=cav_ch,    nqz=cfg["nqz_c"])
        self.declare_gen(ch=qub_ch_D1, nqz=cfg["nqz_q"])
        self.declare_gen(ch=qub_ch_D2, nqz=cfg["nqz_q"])

        # Readout declaration
        readout = self.us2cycles(cfg["readout_length"], ro_ch=cfg["ro_channels"][0])
        for ch in cfg["ro_channels"]:
            self.declare_readout(ch=ch, length=readout,
                                 freq=cfg["cav_freq"], gen_ch=cav_ch)

        # Cavity defaults
        freq_c  = self.freq2reg(cfg["cav_freq"], gen_ch=cav_ch, ro_ch=cfg["ro_channels"][0])
        phase_c = self.deg2reg(cfg["cav_phase"], gen_ch=cav_ch)
        gain_c  = cfg["meas_gain"]
        self.default_pulse_registers(ch=cav_ch, freq=freq_c, phase=phase_c, gain=gain_c)
        self.set_pulse_registers(ch=cav_ch, style="const",
                                 length=self.us2cycles(cfg["meas_window"], gen_ch=cav_ch))

        # Drive frequencies (apply detuning here)
        # NOTE: detuning_MHz is added to BOTH D1 and D2. If you want it on only one drive, change just that one line.
        det = cfg.get("detuning_MHz", 0.0)
        freq_q_D1 = self.freq2reg(cfg["qub_freq_D1"] + det, gen_ch=qub_ch_D1)
        freq_q_D2 = self.freq2reg(cfg["qub_freq_D2"] + det, gen_ch=qub_ch_D2)

        gain_q_D1 = cfg["qub_gain_D1"]
        gain_q_D2 = cfg["qub_gain_D2"]
        self.default_pulse_registers(ch=qub_ch_D1, freq=freq_q_D1, gain=gain_q_D1)
        self.default_pulse_registers(ch=qub_ch_D2, freq=freq_q_D2, gain=gain_q_D2)

        # Gaussian-square definitions (edges from sigma/num_sigma; flat-top 'hold' supplied per-sweep in body)
        sigma_D1     = self.us2cycles(cfg["qub_sigma_D1"], gen_ch=qub_ch_D1)
        sigma_D2     = self.us2cycles(cfg["qub_sigma_D2"], gen_ch=qub_ch_D2)
        num_sigma_D1 = cfg["num_sigma_D1"]
        num_sigma_D2 = cfg["num_sigma_D2"]

        self.add_gauss(ch=qub_ch_D1, name="ex_D1",
                       sigma=sigma_D1,
                       length=int(sigma_D1*num_sigma_D1),
                       maxv=int(self.soccfg['gens'][qub_ch_D1]['maxv']))

        self.add_gauss(ch=qub_ch_D2, name="ex_D2",
                       sigma=sigma_D2,
                       length=int(sigma_D2*num_sigma_D2),
                       maxv=int(self.soccfg['gens'][qub_ch_D2]['maxv']))

        # Base phases (deg -> reg)
        self.base_phase_reg_D1 = self.deg2reg(cfg["qub_phase_D1"], gen_ch=qub_ch_D1)
        self.base_phase_reg_D2 = self.deg2reg(cfg["qub_phase_D2"], gen_ch=qub_ch_D2)

        # Measurement timing constants precomputed
        self.offset    = self.us2cycles(cfg["adc_trig_offset"], gen_ch=cav_ch)
        self.meas_time = self.us2cycles(cfg["meas_time"],       gen_ch=cav_ch)

        # Sync room
        self.synci(200)

    def body(self):
        cfg = self.cfg
        cav_ch     = cfg["cav_channel"]
        qub_ch_D1  = cfg["qub_channel_D1"]
        qub_ch_D2  = cfg["qub_channel_D2"]

        # Shared hold for both pulses (us -> cycles)
        hold_cyc = self.us2cycles(cfg["hold_time_us"], gen_ch=qub_ch_D1)  # same for both (assume same DAC rate)
        # Set flat-top lengths (edges are baked in the arb)
        self.set_pulse_registers(ch=qub_ch_D1, style="flat_top", waveform="ex_D1", length=hold_cyc, phase=self.base_phase_reg_D1)
        self.set_pulse_registers(ch=qub_ch_D2, style="flat_top", waveform="ex_D2", length=hold_cyc, phase=self.base_phase_reg_D2)

        # Start-to-start interval (can be negative). We'll schedule pulses by absolute start times.
        tau_us  = float(cfg.get("tau_us", 0.0))
        tau_cyc = self.us2cycles(tau_us, gen_ch=qub_ch_D1)

        # Optional small gap you sometimes use
        gap_us  = float(cfg.get("gap", 0.0))
        gap_cyc = self.us2cycles(gap_us, gen_ch=qub_ch_D1)

        # Pulse total length (edges + hold)
        sigma_D1     = self.us2cycles(cfg["qub_sigma_D1"], gen_ch=qub_ch_D1)
        sigma_D2     = self.us2cycles(cfg["qub_sigma_D2"], gen_ch=qub_ch_D2)
        num_sigma_D1 = cfg["num_sigma_D1"]
        num_sigma_D2 = cfg["num_sigma_D2"]
        base_phase_deg = {"D1": cfg["qub_phase_D1"], "D2": cfg["qub_phase_D2"]}
        L1 = int(sigma_D1*num_sigma_D1 + hold_cyc)
        L2 = int(sigma_D2*num_sigma_D2 + hold_cyc)

        # Which channel is pulse 1 / pulse 2
        order = cfg["pulse_order"]  # e.g., ('D1','D2'), ('D2','D1'), ('D1','D1'), ...
        ch_map = {"D1": qub_ch_D1, "D2": qub_ch_D2}
        wf_map = {"D1": "ex_D1",   "D2": "ex_D2"}
        len_map = {"D1": L1, "D2": L2}
        ch1 = ch_map[order[0]]
        ch2 = ch_map[order[1]]
        wf1 = wf_map[order[0]]
        wf2 = wf_map[order[1]]
        Lp1 = len_map[order[0]]
        Lp2 = len_map[order[1]]
        
        # convert base phases to regs for the *appropriate* channels
        base1_reg = self.deg2reg(base_phase_deg[order[0]], gen_ch=ch1)
        base2_reg = self.deg2reg(base_phase_deg[order[1]], gen_ch=ch2)

        # Anchor the END of pulse #2 before meas_time, and back-solve for starts with tau
        ends2_later = (tau_cyc + (Lp2 - Lp1)) >= 0
        if ends2_later:
            t2_end   = self.meas_time + self.us2cycles(30, gen_ch=ch2)
            t2_start = t2_end - Lp2
            t1_start = t2_start - tau_cyc - gap_cyc
        else:
            t1_end   = self.meas_time + self.us2cycles(30, gen_ch=ch1)
            t1_start = t1_end - Lp1
            t2_start = t1_start + tau_cyc + gap_cyc

        if min(t1_start, t2_start) < 0:
            raise ValueError("Chosen tau/gap/hold made a pulse start < 0. Reduce |tau| or hold_time_us.")

        # Optional phase reset behavior
        if cfg.get("phase_reset", False):
            self.reset_phase(gen_ch=[ch1, ch2], t=0)

        # Extra phase for pulse #2
        phase2_deg = float(cfg.get("phase2_deg", 0.0))
        phase2_reg = self.deg2reg(phase2_deg, gen_ch=ch2)

        # Trigger ADC
        self.trigger(adcs=self.ro_chs, pins=[0], adc_trig_offset=self.offset)

        # Pulse #1
        # (re-apply style/length to be safe if other calls modify later)
        self.set_pulse_registers(ch=ch1, style="flat_top", waveform=wf1, length=hold_cyc, phase=base1_reg)
        self.pulse(ch=ch1, t=t1_start)

        # Pulse #2 with extra phase
        self.set_pulse_registers(ch=ch2, style="flat_top", waveform=wf2, length=hold_cyc, phase=base2_reg + phase2_reg)
        self.pulse(ch=ch2, t=t2_start)

        # Measurement pulse
        self.pulse(ch=cav_ch, t=self.meas_time)
        self.wait_all()
        self.sync_all(self.us2cycles(cfg["relax_delay"]))


def get_rabi_chevron_settings():

    settings = {}

    settings['scanname']      = 'Rabi_Chevron_two_qubit_channels'
    settings['meas_type']     = 'Rabi_Chevron'
    settings['phase_reset']   = False

    # Cavity
    settings['cav_freq']      = 6e9
    settings['meas_gain']     = 1000

    # D1
    settings['qub_freq_D1']   = 4e9
    settings['qub_gain_D1']   = 1000
    settings['qub_phase_D1']  = 0.0
    settings['qub_sigma_D1']  = 0.02e-6
    settings['num_sigma_D1']  = 4

    # D2
    settings['qub_freq_D2']   = 4e9
    settings['qub_gain_D2']   = 1000
    settings['qub_phase_D2']  = 0.0
    settings['qub_sigma_D2']  = 0.02e-6
    settings['num_sigma_D2']  = 4

    # Fixed scheduling knobs
    settings['pulse_order']     = ('D1','D2')
    settings['tau_us']          = 0.0         # start-to-start interval (can be negative)
    settings['phase2_deg']      = 0.0         # extra phase on pulse #2
    settings['gap']             = 0.0         # optional extra physical gap between starts
    settings['qub_delay_fixed'] = 0.0         # minimum q_pulse->meas end margin

    # 2D sweep definition
    # X: detuning_MHz relative to qub_freq_{D1,D2}
    settings['det_min_MHz']   = -5
    settings['det_max_MHz']   = +5
    settings['det_points']    = 81

    # Y: shared flat-top hold (us)
    settings['hold_min_us']   = 0.02
    settings['hold_max_us']   = 1.00
    settings['hold_points']   = 51

    # Acquisition
    settings['subtract_background'] = False
    settings['reps'] = 1
    settings['soft_avgs'] = 2e3

    return settings


def meas_rabi_chevron(soc, soccfg, instruments, settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse_D1   = exp_globals['qubit_pulse_D1']
    q_pulse_D2   = exp_globals['qubit_pulse_D2']

    # (Optional) LO handling if you actually use LO mixing in your rack
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

    # Base config (MHz for gen regs like your style)
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
        'cav_freq'        : (exp_settings['cav_freq'] - exp_globals['LO_freq']*exp_globals['LO'])/1e6, # MHz

        'nqz_q'           : 2,

        # D1 drives
        'qub_phase_D1'    : q_pulse_D1['qub_phase'],
        'qub_freq_D1'     : (exp_settings['qub_freq_D1'])/1e6, # MHz
        'qub_gain_D1'     : exp_settings['qub_gain_D1'],
        'qub_sigma_D1'    : q_pulse_D1['sigma'],
        'num_sigma_D1'    : q_pulse_D1['num_sigma'],

        # D2 drives
        'qub_phase_D2'    : q_pulse_D2['qub_phase'],
        'qub_freq_D2'     : (exp_settings['qub_freq_D2'])/1e6, # MHz
        'qub_gain_D2'     : exp_settings['qub_gain_D2'],
        'qub_sigma_D2'    : q_pulse_D2['sigma'],
        'num_sigma_D2'    : q_pulse_D2['num_sigma'],

        # Shared hold in µs – will be updated in loop
        'hold_time_us'    : q_pulse_D1['hold_time'],

        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],

        # Scheduling knobs
        'tau_us'          : exp_settings.get('tau_us', 0.0),
        'phase2_deg'      : exp_settings.get('phase2_deg', 0.0),
        'gap'             : exp_settings.get('gap', 0.0),
        'qub_delay_fixed' : exp_globals['qub_delay_fixed'],
        'pulse_order'     : exp_settings['pulse_order'],

        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs'],

        'phase_reset'     : exp_settings['phase_reset'],

        # detuning will be updated in the outer loop
        'detuning_MHz'    : 0.0,
    }

    # Prepare measurement windows
    prog_probe = RabiChevronSequence(soccfg, config)
    meas_start = prog_probe.us2cycles(m_pulse["init_buffer"], ro_ch=0)
    meas_end   = meas_start + prog_probe.us2cycles(m_pulse["meas_window"], ro_ch=0)
    total_samples = prog_probe.us2cycles(config['readout_length'], ro_ch=0)

    # Build sweep axes
    det_list   = np.linspace(exp_settings['det_min_MHz'],  exp_settings['det_max_MHz'],  exp_settings['det_points'])
    hold_list  = np.linspace(exp_settings['hold_min_us'], exp_settings['hold_max_us'], exp_settings['hold_points'])

    Nx = len(det_list)   # columns (x-axis)
    Ny = len(hold_list)  # rows    (y-axis)

    # Allocate results (Ny x Nx)
    amp_map  = np.zeros((Ny, Nx))
    ang_map  = np.zeros((Ny, Nx))
    I_map    = np.zeros((Ny, Nx))
    Q_map    = np.zeros((Ny, Nx))

    # Optional background
    if exp_settings['subtract_background']:
        print('Starting Background Trace')
        bprog = RabiChevronSequence(soccfg, config)
        holder = bprog.acquire(soc, load_pulses=True, progress=False)
        I_back = holder[0][0][0]
        Q_back = holder[1][0][0]
        print('Background Trace Complete')
    else:
        I_back, Q_back = 0.0, 0.0

    # Live plot figure
    fig = plt.figure(1, figsize=(10, 7))
    plt.clf()
    plt.suptitle('Rabi Chevron (live)')
    first_row_done = False
    tstart_global = time.time()

    for iy, hold_us in enumerate(hold_list):
        row_start = time.time()

        # set shared hold for this row
        config['hold_time_us'] = float(hold_us)

        for ix, det_MHz in enumerate(det_list):
            config['detuning_MHz'] = float(det_MHz)

            # Run program at this point
            prog = RabiChevronSequence(soccfg, config)
            holder = prog.acquire(soc, load_pulses=True, progress=False)
            I_sig = holder[0][0][0]
            Q_sig = holder[1][0][0]

            I = I_sig - I_back
            Q = Q_sig - Q_back

            I_map[iy, ix] = I
            Q_map[iy, ix] = Q
            amp_map[iy, ix] = np.sqrt(I**2 + Q**2)
            ang_map[iy, ix] = np.degrees(np.arctan2(Q, I))

        # After finishing this hold row: live colormap refresh
        plt.clf()
        ax = plt.gca()
        # Use your helper if you prefer; here is a direct imshow:
        # extent: x in MHz, y in us
        im = ax.imshow(amp_map[:iy+1, :],
                       origin='lower',
                       aspect='auto',
                       extent=[det_list[0], det_list[-1], hold_list[0], hold_list[iy]])
        plt.colorbar(im, ax=ax, label='Amplitude (arb)')
        ax.set_xlabel('Detuning (MHz)')
        ax.set_ylabel('Hold (µs)')
        stamp = userfuncs.timestamp()
        ax.set_title('Rabi Chevron (live) \n' + filename)
        fig.canvas.draw(); fig.canvas.flush_events()

        # Total-time estimate AFTER the first hold row completes
        if not first_row_done:
            row_stop = time.time()
            # Estimate using rows (Ny total rows)
            estimate_time(row_start, row_stop, Ny)
            first_row_done = True

        # Save rolling progress (optional)
        plt.savefig(os.path.join(saveDir, filename+'_live.png'), dpi=150)
        
        # Save data
        userfuncs.SaveFull(
            saveDir, filename,
            ['det_list', 'hold_list', 'amp_map', 'ang_map', 'I_map', 'Q_map'],
            locals(),
            expsettings=settings,
            instruments=instruments
        )

    tstop_global = time.time()
    print('Elapsed time (s): {:.3f}'.format(tstop_global - tstart_global))

# =============================================================================
#     # Final data saving
#     plt.figure(2, figsize=(11, 7))
#     stamp = userfuncs.timestamp()
#     ttl = f"Rabi Chevron | {stamp}"
#     general_colormap_subplot(
#         np.array(det_list), np.array(hold_list),
#         amp_map,
#         title=ttl,
#         xlabel='Detuning (MHz)',
#         ylabel='Hold (µs)',
#         cbar_label='Amplitude (arb.)'
#     )
#     plt.tight_layout()
#     plt.savefig(os.path.join(saveDir, f"{filename}_Chevron_{stamp}.png"), dpi=200)
# =============================================================================

    # Save data
    userfuncs.SaveFull(
        saveDir, filename,
        ['det_list', 'hold_list', 'amp_map', 'ang_map', 'I_map', 'Q_map'],
        locals(),
        expsettings=settings,
        instruments=instruments
    )

    return prog
