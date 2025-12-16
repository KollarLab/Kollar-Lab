# -*- coding: utf-8 -*-

from qick.averager_program import AveragerProgram

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
import time

from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import estimate_time


class MeasurementSequence(AveragerProgram):
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

        # --- NEW: detuning-aware qubit frequencies ---
        # cfg["qub_freq_D1"], cfg["qub_freq_D2"] are in MHz
        det_MHz = float(cfg.get("detuning_MHz", 0.0))
        freq_q_D1 = self.freq2reg(cfg["qub_freq_D1"] + det_MHz, gen_ch=qub_ch_D1)
        freq_q_D2 = self.freq2reg(cfg["qub_freq_D2"] + det_MHz, gen_ch=qub_ch_D2)

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
            t2_end   = self.meas_time - self.us2cycles(cfg["qub_delay_fixed"], gen_ch=ch2)
            t2_start = t2_end - Lp2
            t1_start = t2_start - tau_cyc
        else:
            t1_end   = self.meas_time - self.us2cycles(cfg["qub_delay_fixed"], gen_ch=ch1)
            t1_start = t1_end - Lp1
            t2_start = t1_start + tau_cyc

        if min(t1_start, t2_start) < 0:
            raise ValueError("Chosen tau/hold made a pulse start < 0. Reduce |tau| or hold_time_us.")

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


#------------------------------------------------------------------------------------
def fm_dependent_phase_offsets(fm):
    """
    Compute phi11_minus_12_deg and phi22_minus_21_deg as functions of modulation frequency fm.

    Relations:
        phi11_minus_12_deg = -1.229  + 360 * (-0.875e-9) * fm
        phi22_minus_21_deg = -5.274  + 360 * ( 0.910e-9) * fm
    """
    phi11_minus_12_deg = -1.229 + 360.0 * (-0.875e-9) * fm
    phi22_minus_21_deg = -5.274 + 360.0 * ( 0.910e-9) * fm
    return phi11_minus_12_deg, phi22_minus_21_deg


def ac_modulation_calculator(
    dbeta1_dV1, phi11_minus_12_deg, dbeta1_dV2,
    dbeta2_dV1, dbeta2_dV2, phi22_minus_21_deg,
    beta1, phi1_prime_deg, beta2, phi2_prime_deg,
):
    """
    Solve for Vpp_1, phi_V1, Vpp_2, phi_V2 ...
    """
    d2r = np.pi / 180.0
    phi11_12 = phi11_minus_12_deg * d2r
    phi22_21 = phi22_minus_21_deg * d2r
    phi1p = phi1_prime_deg * d2r
    phi2p = phi2_prime_deg * d2r

    A = np.array([
        [dbeta1_dV1 * np.exp(1j * phi11_12),      dbeta1_dV2],
        [dbeta2_dV1,                              dbeta2_dV2 * np.exp(1j * phi22_21)]
    ], dtype=complex)

    b = np.array([
        beta1 * np.exp(1j * phi1p),
        beta2 * np.exp(1j * phi2p)
    ], dtype=complex)

    V1_complex, V2_complex = np.linalg.solve(A, b)

    Vpp1 = np.abs(V1_complex)
    Vpp2 = np.abs(V2_complex)
    phi_V1_deg = np.degrees(np.angle(V1_complex)) % 360.0
    phi_V2_deg = np.degrees(np.angle(V2_complex)) % 360.0

    return Vpp1, phi_V1_deg, Vpp2, phi_V2_deg

#-------------------------------------------------------------------------------------
def get_phi_and_detuning_settings():
    """
    Default exp_settings for phi2_prime vs detuning scan
    (two Gaussian-square pulses on D1 & D2 + AC flux modulation).

    These are *experiment-level* knobs; exp_globals still provides:
      - measurement_pulse
      - qubit_pulse_D1 / qubit_pulse_D2 (including hold_time)
      - cav/qubit channels, ro_channels, LO settings, relax_delay, etc.
    """

    settings = {}

    # Bookkeeping / labels
    settings['scanname']    = 'Search_for_EP_and_gyration_phi_and_spec_scan'
    settings['meas_type']   = 'phi_and_spec_scan'
    settings['phase_reset'] = True   # reset DDS phase accumulators at t=0 each shot

    # -------- Cavity --------
    settings['cav_freq']  = 6.0e9     # Hz, can override in driver
    settings['meas_gain'] = 3000

    # -------- Qubit drives (D1 / D2) --------
    # These are the *nominal* qubit frequencies (in Hz) and gains; detuning is added on top.
    settings['qub_freq_D1']  = 4.5e9
    settings['qub_gain_D1']  = 0        # e.g. turn D1 off by default

    settings['qub_freq_D2']  = 4.5e9
    settings['qub_gain_D2']  = 500      # example gain, override in driver

    # The pulse *shapes* (sigma, num_sigma) come from exp_globals['qubit_pulse_D1/D2'],
    # so we only need the phases here:
    settings['qub_phase_D1'] = 0.0
    settings['qub_phase_D2'] = 0.0

    # -------- Scheduling between the two qubit pulses --------
    # pulse_order: which channel is "pulse 1" and "pulse 2"
    settings['pulse_order'] = ('D2', 'D1')  # e.g. ('D1','D2'), ('D1','D1'), etc.
    settings['tau_us']      = 0.0           # start-to-start interval in µs (can be negative)
    settings['phase2_deg']  = 0.0           # extra phase on pulse #2 (deg)

    # -------- AC modulation calibration + targets --------
    # Slopes dβ_j/dV_k (units: your calibration convention)
    settings['dbeta1_dV1'] = 0.0   # TODO: fill from calibration
    settings['dbeta1_dV2'] = 0.0   # TODO
    settings['dbeta2_dV1'] = 0.0   # TODO
    settings['dbeta2_dV2'] = 0.0   # TODO

    # Target modulation amplitudes β1, β2  (Hz or MHz – whatever your convention is)
    settings['beta1'] = 0.0        # TODO
    settings['beta2'] = 0.0        # TODO

    # DC offsets on the two modulation channels (V)
    settings['dc_offset_voltage_1'] = 0.0   # TODO
    settings['dc_offset_voltage_2'] = 0.0   # TODO

    # -------- Fixed fm for this scan --------
    # Single modulation frequency in Hz (used for all rows)
    settings['fm_fixed'] = 1.0e6   # 1 MHz default; override in driver

    # -------- phi2' sweep (rows, y-axis) --------
    # We sweep φ2′ (deg). φ1′ is kept at 0 in the calculator.
    settings['phi2_prime_start']  = -180.0
    settings['phi2_prime_end']    = +180.0
    settings['phi2_prime_points'] = 73

    # -------- Detuning sweep (columns, x-axis) --------
    # Detuning in MHz; actual qubit frequency = qub_freq_D1 + detuning
    settings['det_min_MHz'] = -25.0
    settings['det_max_MHz'] = +25.0
    settings['det_points']  = 51

    # -------- Acquisition --------
    settings['reps']                = 5000
    settings['soft_avgs']           = 1
    settings['subtract_background'] = False

    # Optional debug fields (ignored by core measurement; kept for consistency)
    settings['debug']      = False
    settings['debug_time'] = 5

    return settings



def phi_and_detuning_scan(soc, soccfg, instruments, settings):
    """
    2D scan over phi2_prime_deg and qubit detuning.
    - y-axis: phi2_prime_deg (deg)
    - x-axis: qubit pulse frequency (GHz) = qub_freq_D1_nominal + detuning

    Uses:
      - AC modulation calculator (Vpp1, phi_V1, Vpp2, phi_V2 as a function of φ2′)
      - detuning-aware MeasurementSequence (via cfg['detuning_MHz'])
      - A fixed modulation frequency fm = exp_settings['fm_fixed']

    Time estimates:
      - No per-column prints.
      - After the FIRST row, print estimated TOTAL run time.
      - After EACH subsequent row, update estimated REMAINING time.
    """
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse_D1   = exp_globals['qubit_pulse_D1']
    q_pulse_D2   = exp_globals['qubit_pulse_D2']
    
    Dual_gen = instruments['DCsupply']

    # (Optional) LO handling
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

    # Base config (MHz for gen regs)
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

        # Shared hold in µs
        'hold_time_us'    : q_pulse_D1['hold_time'],

        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],

        # Scheduling knobs
        'tau_us'          : exp_settings.get('tau_us', 0.0),
        'phase2_deg'      : exp_settings.get('phase2_deg', 0.0),
        'qub_delay_fixed' : exp_globals['qub_delay_fixed'],
        'pulse_order'     : exp_settings['pulse_order'],

        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs'],

        'phase_reset'     : exp_settings['phase_reset'],

        # detuning (MHz) will be updated in the inner loop
        'detuning_MHz'    : 0.0,
    }

    # Dummy program to validate timing
    prog_probe = MeasurementSequence(soccfg, config)
    total_samples = prog_probe.us2cycles(config['readout_length'], ro_ch=0)

    # --- Build sweep axes ---

    # phi2' (deg) – rows (y-axis)
    phi2_prime_list = np.linspace(
        exp_settings['phi2_prime_start'],
        exp_settings['phi2_prime_end'],
        exp_settings['phi2_prime_points']
    )

    # detuning (MHz) – columns (x-axis)
    det_list = np.linspace(
        exp_settings['det_min_MHz'],
        exp_settings['det_max_MHz'],
        exp_settings['det_points']
    )

    Nx = len(det_list)        # columns
    Ny = len(phi2_prime_list) # rows

    # x-axis: qubit frequency (GHz) = f_qubit_nominal + detuning
    qub_freq_D1_nom_Hz = exp_settings['qub_freq_D1']  # Hz
    freq_qubit_list_GHz = qub_freq_D1_nom_Hz*1e-9 + det_list*1e-3  # GHz

    # Fixed modulation frequency
    fm_fixed = float(exp_settings['fm_fixed'])  # Hz

    # Precompute fm-dependent offsets once
    phi11_minus_12_deg, phi22_minus_21_deg = fm_dependent_phase_offsets(fm_fixed)

    # Allocate results (Ny x Nx)
    amp_map   = np.zeros((Ny, Nx))
    ang_map   = np.zeros((Ny, Nx))
    I_map     = np.zeros((Ny, Nx))
    Q_map     = np.zeros((Ny, Nx))

    Vpp1_map   = np.zeros((Ny, Nx))
    phi_V1_map = np.zeros((Ny, Nx))
    Vpp2_map   = np.zeros((Ny, Nx))
    phi_V2_map = np.zeros((Ny, Nx))

    # Optional background
    if exp_settings['subtract_background']:
        print('Starting Background Trace')
        bprog = MeasurementSequence(soccfg, config)
        holder = bprog.acquire(soc, load_pulses=True, progress=False)
        I_back = holder[0][0][0]
        Q_back = holder[1][0][0]
        print('Background Trace Complete')
    else:
        I_back, Q_back = 0.0, 0.0

    # ------------------------------------------------------------------------------------

    if Dual_gen.Ch1_output == 0:
        Dual_gen.Ch1_output = 1
    if Dual_gen.Ch2_output == 0:
        Dual_gen.Ch2_output = 1

    dbeta1_dV1 = exp_settings['dbeta1_dV1']
    dbeta1_dV2 = exp_settings['dbeta1_dV2']
    dbeta2_dV1 = exp_settings['dbeta2_dV1']
    dbeta2_dV2 = exp_settings['dbeta2_dV2']
    beta1      = exp_settings['beta1']
    beta2      = exp_settings['beta2']

    # φ1′ is fixed at 0; φ2′ is swept
    phi1_prime_deg = 0.0

    # String to annotate beta and fm in plot titles
    title_params = (r"$\beta_1$={:.3f} MHz, $\beta_2$={:.3f} MHz, "
                    r"$f_m$={:.3f} MHz").format(beta1/1e6, beta2/1e6, fm_fixed/1e6)

    # Live plot figures
    fig_amp = plt.figure(1, figsize=(10, 7))
    plt.clf()
    plt.suptitle('phi2_prime and spec frequency scan (Amplitude, live)')

    fig_phase = plt.figure(2, figsize=(10, 7))
    plt.clf()
    plt.suptitle('phi2_prime and spec frequency scan (Phase, live)')

    # --- Time-estimation bookkeeping (Rabi-chevron style) ---
    tstart_global  = time.time()
    first_row_done = False

    for iy, phi2_prime_deg in enumerate(phi2_prime_list):
        row_start = time.time()

        for ix, det_MHz in enumerate(det_list):
            # AC modulation solution for this (fm_fixed, phi2_prime_deg)
            Vpp1, phi_V1_deg, Vpp2, phi_V2_deg = ac_modulation_calculator(
                dbeta1_dV1, phi11_minus_12_deg,
                dbeta1_dV2, dbeta2_dV1, dbeta2_dV2, phi22_minus_21_deg,
                beta1, phi1_prime_deg, beta2, phi2_prime_deg
            )
            Vpp1_map[iy, ix]   = Vpp1
            phi_V1_map[iy, ix] = phi_V1_deg
            Vpp2_map[iy, ix]   = Vpp2
            phi_V2_map[iy, ix] = phi_V2_deg

            # Program the external generators at fixed fm
            Dual_gen.Ch1_sin_gen(Vpp1, fm_fixed, phase=phi_V1_deg,
                                 offset=exp_settings['dc_offset_voltage_1'])
            Dual_gen.Ch2_sin_gen(Vpp2, fm_fixed, phase=phi_V2_deg,
                                 offset=exp_settings['dc_offset_voltage_2'])
            Dual_gen.phase_sync()
            time.sleep(0.1)

            # Set detuning for this point
            config['detuning_MHz'] = float(det_MHz)

            # Run QICK program
            prog = MeasurementSequence(soccfg, config)
            holder = prog.acquire(soc, load_pulses=True, progress=False)
            I_sig = holder[0][0][0]
            Q_sig = holder[1][0][0]

            I = I_sig - I_back
            Q = Q_sig - Q_back

            I_map[iy, ix]   = I
            Q_map[iy, ix]   = Q
            amp_map[iy, ix] = np.sqrt(I**2 + Q**2)
            ang_map[iy, ix] = np.degrees(np.arctan2(Q, I))

        # --------- Live plot refresh: amplitude ---------
        plt.figure(fig_amp.number)
        plt.clf()
        ax_amp = plt.gca()

        if Ny > 1:
            im_amp = ax_amp.imshow(
                amp_map[:iy+1, :],
                origin='lower',
                aspect='auto',
                extent=[
                    freq_qubit_list_GHz[0],
                    freq_qubit_list_GHz[-1],
                    phi2_prime_list[0],
                    phi2_prime_list[iy]
                ]
            )
            plt.colorbar(im_amp, ax=ax_amp, label='Amplitude (arb)')
            ax_amp.set_xlabel('Frequency (GHz)')
            ax_amp.set_ylabel(r"$\phi_2^\prime$ (deg)")
            ax_amp.set_title(f'phi2_prime and spec frequency scan \n{title_params}\n{filename}')
        else:
            phi2_0 = phi2_prime_list[0]
            ax_amp.plot(freq_qubit_list_GHz, amp_map[0, :], marker='o', linestyle='-')
            ax_amp.set_xlabel('Frequency (GHz)')
            ax_amp.set_ylabel('Amplitude (arb)')
            ax_amp.set_title(
                f'Spec frequency scan at phi2_prime = {phi2_0:.1f}°\n{title_params}\n{filename}'
            )
            ax_amp.grid(True)

        fig_amp.canvas.draw()
        fig_amp.canvas.flush_events()

        # --------- Live plot refresh: phase ---------
        plt.figure(fig_phase.number)
        plt.clf()
        ax_phase = plt.gca()

        if Ny > 1:
            im_phase = ax_phase.imshow(
                ang_map[:iy+1, :],
                origin='lower',
                aspect='auto',
                extent=[
                    freq_qubit_list_GHz[0],
                    freq_qubit_list_GHz[-1],
                    phi2_prime_list[0],
                    phi2_prime_list[iy]
                ]
            )
            plt.colorbar(im_phase, ax=ax_phase, label='Phase (deg)')
            ax_phase.set_xlabel('Frequency (GHz)')
            ax_phase.set_ylabel(r"$\phi_2^\prime$ (deg)")
            ax_phase.set_title(f'phi2_prime and spec frequency scan \n{title_params}\n{filename}')
        else:
            phi2_0 = phi2_prime_list[0]
            ax_phase.plot(freq_qubit_list_GHz, ang_map[0, :], marker='o', linestyle='-')
            ax_phase.set_xlabel('Frequency (GHz)')
            ax_phase.set_ylabel('Phase (deg)')
            ax_phase.set_title(
                f'Spec frequency scan at phi2_prime = {phi2_0:.1f}°\n{title_params}\n{filename}'
            )
            ax_phase.grid(True)

        fig_phase.canvas.draw()
        fig_phase.canvas.flush_events()

        # Save rolling progress
        plt.figure(fig_amp.number)
        plt.savefig(os.path.join(saveDir, filename+'_live_amp.png'), dpi=150)
        plt.figure(fig_phase.number)
        plt.savefig(os.path.join(saveDir, filename+'_live_phase.png'), dpi=150)

        # Save data after each phi2' row
        userfuncs.SaveFull(
            saveDir, filename,
            [
                'phi2_prime_list', 'det_list', 'freq_qubit_list_GHz',
                'amp_map', 'ang_map', 'I_map', 'Q_map',
                'Vpp1_map', 'phi_V1_map', 'Vpp2_map', 'phi_V2_map'
            ],
            locals(),
            expsettings=settings,
            instruments=instruments
        )

        # --------- Time estimate after this row (Rabi-chevron style) ---------
        row_stop = time.time()
        if not first_row_done:
            print(f"Finished first phi2_prime row ({iy+1}/{Ny}).")
            print("Estimated TOTAL run time based on this row:")
            estimate_time(row_start, row_stop, Ny)
            first_row_done = True
        else:
            remaining_rows = Ny - (iy + 1)
            if remaining_rows > 0:
                print(f"Finished phi2_prime row {iy+1}/{Ny}.")
                print("Estimated REMAINING time based on this row:")
                estimate_time(row_start, row_stop, remaining_rows)

    tstop_global = time.time()
    print('Elapsed time (s): {:.3f}'.format(tstop_global - tstart_global))

    # Final save
    userfuncs.SaveFull(
        saveDir, filename,
        [
            'phi2_prime_list', 'det_list', 'freq_qubit_list_GHz',
            'amp_map', 'ang_map', 'I_map', 'Q_map',
            'Vpp1_map', 'phi_V1_map', 'Vpp2_map', 'phi_V2_map'
        ],
        locals(),
        expsettings=settings,
        instruments=instruments
    )

    return prog



