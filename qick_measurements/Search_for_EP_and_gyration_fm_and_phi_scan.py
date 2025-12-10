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

        # Drive frequencies
        freq_q_D1 = self.freq2reg(cfg["qub_freq_D1"], gen_ch=qub_ch_D1)
        freq_q_D2 = self.freq2reg(cfg["qub_freq_D2"], gen_ch=qub_ch_D2)

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


def get_measurement_settings():

    settings = {}

    settings['scanname']      = 'Search_for_EP_and_gyration_fm_and_phi_scan'
    settings['meas_type']     = 'fm_and_phi_scan'
    settings['phase_reset']   = True

    # Cavity
    settings['cav_freq']      = 6e9
    settings['meas_gain']     = 1000

    # D1
    settings['qub_freq_D1']   = 4.5e9
    settings['qub_gain_D1']   = 100
    settings['qub_phase_D1']  = 0.0
    settings['qub_sigma_D1']  = 0.02e-6
    settings['num_sigma_D1']  = 4

    # D2
    settings['qub_freq_D2']   = 4.5e9
    settings['qub_gain_D2']   = 100
    settings['qub_phase_D2']  = 0.0
    settings['qub_sigma_D2']  = 0.02e-6
    settings['num_sigma_D2']  = 4

    # Fixed scheduling knobs
    settings['pulse_order']     = ('D2','D1')
    settings['tau_us']          = 0.0         # start-to-start interval (can be negative)
    settings['phase2_deg']      = 0.0         # extra phase on pulse #2
    settings['qub_delay_fixed'] = 0.0         # minimum q_pulse->meas end margin


    # Acquisition
    settings['subtract_background'] = False
    settings['reps'] = 1
    settings['soft_avgs'] = 2e3

    return settings

# --------------------------------------------------------------------------
def fm_dependent_phase_offsets(fm):
    """
    Compute phi11_minus_12_deg and phi22_minus_21_deg as functions of modulation frequency fm.

    Relations:
        phi11_minus_12_deg = -1.229  + 360 * (-0.875e-9) * fm
        phi22_minus_21_deg = -5.274  + 360 * ( 0.910e-9) * fm

    Parameters
    ----------
    fm : float or array
        Modulation frequency in Hz.

    Returns
    -------
    phi11_minus_12_deg, phi22_minus_21_deg : float or array
        Phase differences in degrees.
    """

    phi11_minus_12_deg = -1.229 + 360.0 * (-0.875e-9) * fm
    phi22_minus_21_deg = -5.274 + 360.0 * ( 0.910e-9) * fm

    return phi11_minus_12_deg, phi22_minus_21_deg
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
def ac_modulation_calculator(
    dbeta1_dV1, phi11_minus_12_deg, dbeta1_dV2,
    dbeta2_dV1, dbeta2_dV2, phi22_minus_21_deg,
    beta1, phi1_prime_deg, beta2, phi2_prime_deg,
):
    """
    Solve for Vpp_1, phi_V1, Vpp_2, phi_V2 in

        [ [db1/dV1 * e^{i(φ11-φ12)},   db1/dV2              ]   [ Vpp1 e^{i φV1} ]   [ β1 e^{i φ1'} ]
          [db2/dV1,                    db2/dV2 * e^{i(φ22-φ21)}] [ Vpp2 e^{i φV2} ] = [ β2 e^{i φ2'} ] ]

    All phase arguments are in degrees. Returned phases are also in degrees (wrapped to [0, 360)).

    Parameters
    ----------
    dbeta1_dV1, dbeta1_dV2, dbeta2_dV1, dbeta2_dV2 : float
        Real calibration slopes dβ_j / dVpp_k.
    phi11_minus_12_deg, phi22_minus_21_deg : float
        Phase offsets (φ11 - φ12) and (φ22 - φ21) in degrees.
    beta1, beta2 : float
        Desired modulation magnitudes β1, β2.
    phi1_prime_deg, phi2_prime_deg : float
        Desired modulation phases φ1', φ2' in degrees.

    Returns
    -------
    Vpp1, phi_V1_deg, Vpp2, phi_V2_deg : floats
        Required peak-to-peak voltages and phases (deg) for the two drives.
    """

    # degrees → radians
    d2r = np.pi / 180.0
    phi11_12 = phi11_minus_12_deg * d2r
    phi22_21 = phi22_minus_21_deg * d2r
    phi1p = phi1_prime_deg * d2r
    phi2p = phi2_prime_deg * d2r

    # Build complex 2×2 calibration matrix
    A = np.array([
        [dbeta1_dV1 * np.exp(1j * phi11_12),      dbeta1_dV2],
        [dbeta2_dV1,                              dbeta2_dV2 * np.exp(1j * phi22_21)]
    ], dtype=complex)

    # Desired complex modulation vector
    b = np.array([
        beta1 * np.exp(1j * phi1p),
        beta2 * np.exp(1j * phi2p)
    ], dtype=complex)

    # Solve A x = b for x = [V1_complex, V2_complex]
    V1_complex, V2_complex = np.linalg.solve(A, b)

    # Convert to amplitude + phase (phase in degrees, wrapped 0–360)
    Vpp1 = np.abs(V1_complex)
    Vpp2 = np.abs(V2_complex)
    phi_V1_deg = np.degrees(np.angle(V1_complex)) % 360.0
    phi_V2_deg = np.degrees(np.angle(V2_complex)) % 360.0

    return Vpp1, phi_V1_deg, Vpp2, phi_V2_deg

#---------------------------------------------------------------------------



def fm_and_phi_scan(soc, soccfg, instruments, settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse_D1   = exp_globals['qubit_pulse_D1']
    q_pulse_D2   = exp_globals['qubit_pulse_D2']
    
    Dual_gen = instruments['DCsupply']

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

        'phase_reset'     : exp_settings['phase_reset']
    }

    # Prepare measurement windows (mainly to get timing aligned; total_samples not used later)
    prog_probe = MeasurementSequence(soccfg, config)
    total_samples = prog_probe.us2cycles(config['readout_length'], ro_ch=0)

    # Build sweep axes
    phi_prime_list = np.linspace(
        exp_settings['phi_prime_start'],
        exp_settings['phi_prime_end'],
        exp_settings['phi_prime_points']
    )
    fm_list = np.linspace(
        exp_settings['fm_start'],
        exp_settings['fm_end'],
        exp_settings['fm_points']
    )

    Nx = len(phi_prime_list)   # columns (x-axis, phi')
    Ny = len(fm_list)          # rows    (y-axis, fm)

    # Allocate results (Ny x Nx)
    amp_map  = np.zeros((Ny, Nx))
    ang_map  = np.zeros((Ny, Nx))
    I_map    = np.zeros((Ny, Nx))
    Q_map    = np.zeros((Ny, Nx))

    # Also save the drive solutions
    Vpp1_map  = np.zeros((Ny, Nx))
    phi_V1_map = np.zeros((Ny, Nx))
    Vpp2_map  = np.zeros((Ny, Nx))
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
    beta1 = exp_settings['beta1']
    beta2 = exp_settings['beta2']

    # Live plot figures
    fig_amp = plt.figure(1, figsize=(10, 7))
    plt.clf()
    plt.suptitle('fm and phi scan (Amplitude, live)')

    fig_phase = plt.figure(2, figsize=(10, 7))
    plt.clf()
    plt.suptitle('fm and phi scan (Phase, live)')

    tstart_global = time.time()

    for iy, fm in enumerate(fm_list):
        row_start = time.time()
        phi11_minus_12_deg, phi22_minus_21_deg = fm_dependent_phase_offsets(fm)

        # For per-row time estimate (after first 10 phi points)
        phi_row_start = time.time()

        for ix, phi_prime in enumerate(phi_prime_list):
            # Print which phi' point we are measuring
            print(
                f"Row {iy+1}/{Ny}, fm = {fm/1e6:.2f} MHz, "
                f"phi index {ix+1}/{Nx}, phi' = {phi_prime:.1f} deg"
            )

            # After first 10 phi points for this row, estimate remaining time for the row
            if ix == 9 and Nx > 10:
                phi_10_stop = time.time()
                print("  Estimated time for this fm row based on first 10 phi points:")
                estimate_time(phi_row_start, phi_10_stop, Nx/10)
                
                if iy == 0:
                    print("  Estimated time for all of the rows:")
                    estimate_time(tstart_global, phi_10_stop, int(Ny*Nx/10))

            Vpp1, phi_V1_deg, Vpp2, phi_V2_deg = ac_modulation_calculator(
                dbeta1_dV1, phi11_minus_12_deg,
                dbeta1_dV2, dbeta2_dV1, dbeta2_dV2, phi22_minus_21_deg,
                beta1, 0.0, beta2, phi_prime
            )
            # Store drive parameters
            Vpp1_map[iy, ix]  = Vpp1
            phi_V1_map[iy, ix] = phi_V1_deg
            Vpp2_map[iy, ix]  = Vpp2
            phi_V2_map[iy, ix] = phi_V2_deg

            # Program the external generator
            Dual_gen.Ch1_sin_gen(Vpp1, fm, phase=phi_V1_deg,
                                 offset=exp_settings['dc_offset_voltage_1'])
            Dual_gen.Ch2_sin_gen(Vpp2, fm, phase=phi_V2_deg,
                                 offset=exp_settings['dc_offset_voltage_2'])
            Dual_gen.phase_sync()
            time.sleep(0.1)

            # Run QICK program at this point
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
            # 2D scan: imshow over fm and phi'
            im_amp = ax_amp.imshow(
                amp_map[:iy+1, :],
                origin='lower',
                aspect='auto',
                extent=[
                    phi_prime_list[0],
                    phi_prime_list[-1],
                    fm_list[0] / 1e6,
                    fm_list[iy] / 1e6
                ]
            )
            plt.colorbar(im_amp, ax=ax_amp, label='Amplitude (arb)')
            ax_amp.set_xlabel(r'$\phi^\prime$ (deg)')
            ax_amp.set_ylabel(r'$f_m$ (MHz)')
            ax_amp.set_title('fm and phi scan \n' + filename)
        else:
            # 1D scan: fixed fm, sweep phi'
            fm_MHz = fm_list[0] / 1e6
            ax_amp.plot(phi_prime_list, amp_map[0, :], marker='o', linestyle='-')
            ax_amp.set_xlabel(r'$\phi^\prime$ (deg)')
            ax_amp.set_ylabel('Amplitude (arb)')
            ax_amp.set_title(f'phi scan at fm = {fm_MHz:.2f} MHz\n' + filename)
            ax_amp.grid(True)

        fig_amp.canvas.draw()
        fig_amp.canvas.flush_events()

        # --------- Live plot refresh: phase (ang_map) ---------
        plt.figure(fig_phase.number)
        plt.clf()
        ax_phase = plt.gca()

        if Ny > 1:
            im_phase = ax_phase.imshow(
                ang_map[:iy+1, :],
                origin='lower',
                aspect='auto',
                extent=[
                    phi_prime_list[0],
                    phi_prime_list[-1],
                    fm_list[0] / 1e6,
                    fm_list[iy] / 1e6
                ]
            )
            plt.colorbar(im_phase, ax=ax_phase, label='Phase (deg)')
            ax_phase.set_xlabel(r'$\phi^\prime$ (deg)')
            ax_phase.set_ylabel(r'$f_m$ (MHz)')
            ax_phase.set_title('fm and phi scan \n' + filename)
        else:
            fm_MHz = fm_list[0] / 1e6
            ax_phase.plot(phi_prime_list, ang_map[0, :], marker='o', linestyle='-')
            ax_phase.set_xlabel(r'$\phi^\prime$ (deg)')
            ax_phase.set_ylabel('Phase (deg)')
            ax_phase.set_title(f'phi scan at fm = {fm_MHz:.2f} MHz\n' + filename)
            ax_phase.grid(True)

        fig_phase.canvas.draw()
        fig_phase.canvas.flush_events()

        # Save rolling progress (optional)
        plt.figure(fig_amp.number)
        plt.savefig(os.path.join(saveDir, filename+'_live_amp.png'), dpi=150)
        plt.figure(fig_phase.number)
        plt.savefig(os.path.join(saveDir, filename+'_live_phase.png'), dpi=150)

        # Save data after each fm row
        userfuncs.SaveFull(
            saveDir, filename,
            [
                'phi_prime_list', 'fm_list',
                'amp_map', 'ang_map', 'I_map', 'Q_map',
                'Vpp1_map', 'phi_V1_map', 'Vpp2_map', 'phi_V2_map'
            ],
            locals(),
            expsettings=settings,
            instruments=instruments
        )

    tstop_global = time.time()
    print('Elapsed time (s): {:.3f}'.format(tstop_global - tstart_global))

    # Final save (same content; last one wins)
    userfuncs.SaveFull(
        saveDir, filename,
        [
            'phi_prime_list', 'fm_list',
            'amp_map', 'ang_map', 'I_map', 'Q_map',
            'Vpp1_map', 'phi_V1_map', 'Vpp2_map', 'phi_V2_map'
        ],
        locals(),
        expsettings=settings,
        instruments=instruments
    )

    return prog

