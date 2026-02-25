# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 15:13:12 2026

@author: KollarLab
"""

from qick.asm_v2 import AveragerProgramV2, QickSweep1D

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
import time

from utility.userfits import fit_T2
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import estimate_time

class MeasurementSequence(AveragerProgramV2):
    def _initialize(self,cfg): 
        ro_ch  = cfg['ro_channel']
        gen_ch = cfg["cav_channel"]
        qub_ch_D1 = cfg["qub_channel_D1"]
        qub_ch_D2 = cfg["qub_channel_D2"]

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=cfg["nqz_c"], mixer_freq=cfg['cav_mixer_freq'],ro_ch=ro_ch)
        
        if cfg.get('qub_mixer_freq_fixed') is not None:
            self.declare_gen(ch=qub_ch_D1, nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq_fixed'])
            self.declare_gen(ch=qub_ch_D2, nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq_fixed'])
            #print(cfg['qub_mixer_freq_fixed'])
            
        else:
            self.declare_gen(ch=qub_ch_D1, nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq'])
            self.declare_gen(ch=qub_ch_D2, nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq'])
        
        # Readout declaration
        self.declare_readout(ch=ro_ch, length=cfg['readout_length'])
        self.add_readoutconfig(ch=ro_ch, name="myro",
                           freq=cfg['cav_freq'],
                           gen_ch=gen_ch,
                           outsel='product')
        
        # Add qubit spec frequency loop
        self.add_loop("freq_sweep", cfg["det_points"])
        
        
        # Cavity pulse
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
                       phrst=1, # testing
                       length=cfg.get("hold_length_D1", 0.0),
                       phase=cfg["qub_phase_D1"],
                       gain=cfg["qub_gain_D1"])
        
        self.add_pulse(ch=qub_ch_D2, name="ex1_D2", ro_ch=ro_ch,
                       style="flat_top", envelope="env_D2",
                       freq=cfg["qub_freq"],
                       phrst=1, # testing
                       length=cfg.get("hold_length_D2", 0.0),
                       phase=cfg["qub_phase_D2"],
                       gain=cfg["qub_gain_D2"])
        
        # pulse #2
        self.add_pulse(ch=qub_ch_D1, name="ex2_D1", ro_ch=ro_ch,
                       style="flat_top", envelope="env_D1",
                       freq=cfg["qub_freq"],
                       phrst=1, # testing
                       length=cfg.get("hold_length_D1", 0.0),
                       phase=cfg.get("qub_phase2_D1", cfg["qub_phase_D1"]), #find the definition of qub_phase2_Dx later
                       gain=cfg.get("qub_gain2_D1", cfg["qub_gain_D1"]))
        
        self.add_pulse(ch=qub_ch_D2, name="ex2_D2", ro_ch=ro_ch,
                       style="flat_top", envelope="env_D2",
                       freq=cfg["qub_freq"],
                       phrst=1, # testing
                       length=cfg.get("hold_length_D2", 0.0),
                       phase=cfg.get("qub_phase2_D2", cfg["qub_phase_D2"]),
                       gain=cfg.get("qub_gain2_D2", cfg["qub_gain_D2"]))


        
# =============================================================================
#         # Phase reset pulses (phrst=1 resets phase accumulator)
#         self.add_pulse(ch=qub_ch_D1, name="phrst_D1", ro_ch=ro_ch,
#                        style="const",
#                        freq=cfg["qub_freq"], phase=0, gain=0,
#                        length=0.010, phrst=1)
#         
#         self.add_pulse(ch=qub_ch_D2, name="phrst_D2", ro_ch=ro_ch,
#                        style="const",
#                        freq=cfg["qub_freq"], phase=0, gain=0,
#                        length=0.010, phrst=1)
# =============================================================================
        
        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)
        
        #self.delay(0.6) # giving some time for pulse to be configured?
        
        
    
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
        
        # Estimate pulse lengths (us) for scheduling
        def pulse_len(which):
            if which == "D1":
                return cfg.get("hold_length_D1", 0.0) + cfg["qub_sigma_D1"]*cfg["num_sigma_D1"]
            else:
                return cfg.get("hold_length_D2", 0.0) + cfg["qub_sigma_D2"]*cfg["num_sigma_D2"]
        
        len1 = pulse_len(order[0])
        len2 = pulse_len(order[1])
        
        # Pulse #2 ends qub_delay_fixed before measurement pulse
        t2_start = meas_time - float(cfg["qub_delay_fixed"]) - len2
        t1_start = t2_start

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

        
        # trigger ADC
        self.trigger(ros=[cfg["ro_channel"]], pins=[0], t=cfg["adc_trig_offset"])

        
        # ---- Play pulses ----
        self.pulse(ch=ch1, name=p1, t=t1_start)
        self.pulse(ch=ch2, name=p2, t=t2_start)
        
        # ---- Measurement ----
        self.pulse(ch=cfg["cav_channel"], name="cav_pulse", t=meas_time)
        self.wait_auto()
        self.delay(cfg["relax_delay"])


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
      
def get_fm_and_detuning_settings():
    
    settings = {}
    
    settings['scanname']    = 'Search_for_EP_and_gyration_fm_and_spec_scan'
    settings['meas_type']   = 'fm_and_spec_scan'
    #settings['phase_reset'] = True   # 

    # -------- Cavity --------
    settings['cav_freq']  = 6.0e9     # Hz, can override in driver
    settings['cav_gain'] = 0
    settings['cav_mixer_detuning'] = 200e6

    # -------- Qubit drives (D1 / D2) --------
    settings['qub_freq']  = 4.5e9
    settings['qub_mixer_detuning'] = -250e6
    settings['qub_mixer_freq_fixed'] = 4.65e9
    settings["qub_gain_D1"]  = 0.15
    settings['hold_time_D1_us'] = 0.1862
    settings["qub_gain_D2"]  = 0.7588
    settings['hold_time_D2_us'] = 0.1846    

    settings["phase2_deg"] = 0

    # -------- Scheduling between the two qubit pulses --------
    # pulse_order: which channel is "pulse 1" and "pulse 2"
    settings['pulse_order'] = ('D2', 'D1')  # e.g. ('D1','D2')
    settings['phase2_deg']  = 0.0           # extra phase on pulse #2 (deg)


    # -------- AC modulation calibration + targets --------
    # Slopes dβ_j/dV_k (units: your calibration convention)
    settings['dbeta1_dV1'] = 0.0   # TODO: fill from calibration
    settings['dbeta1_dV2'] = 0.0   # TODO
    settings['dbeta2_dV1'] = 0.0   # TODO
    settings['dbeta2_dV2'] = 0.0   # TODO

    # Target modulation amplitudes β1, β2
    settings['beta1'] = 0.0        # TODO
    settings['beta2'] = 0.0        # TODO

    # DC offsets on the two modulation channels (V)
    settings['dc_offset_voltage_1'] = 0.0   # TODO
    settings['dc_offset_voltage_2'] = 0.0   # TODO

    # Relative modulation phase φ2′ (φ1′ = 0 by construction in fm_and_detuning_scan)
    settings['phi2_prime_deg'] = 0.0        # deg

    # -------- fm sweep (rows, y-axis) --------
    # Modulation frequency in Hz
    settings['fm_start']  = 0.5e6    # 0.5 MHz
    settings['fm_end']    = 5.0e6    # 5 MHz
    settings['fm_points'] = 41       # number of fm rows

    # -------- Detuning sweep (columns, x-axis) --------
    # Detuning in MHz; actual qubit frequency = qub_freq + detuning
    settings['det_min_MHz'] = -25.0
    settings['det_max_MHz'] = +25.0
    settings['det_points']  = 51

    # -------- Acquisition --------
    settings['reps']   = 5000
    settings['rounds'] = 1
    settings['filter'] = 'all_filter'


    return settings

def fm_and_detuning_scan(soc,soccfg,instruments,settings):
    """
    2D scan over modulation frequency fm and qubit detuning.
    - y-axis: fm (MHz)
    - x-axis: qubit pulse frequency (GHz) = qub_freq_D1_nominal + detuning

    Uses:
      - AC modulation calculator (Vpp1, phi_V1, Vpp2, phi_V2 as a function of fm)
      - detuning-aware MeasurementSequence (via cfg['detuning_MHz'])

    Time estimates:
      - No per-column prints.
      - After the FIRST fm row, print estimated TOTAL run time.
      - After EACH subsequent row, update estimated REMAINING time.
    """
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    
    Dual_gen = instruments['DCsupply']

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

        
        
        'qub_freq'        : QickSweep1D("freq_sweep", exp_settings['det_min_MHz']+exp_settings['qub_freq']/1e6, 
                                        exp_settings['det_max_MHz']+exp_settings['qub_freq']/1e6),
        'det_points'      : exp_settings['det_points'],
        'qub_mixer_freq'  : (exp_settings['qub_freq']+exp_settings['qub_mixer_detuning'])/1e6,
        'qub_mixer_freq_fixed' : exp_settings['qub_mixer_freq_fixed']/1e6,      # useful for fixing the mixer phase
        'qub_gain_D1'     : exp_settings['qub_gain_D1'],
        'qub_gain_D2'     : exp_settings['qub_gain_D2'],
        'qub_gain2_D1'    : exp_settings['qub_gain_D1'], # by default the same, will get updated in the 'gain' mode
        'qub_gain2_D2'    : exp_settings['qub_gain_D2'], # by default the same, will get updated in the 'gain' mode

        
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

        'pulse_order'     : exp_settings['pulse_order'],    # ('D1','D2'), etc.
        
        'phase2_deg'      : exp_settings.get('phase2_deg', 0.0), # used in interval mode (fixed) or overwritten in loop
        'tau_us'          : 0.0,                                 # overwritten in interval loop

        }
    
    # pulse_2 phase block
    order = tuple(config["pulse_order"])

    phi1 = config["qub_phase_D1"] if order[0] == "D1" else config["qub_phase_D2"]
    dphi_fixed = float(exp_settings.get("phase2_deg", 0.0))
    phi2 = phi1 + dphi_fixed

    if order[1] == "D1":
        config["qub_phase2_D1"] = phi2
    else:
        config["qub_phase2_D2"] = phi2
    
    
    # normal channel settings    
    cav_ch = exp_globals['cav_channel']
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


    
    # --- Build sweep axes ---

    # fm (Hz)
    fm_list = np.linspace(
        exp_settings['fm_start'],
        exp_settings['fm_end'],
        exp_settings['fm_points']
    )

    # detuning (MHz)
    det_list = np.linspace(
        exp_settings['det_min_MHz'],
        exp_settings['det_max_MHz'],
        exp_settings['det_points']
    )

    Nx = len(det_list)   # columns (x-axis, detuning)
    Ny = len(fm_list)    # rows    (y-axis, fm)
    
    # x-axis: qubit frequency (GHz) = f_qubit_nominal + detuning
    qub_freq_D1_nom_Hz = exp_settings['qub_freq']  # Hz
    freq_qubit_list_GHz = qub_freq_D1_nom_Hz*1e-9 + det_list*1e-3  # GHz
    
    
    # Allocate results (Ny x Nx)
    amp_map   = np.zeros((Ny, Nx))
    ang_map   = np.zeros((Ny, Nx))
    I_map     = np.zeros((Ny, Nx))
    Q_map     = np.zeros((Ny, Nx))

    Vpp1_vec   = np.zeros(Ny)
    phi_V1_vec = np.zeros(Ny)
    Vpp2_vec   = np.zeros(Ny)
    phi_V2_vec = np.zeros(Ny)
    
    
    
    # Save base phases so we can restore and only rotate pulse #2
    base_phase_D1 = float(config["qub_phase_D1"])
    base_phase_D2 = float(config["qub_phase_D2"])
    
   
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

    # Modulation phases for the AC calculator
    phi1_prime_deg = 0.0
    phi2_prime_deg = float(exp_settings.get('phi2_prime_deg', 0.0))
    
    # String to annotate beta and phi2' in plot titles
    title_params = (r"$\beta_1$={:.3f} MHz, $\beta_2$={:.3f} MHz, "
                    r"$\phi_2^\prime$={:.1f}°").format(beta1/1e6, beta2/1e6, phi2_prime_deg)

    # Live plot figures
    fig_amp = plt.figure(1, figsize=(10, 7))
    plt.clf()
    plt.suptitle('fm and spec frequency scan (Amplitude, live)')

    fig_phase = plt.figure(2, figsize=(10, 7))
    plt.clf()
    plt.suptitle('fm and spec frequency scan (Phase, live)')

    # --- Time-estimation bookkeeping (Rabi-chevron style) ---
    tstart_global  = time.time()
    first_row_done = False
    
    # ----------------------------
    # Main sweep loop
    # ----------------------------
    tstart = time.time()
    
    for iy, fm in enumerate(fm_list):
        row_start = time.time()

        phi11_minus_12_deg, phi22_minus_21_deg = fm_dependent_phase_offsets(fm)

        # AC modulation solution for this fm
        Vpp1, phi_V1_deg, Vpp2, phi_V2_deg = ac_modulation_calculator(
            dbeta1_dV1, phi11_minus_12_deg,
            dbeta1_dV2, dbeta2_dV1, dbeta2_dV2, phi22_minus_21_deg,
            beta1, phi1_prime_deg, beta2, phi2_prime_deg
        )
        
        Vpp1_vec[iy]   = Vpp1
        phi_V1_vec[iy] = phi_V1_deg
        Vpp2_vec[iy]   = Vpp2
        phi_V2_vec[iy] = phi_V2_deg

        # Program the external generators
        Dual_gen.Ch1_sin_gen(Vpp1, fm, phase=phi_V1_deg,
                             offset=exp_settings['dc_offset_voltage_1'])
        Dual_gen.Ch2_sin_gen(Vpp2, fm, phase=phi_V2_deg,
                             offset=exp_settings['dc_offset_voltage_2'])
        Dual_gen.phase_sync()
        time.sleep(0.1)
                
        # Run QICK program
        prog = MeasurementSequence(soccfg, reps=exp_settings["reps"], final_delay=None, final_wait=0, cfg=config)
        iq_list = prog.acquire(soc, rounds=exp_settings["rounds"], load_pulses=True, progress=False)
        
        iq = iq_list[0][0] # (Nx,2)
        I = iq[:, 0]                # (Nx,)
        Q = iq[:, 1]                # (Nx,)
        
        
        # Read I/Q 
        I_map[iy,:] = I
        Q_map[iy,:] = Q
    
        amp_map[iy, :] = np.sqrt(I**2 + Q**2)
        ang_map[iy, :] = np.degrees(np.arctan2(Q, I))
    
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
                    fm_list[0] / 1e6,
                    fm_list[iy] / 1e6
                ]
            )
            plt.colorbar(im_amp, ax=ax_amp, label='Amplitude (arb)')
            ax_amp.set_xlabel('Frequency (GHz)')
            ax_amp.set_ylabel(r'$f_m$ (MHz)')
            ax_amp.set_title(f'fm and spec frequency scan \n{title_params}\n{filename}')
        else:
            fm_MHz = fm_list[0] / 1e6
            ax_amp.plot(freq_qubit_list_GHz, amp_map[0, :], marker='o', linestyle='-')
            ax_amp.set_xlabel('Frequency (GHz)')
            ax_amp.set_ylabel('Amplitude (arb)')
            ax_amp.set_title(f'Spec frequency scan at fm = {fm_MHz:.3f} MHz\n{title_params}\n{filename}')
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
                    fm_list[0] / 1e6,
                    fm_list[iy] / 1e6
                ]
            )
            plt.colorbar(im_phase, ax=ax_phase, label='Phase (deg)')
            ax_phase.set_xlabel('Frequency (GHz)')
            ax_phase.set_ylabel(r'$f_m$ (MHz)')
            ax_phase.set_title(f'fm and spec frequency scan \n{title_params}\n{filename}')
        else:
            fm_MHz = fm_list[0] / 1e6
            ax_phase.plot(freq_qubit_list_GHz, ang_map[0, :], marker='o', linestyle='-')
            ax_phase.set_xlabel('Frequency (GHz)')
            ax_phase.set_ylabel('Phase (deg)')
            ax_phase.set_title(f'Spec frequency scan at fm = {fm_MHz:.3f} MHz\n{title_params}\n{filename}')
            ax_phase.grid(True)

        fig_phase.canvas.draw()
        fig_phase.canvas.flush_events()

        # Save rolling progress
        plt.figure(fig_amp.number)
        plt.savefig(os.path.join(saveDir, filename+'_live_amp.png'), dpi=150)
        plt.figure(fig_phase.number)
        plt.savefig(os.path.join(saveDir, filename+'_live_phase.png'), dpi=150)

        # Save data after each fm row
        userfuncs.SaveFull(
            saveDir, filename,
            [
                'fm_list', 'det_list', 'freq_qubit_list_GHz',
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
            print(f"Finished first fm row ({iy+1}/{Ny}).")
            print("Estimated TOTAL run time based on this row:")
            estimate_time(row_start, row_stop, Ny)
            first_row_done = True
        else:
            remaining_rows = Ny - (iy + 1)
            if remaining_rows > 0:
                print(f"Finished fm row {iy+1}/{Ny}.")
                print("Estimated REMAINING time based on this row:")
                estimate_time(row_start, row_stop, remaining_rows)

    tstop_global = time.time()
    print('Elapsed time (s): {:.3f}'.format(tstop_global - tstart_global))

    # Final save
    userfuncs.SaveFull(
        saveDir, filename,
        [
            'fm_list', 'det_list', 'freq_qubit_list_GHz',
            'amp_map', 'ang_map', 'I_map', 'Q_map',
            'Vpp1_map', 'phi_V1_map', 'Vpp2_map', 'phi_V2_map'
        ],
        locals(),
        expsettings=settings,
        instruments=instruments
    )

    return prog



   
    

        