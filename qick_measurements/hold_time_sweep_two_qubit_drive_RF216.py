# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: kollarlab
"""

from qick.asm_v2 import AveragerProgramV2, QickSweep1D

import numpy as np
import os
import time
import matplotlib.pyplot as plt
import userfuncs
from utility.measurement_helpers import estimate_time
from utility.plotting_tools import simplescan_plot

# --- NEW: for fitting ---
import warnings
try:
    from scipy.optimize import curve_fit
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


      
#Heavily considering getting rid of the initial and post buffers for the speedup classes...
#Don't see the use when we can't acquire_decimated` anyway.

class HoldTimeSweep(AveragerProgramV2):
    def _initialize(self,cfg): 
        ro_ch  = cfg['ro_channel']
        gen_ch = cfg["cav_channel"]
        qub_ch = cfg["qub_channel"]

        # set the nyquist zone
        # Implement an if statement here to catch? 
        self.declare_gen(ch=cfg["cav_channel"], nqz=cfg["nqz_c"], mixer_freq=cfg['cav_mixer_freq'],ro_ch=ro_ch)
        self.declare_gen(ch=cfg["qub_channel"], nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq'])
        
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
        
        
        sigma = cfg["qub_sigma"] 
        num_sigma = cfg["num_sigma"]
        
        self.add_gauss(ch=qub_ch, name='ramp', sigma=sigma,length=sigma*num_sigma)
        
        self.add_pulse(ch=qub_ch, name="qub_pulse", ro_ch=ro_ch,
               style="flat_top",
               envelope="ramp",
               freq=cfg['qub_freq'],
               length= cfg['hold_length'],
               phase=cfg['qub_phase'],
               gain=cfg['qub_gain'],
              )
        
        self.add_pulse(
                ch=qub_ch, name="qub_phrst", ro_ch=ro_ch,
                style="const",
                freq=cfg["qub_freq"],      # doesn't really matter if gain=0
                phase=0,
                gain=0,                   # no output
                length=0.015,              # small but nonzero (us); pick safely > 0
                phrst=1                   # <-- resets phase accumulator
            )
        
        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)
        
        
    
    def _body(self, cfg):
# =============================================================================
#         qub_ch = self.cfg["qub_channel"]
#         self.reset_phase(gen_ch = [qub_ch], t=0)
# =============================================================================
        sigma = float(cfg["qub_sigma"])
        num_sigma = int(cfg["num_sigma"])
        
        pulse_len = float(cfg['hold_length']) + num_sigma*sigma

        offset = cfg["adc_trig_offset"]
        meas_time = self.cfg["meas_time"]
        ex_time = meas_time - cfg['qub_delay'] - pulse_len
        
        if ex_time < 0:
            print("Warning: Time Error, ex_time<0. Desired pulse goes outside the bounds of reality.")
        
        # Optional phase reset behavior
        if cfg.get("phase_reset", False):
            self.pulse(ch=cfg["qub_channel"], name="qub_phrst", t=0)
            
        #Sets off the ADC
        self.trigger(ros=[cfg['ro_channel']],
                    pins=[0],
                    t=offset)
        
        self.pulse(ch=cfg["qub_channel"],name='qub_pulse',t=ex_time)
        self.pulse(ch=cfg["cav_channel"],name='cav_pulse',t=meas_time)
        self.wait_auto()
        self.delay(self.cfg["relax_delay"])
        


def get_hold_time_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'HoldTime'
    
#    settings['cav_freq'] = 1e9
#    settings['cav_gain'] = 1000

#    settings['qub_gain'] = 1000

#    settings['quasi_CW_len'] = 10
    
    #Sweep parameters
#    settings['freq_start']   = 4e9  
#    settings['freq_stop']    = 4.5e9
#    settings['freq_points']  = 6

    #Card settings
#    settings['reps'] = 1
#    settings['soft_avgs'] = 5e3
    
    return settings

# --- NEW: cosine model + robust initial guesses ---
def _cos_model(t, A, omega, phi, C):
    # A*cos(omega*t + phi) + C
    return A*np.cos(omega*t + phi) + C

def _guess_cos_params(t_us, y):
    """
    Make robust initial guesses for A, omega, phi, C.
    t_us is in microseconds; omega will be in rad/us.
    """
    t = np.asarray(t_us)
    y = np.asarray(y)

    # Offset and amplitude guesses
    C0 = np.nanmean(y)
    A0 = 0.5*(np.nanmax(y) - np.nanmin(y))
    if not np.isfinite(A0) or A0 == 0:
        A0 = 1.0

    # Detrend for frequency guess
    yd = y - C0

    # Frequency guess via FFT peak (in cycles/us), then to rad/us
    # Use even spacing assumption (true here)
    if len(t) > 1:
        dt = np.nanmedian(np.diff(t))
        if not np.isfinite(dt) or dt <= 0:
            dt = (t[-1]-t[0])/max(1, (len(t)-1))
    else:
        dt = 1.0

    # Zero-pad a bit for frequency resolution
    n = len(yd)
    nfft = int(2**np.ceil(np.log2(max(8, n*4))))
    # real FFT freqs in cycles/us
    freqs = np.fft.rfftfreq(nfft, d=dt)
    spec = np.abs(np.fft.rfft((yd if np.all(np.isfinite(yd)) else np.zeros_like(yd)), n=nfft))

    # ignore DC for frequency guess
    if len(spec) > 2:
        spec[0] = 0.0

    k = int(np.nanargmax(spec)) if len(spec) > 0 else 1
    f0_cyc_per_us = max(freqs[k], 1e-6)  # avoid zero
    omega0 = 2*np.pi*f0_cyc_per_us      # rad/us

    # Phase guess from first point
    # If y(0) ~ A*cos(phi)+C -> phi ~= arccos((y0 - C)/A) (pick a sign)
    y0 = y[0]
    ratio = np.clip((y0 - C0)/max(A0, 1e-9), -1.0, 1.0)
    phi0 = np.arccos(ratio)  # crude; optimizer will refine

    return A0, omega0, phi0, C0

def fit_cosine(hpts_us, amp):
    """
    Fit A*cos(omega*t + phi) + C to (t, amp).
    Returns dict with params and covariance if available.
    """
    t = np.asarray(hpts_us, float)
    y = np.asarray(amp, float)

    # Clean NaNs/Infs
    m = np.isfinite(t) & np.isfinite(y)
    t = t[m]; y = y[m]

    if len(t) < 4:
        raise RuntimeError("Not enough points to fit.")

    p0 = _guess_cos_params(t, y)

    if _HAS_SCIPY:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, pcov = curve_fit(_cos_model, t, y, p0=p0, maxfev=20000)
        A, omega, phi, C = popt
        # --- enforce A >= 0 convention ---
        if A < 0:
            A = -A
            phi += np.pi
        
        # Wrap phase to [-π, π] for neatness
        phi = (phi + np.pi) % (2*np.pi) - np.pi
        
        return {
            "A": A,
            "omega": omega,
            "phi": phi,
            "C": C,
            "f_MHz": (omega / (2 * np.pi)),            # frequency in MHz
            "T_us": (2 * np.pi / abs(omega)) if omega != 0 else np.inf,
            "popt": np.array([A, omega, phi, C]),
            "pcov": pcov,
            "used_scipy": True,
        }
    
    else:
        # Fallback: linearize w.r.t cos/sin for a first pass, then refine with simple GN steps
        # y = B*cos(ωt) + D*sin(ωt) + C, with ω from guess
        A0, omega0, phi0, C0 = p0
        for _ in range(3):  # refine ω a bit
            X = np.column_stack([np.cos(omega0*t), np.sin(omega0*t), np.ones_like(t)])
            coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            B, D, C_lin = coeffs
            # convert to A/phi
            A_lin = np.hypot(B, D)
            phi_lin = np.arctan2(-D, B)  # because B*cos + D*sin == A*cos(ωt+phi) with phi = atan2(-D,B)
            # crude ω update via phase unwrapping around max slope point
            # (kept simple to avoid heavy code—good enough as fallback)
            omega0 = omega0  # keep initial; the linear pass already does most heavy lifting
            A0, C0, phi0 = A_lin, C_lin, phi_lin
            
        # --- enforce A >= 0 convention ---
        if A0 < 0:
            A0 = -A0
            phi0 += np.pi
        phi0 = (phi0 + np.pi) % (2*np.pi) - np.pi

        return {
            "A": A0, "omega": omega0, "phi": phi0, "C": C0,
            "f_MHz": (omega0/(2*np.pi)),
            "T_us": (2*np.pi/omega0) if omega0!=0 else np.inf,
            "popt": np.array([A0, omega0, phi0, C0]),
            "pcov": None, "used_scipy": False
        }

def hold_time_sweep(soc,soccfg,instruments,settings):

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    
    qub_drive_index = exp_settings['qub_drive_index']
    if qub_drive_index == 'D1':
        qub_ch = exp_globals['qub_channel_1']
        qub_channel = exp_globals['qub_channel_1']['ID']
        q_pulse = exp_globals['qubit_pulse_D1']
        
    elif qub_drive_index == 'D2':
        qub_ch = exp_globals['qub_channel_2']
        qub_channel = exp_globals['qub_channel_2']['ID']
        q_pulse = exp_globals['qubit_pulse_D2']
        
    else:
        print("Wrong qubit drive index, check carefully")
    
    soc.reset_gens()
    
    

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel']['ID'],
        'qub_channel'     : qub_channel,
        'ro_channel'     : exp_globals['ro_channel']['ID'],

        'nqz_c'           : 2,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'cav_gain'        : exp_settings['cav_gain'],
        'cav_freq'        : (exp_settings['cav_freq'])/1e6,
        'cav_mixer_freq'  : (exp_settings['cav_freq'] + exp_settings['cav_mixer_detuning'])/1e6,
        
        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],
        
        'qub_freq'        : exp_settings['qub_freq']/1e6,
        # 'freq_start'      : exp_settings['freq_start']/1e6,
        # 'freq_stop'       : exp_settings['freq_stop']/1e6,
        'qub_gain'        : exp_settings['qub_gain'],
        'qub_mixer_freq'  : (exp_settings['qub_freq']+exp_settings['qub_mixer_detuning'])/1e6,

        'qub_sigma'       : q_pulse['sigma'],
        'qub_delay'       : exp_globals['qub_delay_fixed'],
        'num_sigma'       : q_pulse['num_sigma'],
        'hold_length'        : 1, #Placeholder
        'hold_points'     : exp_settings['hold_points'],

        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        # 'soft_avgs'       : exp_settings['soft_avgs']
        
        'phase_reset'     : exp_settings['phase_reset']
        }

    cav_ch = exp_globals['cav_channel']
    #qub_ch = exp_globals['qub_channel']
    ro_ch  = exp_globals['ro_channel']
    # Set attenuator on DAC.
    soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
    soc.rfb_set_gen_rf(qub_ch['ID'], qub_ch['Atten_1'], qub_ch['Atten_2'])
    # Set attenuator on ADC.
    soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])
    
    if exp_settings['filter'] == 'all_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
        
        soc.rfb_set_gen_filter(config['qub_channel'], fc=config['qub_freq']/1000, ftype='bandpass', bw=qub_ch['BW'])
        
    elif exp_settings['filter'] == 'no_qubit_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
    
        soc.rfb_set_gen_filter(config['qub_channel'], fc=config['qub_freq']/1000, ftype='bypass')
        
    elif exp_settings['filter'] == 'no_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bypass')
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
        
        soc.rfb_set_gen_filter(config['qub_channel'], fc=config['qub_freq']/1000, ftype='bypass')
        
    else:
        print('Please select one option from:')
        print('\'all_filter\', \'no_qubit_filter\', and \'no_filter\'')
        return
        
    prog = HoldTimeSweep(soccfg,reps = exp_settings['reps'], final_delay = None, final_wait = 0, cfg = config)
    rep_period = config['adc_trig_offset'] + config['readout_length'] + config['relax_delay']
    
    
    
    projected_time = exp_settings['reps']*exp_settings['rounds']*config['hold_points']*rep_period/1e6
    print("Projected Time: " + str(projected_time))
    
    tstart = time.time()
    first_it = True
    
    hpts = np.linspace(exp_settings["hold_start"],exp_settings['hold_stop'],exp_settings['hold_points'])
    
    Is = np.zeros(len(hpts))
    Qs = np.zeros(len(hpts))
    powerdat = np.zeros(len(hpts))
    phasedat = np.zeros(len(hpts))
    
    for h in range(0,len(hpts)):
        config['hold_length'] = hpts[h]
        
        prog = HoldTimeSweep(soccfg,reps = exp_settings['reps'], final_delay = None, final_wait = 0, cfg = config)
    
        iq_list = prog.acquire(soc, rounds = exp_settings['rounds'], load_pulses=True, progress=False)
        
        Is[h] = iq_list[0][0][0]
        Qs[h] = iq_list[0][0][1]
        powerdat[h] = np.sqrt(Is[h]**2 + Qs[h]**2)
        phasedat[h] = np.degrees(np.arctan2(Qs[h], Is[h]))
        
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(hpts))
            first_it = False
            
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(hpts, powerdat)
        plt.xlabel('Hold Time (us)')
        plt.ylabel('Amplitude')  
        plt.subplot(122)
        plt.plot(hpts, phasedat)
        plt.xlabel('Hold Time (us)')
        plt.ylabel('Phase')  
        plt.suptitle('Live Rabi data \n'+filename)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'_live.png'), dpi = 150)
     

        userfuncs.SaveFull(saveDir, filename, ['hpts', 'powerdat', 'phasedat'], locals(), expsettings=settings, instruments=instruments)
        
    # --- NEW: Fit A*cos(Ω t + φ) + C to (hpts, powerdat) ---
    try:
        fitres = fit_cosine(hpts, powerdat)  # hpts is in microseconds already
        A, omega, phi, C = fitres["A"], fitres["omega"], fitres["phi"], fitres["C"]
        fMHz, T_us = fitres["f_MHz"], fitres["T_us"]
        
        # --- NEW: find first peak (earliest local maximum within scan) ---
        # Normalize sign of omega for robust timing
        omega_eff = abs(omega)
        phi_eff = phi if omega >= 0 else -phi
        
        tmin = float(np.min(hpts))
        tmax = float(np.max(hpts))
        
        # For A>=0, maxima at cos(...) = +1  -> angle = 2πk
        # For A<0,  maxima at cos(...) = -1  -> angle = π + 2πk
        if A >= 0:
            # t_k = (-phi_eff + 2πk)/omega_eff
            base = 0.0
        else:
            # t_k = (-phi_eff + π + 2πk)/omega_eff
            base = np.pi
        
        # Small helper to get k from a target t
        def k_from_t(t_target):
            return np.ceil((omega_eff*t_target + phi_eff - base) / (2*np.pi))
        
        k0 = int(k_from_t(tmin))
        # Generate a small set of candidate peaks around the start
        ks = np.arange(k0-1, k0+4)  # a few candidates is enough
        t_candidates = ( -phi_eff + (base + 2*np.pi*ks) ) / omega_eff
        
        # Keep only those inside the scan window
        t_candidates = t_candidates[(t_candidates >= tmin) & (t_candidates <= tmax)]
        if len(t_candidates) == 0:
            # Fallback: take argmax of the smooth fit curve near start
            tfit_dense = np.linspace(tmin, tmax, 2001)
            yfit_dense = _cos_model(tfit_dense, A, omega, phi, C)
            t_first_peak = float(tfit_dense[np.argmax(yfit_dense)])
        else:
            t_first_peak = float(np.min(t_candidates))
        
        y_first_peak = float(_cos_model(t_first_peak, A, omega, phi, C))

    
        # High-resolution fit curve for plotting
        tfit = np.linspace(np.min(hpts), np.max(hpts), 1000)
        yfit = _cos_model(tfit, A, omega, phi, C)
    
        # Plot measured amplitudes + fit
        fig_fit = plt.figure(2, figsize=(10,6))
        plt.clf()
        plt.plot(hpts, powerdat, 'o', label='Data')
        plt.plot(tfit, yfit, '-', label='Fit')
        # --- NEW: mark first peak ---
        plt.axvline(t_first_peak, linestyle='--', alpha=0.7, label=f'First peak ≈ {t_first_peak:.4g} μs')
        plt.plot([t_first_peak], [y_first_peak], 's', label='Peak point')
        
        plt.xlabel('Hold Time (us)')
        plt.ylabel('Amplitude')
        title_str = (f'Rabi fit: A*cos(Ω t + φ) + C\n'
                     f'A={A:.3g}, Ω={fMHz:.4g} MHz, φ={phi*180/np.pi:.2g} degrees, C={C:.4g}  |  '
                     f'T={T_us:.4g} μs')
        plt.suptitle(filename)
        plt.title(title_str)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(saveDir, filename+'_fit.png'), dpi=150)
    
        # Optional: print a concise line to console
        print(title_str)
    
        # Save fitted parameters as well
        userfuncs.SaveFull(
            saveDir, filename+'_fit',
            ['A','omega','phi','C','fMHz','T_us'],
            locals(),
            expsettings=settings,
            instruments=instruments
        )
    
    except Exception as e:
        print(f"[WARN] Rabi cosine fit failed: {e}")
    
# =============================================================================
#     hold_times = hpts
#     powerdat = np.sqrt(Is**2 + Qs**2)
#     phasedat = np.arctan2(Qs,Is)*180/np.pi
# =============================================================================
    
    full_data = {} 
    full_data['xaxis'] = hpts
    full_data['mags'] = powerdat 
    full_data['phases'] = phasedat 
    full_data['Is'] = Is 
    full_data['Qs'] = Qs
    

    
# =============================================================================
#     # # ------------------- Plot: one figure, two subplots -------------------
#     fig, (ax1, ax2) = plt.subplots(2, 1, num=1, figsize=(7, 7), sharex=True)
#     fig.clf()
#     fig, (ax1, ax2) = plt.subplots(2, 1, num=1, figsize=(7, 7), sharex=True)
#     fig.suptitle(filename)  # file name as suptitle
#     
#    
#     ax1.plot(full_data['xaxis'], full_data['mags'], label='Data')
#     ax1.set_ylabel('Amplitude')
#     #ax1.legend(loc='best')
#     ax1.grid()
#     
#     # Subplot 2: phase
#     ax2.plot(full_data['xaxis'], full_data['phases'], '.')
#     ax2.set_xlabel('Hold Time (us)')
#     ax2.set_ylabel('Phase (deg)')
#     ax2.grid()
#     
#     fig.tight_layout(rect=[0, 0, 1, 0.95])  # leave space for suptitle
#     plt.savefig(os.path.join(saveDir, filename+'_mag_phase.png'), dpi=150)
#     # # =====================================================================
# 
#     
#     
#     
#     userfuncs.SaveFull(saveDir, filename, ['hold_times','full_data','filename'],
#     locals(), expsettings=settings, instruments={})
# =============================================================================
    


    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))

    data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data}

    return data,prog
