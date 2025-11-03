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

# --- NEW: for fitting ---
import warnings
try:
    from scipy.optimize import curve_fit
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

      
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
        phase_q = self.deg2reg(cfg["qub_phase"], gen_ch=qub_ch)
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
        
        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],
        
# =============================================================================
#         'readout_length'  : m_pulse['init_buffer'] + m_pulse['meas_window'] + m_pulse['post_buffer'],
#         'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'] - m_pulse['init_buffer'],
# =============================================================================


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs']
        }

    hpts = (exp_settings["hold_start"]+exp_settings["hold_step"]*np.arange(exp_settings["hold_points"]))*1e6 #converts s to us

    prog = RabiSequence(soccfg, config)
    
# ============================================================================= old code with limits on readout length and # of ramps
#     meas_start = prog.us2cycles(m_pulse["init_buffer"],ro_ch=0)
#     meas_end = meas_start+prog.us2cycles(m_pulse["meas_window"],ro_ch=0)
#     total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)
# =============================================================================
   

    # powerdat = np.zeros((len(gpts), len(fpts)))
    # phasedat = np.zeros((len(gpts), len(fpts)))
    
    amp_int = np.zeros(len(hpts))
    ang_int = np.zeros(len(hpts))
    
# ============================================================================= also from the old code
#     amps    = np.zeros((len(hpts),total_samples))
#     angles  = np.zeros((len(hpts),total_samples))
# =============================================================================

    #drive_powers_lin = 10**(powers/10) ?
    #drive_amps_lin = np.sqrt(drive_powers_lin)
    
    tstart = time.time()
    first_it = True

    for h in range(0,len(hpts)):
        print("Current Hold Time (us): " + str(np.round(hpts[h], 3)) + ", Max Hold Time (us): " + str(np.round(hpts[-1], 3)))
        config["length"] = hpts[h]
        
        prog = RabiSequence(soccfg, config)

        # Integrated (scalar) I/Q, same as your T2 driver
        holder = prog.acquire(soc, load_pulses=True, progress=False)
        I_sig = holder[0][0][0]   # I, ro_ch 0, avg 0
        Q_sig = holder[1][0][0]   # Q, ro_ch 0, avg 0
        
        amp_int[h] = np.sqrt(I_sig**2 + Q_sig**2)
        ang_int[h] = np.degrees(np.arctan2(Q_sig, I_sig))
        
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(hpts))
            first_it = False

# ============================================================================= from old code
#         total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)
#         
#         prog = RabiSequence(soccfg,config)
#         
#         #Need to assign Iwindow, Qwindow, Ifull, Qfull, xaxis (which should just be timeus)
#         holder = prog.acquire_decimated(soc, load_pulses=True, progress=False, debug=False)
#         I_full = holder[0][0]
#         Q_full = holder[0][1]
#         I_window = I_full[meas_start:meas_end]
#         Q_window = Q_full[meas_start:meas_end]
#        
#         I_final = np.mean(I_window)
#         Q_final = np.mean(Q_window)
# 
#         amps[h] = np.sqrt(I_full**2+Q_full**2)
#         angles[h] = np.arctan2(Q_full,I_full)*180/np.pi
#         amp_int[h] = np.sqrt(I_final**2+Q_final**2)
#         ang_int[h] = np.arctan2(Q_final, I_final)*180/np.pi
# =============================================================================
        
# ============================================================================= old code
#         if first_it:
#             xaxis = np.linspace(0,len(I_full)-1,len(I_full))
#             
#             tstop = time.time()
#             estimate_time(tstart, tstop, len(hpts))
# 
#             for x in range(0,len(xaxis)):
#                 xaxis[x] = prog.cycles2us(xaxis[x],ro_ch=0)
#                 
#             first_it = False
# =============================================================================
            
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
        plt.suptitle('Live Rabi data \n'+filename)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'_live.png'), dpi = 150)
    

        userfuncs.SaveFull(saveDir, filename, ['hpts', 'amp_int', 'ang_int'], locals(), expsettings=settings, instruments=instruments)

    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))
    
    # --- NEW: Fit A*cos(Ω t + φ) + C to (hpts, amp_int) ---
    try:
        fitres = fit_cosine(hpts, amp_int)  # hpts is in microseconds already
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
        plt.plot(hpts, amp_int, 'o', label='Data')
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

        




       
    
    if exp_globals['LO']:
        logen.output = 0
        
    return hpts, amp_int