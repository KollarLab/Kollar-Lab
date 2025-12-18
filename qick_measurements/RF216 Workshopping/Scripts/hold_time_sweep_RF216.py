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
from utility.plotting_tools import simplescan_plot


      
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
        
        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)
        
        
    
    def _body(self, cfg):
# =============================================================================
#         qub_ch = self.cfg["qub_channel"]
#         self.reset_phase(gen_ch = [qub_ch], t=0)
# =============================================================================
        sigma = cfg["qub_sigma"]
        num_sigma = cfg["num_sigma"]
        
        pulse_len = cfg['hold_length'] + int(num_sigma*sigma)

        offset = cfg["adc_trig_offset"]
        meas_time = self.cfg["meas_time"]
        ex_time = meas_time - cfg['qub_delay'] - pulse_len
        
        
        
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

def hold_time_sweep(soc,soccfg,instruments,settings):

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    
    soc.reset_gens()
    
    

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel']['ID'],
        'qub_channel'     : exp_globals['qub_channel']['ID'],
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
        'qub_delay'       : q_pulse['delay'],
        'num_sigma'       : q_pulse['num_sigma'],
        'hold_length'        : 1, #Placeholder
        'hold_points'     : exp_settings['hold_points'],

        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        # 'reps'            : exp_settings['reps'],
        # 'soft_avgs'       : exp_settings['soft_avgs']
        }

    cav_ch = exp_globals['cav_channel']
    qub_ch = exp_globals['qub_channel']
    ro_ch  = exp_globals['ro_channel']
    # Set attenuator on DAC.
    soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
    soc.rfb_set_gen_rf(qub_ch['ID'], qub_ch['Atten_1'], qub_ch['Atten_2'])
    # Set attenuator on ADC.
    soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])
    
    if exp_settings['filter']:
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=1.0)
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=1.0)
        
        soc.rfb_set_gen_filter(config['qub_channel'], fc=config['qub_freq']/1000, ftype='bandpass', bw=1.0)
        

    prog = HoldTimeSweep(soccfg,reps = exp_settings['reps'], final_delay = None, final_wait = 0, cfg = config)
    rep_period = config['adc_trig_offset'] + config['readout_length'] + config['relax_delay']
    
    
    
    projected_time = exp_settings['reps']*exp_settings['rounds']*config['hold_points']*rep_period/1e6
    print("Projected Time: " + str(projected_time))
    
    t_i = time.time()
    
    hpts = np.linspace(exp_settings["hold_start"],exp_settings['hold_stop'],exp_settings['hold_points'])
    
    Is = np.zeros(len(hpts))
    Qs = np.zeros(len(hpts))
    
    for h in range(0,len(hpts)):
        config['hold_length'] = hpts[h]
        
        prog = HoldTimeSweep(soccfg,reps = exp_settings['reps'], final_delay = None, final_wait = 0, cfg = config)
    
        iq_list = prog.acquire(soc, rounds = exp_settings['rounds'], load_pulses=True, progress=False)
        
        Is[h] = iq_list[0][0][0]
        Qs[h] = iq_list[0][0][1]
        
    
    hold_times = hpts
    powerdat = np.sqrt(Is**2 + Qs**2)
    phasedat = np.arctan2(Qs,Is)*180/np.pi
    
    full_data = {} 
    full_data['xaxis'] = hpts
    full_data['mags'] = powerdat 
    full_data['phases'] = phasedat 
    full_data['Is'] = Is 
    full_data['Qs'] = Qs
    
    # # ================= One-figure Lorentzian fit + mag/phase subplots =================
    
    # is_peak = bool(exp_globals.get('hanger', True))  # True => peak, False => dip
    
    # # Data in GHz for x
    # x = holds / 1e9
    # y = powerdat.copy()  # magnitude
    # phase_deg = phasedat # already degrees
    
    # # --- robust baseline guess from edges ---
    # n = max(1, len(x)//10)
    # edge_vals = np.concatenate([y[:n], y[-n:]]) if len(y) >= 2*n else y
    # baseline0 = float(np.median(edge_vals))
    
    # # Peak vs dip handling (use "peakified" data to find extremum)
    # y_work = y if is_peak else (-y)
    # imax = int(np.argmax(y_work))
    # x0_0 = float(x[imax])
    
    # # Initial amplitude (allow negative if dip)
    # A0 = float(y[imax] - baseline0) if is_peak else float(-(y[imax] - baseline0))
    # if A0 == 0:
    #     A0 = (np.max(y) - np.min(y)) * (1.0 if is_peak else -1.0)
    
    # # Rough gamma (HWHM) from crossings; fallback to fraction of span
    # def _gamma_guess(x, y, imax, baseline):
    #     target = baseline + 0.5*(y[imax] - baseline)
    #     left = None
    #     for i in range(imax, 0, -1):
    #         if (y[i-1]-target) * (y[i]-target) <= 0:
    #             t = (target - y[i-1]) / (y[i] - y[i-1] + 1e-18)
    #             left = x[i-1] + t*(x[i]-x[i-1])
    #             break
    #     right = None
    #     for i in range(imax, len(x)-1):
    #         if (y[i]-target) * (y[i+1]-target) <= 0:
    #             t = (target - y[i]) / (y[i+1] - y[i] + 1e-18)
    #             right = x[i] + t*(x[i+1]-x[i])
    #             break
    #     if left is not None and right is not None and right > left:
    #         fwhm = right - left
    #         return float(max(fwhm/2.0, 1e-9))  # gamma = FWHM/2
    #     span = max(1e-6, (np.max(x)-np.min(x)))
    #     return 0.05 * span
    
    # gamma0 = _gamma_guess(x, y, imax, baseline0)
    
    # # Lorentzian model: baseline + A / (1 + ((x - x0)/gamma)^2)
    # def lorentz(xv, baseline, A, x0, gamma):
    #     gamma = np.clip(gamma, 1e-12, None)
    #     return baseline + A / (1.0 + ((xv - x0)/gamma)**2)
    
    # # Fit parameters
    # if curve_fit is not None and len(x) >= 4:
    #     try:
    #         p0 = [baseline0, A0, x0_0, max(gamma0, 1e-6)]
    #         bounds = ([-np.inf, -np.inf, np.min(x)-1.0, 1e-9],
    #                   [ np.inf,  np.inf, np.max(x)+1.0, 1.0])  # gamma up to 1 GHz span
    #         popt, _ = curve_fit(lorentz, x, y, p0=p0, bounds=bounds, maxfev=10000)
    #         baseline_fit, A_fit, x0_fit, gamma_fit = [float(v) for v in popt]
    #     except Exception:
    #         baseline_fit, A_fit, x0_fit, gamma_fit = baseline0, A0, x0_0, max(gamma0, 1e-6)
    # else:
    #     baseline_fit, A_fit, x0_fit, gamma_fit = baseline0, A0, x0_0, max(gamma0, 1e-6)
    
    # # Fitted curve
    # x_fine = np.linspace(np.min(x), np.max(x), 1000)
    # y_fit = lorentz(x_fine, baseline_fit, A_fit, x0_fit, gamma_fit)
    
    # # Reportables
    # f0_GHz   = x0_fit
    # FWHM_MHz = (2.0 * gamma_fit) * 1e3  # GHz -> MHz
    
    # # ------------------- Plot: one figure, two subplots -------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, num=1, figsize=(7, 7), sharex=True)
    fig.clf()
    fig, (ax1, ax2) = plt.subplots(2, 1, num=1, figsize=(7, 7), sharex=True)
    fig.suptitle(filename)  # file name as suptitle
    
    # if exp_settings['fit'] == True:
    
    #     # Subplot 1: magnitude + fit
    #     ax1.plot(x, y, '.', label='Data')
    #     ax1.plot(x_fine, y_fit, '-', label='Lorentzian fit')
    #     ax1.set_ylabel('Amplitude')
    #     peak_or_dip = 'Peak' if is_peak else 'Dip'
    #     ax1.set_title(f'Lorentzian {peak_or_dip} — f0 = {f0_GHz:.6f} GHz, FWHM = {FWHM_MHz:.2f} MHz')
    #     ax1.legend(loc='best')
    #     ax1.grid()
        
    #     # Subplot 2: phase
    #     ax2.plot(x, phase_deg, '.')
    #     ax2.set_xlabel('Frequency (GHz)')
    #     ax2.set_ylabel('Phase (deg)')
    #     ax2.grid()
        
    #     fig.tight_layout(rect=[0, 0, 1, 0.95])  # leave space for suptitle
    #     plt.savefig(os.path.join(saveDir, filename+'_mag_phase_lorentz.png'), dpi=150)
        
    # else:
    # Subplot 1: magnitude + fit
    ax1.plot(full_data['xaxis'], full_data['mags'], label='Data')
    ax1.set_ylabel('Amplitude')
    #ax1.legend(loc='best')
    ax1.grid()
    
    # Subplot 2: phase
    ax2.plot(full_data['xaxis'], full_data['phases'], '.')
    ax2.set_xlabel('Hold Time (us)')
    ax2.set_ylabel('Phase (deg)')
    ax2.grid()
    
    fig.tight_layout(rect=[0, 0, 1, 0.95])  # leave space for suptitle
    plt.savefig(os.path.join(saveDir, filename+'_mag_phase.png'), dpi=150)
    # # =====================================================================

    
    
    
    userfuncs.SaveFull(saveDir, filename, ['hold_times','full_data','filename'],
    locals(), expsettings=settings, instruments={})
    


    t_f    = time.time()
    t_single = t_f - t_i
    
    print("Elapsed Time: " + str(t_single))

    data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data}

    return data,prog
