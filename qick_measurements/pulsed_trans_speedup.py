from qick.averager_program import AveragerProgram

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import estimate_time
import time
from scipy.optimize import curve_fit


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
        #self.reset_phase(gen_ch=self.cfg['cav_channel'],t=0)
        
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
        self.sync_all(self.us2cycles(self.cfg["relax_delay"])) #Syncs to an offset time after the final pulse is sent

def get_trans_settings(): #Default settings dictionary
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'PulsedTrans'
    
    #Sweep parameters
    settings['freq_start']   = 7e9 
    settings['freq_step']    = 100e6
    settings['freq_points']  = 21

    settings['gain_start']  = 500
    settings['gain_step']   = 100
    settings['gain_points'] = 31

    #Card settings
    settings['reps'] = 1
    settings['averages'] = 1e3
    
    return settings

def lorentzian_peak(x, A, f0, gamma, offset):
    """Peak: offset + A * (gamma/2)^2 / ((x-f0)^2 + (gamma/2)^2)"""
    return offset + A * (0.5*gamma)**2 / ((x - f0)**2 + (0.5*gamma)**2)

def lorentzian_dip(x, A, f0, gamma, offset):
    """Dip: offset - A * (gamma/2)^2 / ((x-f0)^2 + (gamma/2)^2)"""
    return offset - A * (0.5*gamma)**2 / ((x - f0)**2 + (0.5*gamma)**2)


def pulsed_trans(soc,soccfg,instruments,settings): #Main measurement function
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    
    #if exp_globals['LO']: # Ignore for now, we no longer use an external LO
    if False:
        logen = instruments['LO']
        
        logen.freq   = exp_globals['LO_freq']
        logen.power  = exp_globals['LO_power']
        logen.output = 1

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel'],
        'ro_channels'     : exp_globals['ro_channels'],

        'nqz_c'           : 1, #irrelevant setting
        'cav_phase'       : m_pulse['cav_phase'], #degrees
        'meas_window'     : m_pulse['meas_window'], #us
        'meas_time'       : m_pulse['meas_pos'], #us
        'meas_gain'       : 500, #Placeholder, gets overwritten in sweep
        'cav_freq'        : 100, #Placeholder
        
        'readout_length'  : m_pulse['meas_window'], #us
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'], #us


        'relax_delay'     : exp_globals['relax_delay'], #us
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs']
        }

    fpts = exp_settings["freq_start"]+exp_settings["freq_step"]*np.arange(exp_settings["freq_points"]) # Defines frequency sweep
    gpts = exp_settings["gain_start"]+exp_settings["gain_step"]*np.arange(exp_settings["gain_points"]) # Defines gain sweep

    phases = (fpts-fpts[0])/1e3 * exp_settings['dphi_df'] # Current solution for the phase wrapping, still need function to extract dphi/df out of a calibration measurement instead of manual
   

    powerdat = np.zeros((len(gpts), len(fpts))) #Initializing empty matrices for data
    phasedat = np.zeros((len(gpts), len(fpts)))
    Is = np.zeros((len(gpts), len(fpts)))
    Qs = np.zeros((len(gpts), len(fpts)))

    #drive_powers_lin = 10**(powers/10) ?
    #drive_amps_lin = np.sqrt(drive_powers_lin)
    
    tstart = time.time()

    for g in range(0,len(gpts)): # Loop through gain values
        print("Current Gain: " + str(gpts[g]) + ", Max Gain: " + str(gpts[-1]))
        
        config["meas_gain"] = gpts[g] 
        
        for f in range(0,len(fpts)): # Loop through freq values
            board_freq = (fpts[f] - exp_globals['LO_freq']*exp_globals['LO'])/1e6
            config["cav_freq"] = board_freq
            config['cav_phase'] = phases[f] # Current solution for the phase wrapping
            prog = CavitySweep(soccfg,config)

            trans_I, trans_Q = prog.acquire(soc,load_pulses=True,progress=False) #Transmission data acquisition occurs here

            mag = np.sqrt(trans_I[0][0]**2 + trans_Q[0][0]**2)
            phase = np.arctan2(trans_Q[0][0], trans_I[0][0])*180/np.pi

            powerdat[g,f] = mag/gpts[g] #bring in the normalization
            phasedat[g,f] = phase
            Is[g,f] = trans_I[0][0]
            Qs[g,f] = trans_Q[0][0]
        
        if g == 0:
            tstop = time.time()
            estimate_time(tstart, tstop, len(gpts))


        full_data = {}
        full_data['xaxis']  = fpts/1e9 # Changing xaxis from Hz to GHz
        full_data['mags']   = powerdat[0:g+1]
        full_data['phases'] = phasedat[0:g+1]
        full_data['Is']     = Is[0:g+1]
        full_data['Qs']     = Qs[0:g+1]


        plot_data = {}
        plot_data['xaxis']  = fpts/1e9
        plot_data['mags']   = powerdat[0:g+1]
        plot_data['phases'] = phasedat[0:g+1]

        single_data = {}
        single_data['xaxis'] = fpts/1e9
        single_data['mag']   = powerdat[g]
        single_data['phase'] = phasedat[g]

        yaxis  = gpts[0:g+1] #- CAV_Attenuation
        labels = ['Freq (GHz)', 'Gain (DAC a.u.)']

        
        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=1, IQdata = False) 
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)


        userfuncs.SaveFull(saveDir, filename, ['gpts','fpts', 'powerdat', 'phasedat','Is','Qs','full_data', 'single_data'],
                             locals(), expsettings=settings, instruments={})
    
    # === Fig. 2: Last-trace transmission magnitude + Lorentzian fit ===
    x = fpts / 1e9                         # GHz
    y = powerdat[-1]                       # normalized magnitude for last gain
    hanger = bool(exp_globals.get('hanger', False))
    
    # Choose model depending on hanger (dip vs. peak)
    model = lorentzian_dip if hanger else lorentzian_peak
    
    # ---- Initial parameter guesses ----
    offset0 = np.median(y)
    idx0 = np.argmin(y) if hanger else np.argmax(y)
    f0_0 = x[idx0]
    step = np.mean(np.diff(x)) if len(x) > 1 else 0.001
    gamma0 = max(5 * step, 0.001)
    A0 = max((offset0 - y[idx0]) if hanger else (y[idx0] - offset0), 1e-6)
    
    p0 = [A0, f0_0, gamma0, offset0]
    bounds = ([0, x.min(), 0, -np.inf], [np.inf, x.max(), np.inf, np.inf])
    
    try:
        popt, pcov = curve_fit(model, x, y, p0=p0, bounds=bounds, maxfev=10000)
    except Exception:
        popt, pcov = p0, None
    
    A, f0, gamma, offset = popt
    FWHM_MHz = (2.0 * gamma) * 1e3
    
    # ---- Plot (Figure 2) ----
    fig2, ax2 = plt.subplots(figsize=(7, 4.5), num=2, clear=True)
    fig2.suptitle(filename, fontsize=11, y=0.98)  # file name as suptitle
    
    ax2.plot(x, y, 'o', ms=4, label='Data (last gain)')
    xf = np.linspace(x.min(), x.max(), 1000)
    ax2.plot(xf, model(xf, *popt), '-', lw=2, label='Lorentzian fit')
    
    ax2.set_xlabel('Frequency (GHz)')
    ax2.set_ylabel('Gain normalized |S21| (a.u.)')
    ax2.legend()
    ax2.grid(alpha=0.3)
    
    kind = 'dip' if hanger else 'peak'
    ax2.set_title(
        f"Lorentzian {kind} fit: "
        f"$f_0$ = {f0:.7f} GHz,  FWHM = {FWHM_MHz:.2f} MHz,  offset = {offset:.4g}",
        fontsize=10
    )
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(saveDir, filename + '_lasttrace_lorentz_fit.png'), dpi=150)




    if exp_globals['LO']:
        pass
        #logen.output = 0
        
    return full_data



