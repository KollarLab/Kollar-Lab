# -*- coding: utf-8 -*-
# %%
"""
Created on Wed Jun 18 16:20:11 2025

@author: jhyang
"""

# Import the QICK drivers and auxiliary libraries
from qick.AveragerProgram import AveragerProgram
#from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import numpy as np
#from scipy.stats import linregress
import userfuncs
import os
import time
from utility.measurement_helpers import estimate_time
from utility.plotting_tools import simplescan_plot

# run initialize_hardware_FPGA for these two
# soc = QickSoc()
# soccfg = soc

class LoopbackProgram(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=1) # nqz zone fixed

        #configure the readout lengths and downconversion frequencies
        self.declare_readout(ch=cfg["ro_channels"][0], length=self.cfg["meas_window"],
                             freq=self.cfg["pulse_freq"], gen_ch=cfg["cav_channel"])

        freq=self.freq2reg(cfg["pulse_freq"], gen_ch=cfg["cav_channel"], 
                                 ro_channel=cfg["ro_channels"][0])
        self.set_pulse_registers(ch=cfg["cav_channel"], style="const", freq=freq,
                                 # converts phase degrees to QICK register val
                                 phase=self.deg2reg(cfg["cav_phase"]), 
                                 gain=cfg["pulse_gain"], 
                                 length=cfg["pulse_length"], mode = "periodic")
        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        self.measure(pulse_ch=self.cfg["cav_channel"], 
             adcs=[self.cfg["ro_channels"]],
             adc_trig_offset=self.cfg["adc_trig_offset"],
             wait=True,
             syncdelay=self.us2cycles(self.cfg["relax_delay"]))

# refactor this to include all settings
# study how settings and global settings are written
def get_CW_trans_settings(): #Default settings dictionary
    settings = {}
    
    settings['scanname'] = 'CW Transmission Scan'
    
    #Sweep parameters
    settings['freq_start']   = 0
    settings['freq_step']    = 2000
    settings['freq_points']  = 1e6
    
    settings['gain_start']     = 4000
    settings['gain_step']      = 2000
    settings['gain_points']    = 1

    #Card settings
    settings['reps'] = 1
    settings['averages'] = 1
    
    # CW settings
    settings['meas_window'] = 900
    settings['initial_phase']   = 0.2 # 0.13 # rad
    
    return settings

def CW_trans(soc, soccfg, instruments, settings):
    # is gain step necessary (like in pulsed_trans)?
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    
    f0_start = settings['freq_start']
    f0_step = settings['freq_step']
    expts = settings['freq_points']
    
    fpts = np.arange(0,expts)*f0_step+f0_start
    soccfg.adcfreq(fpts, gen_ch=exp_globals['cav_channel'], 
                          ro_channel=exp_globals['ro_channels'])
    
    gpts = exp_settings["gain_start"] + exp_settings["gain_step"] * np.arange(exp_settings["gain_points"])
    
    config={"cav_channel":exp_globals['cav_channel'], # --Fixed
            "ro_channels":exp_globals['ro_channels'], # --Fixed
            "relax_delay":m_pulse['relax_delay'], # --Fixed /was 1
            "cav_phase":0, # PLACEHOLDER, loops
            "pulse_style": "const", # --Fixed
            "pulse_length":10, # Fixed, MODE - PERIODIC
            "meas_window":settings['meas_window'], # CW: long readout length
            "pulse_gain":25, # PLACEHOLDER, loops
            "pulse_freq": 75, # PLACEHOLDER, loops [Hz]
            "adc_trig_offset": m_pulse['emp_delay'], # [Clock ticks] /was 100
            "reps":1, 
            "soft_avgs":1,
           }
    
    end_freq = settings["freq_start"] + settings["freq_step"]*settings["freq_points"]
    
    
    slope = settings['initial_phase'] # rad
    phase_slope = (slope) * (180/np.pi) # degrees
    
    # %%
    ##########################################################
    # run this cell again to obtain the phase-corrected output
    ##########################################################
    def freq2phase(fpts):
        return (settings["freq_start"] - fpts) * phase_slope # note the negative sign!
    
    #results=[]
    phase_fpts=freq2phase(fpts)
    phase = iter(phase_fpts)
    
    powerdat = np.zeros((len(gpts), len(fpts)))
    phasedat = np.zeros((len(gpts), len(fpts)))
    Is = np.zeros((len(gpts), len(fpts)))
    Qs = np.zeros((len(gpts), len(fpts)))
    tstart = time.time()
    
    for g in gpts:
        print("Current Gain: " + str(gpts[g]) + ", Max Gain: " + str(gpts[-1]))
        config["pulse_gain"] = gpts[g]
        
        for f in fpts:
            config["pulse_freq"]=int(f/1e6) # convert to MHz
            config["cav_phase"] = next(phase)
            prog =LoopbackProgram(soccfg, config)
            trans_I, trans_Q = prog.acquire(soc,load_pulses=True,progress=False)
            mag = np.sqrt(trans_I[0][0]**2 + trans_Q[0][0]**2)
            phase = np.arctan2(trans_Q[0][0], trans_I[0][0])*180/np.pi

            powerdat[g,f] = mag#/gpts[g] #Taken out normalization for now
            phasedat[g,f] = phase
            Is[g,f] = trans_I[0][0]
            Qs[g,f] = trans_Q[0][0]
            #results.append(prog.acquire(soc))
            
            if g == 0:
                tstop = time.time()
                estimate_time(tstart, tstop, len(gpts))
    #results=np.transpose(results)
    
    #%%
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp
    
    '''
    #######
    DEMO VERSION
    #######
    
    sig = results[0][0][0] + 1j*results[0][0][1]
    avgamp0 = np.abs(sig)
    plt.figure(2)
    plt.plot(fpts, results[0][0][0],label="I value")
    plt.plot(fpts, results[0][0][1],label="Q value")
    plt.plot(fpts, avgamp0,label="Amplitude")
    plt.ylabel("a.u.")
    plt.xlabel("Pulse Freq")
    plt.title("%d-%d MHz Frequency Sweep"%(settings["freq_start"], end_freq))
    plt.legend()
    plt.savefig("images/Freq_sweep_python.pdf", dpi=350)
    
    plt.figure(3)
    plt.plot(fpts, np.arctan2((results[0][0][0]),(results[0][0][1])))
    # about 6 phi / 30 MHz
    
    # if slope is set to 0.2 rad and actual is 0.13 acquired result will be -0.07
    slope += linregress(fpts, np.unwrap
                    (np.arctan2((results[0][0][0]),(results[0][0][1]))))[0]
    print("Slope: ")
    print(slope)
    '''
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
    return full_data
    
    
    
    
    
    
