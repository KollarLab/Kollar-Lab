# -*- coding: utf-8 -*-
# %%
"""
Created on Wed Jun 18 16:20:11 2025

@author: jhyang
"""

# Import the QICK drivers and auxiliary libraries
from qick.AveragerProgram import AveragerProgram
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import userfuncs

# run initialize_hardware_FPGA for these two
#soc = QickSoc()
#soccfg = soc

class LoopbackProgram(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=1)

        #configure the readout lengths and downconversion frequencies
        self.declare_readout(ch=cfg["ro_channel"], length=self.cfg["meas_window"],
                             freq=self.cfg["pulse_freq"], gen_ch=cfg["cav_channel"])

        freq=self.freq2reg(cfg["pulse_freq"], gen_ch=cfg["cav_channel"], 
                                 ro_channel=cfg["ro_channel"])
        self.set_pulse_registers(ch=cfg["cav_channel"], style="const", freq=freq,
                                 # converts phase degrees to QICK register val
                                 phase=soccfg.deg2reg(cfg["cav_phase"]), 
                                 gain=cfg["pulse_gain"], 
                                 length=cfg["pulse_length"], mode = "periodic")
        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        self.measure(pulse_ch=self.cfg["cav_channel"], 
             adcs=[self.cfg["ro_channel"]],
             adc_trig_offset=self.cfg["adc_trig_offset"],
             wait=True,
             syncdelay=self.us2cycles(self.cfg["relax_delay"]))

# refactor this to include all settings
# study how settings and global settings are written
def get_CW_trans_settings(): #Default settings dictionary
    settings = {}
    
    settings['scanname'] = 'CW Transmission Scan'
    
    #Sweep parameters
    settings['freq_start']   = 100 
    settings['freq_step']    = 3
    settings['freq_points']  = 100

    #Card settings
    settings['reps'] = 1
    settings['averages'] = 80
    
    return settings

f0_start=100
f0_step=3
expts=100

f0_v = np.arange(0,expts)*f0_step+f0_start

f0_v = soccfg.adcfreq(f0_v, gen_ch=6, ro_channel=0)

config={"cav_channel":6, # --Fixed for loopback
        "ro_channel":0, # --Fixed for loopback
        "relax_delay":0, # --Fixed /was 1
        "cav_phase":0, # PLACEHOLDER
        "pulse_style": "const", # --Fixed
        "pulse_length":10, # mode - periodic   
        "meas_window":900, # Fixed, long readout length
        "pulse_gain":20, # [DAC units]
        "pulse_freq": 100, # [MHz]
        "adc_trig_offset": 0, # [Clock ticks] /was 100
        "reps":1, 
        "soft_avgs":80,
       }

sweep_cfg={"start":100, "step":3, "expts":100}
sweep_cfg["end"] = sweep_cfg["start"] + sweep_cfg["step"]*sweep_cfg["expts"]
gpts=sweep_cfg["start"] + sweep_cfg["step"]*np.arange(sweep_cfg["expts"])


slope = 0.2 # rad
phase_slope = (slope) * (180/np.pi) # degrees

# %%
# run this cell again to obtain the phase-corrected output
def freq2phase(gpts):
    return (sweep_cfg["start"] - gpts) * phase_slope # note the negative sign!

results=[]
phase_gpts=freq2phase(gpts)
phase = iter(phase_gpts)


for g in tqdm(gpts):
    config["pulse_freq"]=int(g)
    config["cav_phase"] = next(phase)
    #print(config["cav_phase"])
    prog =LoopbackProgram(soccfg, config)
    results.append(prog.acquire(soc))
results=np.transpose(results)

#%%
# refactor this bit to match? (ask kellen) other plot data necessary?
# with savedir timestamp etc
stamp    = userfuncs.timestamp()
saveDir  = userfuncs.saveDir(settings)
filename = exp_settings['scanname'] + '_' + stamp

sig = results[0][0][0] + 1j*results[0][0][1]
avgamp0 = np.abs(sig)
plt.figure(2)
plt.plot(gpts, results[0][0][0],label="I value")
plt.plot(gpts, results[0][0][1],label="Q value")
plt.plot(gpts, avgamp0,label="Amplitude")
plt.ylabel("a.u.")
plt.xlabel("Pulse Freq")
plt.title("%d-%d MHz Frequency Sweep"%(sweep_cfg["start"], sweep_cfg["end"]))
plt.legend()
plt.savefig("images/Freq_sweep_python.pdf", dpi=350)

plt.figure(3)
plt.plot(gpts, np.arctan2((results[0][0][0]),(results[0][0][1])))
# about 6 phi / 30 MHz

# if slope is set to 0.2 rad and actual is 0.13 acquired result will be -0.07
slope += linregress(gpts, np.unwrap
                (np.arctan2((results[0][0][0]),(results[0][0][1]))))[0]
print("Slope: ")
print(slope)






