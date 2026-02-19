# -*- coding: utf-8 -*-
# %%
"""
Created on Wed Jun 18 16:20:11 2025

@author: Kollarlab
"""

from qick.asm_v2 import AveragerProgramV2
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import userfuncs
import os
import time
import logging
from utility.measurement_helpers import estimate_time
from utility.plotting_tools import simplescan_plot

class LoopbackProgramV2(AveragerProgramV2):
    def _initialize(self, cfg):

        ro_ch = cfg['ro_channel']
        gen_ch = cfg['cav_channel']
        self.declare_gen(ch=gen_ch, nqz=cfg['nqz_c'], mixer_freq=cfg['mixer_freq'], ro_ch=ro_ch)
        #for ro_ch in cfg["ro_channels"]:
        self.declare_readout(ch=ro_ch, length=cfg['readout_length'])
        self.add_readoutconfig(ch=ro_ch, name="myro",
                           freq=cfg['cav_freq'],
                           gen_ch=gen_ch,
                           outsel='product')
        self.add_pulse(ch=gen_ch, name="CW_pulse", ro_ch=ro_ch,
                       style="const",
                       freq=cfg['cav_freq'],
                       length= cfg["pulse_length"],
                       phase=cfg['cav_phase'],
                       gain=cfg['meas_gain'],
                       mode='periodic'
                      )
        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)
        
# =============================================================================
#     def initialize(self):
#         cfg=self.cfg   
# 
#         # set the nyquist zone
#         self.declare_gen(ch=cfg["cav_channel"], nqz=1) # nqz zone fixed
# 
#         #configure the readout lengths and downconversion frequencies
#         readout = self.us2cycles(cfg["meas_window"],ro_ch=cfg["ro_channels"][0])
#         self.declare_readout(ch=cfg["ro_channels"][0], length=readout,
#                              freq=self.cfg["pulse_freq"], gen_ch=cfg["cav_channel"])
# 
#         freq=self.freq2reg(cfg["pulse_freq"], gen_ch=cfg["cav_channel"], 
#                                  ro_ch=cfg["ro_channels"][0])
#         self.set_pulse_registers(ch=cfg["cav_channel"], style="const", freq=freq,
#                                  # converts phase degrees to QICK register val
#                                  phase=self.deg2reg(cfg["cav_phase"]), 
#                                  gain=cfg["pulse_gain"], 
#                                  length=self.us2cycles(cfg["pulse_length"],gen_ch=self.cfg["cav_channel"]), mode = "periodic")
#         self.synci(200)  # give processor some time to configure pulses
# =============================================================================

    def _body(self, cfg):

        self.pulse(ch=self.cfg['cav_channel'], name="CW_pulse", t=0.0)
        self.trigger(ros=[self.cfg['ro_channel']], pins=[0], t=self.cfg['adc_trig_offset'])
        self.wait_auto()
        #self.delay(self.cfg["relax_delay"]) # no delay needed here for CW measurements
    
# =============================================================================
#     def body(self):
#         self.measure(pulse_ch=self.cfg["cav_channel"], 
#              adcs=[self.cfg["ro_channels"][0]],
#              adc_trig_offset=self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"]),
#              wait=True,
#              syncdelay=self.us2cycles(self.cfg["relax_delay"],gen_ch=self.cfg["cav_channel"]))
# =============================================================================


def get_CW_trans_settings(): #Default settings dictionary
    settings = {}
    
    settings['scanname'] = 'CW Transmission Scan'
    settings['meas_type'] = 'CWTrans'
    
    #Sweep parameters
    settings['freq_start']   = 4000
    settings['freq_step']    = 5000
    settings['freq_points']  = 1000
    
    settings['gain_start']     = 0.5
    settings['gain_step']      = 0.1
    settings['gain_points']    = 1

    #Card settings
    settings['reps'] = 10
    settings['rounds'] = 1
    
    # CW settings
    settings['readout_length'] = 900
    settings['initial_phase']   = 0.2 # 0.13 # rad
    
    return settings


def CW_trans(soc, soccfg, instruments, settings):
    
    # suppresses the sum buffer overflow warning
    logging.getLogger("qick").setLevel(logging.ERROR)
    
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    
    f0_start = exp_settings['freq_start']
    f0_step = exp_settings['freq_step']
    expts = exp_settings['freq_points']
    
    fpts = np.arange(0,expts)*f0_step+f0_start
    end_freq = exp_settings["freq_start"] + exp_settings["freq_step"]*exp_settings["freq_points"]
    
    if exp_settings['power_sweep_mode'] == 'gain':
        gpts = exp_settings["gain_start"] + (exp_settings["gain_step"] * np.arange(exp_settings["gain_points"]))
        
    
        config={"cav_channel":exp_globals['cav_channel']['ID'], # --Fixed
                "ro_channel":exp_globals['ro_channel']['ID'], # --Fixed
                "relax_delay":exp_globals['relax_delay'], # --Fixed /was 1
                "cav_phase":0, # PLACEHOLDER, loops
                "pulse_style": "const", # --Fixed
                "pulse_length":10, # Fixed, MODE - PERIODIC
                "meas_gain":0.5, # PLACEHOLDER, loops
                "cav_freq":4500, # PLACEHOLDER, loops [MHz]
                "adc_trig_offset": m_pulse['emp_delay'],#+ m_pulse['meas_pos'], # us, maybe in CW the meas_pos can be effectively set to 0?
                "reps":exp_settings['reps'],
                "rounds":exp_settings['rounds'],
                "nqz_c": 2,
                'mixer_freq': (exp_settings["freq_start"] + exp_settings['mixer_detuning'])/1e6, # in MHz
                'readout_length': exp_settings['readout_length'], # this is not the same as pulse_length; pulse_length has a limit, so we keep the readout long
                'cav_atten'       : exp_globals['cav_channel']['Atten_1'] + exp_globals['cav_channel']['Atten_2'],
                'power_sweep_mode': exp_settings['power_sweep_mode']
               }
        
        cav_ch = exp_globals['cav_channel']
        ro_ch  = exp_globals['ro_channel']
        # Set attenuator on DAC.
        soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
        # Set attenuator on ADC.
        soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])
    
        
        powerdat = np.zeros((len(gpts), len(fpts)))
        phasedat = np.zeros((len(gpts), len(fpts)))
        Is = np.zeros((len(gpts), len(fpts)))
        Qs = np.zeros((len(gpts), len(fpts)))
        
        
        
        tstart = time.time()
        
        # gain loop
        for g in range(0,len(gpts)):
            print("Current Gain: " + str(gpts[g]) + ", Max Gain: " + str(gpts[-1]))
            config["meas_gain"] = gpts[g]
            
            # frequency sweep
            for f in range(0,len(fpts)):
                config["cav_freq"]=fpts[f]/1e6 # convert to MHz
                
                board_freq = config['cav_freq']
                
                if exp_settings['filter'] == 'all_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
                    
                elif exp_settings['filter'] == 'cav_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bypass')
                    
                elif exp_settings['filter'] == 'ro_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bypass')
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
                    
                elif exp_settings['filter'] == 'no_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bypass')
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bypass')    
                    
                else:
                    print('Please select one option from:')
                    print('\'all_filter\', \'cav_filter\', \'ro_filter\', and\'no_filter\'')
                    return
                
                if exp_settings['fixed_mixer_frequency'] == False:
                    config['mixer_freq'] = (fpts[f] + exp_settings['mixer_detuning'])/1e6 # in MHz
                #config["cav_phase"] = phase_fpts[f] # currently this feature is not used
                #print(config)
                prog =LoopbackProgramV2(soccfg, reps=exp_settings['reps'], final_delay = None, final_wait=0, cfg=config)
                holder = prog.acquire(soc, rounds=int(config["rounds"]), progress=False) #acquire returns (shape) = (n_ro_ch, n_reads, 2)
                #n_reads is set by how many times adc is triggered in the _body, and the last index is I and Q data
                iq = holder[0]
                soc.reset_gens() #Reset the tProc and run a minimal tProc program that drives all signal generators with 0's. Useful for stopping any periodic or stdysel="last" outputs
                # that may have been driven by a previous program
                
                I_full = iq[:,0]
                Q_full = iq[:,1]
                
                mag = np.sqrt(I_full**2 + Q_full**2)
                phase = np.degrees(np.arctan2(Q_full, I_full))
    
                powerdat[g,f] = mag/gpts[g] #Normalization
                phasedat[g,f] = phase
                Is[g,f] = I_full
                Qs[g,f] = Q_full
                
                    
            if g == 0:
                tstop = time.time()
                estimate_time(tstart, tstop, len(gpts))
                
        #%%
        stamp    = userfuncs.timestamp()
        saveDir  = userfuncs.saveDir(settings)
        filename = exp_settings['scanname'] + '_' + stamp
        
        # if slope is set to 0.2 rad and actual is 0.13 acquired result will be -0.07
        #exp_settings['initial_phase'] += linregress(fpts, np.unwrap
                        #(np.arctan2((results[0][0][0]),(results[0][0][1]))))[0]
    # =============================================================================
    #     print("Slope: ")
    #     print(linregress(fpts, np.unwrap
    #                     (np.arctan2((phase_results[0][0][0]),(phase_results[0][0][1]))))[0])
    # =============================================================================
        
        
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
        identifier = f"Atten={config['cav_atten']} dB"
    
        
        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=1, IQdata = False) 
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
    
    
        userfuncs.SaveFull(saveDir, filename, ['gpts','fpts', 'powerdat', 'phasedat','Is','Qs','full_data', 'single_data'],
                             locals(), expsettings=settings, instruments={})
        return full_data
    
    elif exp_settings['power_sweep_mode'] == 'attenuation':
        apts = exp_settings['atten_start'] + (exp_settings["atten_step"] * np.arange(exp_settings["atten_points"]))
        
        config={"cav_channel":exp_globals['cav_channel']['ID'], # --Fixed
                "ro_channel":exp_globals['ro_channel']['ID'], # --Fixed
                "relax_delay":exp_globals['relax_delay'], # --Fixed /was 1
                "cav_phase":0, # PLACEHOLDER, loops
                "pulse_style": "const", # --Fixed
                "pulse_length":10, # Fixed, MODE - PERIODIC
                "meas_gain":exp_settings['meas_gain'],
                "cav_freq":4500, # PLACEHOLDER, loops [MHz]
                "adc_trig_offset": m_pulse['emp_delay'],#+ m_pulse['meas_pos'], # us, maybe in CW the meas_pos can be effectively set to 0?
                "reps":exp_settings['reps'],
                "rounds":exp_settings['rounds'],
                "nqz_c": 2,
                'mixer_freq': (exp_settings["freq_start"] + exp_settings['mixer_detuning'])/1e6, # in MHz
                'readout_length': exp_settings['readout_length'], # this is not the same as pulse_length; pulse_length has a limit, so we keep the readout long
                'cav_atten'       : 40,#place holder, updated in the loop
                'power_sweep_mode': exp_settings['power_sweep_mode']
               }
        
        cav_ch = exp_globals['cav_channel']
        ro_ch  = exp_globals['ro_channel']
        
        # Set attenuator on ADC.
        soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])
    
        
        powerdat = np.zeros((len(apts), len(fpts)))
        phasedat = np.zeros((len(apts), len(fpts)))
        Is = np.zeros((len(apts), len(fpts)))
        Qs = np.zeros((len(apts), len(fpts)))
        
        
        
        tstart = time.time()
        
        # attenuation loop
        for a in range(0,len(apts)):
            print("Current Attenuation: " + str(apts[a]) + " dB, Final Attenuation: " + str(apts[-1]) + ' dB')
            # Set generator attenuation in sweep
            soc.rfb_set_gen_rf(cav_ch['ID'], apts[a] / 2, apts[a] / 2)
            
            # frequency sweep
            for f in range(0,len(fpts)):
                config["cav_freq"]=fpts[f]/1e6 # convert to MHz
                
                board_freq = config['cav_freq']
                
                if exp_settings['filter'] == 'all_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
                    
                elif exp_settings['filter'] == 'cav_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bypass')
                    
                elif exp_settings['filter'] == 'ro_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bypass')
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
                    
                elif exp_settings['filter'] == 'no_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bypass')
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bypass')    
                    
                else:
                    print('Please select one option from:')
                    print('\'all_filter\', \'cav_filter\', \'ro_filter\', and\'no_filter\'')
                    return
                
                if exp_settings['fixed_mixer_frequency'] == False:
                    config['mixer_freq'] = (fpts[f] + exp_settings['mixer_detuning'])/1e6 # in MHz
                #config["cav_phase"] = phase_fpts[f] # currently this feature is not used
                #print(config)
                prog =LoopbackProgramV2(soccfg, reps=exp_settings['reps'], final_delay = None, final_wait=0, cfg=config)
                holder = prog.acquire(soc, rounds=int(config["rounds"]), progress=False) #acquire returns (shape) = (n_ro_ch, n_reads, 2)
                #n_reads is set by how many times adc is triggered in the _body, and the last index is I and Q data
                iq = holder[0]
                soc.reset_gens() #Reset the tProc and run a minimal tProc program that drives all signal generators with 0's. Useful for stopping any periodic or stdysel="last" outputs
                # that may have been driven by a previous program
                
                I_full = iq[:,0]
                Q_full = iq[:,1]
                
                mag = np.sqrt(I_full**2 + Q_full**2)
                phase = np.degrees(np.arctan2(Q_full, I_full))
    
                powerdat[a,f] = mag
                phasedat[a,f] = phase
                Is[a,f] = I_full
                Qs[a,f] = Q_full
                
                    
            if a == 0:
                tstop = time.time()
                estimate_time(tstart, tstop, len(apts))
                
        #%%
        stamp    = userfuncs.timestamp()
        saveDir  = userfuncs.saveDir(settings)
        filename = exp_settings['scanname'] + '_' + stamp
        
        # if slope is set to 0.2 rad and actual is 0.13 acquired result will be -0.07
        #exp_settings['initial_phase'] += linregress(fpts, np.unwrap
                        #(np.arctan2((results[0][0][0]),(results[0][0][1]))))[0]
    # =============================================================================
    #     print("Slope: ")
    #     print(linregress(fpts, np.unwrap
    #                     (np.arctan2((phase_results[0][0][0]),(phase_results[0][0][1]))))[0])
    # =============================================================================
        
        
        full_data = {}
        full_data['xaxis']  = fpts/1e9 # Changing xaxis from Hz to GHz
        full_data['mags']   = powerdat[0:a+1]
        full_data['phases'] = phasedat[0:a+1]
        full_data['Is']     = Is[0:a+1]
        full_data['Qs']     = Qs[0:a+1]
    
    
        plot_data = {}
        plot_data['xaxis']  = fpts/1e9
        plot_data['mags']   = powerdat[0:a+1]
        plot_data['phases'] = phasedat[0:a+1]
    
        single_data = {}
        single_data['xaxis'] = fpts/1e9
        single_data['mag']   = powerdat[a]
        single_data['phase'] = phasedat[a]
    
        yaxis  = apts[0:a+1] #- CAV_Attenuation
        labels = ['Freq (GHz)', 'Attenuation (dB)']
        identifier = f"Gain={config['meas_gain']}"
    
        
        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=1, IQdata = False) 
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)
    
    
        userfuncs.SaveFull(saveDir, filename, ['apts','fpts', 'powerdat', 'phasedat','Is','Qs','full_data', 'single_data'],
                             locals(), expsettings=settings, instruments={})
        return full_data
        
    else:
        raise ValueError("Please set power_sweep_mode to 'gain' or 'attenuation'")
    
    
    
    
    
    
