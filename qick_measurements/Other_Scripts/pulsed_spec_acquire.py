# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: kollarlab
"""

from qick.averager_program import AveragerProgram

import numpy as np
import os
import matplotlib.pyplot as plt
import time

import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import estimate_time
      
class PulsedSpec(AveragerProgram):
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
        phase_q = self.deg2reg(cfg["qub_phase"], gen_ch=gen_ch)
        gain_q  = cfg["qub_gain"]
        
        self.default_pulse_registers(ch=qub_ch, phase=phase_q, freq=freq_q, gain=gain_q, mode = "oneshot")
        
        sigma = self.us2cycles(cfg["qub_sigma"],gen_ch=qub_ch)
        num_sigma = cfg["num_sigma"]
        
        self.add_gauss(ch=qub_ch, name="ex", sigma=sigma,length=int(sigma*num_sigma))        
        self.set_pulse_registers(ch=qub_ch, style="arb", waveform="ex")

        
        # give processor some time to configure pulses, I believe it sets the offset time 200 cycles into the future?
        # Try varying this and seeing if it moves. Also try putting a synci after the trigger in the body
        self.synci(200)   
    
    def body(self):
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.
        if self.cfg["phase_reset"]:
            self.reset_phase(gen_ch = [self.cfg['cav_channel'], self.cfg['qub_channel']], t=0)
        else:
            pass
        
        sigma = self.us2cycles(self.cfg["qub_sigma"])
        num_sigma = self.cfg["num_sigma"]
        
        offset = self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"])
        meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=self.cfg["cav_channel"])
        ex_time = meas_time - self.us2cycles(self.cfg['qub_delay'],gen_ch=self.cfg["qub_channel"]) - int(num_sigma*sigma)
        #Sets off the ADC
        self.trigger(adcs=self.ro_chs,
                    pins=[0],
                    adc_trig_offset=offset)
        
        #Sends measurement pulse
        self.pulse(ch=self.cfg["qub_channel"],t=ex_time)
        self.pulse(ch=self.cfg["cav_channel"],t=meas_time)
        self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        self.sync_all(self.us2cycles(self.cfg["relax_delay"])) #Syncs to an offset time after the final pulse is sent


def get_spec_settings():
    settings = {}
    
    settings['scanname']    = 'initial_power_scan_q4'
    settings['meas_type']   = 'PulsedSpec'
    settings['phase_reset'] = True
    
    settings['cav_freq']    = 1e9
    settings['cav_gain']    = 1000
    
    #Sweep parameters
    settings['freq_start']  = 4e9  
    settings['freq_step']   = 100e6
    settings['freq_points'] = 6

    settings['gain_start']  = 500
    settings['gain_step']   = 100
    settings['gain_points'] = 11
    #Card settings
    settings['reps']        = 1
    settings['soft_avgs']   = 5e3
    
    return settings

def pulsed_spec(soc,soccfg,instruments,settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    
    soc.reset_gens()

    lo_freq = exp_globals["LO_freq"]
    
    #if exp_globals['LO']:
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
        'qub_freq'        : 4000, #Placeholder
        'qub_gain'        : 4000, #Placeholder
        'qub_sigma'        : q_pulse['sigma'],
        'qub_delay'        : q_pulse['delay'],
        'num_sigma'       : q_pulse['num_sigma'],
        
        'readout_length'  : m_pulse['init_buffer'] + m_pulse['meas_window'] + m_pulse['post_buffer'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'] - m_pulse['init_buffer'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs'],

        'phase_reset'     : exp_settings['phase_reset']
        }

    fpts = exp_settings["freq_start"]+exp_settings["freq_step"]*np.arange(exp_settings["freq_points"])
    gpts = exp_settings["gain_start"]+exp_settings["gain_step"]*np.arange(exp_settings["gain_points"])

    prog = PulsedSpec(soccfg,config)
    meas_start = prog.us2cycles(m_pulse["init_buffer"],ro_ch=0)
    meas_end = meas_start+prog.us2cycles(m_pulse["meas_window"],ro_ch=0)
   

    powerdat = np.zeros((len(gpts), len(fpts)))
    phasedat = np.zeros((len(gpts), len(fpts)))

    #drive_powers_lin = 10**(powers/10) ?
    #drive_amps_lin = np.sqrt(drive_powers_lin)

    tstart = time.time()

    for g in range(0,len(gpts)):
        print("Current Gain: " + str(gpts[g]) + ", Max Gain: " + str(gpts[-1]))
        config["qub_gain"] = gpts[g]

        total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)

        Is  = np.zeros((len(fpts), total_samples))
        Qs  = np.zeros((len(fpts), total_samples))
        
        for f in range(0,len(fpts)):
            board_freq = fpts[f]/1e6
            config["qub_freq"] = board_freq
            prog = PulsedSpec(soccfg,config)
            #Need to assign Iwindow, Qwindow, Ifull, Qfull, xaxis (which should just be timeus)
            holder = prog.acquire(soc, load_pulses=True, progress=False, debug=False)


            I_sig, Q_sig   = [holder[0][0][0], holder[1][0][0]]

            powerdat[g, f] = np.sqrt(I_sig**2 + Q_sig**2)
            phasedat[g, f] = np.arctan2(Q_sig, I_sig)*180/np.pi


            #I_full = holder[0][0]
            #Q_full = holder[0][1]


            #I_window = I_full[meas_start:meas_end]
            #Q_window = Q_full[meas_start:meas_end]
           
            #I_final = np.mean(I_window)
            #Q_final = np.mean(Q_window)

            #Is[f,:]   = I_full
            #Qs[f,:] = Q_full
            #powerdat[g, f] = np.sqrt(I_final**2 + Q_final**2)
            #phasedat[g, f] = np.arctan2(Q_final, I_final)*180/np.pi
        
        if g == 0:
            tstop = time.time()
            estimate_time(tstart, tstop, len(gpts))


            xaxis = []#np.linspace(0,len(I_full)-1,len(I_full))

            for x in range(0,len(xaxis)):
                xaxis = []
                #xaxis[x] = prog.cycles2us(xaxis[x],ro_ch=0)

        full_data = {}
        full_data['xaxis']  = fpts/1e9
        full_data['mags']   = powerdat[0:g+1]
        full_data['phases'] = phasedat[0:g+1]

        #Currently no rescaling is implemented, ask Martin for tips implementing it.
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


        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=1, IQdata = False) #normalized to drive level
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

        full_time = {}
        full_time['xaxis']  = xaxis
        full_time['Is']   = Is
        full_time['Qs'] = Qs

        single_time = {}
        #single_time['xaxis'] = xaxis
        #single_time['I']   = I_full
        #single_time['Q'] = Q_full

        #time_labels = ['Time (us)', 'Freq (GHz)']
        #identifier = 'Gain: {}a.u.'.format(gpts[g])#-CAV_Attenuation)

        #simplescan_plot(full_time, single_time, fpts/1e9, 
        #                'Raw_time_traces\n'+filename, 
        #                time_labels, 
        #                identifier, 
        #                fig_num=2,
        #                IQdata = True)
        #plt.savefig(os.path.join(saveDir, filename+'_RawTimeTraces.png'), dpi = 150)

        userfuncs.SaveFull(saveDir, filename, ['gpts','fpts', 'powerdat', 'phasedat','xaxis','full_data', 'single_data', 'full_time', 'single_time'],
                             locals(), expsettings=settings, instruments={})
    
    tfinish = time.time()
    print(f'Elapsed time: {tfinish-tstart}')

    #if exp_globals['LO']:
    #    logen.output = 0
        
    return full_data