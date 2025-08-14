# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: jhyang

"""

from qick.averager_program import AveragerProgram


import numpy as np
import time
import userfuncs
import os
from utility.measurement_helpers import estimate_time
import logging
import matplotlib.pyplot as plt
import utility.plotting_tools as plots
from utility.userfits_v2 import fit_model

#######################################
# Taken from CW_trans LoopBackProgram #
#######################################
#UNCOMMENT SRS SECTIONS WHEN TESTING! #
#PLOTTING MAY HAVE SOME PROBLEMS!     #
#######################################

# Plays a constant tone, swept under cw_spec_flux()
class CavitySweepFlux(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=1) # nqz zone fixed

        #configure the readout lengths and downconversion frequencies
        readout = self.us2cycles(cfg["cav_meas_window"],ro_ch=cfg["ro_channels"][0])
        self.declare_readout(ch=cfg["ro_channels"][0], length=readout,
                             freq=self.cfg["cav_freq"], gen_ch=cfg["cav_channel"])

        freq=self.freq2reg(cfg["cav_freq"], gen_ch=cfg["cav_channel"], 
                                 ro_ch=cfg["ro_channels"][0])
        self.set_pulse_registers(ch=cfg["cav_channel"], style="const", freq=freq,
                                 # converts phase degrees to QICK register val
                                 phase=self.deg2reg(cfg["cav_phase"]), 
                                 gain=cfg["cav_gain"], 
                                 length=self.us2cycles(cfg["cav_pulse_len"],gen_ch=self.cfg["cav_channel"]), mode = "periodic")
        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        self.measure(pulse_ch=self.cfg["cav_channel"], 
             adcs=[self.cfg["ro_channels"][0]],
             adc_trig_offset=self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"]),
             wait=True,
             syncdelay=self.us2cycles(self.cfg["relax_delay"],gen_ch=self.cfg["cav_channel"]))
        
######################
# Taken from CW_spec #
######################
class CW_spec_flux(AveragerProgram):
    def initialize(self):
        cfg=self.cfg   
        gen_ch = cfg["cav_channel"]
        qub_ch = cfg["qub_channel"]

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=cfg["nqz_c"])
        self.declare_gen(ch=cfg["qub_channel"], nqz=cfg["nqz_q"])
        
        # configure the readout lengths and downconversion frequencies (ensuring it is an available DAC frequency)
        readout = self.us2cycles(cfg["qub_readout_length"],ro_ch=cfg["ro_channels"][0])
        for ch in cfg["ro_channels"]:
            self.declare_readout(ch=ch, length=readout,
                                 freq=self.cfg["cav_freq"], gen_ch=cfg["cav_channel"])

        # Configure cavity DAC
        freq_c  = self.freq2reg(cfg["cav_freq"],gen_ch=gen_ch, ro_ch=cfg["ro_channels"][0])
        phase_c = self.deg2reg(cfg["cav_phase"], gen_ch=gen_ch)
        gain_c  = cfg["cav_gain"]
        self.default_pulse_registers(ch=gen_ch, freq=freq_c, phase=phase_c, gain=gain_c)
        self.set_pulse_registers(ch=gen_ch, style="const", 
                                 length=self.us2cycles(self.cfg["cav_pulse_len"],
                                 gen_ch=gen_ch), mode='periodic')
        
        # Configure qubit DAC
        freq_q  = self.freq2reg(cfg["qub_freq"],gen_ch=qub_ch)
        phase_q = self.deg2reg(cfg["qub_phase"], gen_ch=qub_ch)
        gain_q  = cfg["qub_gain"]
        
        self.default_pulse_registers(ch=qub_ch, freq=freq_q, phase=phase_q, gain=gain_q)
        self.set_pulse_registers(ch=qub_ch, style="const",
                                 length=self.us2cycles(cfg['qub_len'],
                                 gen_ch=qub_ch), mode= 'periodic')
        
        self.synci(200)   
    
    def body(self):
        
        #############################
        # LEGACY CODE FROM QUASI_CW #
        #############################
        
        #self.reset_phase(gen_ch = self.cfg['cav_channel'], t=0)
        #self.reset_phase(gen_ch = self.cfg['qub_channel'], t=0)
        
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.
        #sigma = self.us2cycles(self.cfg["qub_sigma"])
        #num_sigma = self.cfg["num_sigma"]
        #pulse_len = self.us2cycles(self.cfg['qub_len'],gen_ch=self.cfg['qub_channel']) + int(num_sigma*sigma)
        
        offset = self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"]) + self.us2cycles(20,gen_ch = self.cfg["cav_channel"])
        
        #meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=self.cfg["cav_channel"])
        #ex_time = meas_time - self.us2cycles(self.cfg['qub_delay'], gen_ch=self.cfg["qub_channel"]) - pulse_len
        #Sets off the ADC
        
        self.trigger(adcs=self.ro_chs,
                    pins=[0],
                    adc_trig_offset=offset)
        
        # Both measure() and pulse/pulse/wait/sync work.
        # Send pulses and trigger measurement
        #self.measure(pulse_ch=[self.cfg["cav_channel"], self.cfg['qub_channel']], 
             # adcs=[self.cfg["ro_channels"][0]],
             # adc_trig_offset=self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"]),
             # wait=True,
             # syncdelay=self.us2cycles(self.cfg["relax_delay"],gen_ch=self.cfg["cav_channel"]))
        self.pulse(ch=self.cfg["cav_channel"],t=0)
        self.pulse(ch=self.cfg["qub_channel"],t=0)
        self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        self.sync_all(self.us2cycles(self.cfg["relax_delay"])) #Syncs to an offset time after the final pulse is sent


def get_cw_spec_flux_settings():
    settings = {}
    
    settings['scanname'] = 'continuous_power_scan'
    settings['meas_type'] = 'CW_Spec_Flux'
    
    settings['start_voltage']  = 0
    settings['stop_voltage']   = 0.1
    settings['voltage_points'] = 5
    settings['stability'] = False
    
    autoscan = {}
    autoscan['freq_start']     = 4.4e9
    autoscan['freq_stop']      = 4.41e9
    autoscan['freq_points']    = 5
    autoscan['reps']           = 301
    autoscan['soft_avgs']      = 1

    settings['autoscan'] = autoscan
    
    settings['cav_freq'] = 1e9
    settings['cav_phase'] = 0
    settings['cav_gain'] = 1000
    settings['cav_pulse_len'] = 10
    settings['cav_meas_window'] = 900

    settings['qub_gain'] = 0
    settings['qub_len'] = 20
    settings['qub_readout_length'] = 500
    
    #Sweep parameters
    settings['freq_start']   = 4e9  
    settings['freq_step']    = 0.5e9
    settings['freq_points']  = 6

    #Card settings
    settings['reps'] = 1
    settings['soft_avgs'] = 1#5e3
    
    return settings

def cw_spec_flux(soc,soccfg,instruments,settings):
    
    print("Initializing...")
    
    # suppresses sum buffer overflow warning
    logging.getLogger("qick").setLevel(logging.ERROR)

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    
    autoscan_set = exp_settings['autoscan']
    #SRS = instruments['DCsupply']
    soc.reset_gens()
    

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel'],
        'qub_channel'     : exp_globals['qub_channel'],
        'ro_channels'     : exp_globals['ro_channels'],

        'nqz_c'           : 1,
        'cav_phase'       : m_pulse['cav_phase'],
        'cav_pulse_len'   : exp_settings['cav_pulse_len'],
        'cav_meas_window'     : exp_settings['cav_meas_window'],
        # set meas_pos to 0?
        'meas_time'       : m_pulse['meas_pos'],
        #'meas_gain'       : exp_settings['meas_gain'],
        'cav_gain'        : exp_settings['cav_gain'],
        'cav_freq'        : exp_settings['cav_freq']/1e6,
        
        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],
        'freq_start'      : exp_settings['freq_start']/1e6,
        'freq_step'       : exp_settings['freq_step']/1e6,
        'freq_points'     : exp_settings['freq_points'],
        'qub_gain'        : exp_settings['qub_gain'],
        #'qub_sigma'       : q_pulse['sigma'],
        #'qub_delay'       : q_pulse['delay'],
        #'num_sigma'       : q_pulse['num_sigma'],
        'qub_len'         : exp_settings['qub_len'],

        'qub_readout_length'  : exp_settings['qub_readout_length'],
        'adc_trig_offset' : m_pulse['emp_delay'],  #+ m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs']
        }

    #############################
    # Taken from spec_flux_scan #
    #############################
    
    #set voltage sweep
    start_voltage = exp_settings['start_voltage']
    stop_voltage  = exp_settings['stop_voltage']
    voltage_points = exp_settings['voltage_points']
    voltages = np.round(np.linspace(start_voltage, stop_voltage, voltage_points),6)
    max_voltage = 3.5
    if np.max(voltages) > max_voltage:
        raise ValueError('max voltage too large!')
    else:
        settings['voltages'] = voltages

    #SRS.Output = 'On'
    
    #Making array of cavity frequencies for transmission scan (to be looped through later)
    start_freq = autoscan_set['freq_start']
    stop_freq = autoscan_set['freq_stop']
    freq_points = autoscan_set['freq_points']
    trans_fpts = np.linspace(start_freq,stop_freq,freq_points)

    #Dummy arrays for cavity scan
    trans_mags   = np.zeros((voltage_points, autoscan_set['freq_points']))
    trans_phases = np.zeros((voltage_points, autoscan_set['freq_points']))
    
    #Dummy arrays for spec scan
    mags   = np.zeros((voltage_points, exp_settings['freq_points']))
    phases = np.zeros((voltage_points, exp_settings['freq_points']))

    #Defining file variables
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp
    
    #Defining data files
    f0_start = exp_settings['freq_start']
    f0_step = exp_settings['freq_step']
    expts = exp_settings['freq_points']
    spec_fpts = np.arange(0,expts)*f0_step+f0_start
    #powerdat = np.zeros((len(voltages), len(fpts)))
    #phasedat = np.zeros((len(voltages), len(fpts)))
    #Is = np.zeros((len(voltages), len(fpts)))
    #Qs = np.zeros((len(voltages), len(fpts)))
    
    
    
    first_it = True
    for vind in range(len(voltages)):
        
        if exp_settings['stability'] == True: # If running a stability scan, will ramp to starting voltage on the first iteration but then never ramp again.
            if first_it == True:
                print('Stability scan enabled.')
                voltage = voltages[vind]
                print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
                
                #SRS.voltage_ramp(voltage)
                time.sleep(0.1)

                voltages = np.linspace(0,len(voltages)-1,len(voltages))
            voltage = voltages[vind]
        else: # If not, will ramp to appropriate voltage every loop in standard fashion.
            voltage = voltages[vind]
            print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
            
            #SRS.voltage_ramp(voltage)
            time.sleep(0.1)
        
        print('Sweeping transmission.')
        for find in range(len(trans_fpts)):
            tstart = time.time()
            config['cav_freq'] = trans_fpts[find]/1e6 #Update the frequency, board wants it in MHz so converting now
            config['reps']     = autoscan_set['reps']
            config['soft_avgs'] = autoscan_set['soft_avgs']

            trans_prog = CavitySweepFlux(soccfg,config) #Make transmission pulse sequence object

            trans_I, trans_Q = trans_prog.acquire(soc,load_pulses=True,progress=False) #Transmission data acquisition occurs here

            mag = np.sqrt(trans_I[0][0]**2 + trans_Q[0][0]**2)
            phase = np.arctan2(trans_Q[0][0], trans_I[0][0])*180/np.pi

            trans_mags[vind][find] = mag
            trans_phases[vind][find] = phase
            
            projected_time = (time.time() - tstart) * len(trans_fpts)
            
            if find == 0:
                print(f"Projected time for transmission sweep: {projected_time:.1f} seconds. ")
        print("Transmission sweep complete. ")

        hanger = exp_globals['hanger'] #"Fitting" cav freq

        if first_it and exp_settings['stability'] == True: # If running a stability scan, will fix the cavity frequency after first loop.
            first_it = False
            try:
                print("Fitting Lorenzian to Cavity")
                freq_holder = fit_model(trans_fpts, trans_mags[vind], 'lorenz')['center']/1e6

            except:
                print("Fitting Lorenzian Failed, taking extrema...")
                if hanger:
                    freq_holder = trans_fpts[np.argmin(trans_mags[vind])]/1e6

                else:
                    freq_holder = trans_fpts[np.argmax(trans_mags[vind])]/1e6
            config['cav_freq'] = freq_holder
        elif exp_settings['stability'] == False:
            try:
                print("Fitting Lorenzian to Cavity")
                config['cav_freq'] = fit_model(trans_fpts, trans_mags[vind], 'lorenz')['center']/1e6
            except:
                print("Fitting Lorenzian Failed, taking extrema...")
                if hanger:
                    config['cav_freq'] = trans_fpts[np.argmin(trans_mags[vind])]/1e6
                else:
                    config['cav_freq'] = trans_fpts[np.argmax(trans_mags[vind])]/1e6
        else:
            config['cav_freq'] = freq_holder

        print('spec, cav freq: {} GHz'.format(config['cav_freq']/1e3))

        config['reps'] = exp_settings['reps']
        config['soft_avgs'] = exp_settings['soft_avgs']
        
        ######################
        # Taken from CW_spec #
        ######################
        t_start = time.time()
        
        # qubit frequency sweep
        for f in range(0,len(spec_fpts)):
            
            config["qub_freq"]=spec_fpts[f]/1e6 # convert to MHz
            
            prog = CW_spec_flux(soccfg,config) #Make spec pulse sequence object
            avg_di, avg_dq = prog.acquire(soc, progress=False) #Spec data acquisition # REMOVED expt_pts
            soc.reset_gens()
            
            Is = avg_di[0][0] # These ones I copied directly from another script, so I'm more certain of the indexing.
            Qs = avg_dq[0][0]
            #spec_fpts = exp_pts[0]*1e6
            powerdat = np.sqrt(Is**2 + Qs**2)
            phasedat = np.arctan2(Qs,Is)*180/np.pi
            mags[vind] = powerdat
            phases[vind] = phasedat

            #Plots (updates each for loop)
            transdata = {}
            transdata['xaxis'] = trans_fpts/1e9
            transdata['mags'] = trans_mags[0:vind+1,:]
            transdata['phases'] = trans_phases[0:vind+1,:]
            
            specdata = {}
            specdata['xaxis'] = spec_fpts/1e9
            specdata['mags'] = mags[0:vind+1,:]
            specdata['phases'] = phases[0:vind+1,:]
            
            singledata = {}
            singledata['xaxis'] = spec_fpts[f]/1e9
            singledata['mag'] = powerdat
            singledata['phase'] = phasedat
            
            trans_labels = ['Freq (GHz)','Voltage (V)']
            spec_labels  = ['Freq (GHz)','Voltage (V)']
            
            #modify the spec data to subtract the offset in amp and phase
            #and then plot the modified version
            specplotdata = {}
            specplotdata['xaxis']  = specdata['xaxis']
            specplotdata['mags']   = specdata['mags']
            specplotdata['phases'] = specdata['phases']
            
            mat = np.copy(specplotdata['mags'])
            for ind in range(0, mat.shape[0]):
                mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
            specplotdata['mags'] = mat
            
            mat = np.copy(specplotdata['phases'])
            for ind in range(0, mat.shape[0]):
                mat[ind,:]  = mat[ind,:] - np.mean(mat[ind,:])
            specplotdata['phases'] = mat
            
            identifier = 'Cav Gain : ' + str(config['cav_gain'])  + ' au'

            plots.autoscan_plot(transdata, specplotdata, singledata, voltages[0:vind+1], filename, trans_labels, spec_labels, identifier, fig_num = 1)

            #Saving data
            userfuncs.SaveFull(saveDir, filename, ['transdata', 'specdata', 'singledata', 'voltages', 
                                           'filename', 'trans_labels', 'spec_labels'], 
                                           locals(), expsettings=settings, instruments=instruments)
            plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)

            if f == 0:
                t_stop = time.time()
                estimate_time(t_start, t_stop, len(spec_fpts)*len(voltages))

    return transdata,specdata,prog
                
                
    '''
    full_data = {}

    full_data['xaxis']  = fpts/1e9
    full_data['mags']   = powerdat
    full_data['phases'] = phasedat
    full_data['Is']     = Is
    full_data['Qs']     = Qs
    
    plot_data = {}
    plot_data['xaxis']  = fpts/1e9
    plot_data['mags']   = powerdat[0:vind+1]
    plot_data['phases'] = phasedat[0:vind+1]

    single_data = {}
    single_data['xaxis'] = fpts/1e9
    single_data['mag']   = powerdat[vind]
    single_data['phase'] = phasedat[vind]
    
    yaxis  = voltages[0:vind+1]
    labels = ['Freq (GHz)', 'Voltages (V)']
    
    simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier='', fig_num=1, IQdata = False)
    
    userfuncs.SaveFull(saveDir, filename, ['voltages','fpts','full_data','filename'],
    locals(), expsettings=settings, instruments={})

    '''

