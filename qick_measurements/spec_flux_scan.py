# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: kollarlab
"""

from qick.averager_program import QickSweep, AveragerProgram, NDAveragerProgram

import numpy as np
import os
import time
import matplotlib.pyplot as plt
import userfuncs
import utility.plotting_tools as plots
from utility.userfits_v2 import fit_model
      
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

class Quasi_CW(NDAveragerProgram):
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

        # Configure cavity DAC
        freq_c  = self.freq2reg(cfg["cav_freq"],gen_ch=gen_ch, ro_ch=cfg["ro_channels"][0])
        phase_c = self.deg2reg(cfg["cav_phase"], gen_ch=gen_ch)
        gain_c  = cfg["meas_gain"]
        
        self.default_pulse_registers(ch=gen_ch, freq=freq_c, phase=phase_c, gain=gain_c)
        self.set_pulse_registers(ch=gen_ch, style="const", length=self.us2cycles(self.cfg["meas_window"],gen_ch=gen_ch))
        
        # Configure qubit DAC
        freq_q  = self.freq2reg(cfg["freq_start"],gen_ch=qub_ch)
        phase_q = self.deg2reg(cfg["qub_phase"], gen_ch=gen_ch)
        gain_q  = cfg["qub_gain"]
        
        self.default_pulse_registers(ch=qub_ch, phase=phase_q, freq=freq_q, gain=gain_q)
        
        sigma = self.us2cycles(cfg["qub_sigma"],gen_ch=qub_ch) 
        num_sigma = cfg["num_sigma"]
        
        self.add_gauss(ch=qub_ch, name="ex", sigma=sigma,length=int(sigma*num_sigma))        
        self.set_pulse_registers(ch=qub_ch, style="flat_top", waveform="ex",length=self.us2cycles(cfg['quasi_CW_len'],gen_ch=qub_ch))

        ###Start sweep definition
        
        self.qub_r_freq = self.get_gen_reg(qub_ch,"freq")
        
        
        #freq_start_reg = self.freq2reg(cfg["freq_start"],gen_ch=qub_ch)
        #freq_stop_reg = self.freq2reg(cfg["freq_stop"],gen_ch=qub_ch)
        
        self.add_sweep(QickSweep(self, self.qub_r_freq, cfg['freq_start'], cfg['freq_stop'], cfg["freq_points"]))
        
        # give processor some time to configure pulses, I believe it sets the offset time 200 cycles into the future?
        # Try varying this and seeing if it moves. Also try putting a synci after the trigger in the body
        self.synci(200)   
    
    def body(self):
        
        #self.reset_phase(gen_ch = self.cfg['cav_channel'], t=0)
        #self.reset_phase(gen_ch = self.cfg['qub_channel'], t=0)
        
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.
        sigma = self.us2cycles(self.cfg["qub_sigma"])
        num_sigma = self.cfg["num_sigma"]
        
        pulse_len = self.us2cycles(self.cfg['quasi_CW_len'],gen_ch=self.cfg['qub_channel']) + int(num_sigma*sigma)

        offset = self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"])
        meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=self.cfg["cav_channel"])
        ex_time = meas_time - self.us2cycles(self.cfg['qub_delay'],gen_ch=self.cfg["qub_channel"]) - pulse_len
        #Sets off the ADC
        self.trigger(adcs=self.ro_chs,
                    pins=[0],
                    adc_trig_offset=offset)
        
        #Sends measurement pulse
        self.pulse(ch=self.cfg["qub_channel"],t=ex_time)
        self.pulse(ch=self.cfg["cav_channel"],t=meas_time)
        self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        self.sync_all(self.us2cycles(self.cfg["relax_delay"])) #Syncs to an offset time after the final pulse is sent


def get_spec_flux_scan_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'SpecFluxScan'

    settings['start_voltage']  = 0
    settings['stop_voltage']   = 0.1
    settings['voltage_points'] = 5

    settings['meas_gain'] = 10000

    settings['qub_gain']     =  30000#1000 #30000
    settings['quasi_CW_len'] = 150 #us
    settings['freq_start']      = 0
    settings['freq_stop']       = 0
    settings['freq_points']     = 5

    #ADC settings
    settings['reps']      = 301
    settings['soft_avgs']  = 1

    autoscan = {}
    autoscan['freq_start']     = 4.4e9
    autoscan['freq_stop']      = 4.41e9
    autoscan['freq_points']    = 5
    autoscan['reps']           = 301
    autoscan['soft_avgs']      = 1

    settings['autoscan'] = autoscan
    
    return settings

def spec_flux_scan(soc,soccfg,instruments,settings):

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 

    autoscan_set = exp_settings['autoscan']

    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    
    SRS = instruments['DCsupply']


    soc.reset_gens()
    
    lo_freq = exp_globals["LO_freq"] #This is for if we're using a SignalCore for upconversion, in globals this is set to 0 GHz so not a problem
    

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    # config dictionary is what gets passed to the sequence object, this is just a repackaging of the exp_settings and exp_globals
    config = {
        'cav_channel'     : exp_globals['cav_channel'],
        'qub_channel'     : exp_globals['qub_channel'],
        'ro_channels'     : exp_globals['ro_channels'],

        'nqz_c'           : 1,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'meas_gain'       : exp_settings['meas_gain'],
        'cav_freq'        : 0, # Going to sweep through this in for loop, dummy variable for now
        
        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],

        'freq_start'      : exp_settings['freq_start']/1e6, #sweep params for quasi_cw freq, note that board always wants frequency in MHz and time in us
        'freq_stop'       : exp_settings['freq_stop']/1e6,
        'freq_points'     : exp_settings['freq_points'],
        'qub_gain'        : exp_settings['qub_gain'],

        'qub_sigma'       : q_pulse['sigma'],
        'qub_delay'       : q_pulse['delay'],
        'num_sigma'       : q_pulse['num_sigma'],
        'quasi_CW_len'    : exp_settings['quasi_CW_len'],

        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs']
        }


    #rep_period = config['adc_trig_offset'] + config['readout_length'] + config['relax_delay']
    
    
    
    #projected_time = config['reps']*config['soft_avgs']*config['freq_points']*rep_period/1e6
    #print("Projected Time: " + str(projected_time))
    
    t_i = time.time() # Come back and add a time estimation thing
    
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

    SRS.Output = 'On'

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

    first_it = True


    for vind in range(len(voltages)):

        if exp_settings['stability'] == True: # If running a stability scan, will ramp to starting voltage on the first iteration but then never ramp again.
            if first_it == True:
                voltage = voltages[vind]
                print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
                
                SRS.voltage_ramp(voltage)
                time.sleep(0.1)

                voltages = np.linspace(0,len(voltages)-1,len(voltages))
            voltage = voltages[vind]
        else: # If not, will ramp to appropriate voltage every loop in standard fashion.
            voltage = voltages[vind]
            print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
            
            SRS.voltage_ramp(voltage)
            time.sleep(0.1)


        
        print('trans')
        for find in range(len(trans_fpts)):
            config['cav_freq'] = trans_fpts[find]/1e6 #Update the frequency, board wants it in MHz so converting now
            config['reps']     = autoscan_set['reps']
            config['soft_avgs'] = autoscan_set['soft_avgs']

            trans_prog = CavitySweep(soccfg,config) #Make transmission pulse sequence object

            trans_I, trans_Q = trans_prog.acquire(soc,load_pulses=True,progress=False) #Transmission data acquisition occurs here
            #Note to Max: The most likely bug is me screwing up the indexing here on trans_I and trans_Q.
            # trans_I and trans_Q should only have one value but the default is to store it as a list of list, and it's very possible I'm misreading the documentation right now.
            mag = np.sqrt(trans_I[0][0]**2 + trans_Q[0][0]**2)
            phase = np.arctan2(trans_Q[0][0], trans_I[0][0])*180/np.pi

            trans_mags[vind][find] = mag
            trans_phases[vind][find] = phase
        

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

        print('spec, cav freq: {}'.format(config['cav_freq']/1e3))

        config['reps'] = exp_settings['reps']
        config['soft_avgs'] = exp_settings['soft_avgs']
        prog = Quasi_CW(soccfg,config) #Make spec pulse sequence object
        
        exp_pts, avg_di, avg_dq = prog.acquire(soc, load_pulses=True, progress=False) #Spec data acquisition

        Is = avg_di[0][0] # These ones I copied directly from another script, so I'm more certain of the indexing.
        Qs = avg_dq[0][0]
        
        spec_fpts = exp_pts[0]*1e6
        
        powerdat = np.sqrt(Is**2 + Qs**2)
        phasedat = np.arctan(Qs,Is)*180/np.pi

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
        singledata['xaxis'] = spec_fpts/1e9
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
        
        identifier = 'Cav Gain : ' + str(config['meas_gain'])  + ' au'

        plots.autoscan_plot(transdata, specplotdata, singledata, voltages[0:vind+1], filename, trans_labels, spec_labels, identifier, fig_num = 1)

        #Saving data
        userfuncs.SaveFull(saveDir, filename, ['transdata', 'specdata', 'singledata', 'voltages', 
                                       'filename', 'trans_labels', 'spec_labels'], 
                                       locals(), expsettings=settings, instruments=instruments)
        plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)



    return transdata,specdata,prog
