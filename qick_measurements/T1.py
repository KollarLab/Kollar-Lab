# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:16:54 2023

@author: kollarlab
"""

from qick.averager_program import AveragerProgram

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
import time

from utility.userfits import fit_T1
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import estimate_time
import scipy.signal as signal

      
class T1_sequence(AveragerProgram):
    '''
    T1_sequence _summary_

    :param AveragerProgram: _description_
    :type AveragerProgram: _type_
    '''    
    def initialize(self):
        '''
        initialize _summary_
        '''        
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
        
        self.default_pulse_registers(ch=gen_ch, freq=freq_c, phase=phase_c, gain=gain_c)
        self.set_pulse_registers(ch=gen_ch, style="const", length=self.us2cycles(self.cfg["meas_window"],gen_ch=gen_ch))
        
        freq_q  = self.freq2reg(cfg["qub_freq"],gen_ch=qub_ch)
        phase_q = self.deg2reg(cfg["qub_phase"], gen_ch=qub_ch)
        gain_q  = cfg["qub_gain"]
        
        self.default_pulse_registers(ch=qub_ch, phase=phase_q, freq=freq_q, gain=gain_q)
        
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
            self.reset_phase(gen_ch = [self.cfg['qub_channel']], t=0)
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
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.

        #self.reset_phase(gen_ch = self.cfg['cav_channel'], t=0)
        #self.reset_phase(gen_ch = self.cfg['qub_channel'], t=0)

        offset = self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"])
        meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=self.cfg["cav_channel"])
        #Sets off the ADC
        self.trigger(adcs=self.ro_chs,
                    pins=[0],
                    adc_trig_offset=offset)
        
        #Sends measurement pulse
        self.pulse(ch=self.cfg["cav_channel"],t=meas_time)

        self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        self.sync_all(self.us2cycles(self.cfg["relax_delay"]))


def get_T1_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'Tmeas'
    settings['phase_reset'] = True
    
    settings['cav_freq'] = 1e9
    settings['meas_gain'] = 1000
    
    settings['qub_freq'] = 4e9
    settings['qub_gain'] = 1000
    
    #Sweep parameters    
    settings['Tau_min']    = 200e-9
    settings['Tau_max']    = 30e-6
    settings['Tau_points'] = 5
    settings['spacing']    = 'Linear'
    
    settings['T1_guess'] = 10e-6
    #Card settings
    settings['reps'] = 1
    settings['soft_avgs'] = 5e3
    
    return settings

def meas_T1(soc,soccfg,instruments,settings):
    '''
    meas_T1 - This function is the core of the T1 script

    :param soc: _description_
    :type soc: _type_
    :param soccfg: _description_
    :type soccfg: _type_
    :param instruments: _description_
    :type instruments: _type_
    :param settings: _description_
    :type settings: _type_
    '''    

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    
    #if exp_globals['LO']:
    if False:
        logen = instruments['LO']
        
        logen.freq   = exp_globals['LO_freq']
        logen.power  = exp_globals['LO_power']
        logen.output = 1
    #else:
        #print("Warning: No LO has been initialized.")

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
        'meas_gain'       : exp_settings['meas_gain'],
        'cav_freq'        : (exp_settings['cav_freq']-exp_globals['LO_freq']*exp_globals['LO'])/1e6,
        
        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],
        'qub_freq'        : exp_settings['qub_freq']/1e6,
        'qub_gain'         : exp_settings['qub_gain'],
        'qub_sigma'        : q_pulse['sigma'],
        'qub_delay'        : 0, #Placeholder
        'num_sigma'       : q_pulse['num_sigma'],
        
        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['reps'],
        'soft_avgs'       : exp_settings['soft_avgs'],

        'phase_reset'     : exp_settings['phase_reset']
        }


    prog = T1_sequence(soccfg,config)
    meas_start = prog.us2cycles(m_pulse["init_buffer"],ro_ch=0)
    meas_end = meas_start+prog.us2cycles(m_pulse["meas_window"],ro_ch=0)
    total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)
    
    ## Set up array of taus and randomize it
    if exp_settings['spacing']=='Log':
        tau_list = np.logspace(np.log10(exp_settings['Tau_min']), np.log10(exp_settings['Tau_max']), exp_settings['Tau_points'])
    else:
         tau_list = np.linspace(exp_settings['Tau_min'], exp_settings['Tau_max'], exp_settings['Tau_points'])
    taus = np.round(tau_list, 9)
    
    indices = list(range(len(taus)))
    
    amp_int = np.zeros(len(taus))
    ang_int = np.zeros(len(taus))

    amp_orig = np.zeros(len(taus))
    ang_orig = np.zeros(len(taus))


    
    
    first_it = True
    
    if exp_settings['subtract_background']:
        #Acquire background trace
#            qubitgen.freq=3.8e9
#            time.sleep(0.1)
        print('Starting Background Trace')
        bprog = CavitySweep(soccfg,config)
        holder = bprog.acquire(soc, load_pulses=True, progress=False)
        
        I_back = holder[0][0][0]
        Q_back = holder[1][0][0]
        print('Background Trace Complete')
        I_full_b = holder[0][0]
    else:
        I_back, Q_back = 0,0

    tstart = time.time()
    
    for tind in indices:
        
        tau = taus[tind]
        print('Tau: {}'.format(tau))
        
        config['qub_delay'] = tau*1e6
        prog = T1_sequence(soccfg,config)
        holder = prog.acquire(soc, load_pulses=True, progress=False)
        
        ##No background subtraction here!!!!
        I_sig, Q_sig   = [holder[0][0][0], holder[1][0][0]] #<I>, <Q> for signal trace
        

        
#        I_final, Q_final   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back
        
#        amps[tind] = np.sqrt(I_full**2+Q_full**2)
#        angles[tind] = np.arctan2(Q_full,I_full)*180/np.pi
#        amp_int[tind] = np.sqrt(I_final**2+Q_final**2)
#        ang_int[tind] = np.arctan2(Q_final, I_final)*180/np.pi
        
        amp_int[tind] = np.sqrt(I_final**2 + Q_final**2)
        ang_int[tind] = np.arctan2(Q_final, I_final)*180/np.pi
        
        amp_orig[tind] = np.sqrt(I_sig**2+Q_sig**2)
        ang_orig[tind] = np.arctan2(Q_sig, I_sig)*180/np.pi

        
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(taus))
            
                
            first_it = False  
            
        #_______________________________________#
    
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(taus*1e6, amp_int, 'x')
        plt.xlabel('Tau (us)')
        plt.ylabel('Amplitude')  
        plt.subplot(122)
        plt.plot(taus*1e6, ang_int, 'x')
        plt.xlabel('Tau (us)')
        plt.ylabel('Phase')  
        plt.title('Live T1 data (no fit)\n'+filename)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)
    
        # fig2 = plt.figure(2,figsize=(13,8))
        # plt.clf()
    
        # ax = plt.subplot(1,1,1)
        # general_colormap_subplot(ax, xaxis*1e6, taus*1e6, amps, ['Time (us)', 'Tau (us)'], 'Raw data\n'+filename)
    
        # plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
        userfuncs.SaveFull(saveDir, filename, ['taus', 'amp_int'], locals(), expsettings=settings, instruments=instruments)
        
        
        
    if exp_settings['debug']:
        
        
        

        
        config['readout_length'] = 3 # don't mess with this
        config['adc_trig_offset'] = m_pulse['emp_delay'] + m_pulse['meas_pos'] + exp_settings['debug_time']
        
        total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)

        amp_intd = np.zeros(len(taus))
        ang_intd = np.zeros(len(taus))
        ampsd    = np.zeros((len(taus),total_samples))
        anglesd  = np.zeros((len(taus),total_samples))        

        if exp_settings['subtract_background']:
            #Acquire background trace
    #            qubitgen.freq=3.8e9
    #            time.sleep(0.1)
            print('Starting Background Trace')
            bprog = CavitySweep(soccfg,config)
            holder = bprog.acquire_decimated(soc, load_pulses=True, progress=False)
            print('Background Trace Complete')
            I_full_bd = holder[0][0]
            Q_full_bd = holder[0][1]
            I_window_bd = I_full_bd[meas_start:meas_end]
            Q_window_bd = Q_full_bd[meas_start:meas_end]     
        else:
            I_window_bd, Q_window_bd, I_full_bd, Q_full_bd = 0,0,0,0
            
        I_backd, Q_backd = [np.mean(I_window_bd), np.mean(Q_window_bd)]
        
        first_it = True
        
        for tind in indices:
            
            tau = taus[tind]
            print('Tau: {}'.format(tau))
            
            config['qub_delay'] = tau*1e6
            prog = T1_sequence(soccfg,config)
            holder = prog.acquire_decimated(soc, load_pulses=True, progress=False)
            I_fulld = holder[0][0]
            Q_fulld = holder[0][1]
            I_windowd = I_fulld[meas_start:meas_end]
            Q_windowd = Q_fulld[meas_start:meas_end]
            
            ##No background subtraction here!!!!
            I_sigd, Q_sigd   = [np.mean(I_windowd), np.mean(Q_windowd)] #<I>, <Q> for signal trace
            

            
    #        I_final, Q_final   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
            
            I_finald = I_sigd-I_backd #compute <I_net> in the data window
            Q_finald = Q_sigd-Q_backd
            
    #        amps[tind] = np.sqrt(I_full**2+Q_full**2)
    #        angles[tind] = np.arctan2(Q_full,I_full)*180/np.pi
    #        amp_int[tind] = np.sqrt(I_final**2+Q_final**2)
    #        ang_int[tind] = np.arctan2(Q_final, I_final)*180/np.pi
            
            
            ampsd[tind] = np.sqrt((I_fulld-I_full_bd)**2+(Q_fulld-Q_full_bd)**2)
            anglesd[tind] = np.arctan2((Q_fulld-Q_full_bd), (I_fulld-I_full_bd))*180/np.pi
            amp_intd[tind] = np.sqrt(I_finald**2+Q_finald**2)
            ang_intd[tind] = np.arctan2(Q_finald, I_finald)*180/np.pi
            

            
            if first_it:
                tstop = time.time()
                estimate_time(tstart, tstop, len(taus))
                
                xaxis = np.linspace(0,len(I_fulld)-1,len(I_fulld))

                for x in range(0,len(xaxis)):
                    xaxis[x] = prog.cycles2us(xaxis[x],ro_ch=0)
                    
                first_it = False  
                
                
            fig = plt.figure(3, figsize=(13,8))
            plt.clf()
            plt.subplot(121)
            plt.plot(taus*1e6, amp_intd, 'x')
            plt.xlabel('Tau (us)')
            plt.ylabel('Amplitude')  
            plt.subplot(122)
            plt.plot(taus*1e6, ang_intd, 'x')
            plt.xlabel('Tau (us)')
            plt.ylabel('Phase')  
            plt.title('Live T1 data debug (no fit)\n'+filename)
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.savefig(os.path.join(saveDir, filename+'_debug_no_fit.png'), dpi = 150)
        
            fig2 = plt.figure(4,figsize=(13,8))
            plt.clf()
        
            ax = plt.subplot(121)
            general_colormap_subplot(ax, xaxis*1e6, taus*1e6, ampsd, ['Time (us)', 'Tau (us)'], 'Raw data\n'+filename)
            ax = plt.subplot(122)
            general_colormap_subplot(ax, xaxis*1e6, taus*1e6, anglesd, ['Time (us)', 'Tau (us)'], 'Raw data\n'+filename)
            plt.savefig(os.path.join(saveDir, filename+'_debug_fulldata.png'), dpi = 150)
        
    
    t2 = time.time()
    print('Elapsed time: {}'.format(t2-tstart))

    T1_guess = exp_settings['T1_guess']
    amp_guess = max(amp_int)-min(amp_int)
    offset_guess = np.mean(amp_int[-10:])

    fit_guess = [T1_guess, amp_guess, offset_guess]
    T1, amp, offset, fit_xvals, fit_yvals = fit_T1(taus, amp_int, fit_guess)
    fig3 = plt.figure(2)
    plt.clf()
    plt.plot(taus*1e6, amp_int)
    plt.plot(fit_xvals*1e6, fit_yvals)
    plt.title('T1:{}us \n {}'.format(np.round(T1*1e6,3), filename))
    plt.xlabel('Time (us)')
    plt.ylabel('Amplitude')
    fig3.canvas.draw()
    fig3.canvas.flush_events()
    plt.savefig(os.path.join(saveDir, filename+'_fit.png'), dpi=150)

    ang_guess = max(ang_int) - min(ang_int)
    fit_guess_phase = [T1_guess, ang_guess, offset_guess]
    T1, amp, offset, fit_xvals, fit_yvals = fit_T1(taus, ang_int, fit_guess)
    fig4 = plt.figure(3)
    plt.clf()
    plt.plot(taus*1e6, ang_int)
    plt.plot(fit_xvals*1e6, fit_yvals)
    plt.title('Phase T1:{}us \n {}'.format(np.round(T1*1e6,3), filename))
    plt.xlabel('Time (us)')
    plt.ylabel('Phase')
    fig4.canvas.draw()
    fig4.canvas.flush_events()
    plt.savefig(os.path.join(saveDir, filename+'_fit_p.png'), dpi=150)

    userfuncs.SaveFull(saveDir, filename, ['taus','amp_int','ang_int','amp_orig','ang_orig', 'T1', 'amp', 'offset', 'fit_guess'],
                         locals(), expsettings=settings, instruments=instruments)
    
    if exp_globals['LO']:
        pass
        #logen.output = 0

    return T1, taus, amp_int
   
    

        