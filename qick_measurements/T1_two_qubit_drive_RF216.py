# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:16:54 2023

@author: kollarlab
"""

from qick.asm_v2 import AveragerProgramV2, QickSweep1D

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
import time

from utility.userfits import fit_T1
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import estimate_time
import scipy.signal as signal

      
class T1_sequence(AveragerProgramV2):
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
        
        self.add_pulse(
                ch=qub_ch, name="qub_phrst", ro_ch=ro_ch,
                style="const",
                freq=cfg["qub_freq"],      # doesn't really matter if gain=0
                phase=0,
                gain=0,                   # no output
                length=0.015,              # small but nonzero (us); pick safely > 0
                phrst=1                   # <-- resets phase accumulator
            )
        
        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)
        
        
    
    def _body(self, cfg):
# =============================================================================
#         qub_ch = self.cfg["qub_channel"]
#         self.reset_phase(gen_ch = [qub_ch], t=0)
# =============================================================================
        sigma = float(cfg["qub_sigma"])
        num_sigma = int(cfg["num_sigma"])
        
        pulse_len = float(cfg['hold_length']) + num_sigma*sigma
        
        offset = cfg["adc_trig_offset"]
        meas_time = self.cfg["meas_time"]
        ex_time = meas_time - cfg['qub_delay'] - pulse_len
        
        
        # Optional phase reset behavior
        if cfg.get("phase_reset", False):
            self.pulse(ch=cfg["qub_channel"], name="qub_phrst", t=0)
        
        #Sets off the ADC
        self.trigger(ros=[cfg['ro_channel']],
                    pins=[0],
                    t=offset)
        
        self.pulse(ch=cfg["qub_channel"],name='qub_pulse',t=ex_time)
        self.pulse(ch=cfg["cav_channel"],name='cav_pulse',t=meas_time)
        self.wait_auto()
        self.delay(self.cfg["relax_delay"])


def get_T1_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'Tmeas'
    # settings['phase_reset'] = True
    
    # settings['cav_freq'] = 1e9
    # settings['meas_gain'] = 1000
    
    # settings['qub_freq'] = 4e9
    # settings['qub_gain'] = 1000
    
    # #Sweep parameters    
    # settings['Tau_min']    = 200e-9
    # settings['Tau_max']    = 30e-6
    # settings['Tau_points'] = 5
    settings['spacing']    = 'Linear'
    
    # settings['T1_guess'] = 10e-6
    # #Card settings
    # settings['reps'] = 1
    # settings['soft_avgs'] = 5e3
    
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
    
    qub_drive_index = exp_settings['qub_drive_index']
    if qub_drive_index == 'D1':
        qub_ch = exp_globals['qub_channel_1']
        qub_channel = exp_globals['qub_channel_1']['ID']
        q_pulse = exp_globals['qubit_pulse_D1']
        
    elif qub_drive_index == 'D2':
        qub_ch = exp_globals['qub_channel_2']
        qub_channel = exp_globals['qub_channel_2']['ID']
        q_pulse = exp_globals['qubit_pulse_D2']
        
    else:
        print("Wrong qubit drive index, check carefully")
    

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel']['ID'],
        'qub_channel'     : qub_channel,
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
        'qub_gain'        : exp_settings['qub_gain'],
        'qub_mixer_freq'  : (exp_settings['qub_freq']+exp_settings['qub_mixer_detuning'])/1e6,

        'qub_sigma'       : q_pulse['sigma'],
        'qub_delay'       : exp_globals['qub_delay_fixed'],
        'num_sigma'       : q_pulse['num_sigma'],
        'hold_length'     : exp_settings['hold_time'], 

        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'phase_reset'     : exp_settings['phase_reset']
        }

    
    cav_ch = exp_globals['cav_channel']
    #qub_ch = exp_globals['qub_channel']
    ro_ch  = exp_globals['ro_channel']
    # Set attenuator on DAC.
    soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
    soc.rfb_set_gen_rf(qub_ch['ID'], qub_ch['Atten_1'], qub_ch['Atten_2'])
    # Set attenuator on ADC.
    soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])
    
    if exp_settings['filter'] == 'all_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
        
        soc.rfb_set_gen_filter(config['qub_channel'], fc=config['qub_freq']/1000, ftype='bandpass', bw=qub_ch['BW'])
        
    elif exp_settings['filter'] == 'no_qubit_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
    
        soc.rfb_set_gen_filter(config['qub_channel'], fc=config['qub_freq']/1000, ftype='bypass')
        
    elif exp_settings['filter'] == 'no_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bypass')
        soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
        
        soc.rfb_set_gen_filter(config['qub_channel'], fc=config['qub_freq']/1000, ftype='bypass')
        
    else:
        print('Please select one option from:')
        print('\'all_filter\', \'no_qubit_filter\', and \'no_filter\'')
        return
        
    prog = T1_sequence(soccfg,reps = exp_settings['reps'], final_delay = None, final_wait = 0, cfg = config)#Wait this uses reps instead of rounds????

    
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
    
# =============================================================================
#     if exp_settings['subtract_background']:
#         #Acquire background trace
# #            qubitgen.freq=3.8e9
# #            time.sleep(0.1)
#         print('Starting Background Trace')
#         bprog = CavitySweep(soccfg,config)
#         holder = bprog.acquire(soc, load_pulses=True, progress=False)
#         
#         I_back = holder[0][0][0]
#         Q_back = holder[1][0][0]
#         print('Background Trace Complete')
#         I_full_b = holder[0][0]
#     else:
# =============================================================================
    I_back, Q_back = 0,0

    tstart = time.time()
    
    for tind in indices:
        
        tau = taus[tind]
        print('Tau: {}'.format(tau))
        config['qub_delay'] = tau*1e6
        prog = T1_sequence(soccfg,reps = exp_settings['reps'], final_delay = None, final_wait = 0, cfg = config)
        holder = prog.acquire(soc, rounds = exp_settings['rounds'], load_pulses=True, progress=False)
        
        ##No background subtraction here!!!!
        I_sig, Q_sig   = [holder[0][0][0], holder[0][0][1]] #<I>, <Q> for signal trace
        

        
#        I_final, Q_final   = [np.mean(I_window), np.mean(Q_window)] #<I>, <Q> for signal trace
        
        I_final = I_sig-I_back #compute <I_net> in the data window
        Q_final = Q_sig-Q_back
        
#        amps[tind] = np.sqrt(I_full**2+Q_full**2)
#        angles[tind] = np.arctan2(Q_full,I_full)*180/np.pi
#        amp_int[tind] = np.sqrt(I_final**2+Q_final**2)
#        ang_int[tind] = np.arctan2(Q_final, I_final)*180/np.pi
        
        amp_int[tind] = np.sqrt(I_final**2 + Q_final**2)
        ang_int[tind] = np.degrees(np.arctan2(Q_final, I_final))
        
        amp_orig[tind] = np.sqrt(I_sig**2+Q_sig**2)
        ang_orig[tind] = np.degrees(np.arctan2(Q_sig, I_sig))

        
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
        plt.suptitle('Live T1 data (no fit)\n'+filename)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)
    
        # fig2 = plt.figure(2,figsize=(13,8))
        # plt.clf()
    
        # ax = plt.subplot(1,1,1)
        # general_colormap_subplot(ax, xaxis*1e6, taus*1e6, amps, ['Time (us)', 'Tau (us)'], 'Raw data\n'+filename)
    
        # plt.savefig(os.path.join(saveDir, filename+'_fulldata.png'), dpi = 150)
        userfuncs.SaveFull(saveDir, filename, ['taus', 'amp_int'], locals(), expsettings=settings, instruments=instruments)
        
        
        
    
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
   
    

        