# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:16:54 2023

@author: kollarlab
"""

from qick.asm_v2 import AveragerProgramV2

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
import time

from utility.userfits import fit_T2, fit_T1
from utility.plotting_tools import general_colormap_subplot
from utility.measurement_helpers import estimate_time

class T2_sequence(AveragerProgramV2):
    def _initialize(self, cfg):
        ro_ch  = cfg['ro_channel']
        gen_ch = cfg["cav_channel"]
        qub_ch = cfg["qub_channel"]

        # Generators + readout
        self.declare_gen(ch=gen_ch, nqz=cfg["nqz_c"], mixer_freq=cfg['cav_mixer_freq'], ro_ch=ro_ch)
        self.declare_gen(ch=qub_ch, nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq'])
        self.declare_readout(ch=ro_ch, length=cfg['readout_length'])

        self.add_readoutconfig(
            ch=ro_ch, name="myro",
            freq=cfg['cav_freq'],
            gen_ch=gen_ch,
            outsel='product'
        )

        self.add_pulse(
            ch=gen_ch, name="cav_pulse", ro_ch=ro_ch,
            style="const",
            freq=cfg['cav_freq'],
            length=cfg["meas_window"],
            phase=cfg['cav_phase'],
            gain=cfg['cav_gain'],
        )

        # Envelope (all times in us)
        sigma = float(cfg["qub_sigma"])
        ns    = int(cfg["num_sigma"])
        self.add_gauss(ch=qub_ch, name='ramp', sigma=sigma, length=sigma*ns)

        # Pulses: always define all three; in _body we decide whether to play echo
        self.add_pulse(
            ch=qub_ch, name="qub_pulse_1", ro_ch=ro_ch,
            style="flat_top",
            envelope="ramp",
            freq=cfg['qub_freq'],
            length=float(cfg['hold_length']),
            phase=float(cfg.get('qub_pulse_phase_1', cfg.get('qub_phase', 0.0))),
            gain=float(cfg['qub_gain'])/2,
        )

        self.add_pulse(
            ch=qub_ch, name="qub_pulse_echo", ro_ch=ro_ch,
            style="flat_top",
            envelope="ramp",
            freq=cfg['qub_freq'],
            length=float(cfg['hold_length']),
            phase=float(cfg.get('qub_pulse_phase_echo', cfg.get('qub_phase', 0.0))),
            gain=float(cfg['qub_gain']),
        )

        self.add_pulse(
            ch=qub_ch, name="qub_pulse_2", ro_ch=ro_ch,
            style="flat_top",
            envelope="ramp",
            freq=cfg['qub_freq'],
            length=float(cfg['hold_length']),
            phase=float(cfg.get('qub_pulse_phase_2', cfg.get('qub_phase', 0.0))),
            gain=float(cfg['qub_gain'])/2,
        )

        # Phase reset pulse
        self.add_pulse(
            ch=qub_ch, name="qub_phrst", ro_ch=ro_ch,
            style="const",
            freq=cfg["qub_freq"],
            phase=0,
            gain=0,
            length=0.015,   # us
            phrst=1
        )

        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)

    def _body(self, cfg):
        # ----------------------------
        # Time bookkeeping (ALL in us)
        # ----------------------------
        sigma = float(cfg["qub_sigma"])
        ns    = int(cfg["num_sigma"])
        pulse_len = float(cfg['hold_length']) + sigma*ns   # flat_top + gaussian ramp time (us)

        meas_time = float(cfg["meas_time"])                  # us
        t2_start  = meas_time - float(cfg['qub_delay_fixed']) - pulse_len  # start of pulse_2 (us)

        mode = cfg.get("T2_mode", "T2")  # "T2" or "T2_echo"
        tau_us = float(cfg.get("tau_us", 0.0))               # us, defined as END(p1) -> START(p2)

        # pulse_1 start so that end(p1) is tau_us before start(p2)
        # end(p1) = t1_start + pulse_len
        # start(p2) = t2_start
        # => t2_start - (t1_start + pulse_len) = tau_us
        # => t1_start = t2_start - tau_us - pulse_len
        t1_start = t2_start - tau_us - pulse_len

        # Echo: place echo pulse centered in the gap between end(p1) and start(p2)
        # gap = tau_us
        # start_echo = end(p1) + (gap - pulse_len)/2  (if gap < pulse_len, it'll overlap => user should avoid)
        end_p1 = t1_start + pulse_len
        t_echo_start = end_p1 + 0.5*(tau_us - pulse_len)
        

        # Optional phase reset
        if cfg.get("phase_reset", False):
            self.pulse(ch=cfg["qub_channel"], name="qub_phrst", t=0)

        # Trigger ADC
        self.trigger(ros=[cfg['ro_channel']], pins=[0], t=float(cfg["adc_trig_offset"]))

        # Pulses
        self.pulse(ch=cfg["qub_channel"], name='qub_pulse_1', t=t1_start)

        if mode == "T2_echo":
            self.pulse(ch=cfg["qub_channel"], name='qub_pulse_echo', t=t_echo_start)
        elif mode != "T2":
            raise ValueError(f"Unsupported T2_mode: {mode}")

        self.pulse(ch=cfg["qub_channel"], name='qub_pulse_2', t=t2_start)

        # Readout
        self.pulse(ch=cfg["cav_channel"], name='cav_pulse', t=meas_time)
        self.wait_auto()
        self.delay(float(cfg["relax_delay"]))

      


def get_T2_settings():
    settings = {}
    
    settings['scanname'] = 'T2_meas'
    settings['meas_type'] = 'T2_meas'
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
    
    # settings['T2_guess'] = 10e-6
    # settings['detuning'] = 2e6
    # #Card settings
    # settings['reps'] = 1
    # settings['soft_avgs'] = 5e3
    
    return settings

def meas_T2(soc,soccfg,instruments,settings):
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
    
   # if exp_globals['LO']:
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
        'cav_channel'     : exp_globals['cav_channel']['ID'],
        'qub_channel'     : qub_channel,
        'ro_channel'     : exp_globals['ro_channel']['ID'],

        'nqz_c'           : 2,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'cav_gain'       : exp_settings['cav_gain'],
        'cav_freq'        : exp_settings['cav_freq']/1e6,
        'cav_mixer_freq'  : (exp_settings['cav_freq'] + exp_settings['cav_mixer_detuning'])/1e6,
        
        'nqz_q'           : 2,
        'qub_freq'        : exp_settings['qub_freq']/1e6,
        'qub_mixer_freq'  : (exp_settings['qub_freq']+exp_settings['qub_mixer_detuning'])/1e6,
        'qub_gain'        : exp_settings['qub_gain'],
        
        'qub_sigma'       : q_pulse['sigma'],
        'qub_delay_fixed' : exp_globals['qub_delay_fixed'],
        'tau_us'    : 0, #Placeholder
        'num_sigma'       : q_pulse['num_sigma'],
        'hold_length'     : exp_settings['hold_time'],
        
        'T2_mode'         : exp_settings.get('T2_mode', 'T2'),   # "T2" or "T2_echo"
        'synthetic_detuning': exp_settings.get('synthetic_detuning', 0.0)/1e6,
        
        # default phases (degrees)
        'qub_phase'           : q_pulse['qub_phase'],   # base phase
        'qub_pulse_phase_1'   : q_pulse['qub_phase'],
        'qub_pulse_phase_echo': q_pulse['qub_phase'],
        'qub_pulse_phase_2'   : q_pulse['qub_phase'],

        
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


#    prog = T2_sequence(soccfg,config)
    # meas_start = prog.us2cycles(m_pulse["init_buffer"],ro_ch=0)
    # meas_end = meas_start+prog.us2cycles(m_pulse["meas_window"],ro_ch=0)
    # total_samples = prog.us2cycles(config['readout_length'],ro_ch=0)
    
    ## Set up array of taus and randomize it
    if exp_settings['spacing'] == 'Log':
        tau_list = np.logspace(np.log10(exp_settings['Tau_min']), np.log10(exp_settings['Tau_max']), exp_settings['Tau_points'])
    else:
         tau_list = np.linspace(exp_settings['Tau_min'], exp_settings['Tau_max'], exp_settings['Tau_points'])
    taus = np.round(tau_list, 9)
    
    indices = list(range(len(taus)))
#    np.random.shuffle(indices)
    
    amp_int = np.zeros(len(taus))
    amp_orig = np.zeros(len(taus))

    
    ang_int = np.zeros(len(taus))
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
#         print('Background Trace Complete')
#         I_back = holder[0][0][0]
#         Q_back = holder[1][0][0]
#     else:
# =============================================================================
    I_back, Q_back = 0,0
    
    tstart = time.time()
    
    for tind in indices:
        
        tau = taus[tind]
        print('Tau: {} us'.format(tau))
        
        config['tau_us'] = float(tau)
        
        # For echo mode: require tau_us >= pulse_len so echo can fit in the gap.
        if config['T2_mode'] == "T2_echo":
            pulse_len_us = float(config['hold_length']) + float(config['qub_sigma']) * int(config['num_sigma'])
            if float(tau) < pulse_len_us:
                raise ValueError(f"tau_us={tau:.6g} < pulse_len_us={pulse_len_us:.6g}. Increase tau_min.")


        base_phi = float(config['qub_phase'])  # degrees
        
        if config['T2_mode'] == "T2":
            # synthetic phase accumulation: phi2 = phi1 + 360 * f_syn * tau
            f_syn = float(config.get('synthetic_detuning', 0.0))  # MHz
            dphi = 360.0 * f_syn * float(tau)               # degrees
        
            config['qub_pulse_phase_1'] = base_phi
            config['qub_pulse_phase_2'] = base_phi + dphi
        
            # echo pulse isn't played in T2 mode; phase value doesn't matter, but keep sane:
            config['qub_pulse_phase_echo'] = base_phi
        
        elif config['T2_mode'] == "T2_echo":
            # no synthetic detuning in echo mode: all three pulses same phase
            config['qub_pulse_phase_1'] = base_phi
            config['qub_pulse_phase_echo'] = base_phi
            config['qub_pulse_phase_2'] = base_phi
        
        else:
            raise ValueError(f"Unsupported T2_mode: {config['T2_mode']}")

        prog = T2_sequence(soccfg,reps = exp_settings['reps'], final_delay = None, final_wait = 0, cfg = config)
        holder = prog.acquire(soc, rounds = exp_settings['rounds'], load_pulses=True, progress=False)
        
        
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

        amp_orig[tind] = np.sqrt(I_sig**2 + Q_sig**2)
        ang_orig[tind] = np.degrees(np.arctan2(Q_sig, I_sig))
        
        if first_it:
            tstop = time.time()
            estimate_time(tstart, tstop, len(taus))
            
                
            first_it = False  
            
        #_______________________________________#
    
        fig = plt.figure(1, figsize=(13,8))
        plt.clf()
        plt.subplot(121)
        plt.plot(taus, amp_int)
        plt.suptitle(f'{filename}\nLive T2 data (no fit)')
        plt.xlabel('Tau (us)')
        plt.ylabel('Amplitude')
        plt.subplot(122)
        plt.plot(taus, ang_int)
        plt.xlabel('Tau (us)')
        plt.ylabel('Phase')
        fig.canvas.draw()
        fig.canvas.flush_events()

        # plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)
        # userfuncs.SaveFull(saveDir, filename, ['taus', 'amp_int', 'ang_int'], locals(), expsettings=settings, instruments=instruments)

    #last save at the end    
    plt.savefig(os.path.join(saveDir, filename+'_no_fit.png'), dpi = 150)
    userfuncs.SaveFull(saveDir, filename, ['taus', 'amp_int', 'ang_int'], locals(), expsettings=settings, instruments=instruments)
    
    


    t2 = time.time()
    
    print('Elapsed time: {}'.format(t2-tstart))
    
    # -------------------------
    # Fit + final plot
    # -------------------------
    if config['T2_mode'] == 'T2':
        # ----- T2 fit (oscillatory) -----
        T2_guess_s   = exp_settings['T2_guess_us'] * 1e-6  # s
        amp_guess    = np.max(amp_int) - np.min(amp_int)
        offset_guess = np.mean(amp_int[-10:])
    
        taus_s = taus * 1e-6  # us -> s
    
        freq_guess_hz = exp_settings.get('synthetic_detuning', 0.0)  # Hz
        phi_guess     = 0.0
    
        fit_guess = [T2_guess_s, amp_guess, offset_guess, freq_guess_hz, phi_guess]
        T2_s, amp, offset, freq_hz, phi, fit_x_s, fit_y = fit_T2(taus_s, amp_int, fit_guess)
    
        fig3 = plt.figure(3)
        plt.clf()
        plt.plot(taus, amp_int)
        plt.plot(fit_x_s * 1e6, fit_y)
        plt.title(f'T2: {T2_s*1e6:.3f} us, freq: {freq_hz/1e6:.3f} MHz\n{filename}')
        plt.xlabel('Time (us)')
        plt.ylabel('Amplitude')
        fig3.canvas.draw()
        fig3.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename + '_fit.png'), dpi=150)
    
        # unify outputs/fields for saving + return
        T2 = T2_s
        freq = freq_hz
    
    elif config['T2_mode'] == 'T2_echo':
        # ----- T2_echo fit (no oscillations) -----
        T2_guess_s   = exp_settings['T2_guess_us'] * 1e-6  # s
        amp_guess    = np.max(amp_int) - np.min(amp_int)
        offset_guess = np.mean(amp_int[-10:])
    
        taus_s = taus * 1e-6  # us -> s
        fit_guess = [T2_guess_s, amp_guess, offset_guess]
    
        T2_s, amp, offset, fit_x_s, fit_y = fit_T1(taus_s, amp_int, fit_guess)
    
        fig3 = plt.figure(3)
        plt.clf()
        plt.plot(taus, amp_int)
        plt.plot(fit_x_s * 1e6, fit_y)
        plt.title(f'T2_echo: {T2_s*1e6:.3f} us\n{filename}')
        plt.xlabel('Time (us)')
        plt.ylabel('Amplitude')
        fig3.canvas.draw()
        fig3.canvas.flush_events()
        plt.savefig(os.path.join(saveDir, filename + '_fit.png'), dpi=150)
    
        # unify outputs/fields for saving + return
        T2 = T2_s
        freq = 0.0
        phi = 0.0  # for SaveFull consistency
    
    else:
        raise ValueError(f"Unsupported T2_mode: {config['T2_mode']}")
    
    # -------------------------
    # Save + return
    # -------------------------
    userfuncs.SaveFull(
        saveDir, filename,
        ['taus','ang_int','amp_int','amp_orig','ang_orig','amp','offset','freq','phi','fit_guess'],
        locals(),
        expsettings=settings,
        instruments=instruments
    )

    return T2, freq, taus, amp_int

    


   
    

        