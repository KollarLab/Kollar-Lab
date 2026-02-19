# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: kollarlab
"""

from qick.asm_v2 import AveragerProgramV2, QickSweep1D


import numpy as np
import time
import userfuncs
from utility.measurement_helpers import estimate_time
import logging
from utility.plotting_tools import simplescan_plot


class CW_spec(AveragerProgramV2):
    def _initialize(self, cfg):
        
        ro_ch = cfg['ro_channel']
        gen_ch = cfg['cav_channel']
        qub_ch = cfg["qub_channel"] # These are just channel ID, not the whole dict.
        
        #self.add_loop("freq_sweep", cfg["freq_points"]) # if you want the QickSweep1D thing, turn this back

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=cfg["nqz_c"], mixer_freq=cfg['cav_mixer_freq'],ro_ch = ro_ch)
        self.declare_gen(ch=cfg["qub_channel"], nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq'])
        
        self.declare_readout(ch=ro_ch, length=cfg['readout_length'])
        
        
        self.add_readoutconfig(ch=ro_ch, name="myro",
                           freq=cfg['cav_freq'],
                           gen_ch=gen_ch,
                           outsel='product')
        
        # configure cavity pulse
        self.add_pulse(ch=gen_ch, name="cav_pulse", ro_ch=ro_ch,
                       style="const",
                       freq=cfg['cav_freq'],
                       length= cfg["cav_pulse_len"],
                       phase=cfg['cav_phase'],
                       gain=cfg['meas_gain'],
                       mode='periodic'
                      )
        
        # configure qubit pulse
        self.add_pulse(ch=qub_ch, name="qub_pulse", ro_ch=ro_ch,
               style="const",
               freq=cfg['qub_freq'],
               length= cfg['qub_pulse_len'],
               phase=cfg['qub_phase'],
               gain=cfg['qub_gain'],
               mode='periodic'
              )
        
        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)

    
    def _body(self, cfg):
        
        #Sets off the ADC
        self.trigger(ros=[cfg['ro_channel']],
                    pins=[0],
                    t=self.cfg['adc_trig_offset'])
        
        self.pulse(ch=cfg["qub_channel"],name='qub_pulse',t=0.0)
        self.pulse(ch=cfg["cav_channel"],name='cav_pulse',t=0.0)
        #self.wait_auto()
        # DO NOT rely on wait_auto() here because wait_auto() often cannot advance the timeline the way you think during a sweep loop
        # since for periodic pulses, the generator is considered “running continuously” (it doesn’t have a finite end time)
        self.delay_auto(cfg['adc_trig_offset'] + cfg['readout_length'] + 4.0)
        #self.delay(self.cfg["relax_delay"])
        


def get_cw_spec_settings():
    settings = {}
    
    settings['scanname'] = 'continuous_power_scan'
    settings['meas_type'] = 'CWSpec'
    
    settings['cav_freq'] = 1e9
    settings['meas_gain'] = 1000
    settings['cav_pulse_len'] = 10
    settings['readout_length'] = 900

    settings['qub_gain_start']     = 4000
    settings['qub_gain_step']      = 2000
    settings['qub_gain_points']    = 1
    settings['qub_pulse_len'] = 20
    
    #Sweep parameters
    settings['freq_start']   = 4e9  
    settings['freq_step']    = 0.5e9
    settings['freq_points']  = 6
    settings['freq_stop'] = settings['freq_start'] + settings['freq_step']*(settings['freq_points']-1)


    #Card settings
    settings['reps'] = 1
    settings['rounds'] = 1#5e3
    
    settings['filter'] = 'all_filter'
    settings['cav_mixer_detuning'] = 200e6  # Hz
    settings['qub_mixer_detuning'] = 200e6  # Hz

    
    return settings

def cw_spec(soc,soccfg,instruments,settings):
    
    # suppresses sum buffer overflow warning
    logging.getLogger("qick").setLevel(logging.ERROR)

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
    
    if exp_settings['power_sweep_mode'] == 'gain':

        config = {
            'cav_channel'     : exp_globals['cav_channel']['ID'],
            'qub_channel'     : qub_channel,
            'ro_channel'     : exp_globals['ro_channel']['ID'],
    
            'nqz_c'           : 2,
            'cav_phase'       : m_pulse['cav_phase'],
            'cav_pulse_len'   : exp_settings['cav_pulse_len'],
            'readout_length'  : exp_settings['readout_length'],
            'meas_gain'       : exp_settings['meas_gain'],
            'cav_freq'        : exp_settings['cav_freq']/1e6,
            'cav_mixer_freq'  : (exp_settings['cav_freq'] + exp_settings['cav_mixer_detuning'])/1e6,
            
            'nqz_q'           : 2,
            'qub_phase'       : q_pulse['qub_phase'],
    
            #'qub_freq'        : QickSweep1D("freq_sweep", exp_settings['freq_start']/1e6, exp_settings['freq_stop']/1e6), # if you want to use QickSweep1D, bring this back
            'qub_freq'        : 0, #placeholder, freq loop
            'freq_points'     : exp_settings['freq_points'],
            'qub_gain'        : 0, # PLACEHOLDER, gain loop
            #'qub_mixer_freq'  : (exp_settings['freq_start']+exp_settings['qub_mixer_detuning'])/1e6, # bring this back if using QickSweep1D
            'qub_mixer_freq'  : 0, #placeholder
    
            'qub_pulse_len'   : exp_settings['qub_pulse_len'],
            'adc_trig_offset' : m_pulse['emp_delay'],# + m_pulse['meas_pos'] # similar to CW_trans_RF216, this meas_pos can be effectively set to zero
            'reps'            : exp_settings['reps'],
            'rounds'          : exp_settings['rounds'],
            
            'cav_atten'       : exp_globals['cav_channel']['Atten_1'] + exp_globals['cav_channel']['Atten_2'],
            'qub_atten'       : qub_ch['Atten_1'] + qub_ch['Atten_2']
            }
        
        cav_ch = exp_globals['cav_channel']
        #qub_ch = exp_globals['qub_channel'] defined eralier in the script
        ro_ch  = exp_globals['ro_channel']
        # Set attenuator on DAC.
        soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
        soc.rfb_set_gen_rf(qub_ch['ID'], qub_ch['Atten_1'], qub_ch['Atten_2'])
        # Set attenuator on ADC.
        soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])
        
# ============================================================================= #bring this back if using QickSweep1D
#         if exp_settings['filter'] == 'all_filter':
#             soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
#             soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
#             
#             center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
#             soc.rfb_set_gen_filter(config['qub_channel'], fc=center_freq, ftype='bandpass', bw=qub_ch['BW'])
#             
#         elif exp_settings['filter'] == 'no_ro_filter':
#             soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
#             soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
#             
#             center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
#             soc.rfb_set_gen_filter(config['qub_channel'], fc=center_freq, ftype='bypass')
#             
#         elif exp_settings['filter'] == 'no_qubit_filter':
#             soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
#             soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
#             
#             center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
#             soc.rfb_set_gen_filter(config['qub_channel'], fc=center_freq, ftype='bypass')
#             
#         elif exp_settings['filter'] == 'no_filter':
#             soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bypass')
#             soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
#             
#             center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
#             soc.rfb_set_gen_filter(config['qub_channel'], fc=center_freq, ftype='bypass')
#             
#         else:
#             print('Please select one option from:')
#             print('\'all_filter\', \'no_qubit_filter\', and \'no_filter\'')
#             return
# =============================================================================
    
        
        stamp    = userfuncs.timestamp()
        saveDir  = userfuncs.saveDir(settings)
        filename = exp_settings['scanname'] + '_' + stamp
        
        f0_start = exp_settings['freq_start']
        f0_step = exp_settings['freq_step']
        expts = exp_settings['freq_points']
        fpts = np.arange(0,expts)*f0_step+f0_start #in Hz
        
        g_start = exp_settings['qub_gain_start']
        g_step = exp_settings['qub_gain_step']
        g_expts = exp_settings['qub_gain_points']
        gpts = exp_settings["qub_gain_start"] + (exp_settings["qub_gain_step"] * np.arange(exp_settings["qub_gain_points"]))
        
        powerdat = np.zeros((len(gpts), len(fpts)))
        phasedat = np.zeros((len(gpts), len(fpts)))
        Is = np.zeros((len(gpts), len(fpts)))
        Qs = np.zeros((len(gpts), len(fpts)))
    
        
        t_start = time.time()
        
        # qubit gain sweep
        for g in range(0,len(gpts)):
            print("Current Qubit Gain: " + str(gpts[g]) + ", Max Gain: " + str(gpts[-1]))
            config["qub_gain"] = gpts[g]
            
            for f in range(0,len(fpts)): ###delete this loop if you want to use QickSweep1D
                config['qub_freq'] = fpts[f]/1e6
                config['qub_mixer_freq'] = fpts[f]/1e6 + exp_settings['qub_mixer_detuning']/1e6
            
                if exp_settings['filter'] == 'all_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
                    
                    center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
                    soc.rfb_set_gen_filter(config['qub_channel'], fc=fpts[f]/1e9, ftype='bandpass', bw=qub_ch['BW'])
                    
                elif exp_settings['filter'] == 'no_ro_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
                    
                    center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
                    soc.rfb_set_gen_filter(config['qub_channel'], fc=fpts[f]/1e9, ftype='bypass')
                    
                elif exp_settings['filter'] == 'no_qubit_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
                    
                    center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
                    soc.rfb_set_gen_filter(config['qub_channel'], fc=fpts[f]/1e9, ftype='bypass')
                    
                elif exp_settings['filter'] == 'no_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bypass')
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
                    
                    center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
                    soc.rfb_set_gen_filter(config['qub_channel'], fc=fpts[f]/1e9, ftype='bypass')
                    
                else:
                    print('Please select one option from:')
                    print('\'all_filter\', \'no_qubit_filter\', and \'no_filter\'')
                    return
                
                prog = CW_spec(soccfg,reps = exp_settings['reps'], final_delay = None, final_wait = 0, cfg = config)
                holder = prog.acquire(soc, rounds = exp_settings['rounds'], progress=False) # shape (1, 1, 2)
                
        # =============================================================================
        #         print(type(iq_list), np.shape(iq_list)) 
        #         print("iq_list[0] type/shape:", type(iq_list[0]), np.shape(iq_list[0]))
        #         try:
        #             print("iq_list[0][0] type/shape:", type(iq_list[0][0]), np.shape(iq_list[0][0]))
        #         except Exception as e:
        #             print("iq_list[0][0] indexing failed:", e)
        # =============================================================================
                # debugging
                
                iq = holder[0] 
                soc.reset_gens()
        
                I_full = iq[:,0]
                Q_full = iq[:,1]
                
                # debugging
        # =============================================================================
        #         print("I first/last:", I_fulls[0], I_fulls[-1])
        #         print("nonzero count:", np.count_nonzero(np.abs(I_fulls) + np.abs(Q_fulls)))
        # =============================================================================
        
                
                #print (np.shape(I_fulls))
                
                powerdat[g,f] = np.sqrt(I_full**2 + Q_full**2)
                phasedat[g,f] = np.degrees(np.arctan2(Q_full,I_full))
                Is[g,f] = I_full
                Qs[g,f] = Q_full
                
            if g == 0:
                t_stop = time.time()
                estimate_time(t_start, t_stop, len(gpts))
    
        full_data = {}
    
        full_data['xaxis']  = fpts/1e9
        full_data['mags']   = powerdat
        full_data['phases'] = phasedat
        full_data['Is']     = Is
        full_data['Qs']     = Qs
        full_data['cav_atten'] = config['cav_atten']
        full_data['qub_atten'] = config['qub_atten']
        
        plot_data = {}
        plot_data['xaxis']  = fpts/1e9
        plot_data['mags']   = powerdat[0:g+1]
        plot_data['phases'] = phasedat[0:g+1]
    
        single_data = {}
        single_data['xaxis'] = fpts/1e9
        single_data['mag']   = powerdat[g]
        single_data['phase'] = phasedat[g]
        
        yaxis  = gpts[0:g+1]
        labels = ['Freq (GHz)', 'Gain (DAC a.u.)']
        identifier = f"Cav_Atten={config['cav_atten']} dB, Qub_Atten={config['qub_atten']} dB"
        
        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=1, IQdata = False)
        
        userfuncs.SaveFull(saveDir, filename, ['gpts','fpts','full_data','filename'],
        locals(), expsettings=settings, instruments={})
    
        data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data}
    
        return data,prog
    
    
    elif exp_settings['power_sweep_mode'] == 'attenuation':
        config = {
            'cav_channel'     : exp_globals['cav_channel']['ID'],
            'qub_channel'     : qub_channel,
            'ro_channel'     : exp_globals['ro_channel']['ID'],
    
            'nqz_c'           : 2,
            'cav_phase'       : m_pulse['cav_phase'],
            'cav_pulse_len'   : exp_settings['cav_pulse_len'],
            'readout_length'  : exp_settings['readout_length'],
            'meas_gain'       : exp_settings['meas_gain'],
            'cav_freq'        : exp_settings['cav_freq']/1e6,
            'cav_mixer_freq'  : (exp_settings['cav_freq'] + exp_settings['cav_mixer_detuning'])/1e6,
            
            'nqz_q'           : 2,
            'qub_phase'       : q_pulse['qub_phase'],
    
            #'qub_freq'        : QickSweep1D("freq_sweep", exp_settings['freq_start']/1e6, exp_settings['freq_stop']/1e6),
            'qub_freq'        : 0, #placeholder, freq loop
            'freq_points'     : exp_settings['freq_points'],
            'qub_gain'        : exp_settings['qub_gain'],
            #'qub_mixer_freq'  : (exp_settings['freq_start']+exp_settings['qub_mixer_detuning'])/1e6, # bring this back if using QickSweep1D
            'qub_mixer_freq'  : 0, #placeholder
    
            'qub_pulse_len'   : exp_settings['qub_pulse_len'],
            'adc_trig_offset' : m_pulse['emp_delay'],# + m_pulse['meas_pos'] # similar to CW_trans_RF216, this meas_pos can be effectively set to zero
            'reps'            : exp_settings['reps'],
            'rounds'          : exp_settings['rounds'],
            
            'cav_atten'       : exp_globals['cav_channel']['Atten_1'] + exp_globals['cav_channel']['Atten_2'],
            'qub_atten'       : 0, #placeholder, updated in the loop
            }
        
        cav_ch = exp_globals['cav_channel']
        #qub_ch = exp_globals['qub_channel'] defined eralier in the script
        ro_ch  = exp_globals['ro_channel']
        # Set attenuator on DAC.
        soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
        
        # Set attenuator on ADC.
        soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])
        
        
# =============================================================================
#         if exp_settings['filter'] == 'all_filter':
#             soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
#             soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
#             
#             center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
#             soc.rfb_set_gen_filter(config['qub_channel'], fc=center_freq, ftype='bandpass', bw=qub_ch['BW'])
#             
#         elif exp_settings['filter'] == 'no_ro_filter':
#             soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
#             soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
#             
#             center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
#             soc.rfb_set_gen_filter(config['qub_channel'], fc=center_freq, ftype='bypass')
#             
#         elif exp_settings['filter'] == 'no_qubit_filter':
#             soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
#             soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
#             
#             center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
#             soc.rfb_set_gen_filter(config['qub_channel'], fc=center_freq, ftype='bypass')
#             
#         elif exp_settings['filter'] == 'no_filter':
#             soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bypass')
#             soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
#             
#             center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
#             soc.rfb_set_gen_filter(config['qub_channel'], fc=center_freq, ftype='bypass')
#             
#         else:
#             print('Please select one option from:')
#             print('\'all_filter\', \'no_qubit_filter\', and \'no_filter\'')
#             return
# =============================================================================
    
        
        stamp    = userfuncs.timestamp()
        saveDir  = userfuncs.saveDir(settings)
        filename = exp_settings['scanname'] + '_' + stamp
        
        f0_start = exp_settings['freq_start']
        f0_step = exp_settings['freq_step']
        expts = exp_settings['freq_points']
        fpts = np.arange(0,expts)*f0_step+f0_start
        
        a_start = exp_settings['qub_atten_start']
        a_step = exp_settings['qub_atten_step']
        a_expts = exp_settings['qub_atten_points']
        apts = exp_settings["qub_atten_start"] + (exp_settings["qub_atten_step"] * np.arange(exp_settings["qub_atten_points"]))
        
        powerdat = np.zeros((len(apts), len(fpts)))
        phasedat = np.zeros((len(apts), len(fpts)))
        Is = np.zeros((len(apts), len(fpts)))
        Qs = np.zeros((len(apts), len(fpts)))
    
        
        t_start = time.time()
        
        # qubit attenuation sweep
        for a in range(0,len(apts)):
            print("Current Qubit Attenuation: " + str(apts[a]) + " dB, Final Attenuation: " + str(apts[-1]) + ' dB')
            soc.rfb_set_gen_rf(qub_ch['ID'], apts[a]/2, apts[a]/2)
            
            for f in range(0,len(fpts)): ###delete this loop if you want to use QickSweep1D
                config['qub_freq'] = fpts[f]/1e6
                config['qub_mixer_freq'] = fpts[f]/1e6 + exp_settings['qub_mixer_detuning']/1e6
                
                if exp_settings['filter'] == 'all_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
                    
                    center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
                    soc.rfb_set_gen_filter(config['qub_channel'], fc=fpts[f]/1e9, ftype='bandpass', bw=qub_ch['BW'])
                    
                elif exp_settings['filter'] == 'no_ro_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
                    
                    center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
                    soc.rfb_set_gen_filter(config['qub_channel'], fc=fpts[f]/1e9, ftype='bypass')
                    
                elif exp_settings['filter'] == 'no_qubit_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])
                    
                    center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
                    soc.rfb_set_gen_filter(config['qub_channel'], fc=fpts[f]/1e9, ftype='bypass')
                    
                elif exp_settings['filter'] == 'no_filter':
                    soc.rfb_set_gen_filter(config['cav_channel'], fc=config['cav_freq']/1000, ftype='bypass')
                    soc.rfb_set_ro_filter(config['ro_channel'], fc=config['cav_freq']/1000, ftype='bypass')
                    
                    center_freq = (exp_settings['freq_start']+exp_settings['freq_stop'])/2e9
                    soc.rfb_set_gen_filter(config['qub_channel'], fc=fpts[f]/1e9, ftype='bypass')
                    
                else:
                    print('Please select one option from:')
                    print('\'all_filter\', \'no_qubit_filter\', and \'no_filter\'')
                    return
            
                prog = CW_spec(soccfg,reps = exp_settings['reps'], final_delay = None, final_wait = 0, cfg = config)
                holder = prog.acquire(soc, rounds = exp_settings['rounds'], progress=False) # shape (1, 1, 2)
                
        # =============================================================================
        #         print(type(iq_list), np.shape(iq_list)) 
        #         print("iq_list[0] type/shape:", type(iq_list[0]), np.shape(iq_list[0]))
        #         try:
        #             print("iq_list[0][0] type/shape:", type(iq_list[0][0]), np.shape(iq_list[0][0]))
        #         except Exception as e:
        #             print("iq_list[0][0] indexing failed:", e)
        # =============================================================================
                # debugging
                
                iq = holder[0] 
                soc.reset_gens()
        
                I_full = iq[:,0] 
                Q_full = iq[:,1]
                
                # debugging
        # =============================================================================
        #         print("I first/last:", I_fulls[0], I_fulls[-1])
        #         print("nonzero count:", np.count_nonzero(np.abs(I_fulls) + np.abs(Q_fulls)))
        # =============================================================================
        
                
                #print (np.shape(I_fulls))
                
                powerdat[a,f] = np.sqrt(I_full**2 + Q_full**2)
                phasedat[a,f] = np.degrees(np.arctan2(Q_full,I_full))
                Is[a,f] = I_full
                Qs[a,f] = Q_full
            
            if a == 0:
                t_stop = time.time()
                estimate_time(t_start, t_stop, len(apts))
    
        full_data = {}
    
        full_data['xaxis']  = fpts/1e9
        full_data['mags']   = powerdat
        full_data['phases'] = phasedat
        full_data['Is']     = Is
        full_data['Qs']     = Qs
        full_data['cav_atten'] = config['cav_atten']
        full_data['qub_gain'] = config['qub_gain']
        
        plot_data = {}
        plot_data['xaxis']  = fpts/1e9
        plot_data['mags']   = powerdat[0:a+1]
        plot_data['phases'] = phasedat[0:a+1]
    
        single_data = {}
        single_data['xaxis'] = fpts/1e9
        single_data['mag']   = powerdat[a]
        single_data['phase'] = phasedat[a]
        
        yaxis  = apts[0:a+1]
        labels = ['Freq (GHz)', 'Attenuation (dB)']
        identifier = f"Cav_Atten={config['cav_atten']} dB, Qub_Gain={config['qub_gain']}"
        
        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=1, IQdata = False)
        
        userfuncs.SaveFull(saveDir, filename, ['apts','fpts','full_data','filename'],
        locals(), expsettings=settings, instruments={})
    
        data = {'saveDir': saveDir, 'filename': filename, 'full_data': full_data}
    
        return data,prog