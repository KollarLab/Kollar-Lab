# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:52:24 2022

@author: kollarlab

Modified by: jhyang
"""

import os
from qick.asm_v2 import AveragerProgramV2, QickSweep1D


import numpy as np
import time
import userfuncs
import matplotlib.pyplot as plt
from utility.measurement_helpers import estimate_time
import logging
from utility.plotting_tools import simplescan_plot
import utility.plotting_tools as plots
from utility.userfits_v2 import fit_model

class CavitySweep(AveragerProgramV2):
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
        self.add_pulse(ch=gen_ch, name="mypulse", ro_ch=ro_ch,
                       style="const",
                       freq=cfg['cav_freq'],
                       length= cfg["meas_window"],
                       phase=cfg['cav_phase'],
                       gain=cfg['meas_gain'],
                      )
        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)
        
    def _body(self, cfg):
        
        self.pulse(ch=self.cfg['cav_channel'], name="mypulse", t=self.cfg["meas_time"])
        self.trigger(ros=[self.cfg['ro_channel']], pins=[0], t=self.cfg['adc_trig_offset'])
        self.wait_auto()
        self.delay(self.cfg["relax_delay"])


class Quasi_CW(AveragerProgramV2):
    def _initialize(self,cfg): 
        ro_ch  = cfg['ro_channel']
        gen_ch = cfg["cav_channel"]
        qub_ch = cfg["qub_channel"] # These are just channel ID, not the whole dict.

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"], nqz=cfg["nqz_c"], mixer_freq=cfg['cav_mixer_freq'],ro_ch = ro_ch)
        
        if cfg.get('qub_mixer_freq_fixed') is not None:
            self.declare_gen(ch=cfg["qub_channel"], nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq_fixed'])
            
        else:
            self.declare_gen(ch=cfg["qub_channel"], nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq'])
        
        self.declare_readout(ch=ro_ch, length=cfg['readout_length'])
        
        self.add_loop("freq_sweep", cfg["freq_points"])
        
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
               length= cfg['quasi_CW_len'],
               phase=cfg['qub_phase'],
               gain=cfg['qub_gain'],
               phrst=1
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
        
        # configure the readout lengths and downconversion frequencies (ensuring it is an available DAC frequency)
        
        ###Start sweep definition
        
        # self.qub_r_freq = self.get_gen_reg(qub_ch,"freq")
        
        
        # #freq_start_reg = self.freq2reg(cfg["freq_start"],gen_ch=qub_ch)
        # #freq_stop_reg = self.freq2reg(cfg["freq_stop"],gen_ch=qub_ch)
        
        # self.add_sweep(QickSweep(self, self.qub_r_freq, cfg['freq_start'], cfg['freq_stop'], cfg["freq_points"]))
        
    
    def _body(self, cfg):
# =============================================================================
#         qub_ch = self.cfg["qub_channel"]
#         self.reset_phase(gen_ch = [qub_ch], t=0)
# =============================================================================
        
        #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.
        sigma = cfg["qub_sigma"]
        num_sigma = cfg["num_sigma"]
        
        pulse_len = cfg['quasi_CW_len'] + int(num_sigma*sigma)

        offset = cfg["adc_trig_offset"]
        meas_time = self.cfg["meas_time"]
        ex_time = meas_time - cfg['qub_delay'] - pulse_len
        
        
# =============================================================================
#         # Optional phase reset behavior
#         if cfg.get("phase_reset", False):
#             self.pulse(ch=cfg["qub_channel"], name="qub_phrst", t=0)
# =============================================================================
        
        #Sets off the ADC
        self.trigger(ros=[cfg['ro_channel']],
                    pins=[0],
                    t=offset)
        
        self.pulse(ch=cfg["qub_channel"],name='qub_pulse',t=ex_time)
        self.pulse(ch=cfg["cav_channel"],name='cav_pulse',t=meas_time)
        self.wait_auto()
        self.delay(self.cfg["relax_delay"])

def get_quasi_cw_flux_settings():
    fullsettings = {}
    settings = {}
    autoscan_settings = {}

    settings['scanname'] = 'QCW_flux_scan'
    settings['meas_type'] = 'QuasiCW_Flux'

    settings['qub_drive_index'] = 'D1'
    settings['cav_freq'] = 6.325750e9#6.10208e9 
    settings['cav_mixer_detuning'] = 200e6
    settings['qub_mixer_detuning'] = -250e6 # currently not used
    settings['cav_gain'] = 0.7

    settings['qub_mixer_freq_fixed'] = 4.65e9 # set to be fixed for stable phase of the mixer

    settings['fit'] = False

    settings['qub_gain']     = 0.015
    settings['quasi_CW_len'] = 30 #us
    settings['phase_reset'] = False # does not do anything now
    settings['filter'] = 'all_filter'

    #Sweep Parameters
    settings['freq_start']      = 4.4e9-50e6/2
    settings['freq_stop']       = 4.4e9+50e6/2
    settings['freq_points']     = 141

    #ADC settings
    settings['reps']      = 3000
    settings['rounds']  =1

    
    #Cavity Sweep parameters
    autoscan_settings['freq_start']   = 7e9 
    autoscan_settings['freq_step']    = 100e6
    autoscan_settings['freq_points']  = 21
    autoscan_settings['freq_stop']       = 7e9 + 100e6 * 21

    autoscan_settings['gain']  = 0.7

    #Card settings
    autoscan_settings['reps'] = 1e3
    autoscan_settings['rounds'] = 1

    fullsettings['spec'] = settings
    fullsettings['autoscan'] = autoscan_settings

    return fullsettings

def quasi_cw_flux(soc,soccfg,instruments,settings):
    
    # suppresses sum buffer overflow warning
    logging.getLogger("qick").setLevel(logging.ERROR)

    SRS = instruments['DCsupply']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    spec_set = exp_settings['spec']
    autoscan_set = exp_settings['autoscan']

    qub_drive_index = exp_settings['spec']['qub_drive_index']
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
        
    
    soc.reset_gens()
    
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['spec']['scanname'] + '_' + stamp

    config_trans = {
        'cav_channel'     : exp_globals['cav_channel']['ID'],
        'ro_channel'     : exp_globals['ro_channel']['ID'],

        'nqz_c'           : 2,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'meas_gain'       : exp_settings['autoscan']['gain'],
        'cav_atten'       : exp_globals['cav_channel']['Atten_1'] + exp_globals['cav_channel']['Atten_2'],
        'cav_freq'        : 6000, #Placeholder, MHz
        'mixer_freq'      : 6000, #Placeholder, MHz 
        
        'ramp_len'        :5e-3, #Placeholder, us
        'flat_len'        :0.50, #Placeholder, us
        
        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay']
        }

    config_qubit = {
        'cav_channel'     : exp_globals['cav_channel']['ID'],
        'qub_channel'     : qub_channel, # just the channel ID here
        'ro_channel'     : exp_globals['ro_channel']['ID'],

        'nqz_c'           : 2,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'cav_gain'       : exp_settings['autoscan']['gain'],
        'cav_freq'        : (exp_settings['spec']['cav_freq'])/1e6,
        'cav_mixer_freq'  : (exp_settings['spec']['cav_freq'] + exp_settings['spec']['cav_mixer_detuning'])/1e6,

        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],

        'qub_freq'        : QickSweep1D("freq_sweep", exp_settings['spec']['freq_start']/1e6, exp_settings['spec']['freq_stop']/1e6),
        # 'freq_start'      : exp_settings['freq_start']/1e6,
        # 'freq_stop'       : exp_settings['freq_stop']/1e6,
        'freq_points'     : exp_settings['spec']['freq_points'],
        'qub_gain'        : exp_settings['spec']['qub_gain'],
        'qub_mixer_freq'  : (exp_settings['spec']['freq_start']+exp_settings['spec']['qub_mixer_detuning'])/1e6,
        'qub_mixer_freq_fixed' : exp_settings['spec']['qub_mixer_freq_fixed']/1e6,      # this is useful for maintaining a fixed phase offset between two channels

        'qub_sigma'       : q_pulse['sigma'],
        'qub_delay'       : exp_globals['qub_delay_fixed'],
        'num_sigma'       : q_pulse['num_sigma'],
        'quasi_CW_len'    : exp_settings['spec']['quasi_CW_len'],

        'readout_length'  : m_pulse['meas_window'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : exp_settings['spec']['reps'],
        # 'soft_avgs'       : exp_settings['soft_avgs']
        
        'phase_reset'     : exp_settings['spec']['phase_reset']
        }

    
    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = spec_set['scanname'] + '_' + stamp


    #set voltage sweep
    start_voltage = spec_set['start_voltage']
    stop_voltage  = spec_set['stop_voltage']
    voltage_points = spec_set['voltage_points']
    voltages = np.round(np.linspace(start_voltage, stop_voltage, voltage_points),6)
    max_voltage = 10#3.5
    if np.max(voltages) > max_voltage:
        raise ValueError('max voltage too large!')
    else:
        settings['voltages'] = voltages

    SRS.Output = 'On'
    
    f_start_trans = autoscan_set['freq_start']
    f_stop_trans  = autoscan_set['freq_stop']
    expts_trans   = autoscan_set['freq_points']

    fpts_trans = np.linspace(f_start_trans,f_stop_trans,expts_trans)

    #fpts_trans = np.arange(0,expts_trans)*f_step_trans+f_start_trans

    trans_mags   = np.zeros((voltage_points, autoscan_set['freq_points']))
    trans_phases = np.zeros((voltage_points, autoscan_set['freq_points']))

    raw_trans = {}
    raw_trans['Is'] = np.zeros((voltage_points, autoscan_set['freq_points']))
    raw_trans['Qs'] = np.zeros((voltage_points, autoscan_set['freq_points']))   

###################

    mags = np.zeros((voltage_points, spec_set['freq_points']))
    phases = np.zeros((voltage_points, spec_set['freq_points']))

    raw_spec = {}
    raw_spec['Is'] = np.zeros((voltage_points, spec_set['freq_points']))
    raw_spec['Qs'] = np.zeros((voltage_points, spec_set['freq_points']))
    
    t_start = time.time()
    
    identifier = 'Cav Gain : ' + str(autoscan_set['gain'])  + ' au'

    cav_ch = exp_globals['cav_channel']
    ro_ch  = exp_globals['ro_channel']
    # Set attenuator on DAC.
    soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
    # Set attenuator on ADC.
    soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])

    # voltage sweep
    for vind in range(len(voltages)):
        voltage = voltages[vind]
        print('Voltage: {}, final voltage: {}'.format(voltage, voltages[-1]))
        
        SRS.voltage_ramp(voltage)
        time.sleep(0.1)

        print('trans')

        for tfind in range(len(fpts_trans)):

            config_trans['cav_freq'] = fpts_trans[tfind]/1e6
            config_trans['mixer_freq'] = config_trans['cav_freq'] + autoscan_set['mixer_detuning']/1e6
            
            prog = CavitySweep(soccfg,reps = exp_settings['autoscan']['reps'], final_delay = None, final_wait = 0, cfg =config_trans)
            setfilters_trans(soc, config_trans, exp_settings, exp_globals, config_trans['cav_freq'])

            holder = prog.acquire(soc,rounds = exp_settings['autoscan']['rounds'], load_pulses=True, progress=False)
            #print(holder)
            iq = holder[0][0]
            trans_I = iq[0] 
            trans_Q = iq[1]
            soc.reset_gens()

            mag = np.sqrt(trans_I**2 + trans_Q**2)
            phase = np.arctan2(trans_Q, trans_I)*180/np.pi

            trans_mags[vind,tfind] = mag
            trans_phases[vind,tfind] = phase
            raw_trans['Is'][vind,tfind] = trans_I
            raw_trans['Qs'][vind,tfind] = trans_Q

        hanger = exp_globals['hanger'] #"Fitting" cav freq

        if autoscan_set['fit_cav']:
            try:
                print("Fitting Lorenzian to Cavity")
                config_qubit['cav_freq'] = fit_model(fpts_trans, trans_mags[vind], 'lorenz')['center']/1e6
            except:
                print("Fitting Lorenzian Failed, taking extrema...")
                if hanger:
                    config_qubit['cav_freq'] = fpts_trans[np.argmin(trans_mags[vind])]/1e6
                else:
                    config_qubit['cav_freq'] = fpts_trans[np.argmax(trans_mags[vind])]/1e6
        else:
            if hanger:
                    config_qubit['cav_freq'] = fpts_trans[np.argmin(trans_mags[vind])]/1e6
            else:
                config_qubit['cav_freq'] = fpts_trans[np.argmax(trans_mags[vind])]/1e6

        config_qubit['cav_mixer_freq'] = config_qubit['cav_freq'] + autoscan_set['mixer_detuning']/1e6

        print('spec, cav freq: {}'.format(config_qubit['cav_freq']))

        # qubit frequency sweep
        cav_ch = exp_globals['cav_channel']
        #qub_ch = exp_globals['qub_channel'] defined eralier in the script
        ro_ch  = exp_globals['ro_channel']
        # Set attenuator on DAC.
        soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
        soc.rfb_set_gen_rf(qub_ch['ID'], qub_ch['Atten_1'], qub_ch['Atten_2'])
        # Set attenuator on ADC.
        soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])
        
        if exp_settings['spec']['filter'] == 'all_filter':
            soc.rfb_set_gen_filter(config_qubit['cav_channel'], fc=config_qubit['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
            soc.rfb_set_ro_filter(config_qubit['ro_channel'], fc=config_qubit['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])

            center_freq = (exp_settings['spec']['freq_start']+exp_settings['spec']['freq_stop'])/2e9
            soc.rfb_set_gen_filter(config_qubit['qub_channel'], fc=center_freq, ftype='bandpass', bw=qub_ch['BW'])

        elif exp_settings['spec']['filter'] == 'no_ro_filter':
            soc.rfb_set_gen_filter(config_qubit['cav_channel'], fc=config_qubit['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
            soc.rfb_set_ro_filter(config_qubit['ro_channel'], fc=config_qubit['cav_freq']/1000, ftype='bypass')

            center_freq = (exp_settings['spec']['freq_start']+exp_settings['spec']['freq_stop'])/2e9
            soc.rfb_set_gen_filter(config_qubit['qub_channel'], fc=center_freq, ftype='bypass')

        elif exp_settings['spec']['filter'] == 'no_qubit_filter':
            soc.rfb_set_gen_filter(config_qubit['cav_channel'], fc=config_qubit['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
            soc.rfb_set_ro_filter(config_qubit['ro_channel'], fc=config_qubit['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])

            center_freq = (exp_settings['spec']['freq_start']+exp_settings['spec']['freq_stop'])/2e9
            soc.rfb_set_gen_filter(config_qubit['qub_channel'], fc=center_freq, ftype='bypass')

        elif exp_settings['spec']['filter'] == 'no_filter':
            soc.rfb_set_gen_filter(config_qubit['cav_channel'], fc=config_qubit['cav_freq']/1000, ftype='bypass')
            soc.rfb_set_ro_filter(config_qubit['ro_channel'], fc=config_qubit['cav_freq']/1000, ftype='bypass')

            center_freq = (exp_settings['spec']['freq_start']+exp_settings['spec']['freq_stop'])/2e9
            soc.rfb_set_gen_filter(config_qubit['qub_channel'], fc=center_freq, ftype='bypass')

        else:
            print('Please select one option from:')
            print('\'all_filter\', \'no_qubit_filter\', and \'no_filter\'')
            return
            
        prog = Quasi_CW(soccfg,reps = exp_settings['spec']['reps'], final_delay = None, final_wait = 0, cfg = config_qubit)
        rep_period = config_qubit['adc_trig_offset'] + config_qubit['readout_length'] + config_qubit['relax_delay']


        projected_time = exp_settings['spec']['reps']*exp_settings['spec']['rounds']*rep_period/1e6
        print("Projected Time: " + str(projected_time))
        
        
        iq_list = prog.acquire(soc, rounds = exp_settings['spec']['rounds'], load_pulses=True, progress=False)
        #iq_list = prog.acquire(soc, reps = exp_settings['reps'], load_pulses=True, progress=False)
        
        iq = iq_list[0][0]
        
        #print(iq.shape)
        
        
        Is = iq[:,0]
        Qs = iq[:,1]

        freqs = np.linspace(exp_settings['spec']['freq_start'],exp_settings['spec']['freq_stop'],exp_settings['spec']['freq_points'])

        mag = np.sqrt(Is**2 + Qs**2)
        phase = np.degrees(np.arctan2(Qs,Is))
            
        mags[vind,:] = mag
        phases[vind,:] = phase
        raw_spec['Is'][vind,:] = Is
        raw_spec['Qs'][vind,:] = Qs
            
        if vind == 0:
            t_stop = time.time()
            estimate_time(t_start, t_stop, len(voltages))
    
        transdata = {}
        transdata['xaxis'] = fpts_trans
        transdata['mags'] = trans_mags[0:vind+1,:]
        transdata['phases'] = trans_phases[0:vind+1,:]

        specdata = {}
        specdata['xaxis'] = freqs/1e9 # Convert to GHz
        specdata['mags'] = mags[0:vind+1,:]
        specdata['phases'] = phases[0:vind+1,:]

        singledata = {}
        singledata['xaxis'] = freqs/1e9
        singledata['mag']   = specdata['mags'][vind]
        singledata['phase'] = specdata['phases'][vind]

        trans_labels = ['Freq (GHz)','Voltage (V)']
        spec_labels  = ['Freq (GHz)','Voltage (V)']
        
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
        
        plots.autoscan_plot(transdata, specplotdata, singledata, voltages[0:vind+1], filename, trans_labels, spec_labels, identifier, fig_num = 1)
        userfuncs.SaveFull(saveDir, filename, ['transdata', 'raw_trans', 'specdata', 'raw_spec', 'singledata', 'voltages', 
                                       'filename', 'trans_labels', 'spec_labels'], 
                                       locals(), expsettings=settings, instruments=instruments)
        plt.savefig(os.path.join(saveDir, filename+'.png'), dpi = 150)    
    

    data = {'saveDir': saveDir, 'filename': filename, 'specdata': specdata}

    return data,prog

    
def setfilters_trans(soc, config, exp_settings, exp_globals, board_freq):
    if exp_settings['spec']['filter'] == 'all_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])

    elif exp_settings['spec']['filter'] == 'cav_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bypass')

    elif exp_settings['spec']['filter'] == 'ro_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bypass')
        soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])

    elif exp_settings['spec']['filter'] == 'no_filter':
        soc.rfb_set_gen_filter(config['cav_channel'], fc=board_freq/1000, ftype='bypass')
        soc.rfb_set_ro_filter(config['ro_channel'], fc=board_freq/1000, ftype='bypass')

    else:
        print('Please select one option from:')
        print('\'all_filter\', \'cav_filter\', \'ro_filter\', and\'no_filter\'')
        return