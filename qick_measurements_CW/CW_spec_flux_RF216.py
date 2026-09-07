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

class CW_trans(AveragerProgramV2):
    def _initialize(self, cfg):

        ro_ch = cfg['ro_channel']['ID']
        gen_ch = cfg['cav_channel']['ID']
        self.declare_gen(ch=gen_ch, nqz=cfg['nqz_c'], mixer_freq=cfg['mixer_freq'], ro_ch=ro_ch)
        #for ro_ch in cfg["ro_channels"]:
        self.declare_readout(ch=ro_ch, length=cfg['readout_length'])
        self.add_readoutconfig(ch=ro_ch, name="myro",
                           freq=cfg['pulse_freq'],
                           gen_ch=gen_ch,
                           outsel='product')
        self.add_pulse(ch=gen_ch, name="CW_pulse", ro_ch=ro_ch,
                       style="const",
                       freq=cfg['pulse_freq'],
                       length= cfg["pulse_length"],
                       phase=cfg['cav_phase'],
                       gain=cfg['pulse_gain'],
                       mode='periodic'
                      )
        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)
        
    def _body(self, cfg):

        self.pulse(ch=self.cfg['cav_channel']['ID'], name="CW_pulse", t=0.0)
        self.trigger(ros=[self.cfg['ro_channel']['ID']], pins=[0], t=self.cfg['adc_trig_offset'])
        self.delay_auto(t=4, gens=False, ros=True) #Set the reference time to the end of the last pulse/readout, plus the specified value.
                                                     #You can select whether this accounts for pulses, readout windows, or both.
        #self.delay(self.cfg["relax_delay"]) # no delay needed here for CW measurements


class CW_spec(AveragerProgramV2):
    def _initialize(self, cfg):

        ro_ch = cfg['ro_channel']['ID']
        gen_ch = cfg['cav_channel']['ID']
        qub_ch = cfg["qub_channel"]["ID"] # These are just channel ID, not the whole dict.

        #self.add_loop("freq_sweep", cfg["freq_points"]) # if you want the QickSweep1D thing, turn this back

        # set the nyquist zone
        self.declare_gen(ch=cfg["cav_channel"]["ID"], nqz=cfg["nqz_c"], mixer_freq=cfg['cav_mixer_freq'],ro_ch = ro_ch)
        self.declare_gen(ch=cfg["qub_channel"]["ID"], nqz=cfg["nqz_q"], mixer_freq=cfg['qub_mixer_freq'])

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
                       gain=cfg['cav_gain'],
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
        self.trigger(ros=[cfg['ro_channel']['ID']],
                    pins=[0],
                    t=self.cfg['adc_trig_offset'])
        
        self.pulse(ch=cfg["qub_channel"]['ID'],name='qub_pulse',t=0.0)
        self.pulse(ch=cfg["cav_channel"]['ID'],name='cav_pulse',t=0.0)
        #self.wait_auto()
        # DO NOT rely on wait_auto() here because wait_auto() often cannot advance the timeline the way you think during a sweep loop
        # since for periodic pulses, the generator is considered “running continuously” (it doesn’t have a finite end time)
        self.delay_auto(t=4, gens=False, ros=True) #Set the reference time to the end of the last pulse/readout, plus the specified value.
                                                     #You can select whether this accounts for pulses, readout windows, or both.
        #self.delay(self.cfg["relax_delay"])

def get_cw_spec_flux_settings():
    fullsettings = {}
    settings = {}
    autoscan_settings = {}

    settings['scanname'] = 'continuous_power_scan'
    settings['meas_type'] = 'CW_Spec_Flux'
    
    settings['cav_gain'] = 0.6
    settings['meas_window'] = 900

    settings['qub_gain'] = 0
    
    #Freq sweep parameters
    settings['freq_start']   = 4e9  
    settings['freq_stop']    = 4.5e9
    settings['freq_points']  = 6

    #Voltage sweep parameters
    settings['start_voltage']   = 0
    settings['stop_voltage']    = 1
    settings['voltage_points']  = 2

    #Card settings
    settings['reps'] = 1
    settings['soft_avgs'] = 1#5e3

    autoscan_settings['reps'] = 1
    autoscan_settings['soft_avgs'] = 1 

    autoscan_settings['freq_start'] = 5e6
    autoscan_settings['freq_stop']  = 5.5e6
    autoscan_settings['freq_points'] = 6

    fullsettings['spec'] = settings
    fullsettings['autoscan'] = autoscan_settings

    return fullsettings

def cw_spec_flux(soc,soccfg,instruments,settings):
    
    # suppresses sum buffer overflow warning
    logging.getLogger("qick").setLevel(logging.ERROR)

    SRS = instruments['DCsupply']

    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']

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

    m_pulse      = exp_globals['measurement_pulse']
    q_pulse      = exp_globals['qubit_pulse']
    

    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = spec_set['scanname'] + '_' + stamp
    config_trans = {
        'cav_channel'     : exp_globals['cav_channel'],
        'ro_channel'     : exp_globals['ro_channel'],
        'relax_delay'     : exp_globals['relax_delay'],
        'cav_atten'       : exp_globals['cav_channel']['Atten_1'] + exp_globals['cav_channel']['Atten_2'],
        'nqz_c'           : 2,
        'cav_phase'       : m_pulse['cav_phase'],
        'pulse_length'    : 1,
        'meas_window'     : autoscan_set['readout_length'],
        'pulse_gain'      : autoscan_set['cav_gain'],
        'pulse_freq'      : 0, #Placeholder
        'adc_trig_offset' : m_pulse['emp_delay'],
        'reps'            : autoscan_set['reps'],
        'soft_avgs'       : autoscan_set['soft_avgs'],
        'readout_length'  : autoscan_set['readout_length']
    }

    config_qubit = {
        'cav_channel'     : exp_globals['cav_channel'],
        'qub_channel'     : qub_ch,
        'ro_channel'     : exp_globals['ro_channel'],

        'nqz_c'           : 2,
        'cav_phase'       : m_pulse['cav_phase'],
        'cav_pulse_len'   : 1,

        'meas_time'       : m_pulse['meas_pos'],
        'cav_gain'        : autoscan_set['cav_gain'],
        'cav_freq'        : 0, #exp_settings['cav_freq']/1e6,
        
        'nqz_q'           : 2,
        'qub_phase'       : q_pulse['qub_phase'],
        'qub_gain'        : spec_set['qub_gain'], 

        #'qub_sigma'       : q_pulse['sigma'],
        #'qub_delay'       : q_pulse['delay'],
        #'num_sigma'       : q_pulse['num_sigma'],
        'qub_len'         : 1,

        'readout_length'  : autoscan_set['readout_length'],
        'qub_pulse_len'   : spec_set['qub_pulse_len'],
        'adc_trig_offset' : m_pulse['emp_delay'],  #+ m_pulse['meas_pos'],


        'relax_delay'     : exp_globals['relax_delay'],
        'reps'            : spec_set['reps'],
        'soft_avgs'       : spec_set['soft_avgs']
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

    f0_start = spec_set['freq_start']
    f0_stop = spec_set['freq_stop']
    expts = spec_set['freq_points']
    
    fpts = np.linspace(f0_start,f0_stop,expts) #np.arange(0,expts)*f0_step+f0_start


    mags = np.zeros((voltage_points, spec_set['freq_points']))
    phases = np.zeros((voltage_points, spec_set['freq_points']))

    raw_spec = {}
    raw_spec['Is'] = np.zeros((voltage_points, spec_set['freq_points']))
    raw_spec['Qs'] = np.zeros((voltage_points, spec_set['freq_points']))
    
    t_start = time.time()
    
    identifier = 'Cav Gain : ' + str(autoscan_set['cav_gain'])  + ' au'

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

            config_trans['pulse_freq'] = fpts_trans[tfind]/1e6
            config_trans['mixer_freq'] = config_trans['pulse_freq'] + autoscan_set['cav_mixer_detuning']/1e6
            setfilters_trans(soc, config_trans, exp_settings, exp_globals, config_trans['pulse_freq'])
            prog = CW_trans(soccfg,reps = exp_settings['autoscan']['reps'], final_delay = None, final_wait = 0, cfg =config_trans)

            holder = prog.acquire(soc,progress=False)
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

        config_qubit['cav_mixer_freq'] = config_qubit['cav_freq'] + autoscan_set['cav_mixer_detuning']/1e6

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

        for find in range(0,len(fpts)):

            config_qubit['qub_freq'] = fpts[find]/1e6 # convert to MHz
            config_qubit['qub_mixer_freq'] = config_qubit['qub_freq'] + exp_settings['spec']['qub_mixer_detuning']/1e6
            #print(config_qubit)

            setfilters_spec(soc, config_qubit, exp_settings, exp_globals, fpts, find, qub_ch)

            
            prog = CW_spec(soccfg,reps = exp_settings['spec']['reps'], final_delay = None, final_wait = 0, cfg = config_qubit)
            holder = prog.acquire(soc,progress=False)
            iq = holder[0][0]
            I = iq[0] 
            Q = iq[1]
            soc.reset_gens()
            
            mag = np.sqrt(I**2 + Q**2)
            phase = np.arctan2(Q, I)*180/np.pi
            
            mags[vind,find] = mag
            phases[vind,find] = phase
            raw_spec['Is'][vind,find] = I
            raw_spec['Qs'][vind,find] = Q
            
        if vind == 0:
            t_stop = time.time()
            estimate_time(t_start, t_stop, len(voltages))
    
        transdata = {}
        transdata['xaxis'] = fpts_trans
        transdata['mags'] = trans_mags[0:vind+1,:]
        transdata['phases'] = trans_phases[0:vind+1,:]

        specdata = {}
        specdata['xaxis'] = fpts/1e9 # Convert to GHz
        specdata['mags'] = mags[0:vind+1,:]
        specdata['phases'] = phases[0:vind+1,:]

        singledata = {}
        singledata['xaxis'] = fpts/1e9
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

def setfilters_spec(soc, config, exp_settings, exp_globals, fpts, f, qub_ch):
    if exp_settings['spec']['filter'] == 'all_filter':
        soc.rfb_set_gen_filter(config['cav_channel']['ID'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel']['ID'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])

        soc.rfb_set_gen_filter(config['qub_channel']['ID'], fc=fpts[f]/1e9, ftype='bandpass', bw=qub_ch['BW'])

    elif exp_settings['spec']['filter'] == 'no_ro_filter':
        soc.rfb_set_gen_filter(config['cav_channel']['ID'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel']['ID'], fc=config['cav_freq']/1000, ftype='bypass')

        soc.rfb_set_gen_filter(config['qub_channel']['ID'], fc=fpts[f]/1e9, ftype='bypass')

    elif exp_settings['spec']['filter'] == 'no_qubit_filter':
        soc.rfb_set_gen_filter(config['cav_channel']['ID'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel']['ID'], fc=config['cav_freq']/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])

        soc.rfb_set_gen_filter(config['qub_channel']['ID'], fc=fpts[f]/1e9, ftype='bypass')

    elif exp_settings['spec']['filter'] == 'no_filter':
        soc.rfb_set_gen_filter(config['cav_channel']['ID'], fc=config['cav_freq']/1000, ftype='bypass')
        soc.rfb_set_ro_filter(config['ro_channel']['ID'], fc=config['cav_freq']/1000, ftype='bypass')

        soc.rfb_set_gen_filter(config['qub_channel']['ID'], fc=fpts[f]/1e9, ftype='bypass')

    else:
        print('Please select one option from:')
        print('\'all_filter\', \'no_qubit_filter\', and \'no_filter\'')
        return
    
def setfilters_trans(soc, config, exp_settings, exp_globals, board_freq):
    if exp_settings['spec']['filter'] == 'all_filter':
        soc.rfb_set_gen_filter(config['cav_channel']['ID'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel']['ID'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])

    elif exp_settings['spec']['filter'] == 'cav_filter':
        soc.rfb_set_gen_filter(config['cav_channel']['ID'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['cav_channel']['BW'])
        soc.rfb_set_ro_filter(config['ro_channel']['ID'], fc=board_freq/1000, ftype='bypass')

    elif exp_settings['spec']['filter'] == 'ro_filter':
        soc.rfb_set_gen_filter(config['cav_channel']['ID'], fc=board_freq/1000, ftype='bypass')
        soc.rfb_set_ro_filter(config['ro_channel']['ID'], fc=board_freq/1000, ftype='bandpass', bw=exp_globals['ro_channel']['BW'])

    elif exp_settings['spec']['filter'] == 'no_filter':
        soc.rfb_set_gen_filter(config['cav_channel']['ID'], fc=board_freq/1000, ftype='bypass')
        soc.rfb_set_ro_filter(config['ro_channel']['ID'], fc=board_freq/1000, ftype='bypass')

    else:
        print('Please select one option from:')
        print('\'all_filter\', \'cav_filter\', \'ro_filter\', and\'no_filter\'')
        return