from qick.asm_v2 import AveragerProgramV2

import numpy as np
import os
import matplotlib.pyplot as plt
import userfuncs
from utility.plotting_tools import simplescan_plot
from utility.measurement_helpers import estimate_time
import time

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
        # #self.reset_phase(gen_ch=self.cfg['cav_channel'],t=0)
        
        # #The body sets the pulse sequence, it runs through it a number of times specified by "reps" and takes averages
        # #specified by "soft_averages." Both are required if you wish to acquire_decimated, only "reps" is otherwise.
        # offset = self.us2cycles(self.cfg["adc_trig_offset"],gen_ch=self.cfg["cav_channel"])
        # meas_time = self.us2cycles(self.cfg["meas_time"],gen_ch=self.cfg["cav_channel"])
        # #Sets off the ADC
        # self.trigger(adcs=self.ro_chs,
        #             pins=[0],
        #             adc_trig_offset=offset)
        
        # #Sends measurement pulse
        # self.pulse(ch=self.cfg["cav_channel"],t=meas_time)
        # self.wait_all() #Tells TProc to wait until pulses are complete before sending out the next command
        # self.sync_all(self.us2cycles(self.cfg["relax_delay"])) #Syncs to an offset time after the final pulse is sent

        self.pulse(ch=self.cfg['cav_channel'], name="mypulse", t=self.cfg["meas_time"])
        self.trigger(ros=[self.cfg['ro_channel']], pins=[0], t=self.cfg['adc_trig_offset'])
        self.wait_auto()
        self.delay(self.cfg["relax_delay"])

def get_trans_settings():
    settings = {}
    
    settings['scanname'] = 'initial_power_scan_q4'
    settings['meas_type'] = 'PulsedTrans'
    
    #Sweep parameters
    settings['freq_start']   = 7e9 
    settings['freq_step']    = 100e6
    settings['freq_points']  = 21

    settings['gain_start']  = 500
    settings['gain_step']   = 100
    settings['gain_points'] = 31

    #Card settings
    settings['reps'] = 1
    settings['averages'] = 1e3
    
    return settings

def pulsed_trans(soc,soccfg,instruments,settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings'] 
    m_pulse      = exp_globals['measurement_pulse']
    
    # print(m_pulse)
    



    stamp    = userfuncs.timestamp()
    saveDir  = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel'     : exp_globals['cav_channel']['ID'],
        'ro_channel'     : exp_globals['ro_channel']['ID'],

        'nqz_c'           : 2,
        'cav_phase'       : m_pulse['cav_phase'],
        'meas_window'     : m_pulse['meas_window'],
        'meas_time'       : m_pulse['meas_pos'],
        'meas_gain'       : 1.0, #Placeholder, gets overwritten in sweep
        'cav_atten'       : exp_globals['cav_channel']['Atten_1'] + exp_globals['cav_channel']['Atten_2'],
        'cav_freq'        : 6000, #Placeholder, MHz
        'mixer_freq'      : 6000, #Placeholder, MHz 
        
        'ramp_len'        :5e-3, #Placeholder, us
        'flat_len'        :0.50, #Placeholder, us
        
        'readout_length'  : m_pulse['init_buffer'] + m_pulse['meas_window'] + m_pulse['post_buffer'],
        'adc_trig_offset' : m_pulse['emp_delay'] + m_pulse['meas_pos'] - m_pulse['init_buffer'],


        'relax_delay'     : exp_globals['relax_delay']
        }

    cav_ch = exp_globals['cav_channel']
    ro_ch  = exp_globals['ro_channel']
    # Set attenuator on DAC.
    soc.rfb_set_gen_rf(cav_ch['ID'], cav_ch['Atten_1'], cav_ch['Atten_2'])
    # Set attenuator on ADC.
    soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])


    fpts = exp_settings["freq_start"]+exp_settings["freq_step"]*np.arange(exp_settings["freq_points"])
    gpts = exp_settings["gain_start"]+exp_settings["gain_step"]*np.arange(exp_settings["gain_points"])

    prog = CavitySweep(soccfg, reps=exp_settings['reps'], final_delay = None, final_wait=0, cfg=config)
   

    powerdat = np.zeros((len(gpts), len(fpts)))
    normdat  = np.zeros((len(gpts), len(fpts)))
    phasedat = np.zeros((len(gpts), len(fpts)))
    Is = np.zeros((len(gpts), len(fpts)))
    Qs = np.zeros((len(gpts), len(fpts)))
    
    
    tstart = time.time()
    
    # Here we try to fix the mixer frequency to see if it makes phase data better
    config["mixer_freq"] = (exp_settings["freq_start"] + exp_settings['mixer_detuning'])/1e6


    for g in range(0,len(gpts)):
        print("Current Gain: " + str(gpts[g]) + ", Max Gain: " + str(gpts[-1]))
        
        config["meas_gain"] = gpts[g]
        
        for f in range(0,len(fpts)):
            board_freq = fpts[f]/1e6
            
            config["cav_freq"] = board_freq
            #config["mixer_freq"] = board_freq + exp_settings['mixer_detuning']/1e6
            # the above line is only useful when the span is very wide such that the IF is out of range
            # but it mess up the phase data
            prog = CavitySweep(soccfg, reps=exp_settings['reps'], final_delay = None, final_wait=0, cfg=config)
            #Need to assign Iwindow, Qwindow, Ifull, Qfull, xaxis (which should just be timeus)
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
            
            holder = prog.acquire(soc, reps = exp_settings['reps'], load_pulses=True, progress=False)
            #print(holder)
            iq = holder[0]
            
            I_full = iq[:,0]
            Q_full = iq[:,1]
            
            
            mag = np.sqrt(I_full**2 + Q_full**2)
            phase = np.arctan2(Q_full, I_full)*180/np.pi

            powerdat[g,f] = mag/gpts[g] #bring in the normalization
            phasedat[g,f] = phase
            Is[g,f] = I_full
            Qs[g,f] = Q_full
            
        
        if g == 0:
            tstop = time.time()
            estimate_time(tstart, tstop, len(gpts))

            xaxis = np.linspace(0,len(I_full)-1,len(I_full))

            for x in range(0,len(xaxis)):
                xaxis[x] = prog.cycles2us(xaxis[x],ro_ch=0)

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
        identifier = f"Atten={config['cav_atten']} dB"

        
        simplescan_plot(plot_data, single_data, yaxis, filename, labels, identifier=identifier, fig_num=1, IQdata = False) #normalized to drive level
        plt.savefig(os.path.join(saveDir, filename+'_fullColorPlot.png'), dpi = 150)

        userfuncs.SaveFull(saveDir, filename, ['gpts','fpts', 'powerdat', 'phasedat','xaxis','Is', 'Qs', 'full_data', 'single_data'],
                             locals(), expsettings=settings, instruments={})
    
    # if exp_globals['LO']:
    #     pass
    #     #logen.output = 0
        
    return full_data



