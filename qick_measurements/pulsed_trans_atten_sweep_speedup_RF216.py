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
        self.declare_readout(ch=ro_ch, length=cfg['readout_length'])

        self.add_readoutconfig(
            ch=ro_ch, name="myro",
            freq=cfg['cav_freq'],
            gen_ch=gen_ch,
            outsel='product'
        )

        self.add_pulse(
            ch=gen_ch, name="mypulse", ro_ch=ro_ch,
            style="const",
            freq=cfg['cav_freq'],
            length=cfg["meas_window"],
            phase=cfg['cav_phase'],
            gain=cfg['meas_gain'],
        )

        self.send_readoutconfig(ch=ro_ch, name="myro", t=0)

    def _body(self, cfg):
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

def pulsed_trans(soc, soccfg, instruments, settings):
    exp_globals = settings['exp_globals']
    exp_settings = settings['exp_settings']
    m_pulse = exp_globals['measurement_pulse']

    stamp = userfuncs.timestamp()
    saveDir = userfuncs.saveDir(settings)
    filename = exp_settings['scanname'] + '_' + stamp

    config = {
        'cav_channel': exp_globals['cav_channel']['ID'],
        'ro_channel': exp_globals['ro_channel']['ID'],

        'nqz_c': 2,
        'cav_phase': m_pulse['cav_phase'],
        'meas_window': m_pulse['meas_window'],
        'meas_time': m_pulse['meas_pos'],
        'meas_gain': exp_settings['meas_gain'],
        'meas_atten': 0.0,  # overwritten in sweep

        'cav_freq': 6000,   # placeholder, MHz
        'mixer_freq': 6000, # placeholder, MHz

        'readout_length': m_pulse['meas_window'],
        'adc_trig_offset': m_pulse['emp_delay'] + m_pulse['meas_pos'],
        'relax_delay': exp_globals['relax_delay']
    }

    cav_ch = exp_globals['cav_channel']
    ro_ch = exp_globals['ro_channel']

    # ADC attenuation
    soc.rfb_set_ro_rf(ro_ch['ID'], ro_ch['Atten'])

    # Sweep axes
    fpts = exp_settings["freq_start"] + exp_settings["freq_step"] * np.arange(exp_settings["freq_points"])
    apts = exp_settings["atten_start"] + exp_settings["atten_step"] * np.arange(exp_settings["atten_points"])

    powerdat = np.zeros((len(apts), len(fpts)))
    phasedat = np.zeros((len(apts), len(fpts)))
    Is = np.zeros((len(apts), len(fpts)))
    Qs = np.zeros((len(apts), len(fpts)))


    tstart = time.time()
    # Here we try to fix the mixer frequency to see if it makes phase data better
    f_center_MHz = (np.mean(fpts)) / 1e6
    config["mixer_freq"] = f_center_MHz + exp_settings['mixer_detuning']/1e6

    for a in range(len(apts)):
        print(f"Current Attenuation: {apts[a]}, Max Attenuation: {apts[-1]}")
        config["meas_atten"] = apts[a]

        # Set generator attenuation in sweep
        soc.rfb_set_gen_rf(cav_ch['ID'], apts[a] / 2, apts[a] / 2)

        for f in range(len(fpts)):
            board_freq = fpts[f] / 1e6  # MHz

            config["cav_freq"] = board_freq
            #config["mixer_freq"] = board_freq + exp_settings['mixer_detuning'] / 1e6

            prog = CavitySweep(soccfg, reps=exp_settings['reps'], final_delay=None, final_wait=0, cfg=config)

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

            holder = prog.acquire(soc, rounds = exp_settings['rounds'], load_pulses=True, progress=False)
            
# =============================================================================
# #-----------------------------------------------debugging---------------------------------------------------------            
#             print("type(holder):", type(holder))
#             try:
#                 print("len(holder):", len(holder))
#             except Exception as e:
#                 print("holder has no len:", e)
#             
#             # common pattern: holder[0] is IQ
#             iq = holder[0]
#             print("type(iq):", type(iq))
#             print("iq shape:", getattr(iq, "shape", None))
#             print("iq dtype:", getattr(iq, "dtype", None))
#             
#             # print a few rows safely
#             print("iq preview:", iq[:min(5, len(iq))])
# #-----------------------------------------------debugging--------------------------------------------------------- 
# 
# =============================================================================
            iq = holder[0]          # typically shape (1, 2)
            I = iq[0, 0]
            Q = iq[0, 1]
            
            powerdat[a, f] = np.sqrt(I**2 + Q**2)
            phasedat[a, f] = np.arctan2(Q, I) * 180/np.pi
            
            Is[a, f] = I
            Qs[a, f] = Q

        # After first attenuation point, estimate total time
        if a == 0:
            tstop = time.time()
            estimate_time(tstart, tstop, len(apts))

        # Pack data for plots/saving
        full_data = {
            'xaxis': fpts / 1e9,
            'mags': powerdat[0:a+1],
            'phases': phasedat[0:a+1],
        }

        plot_data = {
            'xaxis': fpts / 1e9,
            'mags': powerdat[0:a+1],
            'phases': phasedat[0:a+1],
        }

        single_data = {
            'xaxis': fpts / 1e9,
            'mag': powerdat[a],
            'phase': phasedat[a],
        }

        yaxis = apts[0:a+1]
        labels = ['Freq (GHz)', 'Attenuation (dB)']
        identifier = f"Gain={config['meas_gain']} a.u."

        simplescan_plot(plot_data, single_data, yaxis, filename, labels,
                        identifier=identifier, fig_num=1, IQdata=False)
        plt.savefig(os.path.join(saveDir, filename + '_fullColorPlot.png'), dpi=150)

        userfuncs.SaveFull(
            saveDir, filename,
            ['apts', 'fpts', 'powerdat', 'phasedat', 'Is', 'Qs','full_data', 'single_data'],
            locals(), expsettings=settings, instruments={}
        )

    return full_data
