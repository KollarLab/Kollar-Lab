# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 15:08:04 2025

@author: BF2-meas
"""

import numpy as np
import matplotlib.pyplot as plt
from qick.asm_v2 import AveragerProgramV2


GEN_CH=5 # DAC 
RO_CH=0  # ADC

# print("read bias:", soc.rfb_get_bias(0))
# print("set bias:", soc.rfb_set_bias(0, 0.0))
# print("read bias:", soc.rfb_get_bias(0))

class LoopbackProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        ro_ch = cfg['ro_ch']
        gen_ch = cfg['gen_ch']
        self.declare_gen(ch=gen_ch, nqz=cfg['nqz'], mixer_freq=cfg['mixer_freq'])#, ro_ch=ro_ch)
        self.declare_readout(ch=ro_ch, length=cfg['ro_len'])
        self.add_readoutconfig(ch=ro_ch, name="myro",
                               freq=cfg['freq'],
                               gen_ch=gen_ch,
                               outsel='product')
        self.add_cosine(ch=gen_ch, name="ramp", length=cfg['ramp_len'], even_length=True)
        self.add_pulse(ch=gen_ch, name="mypulse", ro_ch=ro_ch,
                       # style="const",
                       style="flat_top", 
                       envelope="ramp", 
                       freq=cfg['freq'],
                       length=cfg['flat_len'],
                       phase=cfg['phase'],
                       gain=cfg['gain'],
                      )
        self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)
    def _body(self, cfg):
        self.delay_auto()
        self.pulse(ch=cfg['gen_ch'], name="mypulse", t=0)
        self.trigger(ros=[cfg['ro_ch']], pins=[0], t=cfg['trig_time'], mr=True)
        self.wait_auto()
        # self.delay(self.cfg["relax_delay"])
        #self.sync_all(self.us2cycles(self.cfg["relax_delay"]))

        

config = {'gen_ch': GEN_CH,
      'ro_ch': RO_CH,
      'mixer_freq': 4000,
      'freq': 4200,
      'nqz': 1,
      'trig_time': 0.0,
      'ro_len': 3.0,
      'flat_len': 2.0,
      'ramp_len': 1.0,
      'phase': 90,
      'gain': 1.0
     }

# config = {'gen_ch': GEN_CH,
#           'ro_ch': RO_CH,
#           'mixer_freq': 4200,
#           'freq': 4500,
#           'nqz': 1,
#           'trig_time': 0.0,
#           'ro_len': 3.0,
#           'flat_len': 1.0,
#           'ramp_len': 1.0,
#           'phase': 0,
#           'gain': 1.0,
#           'relax_delay' : 0.5
#          }

prog = LoopbackProgram(soccfg, reps=1, final_delay = None, final_wait=0, cfg=config)

freq = config['freq']
soc.rfb_set_gen_filter(config['gen_ch'], fc=freq/1000, ftype='bandpass', bw=1.0)
soc.rfb_set_ro_filter(config['ro_ch'], fc=freq/1000, ftype='bandpass', bw=1.0)


# Set attenuator on DAC.
soc.rfb_set_gen_rf(config['gen_ch'], 0, 10)
# Set attenuator on ADC.
soc.rfb_set_ro_rf(config['ro_ch'], 25)

iq_list = prog.acquire_decimated(soc, rounds=10,progress=False)


t = prog.get_time_axis(ro_index=0)
iq = iq_list[0]

plt.figure(1)
plt.clf()
plt.plot(t, iq[:,0], label="I value")
plt.plot(t, iq[:,1], label="Q value")
plt.legend()
plt.ylabel("amplitude [ADU]")
plt.xlabel("time [us]");



