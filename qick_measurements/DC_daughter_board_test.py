# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 13:51:24 2026

@author: KollarLab
"""

from qick.asm_v2 import AveragerProgramV2, QickSweep1D


import numpy as np
import time
import userfuncs
from utility.measurement_helpers import estimate_time
import logging
from utility.plotting_tools import simplescan_plot
import matplotlib.pyplot as plt

class DC_Board_Output(AveragerProgramV2):
    def _initialize(self, cfg):
        
        gen_ch = cfg['gen_ch']
        self.declare_gen(ch=gen_ch, nqz=1)
        
        self.add_pulse(ch=gen_ch, name="dc_pulse",
                       style="const", 
                       freq=0, # MHz
                       length=10, # us
                       phase=0,
                       gain=cfg['gain'],
                       mode = "periodic"#"oneshot"
                      )
        
        self.add_pulse(ch=gen_ch, name="ac_pulse",
                       style='const',
                       freq=20, # MHz
                       length=20, # us
                       phase=0,
                       phrst=1,
                       gain=cfg['gain'],
                       mode = 'periodic'
                       )
        
    def _body(self, cfg):
        #self.pulse(ch=cfg["gen_ch"],name='dc_pulse',t=0.0)
        self.pulse(ch=cfg["gen_ch"],name='ac_pulse',t=0.0)
        self.delay_auto(t=4.0, gens=False, ros=True)
        
#soc.rfb_set_dac_dc(0)
        
config = {
    'gen_ch' : 3,
    'ro_ch'  : 0,
    'nqz'    : 1,
    'gain'   : 0.006 
    }

prog = DC_Board_Output(soccfg, reps=1, final_delay=1, cfg=config)
prog.run_rounds(soc, rounds = 1, load_envelopes=False, progress=False)

#soc.reset_gens()


# =============================================================================
# class ScopeTest(AveragerProgramV2):
#     def _initialize(self, cfg):
#         ro_ch = cfg['ro_ch']
#         gen_ch = cfg['gen_ch']
#         self.declare_gen(ch=gen_ch, nqz=cfg['nqz'], ro_ch=ro_ch)
#         self.declare_readout(ch=ro_ch, length=cfg['ro_len'])
#         self.add_readoutconfig(ch=ro_ch, name="myro",
#                                freq=cfg['freq'],
#                                gen_ch=gen_ch,
#                                outsel='product')
#         self.add_cosine(ch=gen_ch, name="ramp", length=cfg['ramp_len'], even_length=True)
#         self.add_pulse(ch=gen_ch, name="mypulse", ro_ch=ro_ch,
#                        #style="const",
#                        style="flat_top", 
#                        envelope="ramp", 
#                        freq=cfg['freq'],
#                        length=cfg['flat_len'],
#                        phase=cfg['phase'],
#                        gain=cfg['gain'],
#                        #mode = "periodic"
#                       )
#         self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)
#     def _body(self, cfg):
#         self.delay_auto()
#         self.pulse(ch=cfg['gen_ch'], name="mypulse", t=0.0)
#         #self.trigger(ros=[cfg['ro_ch']], pins=[0], t=cfg['trig_time'], mr=True)
# 
# 
# config = {'gen_ch': 0,
#           'ro_ch': 0,
#           'freq': 0,
#           'nqz': 1,
#           'trig_time': 0.0,
#           'ro_len': 3.0,
#           'flat_len': 2.0,
#           'ramp_len': 4.0,
#           'phase': 0,
#           'gain': 0.4
#          }
# 
# prog = ScopeTest(soccfg, reps=1, final_delay=0.5, cfg=config)
# =============================================================================
    
    
    
    
    
    
        
        