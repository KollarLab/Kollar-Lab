# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 16:21:28 2024

@author: KollarLab
"""

#Needs to have Signal Core software downloaded to function!

import sc5511a_JY
from sc5511a_JY import *

SC = SignalCore_SC5511A(name='SC1',serial_number='10004712')

#Think it's worth looking through the holzworth code. That has two subclasses, one for the RF channels and one for reference clock. 
#To talk with the channel, we go "holz.ch1.freq" and to talk with the reference we go "holz.ref.mode"
# SC.settings (would give settings for whole)

# SC.set_open works but then setting RF1 output fails with memory access violation.
#SC.set_open(True)
SC.standby
SC.standby = True
SC.standby
SC.rf2_standby
SC.output
SC.output = True
SC.rf_mode = 1 # 0 Single 1 List/Sweep

SC.start_freq = 321e6
SC.stop_freq = 123
SC.step_freq
SC.list_cycle_count = 17

SC.freq = 6.17e9 # SC.ch1.freq = 5e9
SC.freq
SC.level  = 4 # SC.ch1.power = 13 #SET BACK TO THREE
SC.rf2_freq = 1234

SC.phase = 60
print(SC.phase)



ext_ref_freq = 0 #(selects input as 10 MHz or 100 MHz) We'll do 10 MHz since that's our RB clock
ext_direct_clk = 0 # Chooses whether to directly clock synthesizer with external ref, think this only works if we're using a 100 MHz ref
select_high = 1 #(selects 10 MHz or 100 MHz) but no clue what it does, maybe internal clock freq? Example sets it to 100 MHz, let's keep it that way
lock_external = 1 #tells it to lock internal clock to external ref.

SC.set_clock_reference(ext_ref_freq, ext_direct_clk, select_high, lock_external)
print(SC.get_rf_parameters())
SC._close()

#Example ref clock functions
# SC.ref.mode = 'Ext'
# SC.ref.freq = 10e6 or 100e6




#SC._close() #This doesn't work
#SC._dll.sc5511a_close_device(SC._handle) #This does work




#Example property for SC.ch1
# @property
# def freq(self):
#     val = float(self.get_frequency())
#     return val
# @freq.setter(x)
# def freq(self, val):
#     self.set_frequency(val)

#The 