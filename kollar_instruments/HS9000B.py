# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 18:04:20 2020

@author: Kollarlab
"""
from userfuncs import freeze
import pyvisa

@freeze
class HS9000B():
    def __init__(self, address):
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        self.inst.baud_rate = 115200

        self.ch1 = channel(self.inst, 1)
        self.ch2 = channel(self.inst, 2)
    
@freeze
class channel():
    def __init__(self, inst, number):
        self.ID = number
        self.inst = inst

    @property 
    def freq(self):
        return self.inst.query(':CH{}:FREQ?'.format(self.ID))
    @freq.setter
    def freq(self, value):
        self.inst.write(':CH{}:FREQ:{}'.format(self.ID, value))
    
    @property
    def output(self):
        return self.inst.query(':CH{}:PWR:RF?'.format(self.ID))
    @output.setter
    def output(self, state):
        self.inst.write(':CH{}:PWR:RF:{}'.format(self.ID, value))
    
    @property
    def power(self):
        return self.inst.query(':CH{}:PWR?'.format(self.ID))
    @power.setter
    def power(self, value):
        self.inst.write(':CH{}:PWR:{}'.format(self.ID, value))
    
    @property
    def phase(self):
        return self.inst.query(':CH{}:PHASE?'.format(self.ID))
    @phase.setter
    def phase(self, value):
        self.inst.write(':CH{}:PHASe:{}'.format(self.ID, value))

    @property
    def ext_mod_enable(self):
        self.inst.query(':CH{}:MOD?'.format(self.ID))
    
    @property
    def ext_mod(self):
        return self.inst.query(':CH{}:MOD:MODE?'.format(self.ID))
    @ext_mod.setter
    def ext_mod(self, mode):
        self.inst.write(':CH{}:MOD:MODE:{}'.format(self.ID, mode))