# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 18:04:20 2020

@author: Kollarlab
"""
import pyvisa

class HS9000B():
    def __init__(self, address):
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        self.inst.baud_rate = 115200
        self.inst.write_termination = '\n'
        self.inst.read_termination = '\n'

        self.ch1 = channel(self.inst, 1)
        self.ch2 = channel(self.inst, 2)
        self.ref = reference(self.inst)

    @property
    def settings(self):
        fullsettings = {}
        fullsettings['ch1'] = self.ch1.settings
        fullsettings['ch2'] = self.ch2.settings
        fullsettings['ref'] = self.ref.settings
    
        return fullsettings

    def reset(self):
        self.ch1.reset()
        self.ch2.reset()

    def close(self):
        self.inst.close()

class reference(object):
    def __init__(self, inst):
        self.inst = inst
        self._source = 'INT'
        self._frequency = '100MHz'
    
    @property
    def source(self):
        print('reading ref')
        [mode, freq] = self.inst.query(':REF:STATUS?').split()
        return mode
    @source.setter
    def source(self, mode):
        print('setting mode')
        self.inst.query(':REF:{}:{}'.format(mode, self._frequency))
        self._source = mode
    
    @property
    def frequency(self):
        print('reading ref(freq)')
        [mode, freq] = self.inst.query(':REF:STATUS?').split()
        return freq
    @frequency.setter
    def frequency(self, value):
        print('{},{}'.format(self._source, value))
        self.inst.query(':REF:{}:{}'.format(self._source, value))
        self._frequency  = value
    
    @property
    def settings(self):
        fullsettings = {}
        fullsettings['source'] = self.source
        fullsettings['frequency'] = self.frequency
        return fullsettings

class channel(object):
    def __init__(self, inst, number):
        #super().__setattr__('init', True)
        self.init = True

        self.ID = number
        self.inst = inst

        self.channel = ':CH{}'.format(self.ID)
        self.props = {}
        self.props['freq'] = ':FREQ'
        self.props['phase'] = ':PHASE'
        self.props['output'] = ':PWR:RF'
        self.props['power'] = ':PWR'
        self.props['ext_mod_enable'] = ':MOD'
        self.props['ext_mod'] = ':MOD:MODE'

        #super.__setattr__(self, 'init', False)
        self.init = False

    @property
    def settings(self):
        fullsettings = {}
        for setting in self.props:
            fullsettings[setting] = self.__getattr__(setting)
        return fullsettings
    
    def reset(self):
        self.inst.query(':CH{}*RST'.format(self.ID))
        self.output = 'OFF'

    def __getattr__(self, name):
        initializing = self.__dict__['init']
        full_dict = self.__dict__
        if initializing or name in full_dict.keys():
            super().__getattribute__(self, name)
        else:
            prop_dict = self.__dict__['props']
            channel = self.__dict__['channel']
            if name in prop_dict.keys():
                val =  self.inst.query(channel+prop_dict[name]+'?')
                try:
                    return eval(val)
                except:
                    return val
            else:
                print('invalid command')
    
    def __setattr__(self, name, value):
        try:
            initializing = self.__dict__['init']
        except:
            initializing = True
        full_dict = self.__dict__
        if initializing or name in full_dict.keys():
            super().__setattr__(name, value)
        else:
            prop_dict = self.__dict__['props']
            channel = self.__dict__['channel']
            if name in prop_dict.keys():
                self.inst.query(channel+prop_dict[name]+':{}'.format(value))
            else:
                print('setting an invalid command')