# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 18:04:20 2020

@author: Kollarlab
"""
import pyvisa

class HS9000B():
    '''
    Holzworth class, currently hard coded to 2 channels.
    Attributes:
        settings: property that returns the settings of the submodules (channels and reference)
        ch1, ch2: channel objects with frequency, phase, power, output control
        ref: reference objet that controls the reference clock source (INT/EXT) and frequency (10/100MHz)
    '''
    def __init__(self, address):
        '''
        Constructor for Holzworth
        Arguments:
            address (string): pyvisa resource string (e.g. 'ASRL3::INSTR')

        '''
        rm = pyvisa.ResourceManager()
        self.inst = rm.open_resource(address)
        #We have to explictly set the baud rate and termination because 
        #they don't match the pyvisa defaults
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
    '''
    Channel object for Holzworth. Very similar to SCPIinst object but with tweaked syntax
    because Holzworth doesn't follow the SCPI standard...
    Behavior is pretty much the same: a dictionary of commands is defined and the keys define 
    properties of the object that are accessed (e.g,: channel.freq queries the instrument and
    channel.freq = '5GHz' sets the instrument frequency to 5 GHz)
    All interactions with the instrument happen as queries because that's just how it is
    not sure why.
    Attributes:
        ID: channel number
        inst: pyvisa handle
        props: dictionary of commands that can be set/ read from the instrument
    '''
    def __init__(self, inst, number):

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
        '''
        Overwriting the get attr method so that we can have the instrument
        behave like a normal python object without having to write properties
        or getters/setters for every single setting
        Avoids infinite loops by invoking the standard getattr for the init 
        phase
        '''
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
        '''
        Overwriting the set attribute method
        It will behave 'normally' for the defined properties (ID etc) but
        execute a command on the instrument for the standard operations (freq, phase
        etc.)
        The code at the beginning avoids infinite recursion for the very first setting (init)
        '''
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