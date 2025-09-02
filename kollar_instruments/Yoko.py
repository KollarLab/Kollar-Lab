# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 08:15:08 2024

@author: gchis
"""
from kollar_instruments.SCPIinst import SCPIinst
#from SCPIinst import SCPIinst
from bidict import bidict
import numpy as np
import time

class Yoko(SCPIinst):
    
    errcmds           = {}
    errcmds['error']  = ':SYSTem:ERRor?'
    
    commandlist = {}
    commandlist['core']   = {}

    
    core = {}

    core['mode']    = ':SOURce:FUNCtion' #'CURRent or VOLTage (case insensitive)
    core['output']  = ':OUTPut' #0 or 1, ON or OFF
    
    voltage = {}
    
    voltage['level']       = ':SOURce:LEVEL' #V
    voltage['range']       = ':SOURce:RANGe' #V
    
    current = {}
    
    current['level']       = ':SOURce:LEVEL' #mA
    current['range']       = ':SOURce:RANGe' #A
    
    
    commandlist['core']    = core
    core['voltage']        = voltage
    core['current']        = current
    
    commandlist['voltage'] = voltage
    commandlist['current'] = current
    
    settinglist = {}
    settinglist['mode']  = core['mode']
    settinglist['level'] = voltage['level']
    settinglist['range'] = voltage['range']
    
    '''
    Reminder: Since level and range are shared across voltage and current, we have to decide whether we
    want the settings to be querried from the machine directly (which would involved switching modes and
    resetting the level) or taken from the global variables we are using (crange, vrange).
    '''
    
    
    
    
    '''Here, vrange, and crange are the ranges set at init. They get copied to the class variables to be 
    used elsewhere for safety checks. Here, input voltage and current ranges in A and V.
    Allowed crange inputs: .001, .01, .1, .2
    Allowed vrange inputs: .01, .1, .2, 1, 10, 30
    Range limits are remembered when switching modes.'''
    
    def __init__(self, address, reset=True, mode='CURR', crange = .1, vrange = 10, Output = 1): 
        
        self.init = True
        super().__init__(address, self.commandlist, self.errcmds, reset = reset, baud_rate = 9600)
        
        global curr_range
        global volt_range
        
        
        self.inst.read_termination = '\n'
        self.inst.write_termination = '\r'
        
        self.mode = 'CURR'
        curr_range = crange

        self.current.range = curr_range
        #print(float(self.inst.query('SOUR:RANG?')))

        if curr_range == float(self.inst.query('SOUR:RANG?')):    
            print('Current Range is set to',float(self.inst.query('SOUR:RANG?'))*1000,'mA')
        else:
            print("Curr Error. Current range set to",float(self.inst.query('SOUR:RANG?'))*1000,'mA')
        
        self.mode = 'VOLT'
        volt_range = vrange
        self.voltage.range = volt_range
        if volt_range == float(self.inst.query('SOUR:RANG?')):
            print('Voltage Range is set to',float(self.inst.query('SOUR:RANG?')),'V')
        else:
            print("Volt Error. Voltage range is set to",float(self.inst.query('SOUR:RANG?')),'V')
        
        
        self.mode = mode

        
        
    def set_voltage(self,value):
        
        #print(volt_range)
        read_mode = self.inst.query('SOUR:FUNC?')
        
        if self.mode == "VOLT" and read_mode == "VOLT":
            if value <= volt_range:
                self.voltage.level = value
                return True
            else:
                print('Input,',value,'V, is above limit,',volt_range,'V.')
                return False
        elif self.mode != read_mode:
            print("Machine function,",read_mode,", is not the same as the function in the code,",self.mode,". There's probably an issue with the code.")
            return False
        else:
            print("Machine is currently set to",self.mode,"mode!")
            return False
        
    def set_current(self,value): #input desired value in A
        
        #print(curr_range)
        read_mode = self.inst.query('SOUR:FUNC?')
        
        if self.mode == "CURR" and read_mode == "CURR": #ensuring that the mode is correctly set
            if value <= curr_range:
                self.current.level = value
                return True
            else:
                print('Input,',value,'A, is above limit,',curr_range,'A.')
                return False

        elif self.mode != read_mode:
            print("Machine function,",read_mode,", is not the same as the function in the code,",self.mode,". There's probably an issue with the code.")
            return False
        
        else:
            print("Machine is currently set to",self.mode,"mode!")
            return False

    def change_function(self,value): #changing mode always resets the level to 0
        read_mode = self.inst.query('SOUR:FUNC?')
        if self.mode != value and read_mode != value:
            if value.casefold() == 'curr':
                self.voltage_ramp(0)
            else:
                self.current_ramp(0)
            
            self.mode = value
        elif self.mode != read_mode:
            print("Machine function,",read_mode,", is not the same as the function in the code,",self.mode,". There's probably an issue with the code.")
        
        else:
            print("Machine is currently set to",self.mode,"mode!")
            
            
            
    def settings(self):
        fullsettings = {}
        for setting in list(self.settinglist.values()):
            print(setting)
            fullsettings[setting] = self.inst.query(setting+'?')
            #print(fullsettings)
        return fullsettings
        
    def voltage_ramp(self, newV, step_size = 0.00005, step_time = 0.001): #change conditional order
        
        if newV > volt_range:
            print('Final voltage,',newV,'is above limit,',volt_range)
            return
        deltaV = newV - self.voltage.level
        numSteps = max(2, int(np.abs(np.ceil(deltaV/step_size))))
        vsteps = np.linspace(self.voltage.level, newV, numSteps)
        for vstep in vsteps:
            if self.set_voltage(np.round(vstep,6)) == False:
                break
            time.sleep(step_time)
        return
    
    def current_ramp(self, newC, step_size = 0.00005, step_time = 0.001):
        
        if newC > curr_range:
            print('Final voltage,',newC,'is above limit,',curr_range)
            return
        deltaC = newC - self.current.level
        numSteps = max(2, int(np.abs(np.ceil(deltaC/step_size))))
        vsteps = np.linspace(self.current.level, newC, numSteps)
        for vstep in vsteps:
            if self.set_current(np.round(vstep,6)) == False:
                break
            time.sleep(step_time)
        return
    
    def set_limit_curr(self, clim):
        global curr_range
        curr_range = clim
        if self.mode == 'VOLT':
            print("Machine is currently set to",self.mode,"mode!")
            return
        self.current.range = clim
        if curr_range == float(self.inst.query('SOUR:RANG?')):    
            print('Current Range is set to',float(self.inst.query('SOUR:RANG?'))*1000,'mA')
        else:
            print("Curr Error. Current range set to",float(self.inst.query('SOUR:RANG?'))*1000,'mA')    
    
    def set_limit_volt(self, vlim):
        global volt_range
        volt_range = vlim
        if self.mode == 'CURR':
            print("Machine is currently set to",self.mode,"mode!")
            return
        self.voltage.range = vlim
        if volt_range == float(self.inst.query('SOUR:RANG?')):
            print('Voltage Range is set to',float(self.inst.query('SOUR:RANG?')),'V')
        else:
            print("Volt Error. Voltage range is set to",float(self.inst.query('SOUR:RANG?')),'V')
        
    
    
    
    
    
    
    
    '''fullsettings = {}
        #for setting in self.commandlist:
#            fullsettings[setting] = self.__getattr__(setting)
        #    fullsettings[setting] =getattr(self,setting)
        fullsetings[self.commandlist.voltage] = getattr(self,commandlist.voltage)
        return fullsettings
    '''
    '''errcmds           = {}
    errcmds['error']  = ':SYSTem:ERRor?'
    
    #commandlist = {}
    #commandlist['core']   = {}

    
    core = {}
    #0 or 1
    #'CURRent or VOLTage (case insensitive)
    #core['voltage'] = 'SOURce:CURRent:LEVEL'
    
    #voltage = {}
    #voltage['level']        = ':SOURce:LEVel'
    

    #current = {}   
    #current['level']        = 'SOURce:LEVel'
    
    #mode = {}
    #mode['set'] = ':SOURce:FUNCtion'
    
    #commandlist['core']    = core
    #commandlist['voltage'] = voltage
    #commandlist['current'] = current
    #commandlist['mode'] = mode
    
    #commandlist['level']    = ':SOURce:LEVel'
    #commandlist['mode']     = 'SOURce:FUNCtion'
    '''
    
    
    
    
    
    
    
    
    
    
    