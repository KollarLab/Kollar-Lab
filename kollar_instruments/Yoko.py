from .SCPIinst import SCPIinst
import numpy as np

class Yoko(SCPIinst):
    errcmds           = {}
    errcmds['error']  = ':SYSTem:ERRor?'
    
    commandlist = {}
    commandlist['core'] = {}

    core = {}
    core['Output']  = ':OUTPut' #0 or 1, ON or OFF
    core['mode']    = ':SOURce:FUNCtion' #'CURRent or VOLTage (case insensitive)
    
    voltage = {}
    voltage['level'] = ':SOURce:LEVEL' #actual output level [V]
    voltage['range'] = ':SOURce:RANGe' #max limit [V]

    current = {}
    current['level'] = ':SOURce:LEVEL' #actual output level [A]
    current['range'] = ':SOURce:RANGe' #max limit [A]

    commandlist['core']    = core
    commandlist['voltage'] = voltage
    commandlist['current'] = current

        
    '''Here, vrange, and crange are the ranges set at init. They get copied to the class variables to be 
    used elsewhere for safety checks. Here, input voltage and current ranges in A and V.
    Allowed crange inputs: .001, .01, .1, .2
    Allowed vrange inputs: .01, .1, .2, 1, 10, 30
    Range limits are remembered when switching modes.'''

    def __init__(self, address, reset=True, mode='CURR', volt_range=10, curr_range=0.01, Output=1):
        self.instrument_type = 'Yokogawags200'
        
        super().__init__(address, self.commandlist, self.errcmds, reset = reset, baud_rate = 9600)
        self.inst.read_termination = '\n'
        self.inst.write_termination = '\r'

        self.mode = mode

        self.current.range = curr_range
        self.voltage.range = volt_range
    
    def set_voltage(self,value): #0.13s operation
        if self.mode != 'VOLT':
            raise ValueError("mode is not voltage")
        self.voltage.level = value
    
    def set_current(self,value, force=False): #0.13s operation
        if self.mode != 'CURR':
            raise ValueError("mode is not current")
        if value>0.01:
            print('Caution, current units are A not mA. Set "force" flag to true')
            if not force:
                return
        self.current.level = value
    
    def current_ramp(self,value, step_size=0.0001):
        init_val = self.current.level
        delta = abs(value-init_val)
        ramp = np.linspace(init_val, value, max(2,int(delta/step_size)+1))
        for amp in ramp:
            self.set_current(amp)
    
    def voltage_ramp(self,value,step_size=0.001):
        init_val = self.voltage.level
        delta = abs(value-init_val)
        ramp = np.linspace(init_val, value, max(2,int(delta/step_size)+1))
        for volt in ramp:
            self.set_voltage(volt)