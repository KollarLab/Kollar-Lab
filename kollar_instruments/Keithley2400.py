from .SCPIinst import SCPIinst
import numpy as np

#currently, we only want to use keithley as a stable source, so only set_voltage and set_current should be important
#If trying to measure something, then a flashing 'cmpl' or value to the right of it means we've hit a range limit, and the range needs to be changed to measure what we have

class Keithley2400(SCPIinst):
    errcmds           = {}
    errcmds['error']  = 'STAT:QUE?'
    
    commandlist = {}
    commandlist['core']   = {}

    
    core = {}
    core['Output']  = 'OUTPut:STATe' #0 or 1
    core['mode']    = 'SOURce:FUNCtion:MODE' #'CURRent or VOLTage (case insensitive)
    #core['voltage'] = 'SOURce:CURRent:LEVEL'
    
    voltage = {}
    voltage['mode']         = 'SOURce:VOLTage:MODE' #FIXed should be good for all purposes
    voltage['level']        = 'SOURce:VOLTage:LEVEL'
    voltage['output_range'] = 'SOURce:VOLTage:Range'
    voltage['meas_range']   = 'SENSe:VOLTage:Range'
    voltage['cmpl']         = 'SENSe:VOLTage:PROTection'
    voltage['measurement']  = 'MEASure:VOLTage'

    current = {}
    current['mode']         = 'SOURce:CURRent:MODE' #FIXed should be good for all purposes
    current['level']        = 'SOURce:CURRent:LEVEL'
    current['output_range'] = 'SOURce:CURRent:Range'
    current['meas_range']   = 'SENSe:CURRent:Range'
    current['cmpl']         = 'SENSe:CURRent:PROTection'
    current['measurement']  = 'MEASure:CURRent'

    commandlist['core']    = core
    commandlist['voltage'] = voltage
    commandlist['current'] = current

    def __init__(self, address, reset=True, mode='current', volt_range=10, Output=1):
        self.instrument_type = 'keithley2400'
        
        super().__init__(address, self.commandlist, self.errcmds, reset = reset, baud_rate = 9600)
        self.inst.read_termination = '\r'
        self.inst.write_termination = '\n'
        #self.mode=mode
        #if mode == 'current':
        #    self.set_voltage_meas_range(volt_range)
        #self.Output=1
    
    def set_voltage(self,value): #0.13s operation
        self.voltage.output_range = value
        self.voltage.level = value
    
    def set_current(self,value, force=False, verbose=False): #0.13s operation
        if value>0.01:
            print('Caution, current units are A not mA. Set "force" flag to true')
            if not force:
                return
        self.current.output_range = value
        self.current.level = value
        if verbose:
            print(self.current.measurement[1])
    
    def set_voltage_meas_range(self,value): #0.13s operation
        #the user should set the range to encompass the expected measurment value
        self.voltage.cmpl = value
        self.voltage.meas_range = value

    def set_current_meas_range(self,value): #0.13s operation
        #the user should set the range to encompass the expected measurment value
        self.current.cmpl = value
        self.current.meas_range = value
    
    def measure(self): #0.2s operation
        if self.mode == 'VOLT':
            return self.current.measurement[1]
        elif self.mode == 'CURR':
            return self.voltage.measurement[0]
        else:
            raise('mode is not current or voltage')
    
    def ramp_current(self,value, step_size=0.0002):
        init_val = self.current.level
        delta = abs(value-init_val)
        ramp = np.linspace(init_val, value, int(delta/step_size)+2)
        for amp in ramp:
            self.set_current(amp)
        self.measure()
    
    def ramp_voltage(self,value):
        npoints=10
        for i in range(0,npoints):
            self.set_voltage(value*i/npoints)




