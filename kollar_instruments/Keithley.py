from .SCPIinst import SCPIinst

#currently, we only want to use keithley as a stable source, so only set_voltage and set_current should be important
#If trying to measure something, then a flashing 'cmpl' or value to the right of it means we've hit a range limit, and the range needs to be changed to measure what we have

class keithley(SCPIinst):
    errcmds           = {}
    errcmds['error']  = 'STAT:QUE?'
    
    commandlist = {}
    commandlist['core']   = {}

    
    core = {}
    core['Output']  = 'OUTPut:STATe' #0 or 1
    core['mode']    = 'SOURce:FUNCtion:MODE' #'CURRent or VOLTage
    core['voltage'] = 'SOURce:CURRent:LEVEL'
    
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

    def __init__(self, address):
        self.instrument_type = 'keithley'
        
        super().__init__(address, self.commandlist, self.errcmds)
    
    def set_voltage(self,value):
        self.mode = 'VOLTage'
        self.voltage.mode = 'FIXed'
        self.voltage.output_range = value
        self.voltage.level = value
    
    def set_current(self,value):
        self.mode = 'CURRent'
        self.current.mode = 'FIXed'
        self.current.output_range = value
        self.current.level = value
    
    def set_voltage_range(self,value):
        #given an expected measurment value, sets the measurment range accordingly
        self.voltage.meas_range = value*2
        self.voltage.cmpl = value*2
    
    def set_current_range(self,value):
        #given an expected measurment value, sets the measurment range accordingly
        self.current.meas_range = value*2
        self.current.cmpl = value*2

