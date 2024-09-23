from .SCPIinst import SCPIinst
import numpy as np
import time

class Keysight33622A(SCPIinst):
    errcmds = {}
    errcmds['error'] = 'SYST:ERR?'

    commandlist = {}
    commandlist['core'] = {}

    core = {}

    core['phase_unit'] = 'UNIT:ANGLe' #could be DEGree, RADian, SECond, DEFault

    core['Ch1_waveform'] = 'SOURce1:FUNCtion' #could be SINusoid, DC, or others
    core['Ch1_frequency'] = 'SOURce1:FREQuency' 
    core['Ch1_voltage'] = 'SOURce1:VOLTage' #by default the unit is Vpp
    core['Ch1_phase'] = 'SOURce1:PHASe' #by default the unit is DEGree
    core['Ch1_offset'] = 'SOURce1:VOLTage:OFFSet'  #by default the unit is V
    core['Ch1_output'] = 'OUTPut1' #takes ON/1, OFF/0

    core['Ch2_waveform'] = 'SOURce2:FUNCtion'
    core['Ch2_frequency'] = 'SOURce2:FREQuency'
    core['Ch2_voltage'] = 'SOURce2:VOLTage'
    core['Ch2_phase'] = 'SOURce2:PHASe'
    core['Ch2_offset'] = 'SOURce2:VOLTage:OFFSet' 
    core['Ch2_output'] = 'OUTPut2'

    commandlist['core'] = core

    def __init__(self, address, reset = False):
        '''
        Initialize connection with instrument, reset it if specified, and clear all errors
        '''
        self.instrument_type = 'Keysight33622A'

        super().__init__(address, self.commandlist, self.errcmds, reset)
        self.Ch1_output = 0
        self.Ch2_output = 0

        if reset:
            self.reset()

    def reset(self):
        self.inst.write('*RST; *CLS')
        self.inst.write('OUTP1 OFF')
        self.inst.write('OUTP2 OFF')

    def Ch1_dc_voltage_ramp(self, newV, step_size = 0.005, step_time = 0.001):
        self.Ch1_waveform = 'DC'
        #self.Ch1_output = 1
        deltaV = newV - self.Ch1_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV/step_size))))
        vsteps = np.linspace(self.Ch1_offset, newV, numSteps)
        for vstep in vsteps:
            self.Ch1_offset = np.round(vstep,4)
            time.sleep(step_time)

    def Ch2_dc_voltage_ramp(self, newV, step_size = 0.005, step_time = 0.001):
        self.Ch2_waveform = 'DC'
        #self.Ch2_output = 1
        deltaV = newV - self.Ch2_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV/step_size))))
        vsteps = np.linspace(self.Ch2_offset, newV, numSteps)
        for vstep in vsteps:
            self.Ch2_offset = np.round(vstep,4)
            time.sleep(step_time)

    def Ch1_sin_gen(self, amplitude, frequency, phase=0, offset=0, step_size = 0.005, step_time = 0.001):
        self.Ch1_waveform = 'SIN'
        self.Ch1_frequency = frequency
        self.Ch1_voltage = str(amplitude) + ' Vpp'
        self.Ch1_phase = phase
        #self.Ch1_output = 1

        deltaV = offset - self.Ch1_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV/step_size))))
        vsteps = np.linspace(self.Ch1_offset, offset, numSteps)
        for vstep in vsteps:
            self.Ch1_offset = np.round(vstep,4)
            time.sleep(step_time)

    def Ch2_sin_gen(self, amplitude, frequency, phase=0, offset=0, step_size = 0.005, step_time = 0.001):
        self.Ch2_waveform = 'SIN'
        self.Ch2_frequency = frequency
        self.Ch2_voltage = str(amplitude) + ' Vpp'
        self.Ch2_phase = phase
        #self.Ch2_output = 1

        deltaV = offset - self.Ch2_offset
        numSteps = max(2, int(np.abs(np.ceil(deltaV/step_size))))
        vsteps = np.linspace(self.Ch2_offset, offset, numSteps)
        for vstep in vsteps:
            self.Ch2_offset = np.round(vstep,4)
            time.sleep(step_time)

    def phase_sync(self):
        self.inst.write('SOURce1:PHASe:SYNChronize') #source1 or source2 mean nothing here, simply give two channels a common reference point