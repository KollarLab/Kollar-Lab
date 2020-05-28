import pyvisa
from userfuncs import freeze

@freeze
class Generator():

    def __init__(self, address):
        rm=pyvisa.ResourceManager()
        self.inst=rm.open_resource(address)
        self.inst.write('*RST')

    @property
    def waveform(self):
        wtype = self.inst.query(':FUNC?')
        return wtype
    @waveform.setter
    def waveform(self, wtype='SIN'):
        self.inst.write(':FUNC {}'.format(wtype))
    
    @property
    def frequency(self):
        freq = self.inst.query_ascii_values(':FREQ?')
        return freq
    @frequency.setter
    def frequency(self, freq):
        self.inst.write(':FREQ {}'.format(freq))

    @property
    def volts(self):
        volt = self.inst.query_ascii_values(':VOLT?')
        return volt
    @volts.setter
    def volts(self, volt):
        self.inst.write(':VOLT {}'.format(volt))
    
    @property
    def offset(self):
        off = self.inst.query_ascii_values('VOLT:OFFS?')
        return off
    @offset.setter
    def offset(self, off):
        self.inst.write('VOLT:OFFS {}'.format(off))
    
    @property
    def reference(self):
        ref = self.inst.query(':ROSC:SOUR?')
        return ref
    @reference.setter
    def reference(self, ref):
        self.inst.write('ROSC:SOUR {}')

    @property
    def pulse_width(self):
        width = self.inst.query(':FUNC:PULS:WIDT?')
        return width
    @pulse_width.setter
    def pulse_width(self, width):
        self.inst.query_ascii_values(':FUNC:PULS:WIDT {}'.format(width))

    @property
    def output(self):
        out = self.inst.query(':OUT?')
        return out
    @output.setter
    def output(self, state):
        self.inst.write(':OUT {}'.format(state))
