import pyvisa
from userfuncs import freeze

@freeze
class Generator():

    def __init__(self, address):
        rm=pyvisa.ResourceManager()
        self.inst=rm.open_resource(address)
        self.inst.write('*RST')

    def get_settings(self):
        settings = {}
        for setting in dir(self):
            if not setting.startswith('__') and not callable(getattr(self,setting)) and setting!='inst':
                settings[setting] = getattr(self,setting)
        return settings

    def set_settings(self, settings):
        for key in settings.keys():
            try:
                setattr(self, key, settings[key])
            except:
                print('Unknown attribute {}, ignoring'.format(key))
            
    @property
    def waveform(self):
        wtype = self.inst.query(':FUNC?')
        return wtype
    @waveform.setter
    def waveform(self, wtype):
        funclist = ['ARB', 'NOIS','PRBS', 'PULS', 'RAMP', 'SIN', 'SQU']
        if wtype not in funclist:
            print('Invalid option, valid options are: {}'.format(funclist))
        self.inst.write(':FUNC {}'.format(wtype))
    
    @property
    def frequency(self):
        freq = self.inst.query_ascii_values(':FREQ?')
        return freq[0]
    @frequency.setter
    def frequency(self, freq):
        self.inst.write(':FREQ {}'.format(freq))

    @property
    def volts(self):
        volt = self.inst.query_ascii_values(':VOLT?')
        return volt[0]
    @volts.setter
    def volts(self, volt):
        self.inst.write(':VOLT {}'.format(volt))
    
    @property
    def offset(self):
        off = self.inst.query_ascii_values('VOLT:OFFS?')
        return off[0]
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
        width = self.inst.query_ascii_values(':FUNC:PULS:WIDT?')
        return width[0]
    @pulse_width.setter
    def pulse_width(self, width):
        self.inst.write(':FUNC:PULS:WIDT {}'.format(width))

    @property
    def output(self):
        out = self.inst.query_ascii_values(':OUTP?')
        return out[0]
    @output.setter
    def output(self, state):
        self.inst.write(':OUTP {}'.format(state))
