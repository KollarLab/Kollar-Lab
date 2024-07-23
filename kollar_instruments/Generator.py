from userfuncs import freeze
from .SCPIinst import SCPIinst_c

class Keysight33500B(SCPIinst_c):
    errorcmd    = 'SYST:ERR?'
    commandlist = {}
    commandlist['core']  = {}
    commandlist['Pulse'] = {}
    
    core = {}
    core['Freq']     = 'FREQuency'
    core['Volt']     = 'VOLTage'
    core['Offset']   = 'VOLTage:OFFSet'
    core['Waveform'] = 'FUNCtion'
    core['Output']   = 'OUTPut'
    core['Duty_cycle'] = 'FUNCtion:SQUare:DCYCLe'
    
    pulse = {}
    pulse['Width'] = 'FUNCtion:PULSe:WIDTh'
    
    ref = {}
    ref['Source'] = 'ROSCillator:SOURce'
    
    commandlist['core']  = core
    commandlist['Pulse'] = pulse
    commandlist['Ref']   = ref
    
    def __init__(self, address):
        super().__init__(address, self.commandlist, self.errorcmd)