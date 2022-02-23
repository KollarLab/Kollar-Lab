from .SCPIinst import SCPIinst

class SGS100A(SCPIinst):
    errcmds           = {}
    errcmds['error']  = 'SYST:ERR?'
    errcmds['serror'] = 'SYST:SERR?'
    
    commandlist = {}
    commandlist['core']   = {}
    commandlist['IQ']     = {}
    commandlist['Ref']    = {}
    commandlist['LO']     = {}
    commandlist['RefOut'] = {}
    
    core = {}
    core['Output'] = 'OUTPut:STATe'
    core['Power']  = 'SOURce:POWer'
    core['Phase']  = 'SOURce:PHASe'
    core['Freq']   = 'SOURce:FREQuency'
    core['Offset'] = 'SOURce:FREQuency:OFFSet'
    
    IQ = {}
    IQ['Mod']   = 'SOURce:IQ:STATe'
    IQ['Imp']   = 'SOURce:IQ:IMPairment'
    IQ['Ileak'] = 'SOURce:IQ:IMPairment:LEAKage:I'
    IQ['Qleak'] = 'SOURce:IQ:IMPairment:LEAKage:Q'
    
    Pulse = {}
    Pulse['Mod']      = 'SOURce:PULM:STATe'
    Pulse['Source']   = 'SOURce:PULM:SOURce'
    Pulse['Polarity'] = 'SOURce:PULM:POLarity'
    Pulse['Trig_out'] = 'CONNector:TRIGger:OMODe'
    
    Ref = {}
    Ref['Source']    = 'SOURce:ROSCillator:SOURce'
    Ref['Frequency'] = 'SOURce:ROSCillator:EXTernal:FREQuency'
    
    LO = {}
    LO['Source'] = 'SOURce:LOSCillator:SOURce'
    
    RefOut = {}
    RefOut['Source']    = 'CONNector:REFLo:OUTPut'
    RefOut['Frequency'] = 'SOURce:ROSCillator:OUTPut:FREQuency'
    
    commandlist['core']   = core
    commandlist['IQ']     = IQ
    commandlist['Pulse']  = Pulse
    commandlist['Ref']    = Ref
    commandlist['LO']     = LO
    commandlist['RefOut'] = RefOut

    def __init__(self, address):
        self.instrument_type = 'SGS'
        
        super().__init__(address, self.commandlist, self.errcmds)
    
    def enable_IQ(self, Ileak=0, Qleak=0):
        self.IQ.Mod = 'On'
        self.IQ.Imp = 'On'
        self.IQ.Ileak = Ileak
        self.IQ.Qleak = Qleak
    
    def disable_IQ(self):
        self.IQ.Mod = 'Off'
        
    def enable_pulse(self, source='EXT', polarity='NORM', trig_conn='PEMSource'):
        self.Pulse.Mod    = 'On'
        self.Pulse.Source = source
        self.Pulse.Polarity = polarity
        self.Pulse.Trig_out = trig_conn
    
    def disable_pulse(self):
        self.Pulse.Mod = 'Off'
    