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
    commandlist['Ref']    = Ref
    commandlist['LO']     = LO
    commandlist['RefOut'] = RefOut

    def __init__(self, address):
        self.instrument_type = 'SGS'
        
        super().__init__(address, self.commandlist, self.errcmds) 