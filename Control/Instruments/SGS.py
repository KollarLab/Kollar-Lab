from Instruments.SCPIinst import SCPIinst

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
    core['Output']  = 'OUTPut:STATe'
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
        super().__init__(address, self.commandlist, self.errcmds) 
       
#class RFgen():
#
#    def __init__(self, address):
##        rm=visa.ResourceManager()
#        rm=pyvisa.ResourceManager()
#        self.inst=rm.open_resource(address)
#        self.inst.write('*RST')
#
#    @property
#    def settings(self):
#        return {}
#    @settings.setter
#    def settings(self, fullsettings):
#        print('Tried to set SGS settings')
#
#    def set_Amp_V(self,amp):
#        self.inst.write('SOURce:POWer {} V'.format(amp))
#        
#    def set_Amp(self,amp):
#        #set power in dB
#        self.inst.write('SOURce:POWer {}'.format(amp))
#
#    def set_Freq(self,freq):
#        self.inst.write('SOURce:FREQuency:CW {} GHz'.format(freq))
#        
#    def set_Offset(self, offset):
#        self.inst.write('SOURce:FREQuency:OFFSet {} MHz'.format(offset))
#        
#    def set_External_Reference(self, freq=10):
#        self.inst.write(':SOURce:ROSCillator:SOURce EXTernal')
#        self.inst.write(':SOURce:ROSCillator:EXTernal:FREQuency {} MHz'.format(freq))
#        
#    def set_Internal_Reference(self):
#        self.inst.write(':SOURce:ROSCillator:SOURce Internal')
#        
#    def set_External_LO(self):
#        self.inst.write(':SOURce:LOSCillator:SOURce EXTernal')
#    
#    def set_Internal_LO(self):
#        self.inst.write(':SOURce:LOSCillator:SOURce INTernal')
#    
#    def set_RefLO_output(self, output = 'Ref', freq = 10):
#        self.inst.write(':CONNector:REFLo:OUTPut {}'.format(output))
#        if output == 'Ref':
#            self.inst.write(':SOURce:ROSCillator:OUTPut:FREQuency {} MHz'.format(freq))
#        
#    def set_Phase(self, phase):
#        #degrees
#        self.inst.write(':SOURce:PHASe {}'.format(phase))
#        
#        
#    def power_On(self):
#        self.inst.write(':OUTPut:STATe ON')
#        
#    def power_Off(self):
#        self.inst.write(':OUTPut:STATe OFF')
#        
#        
#    def mod_On(self):
#        self.inst.write(':SOURce:IQ:STATe ON')
#        
#    def mod_Off(self):
#        self.inst.write(':SOURce:IQ:STATe OFF')
#        
#    def imp_Off(self):
#        self.inst.write(':SOURce:IQ:IMPairment OFF')
#        
#    def imp_On(self):
#        self.inst.write(':SOURce:IQ:IMPairment ON')
#        
#    def leak_I(self, ileak):
#        self.inst.write(':SOURce:IQ:IMPairment:LEAKage:I {}'.format(ileak))
#        
#    def leak_Q(self, qleak):
#        self.inst.write(':SOURce:IQ:IMPairment:LEAKage:Q {}'.format(qleak))     
#
#    def send_cmd(self,cmd):
#        self.inst.write(cmd)
#
#    def query(self,cmd):
#        print(self.inst.query(cmd))
#
##    def run(self):
##        self.inst.write('OUTPut:STATe ON')