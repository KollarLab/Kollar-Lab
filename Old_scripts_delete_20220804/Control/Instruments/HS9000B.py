# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 18:04:20 2020

@author: Kollarlab
"""

from Instruments.SCPIinst import SCPIinst

class HS9000B(SCPIinst):
    errcmds           = {}
    errcmds['error']  = 'SYST:ERR?'
    errcmds['serror'] = 'SYST:SERR?'
    
    commandlist = {}
    commandlist['core']   = {}
    commandlist['Mod']     = {}
    commandlist['Ref']    = {}
    commandlist['RefOut'] = {}
    
    core = {}
    core['Output']  = 'CH1:PWR:RF'
    core['Power']  = 'CH1:PWR'
    core['Freq']   = 'CH1:FREQ'
    core['Phase'] = 'CH1:PHASE'
    
    Mod = {}
    Mod['Mod']   = 'SOURce:IQ:STATe'
    Mod['Imp']   = 'SOURce:IQ:IMPairment'
    
    Ref = {}
    Ref['Source']    = 'SOURce:ROSCillator:SOURce'
    Ref['Frequency'] = 'SOURce:ROSCillator:EXTernal:FREQuency'
    
    RefOut = {}
    RefOut['Source']    = 'CONNector:REFLo:OUTPut'
    RefOut['Frequency'] = 'SOURce:ROSCillator:OUTPut:FREQuency'
    
    commandlist['core']   = core
    commandlist['Ref']    = Ref
    commandlist['RefOut'] = RefOut

    def __init__(self, address):
        super().__init__(address, self.commandlist, self.errcmds) 