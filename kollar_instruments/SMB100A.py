# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 12:10:12 2020

@author: Kollarlab
"""

from .SCPIinst import SCPIinst

class SMB100A(SCPIinst):
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
    
    Ref = {}
    Ref['Source']    = 'SOURce:ROSCillator:SOURce'
    Ref['Frequency'] = 'SOURce:ROSCillator:EXTernal:FREQuency'
    
    commandlist['core']   = core
    commandlist['Ref']    = Ref

    def __init__(self, address):
        super().__init__(address, self.commandlist, self.errcmds) 