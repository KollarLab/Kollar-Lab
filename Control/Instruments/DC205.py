# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:10:29 2020

@author: Kollarlab
"""

from Instruments.SCPIinst import SCPIinst

class DC205(SCPIinst):
    commandlist = {}
    commandlist['core'] = {}
    
    core = {}
    core['Volt'] = 'VOLT'
    core['Output'] = 'SOUT'
    core['Range'] = 'RNGE'
    core['Ground'] = 'ISOL'
    
    commandlist['core'] = core
    
    def __init__(self, address):
        super().__init__(address, self.commandlist)