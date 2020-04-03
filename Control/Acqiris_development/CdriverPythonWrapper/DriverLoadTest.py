# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 13:42:42 2020

@author: Kollarlab
"""
#from ctypes import *
import ctypes
import comtypes
import os
import time
#import subprocess

#import re
import scipy
import pylab
#import tarfile
#import struct
#import glob
import numpy
import time
import sys


#%%

from comtypes.client import GetModule
from comtypes.client import CreateObject

if not hasattr(sys, "frozen"):
    GetModule("C:\Program Files\IVI Foundation\IVI\Bin\AqMD2_64.dll")
    
#    os.chdir(r'C:\Program Files\IVI Foundation\IVI\Microsoft.NET\Framework64\v4.0.30319\Acqiris.AqMD3 3.4.2029')
#    GetModule("Acqiris.AqMD3.Fx40.dll")
    
    ###junk
#    GetModule("C:\Program Files\IVI Foundation\IVI\Bin\AqMD3_64.dll")
#    GetModule("C:\Program Files (x86)\IVI Foundation\IVI\Bin\AqMD3.dll")
    
    
#from comtypes.gen import AqMD3Lib
from comtypes.gen import AqMD2Lib

#myCard = CreateObject('AqMD2.AqMD2')
myCard = CreateObject('AqMD2.AqMD2')

hardwareAddress = 'PXI23::0::0::INSTR'
#initOptions = 'Simulate=True'
initOptions ='Simulate=True, DriverSetup= CAL=0, Trace=false, model = U5309A'

myCard.Initialize(hardwareAddress, True, True, initOptions )


#%%
#IVIFounationPath = os.path.join(r'C:\Program Files\IVI Foundation\IVI\Microsoft.NET\Framework64\v2.0.50727','IviFoundationSharedComponents 1.4.1')
IVIFounationPath = "C:\\Program Files\\IVI Foundation\\IVI\\Microsoft.NET\\Framework64\\v2.0.50727\\IviFoundationSharedComponents 1.4.1\\"
IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
if not IVIFounationPath in sys.path:
    sys.path.append(IVIFounationPath)
if not IVIbinPath in sys.path:
    sys.path.append(IVIbinPath)

netdllfolder = r'C:\Program Files\IVI Foundation\IVI\Microsoft.NET\Framework64\v4.0.30319\Acqiris.AqMD3 3.4.2029'
netdllpath = os.path.join(netdllfolder, 'Acqiris.AqMD3.Fx40')

netdllpath2 = os.path.join(netdllfolder, 'Acqiris.AqMD3Version')


if not netdllpath in sys.path:
    sys.path.append(netdllpath)
if not netdllpath2 in sys.path:
    sys.path.append(netdllpath2)
if not netdllfolder in sys.path:
    sys.path.append(netdllfolder)

#import Ivi
import clr   
#clr.FindAssembly(netdllpath) 
clr.FindAssembly('Acqiris.AqMD3.Fx40') 
#temp = clr.AddReference('Acqiris.AqMD3')
temp = clr.AddReference(netdllpath)
#temp = clr.AddReference(netdllpath2)

import Acqiris.AqMD3
#from Acqiris.AqMD3 import AqMD3
#from Acqiris.AqMD3 import AqMD3ChannelCollection

hardwareAddress = 'PXI23::0::0::INSTR'
#hardwareAddress = 'PXI21::0::0::INSTR' #junk address
idQuery = True
#idQuery = False
reset   = True
#reset   = False
initOptions ='Simulate=False,  DriverSetup= model = SA220P'
#initOptions ='Simulate=True,  DriverSetup= model = SA220P'
#initOptions ='Simulate=True,  DriverSetup= model = U5309A'

card = Acqiris.AqMD3.AqMD3(hardwareAddress, idQuery, reset, initOptions)
#card = Acqiris.AqMD3.AqMD3(hardwareAddress, idQuery, reset)

#%%
#cdllpath = r'C:\Program Files\IVI Foundation\IVI\Bin\AqMD3_64.dll'
#cdllpath2 = r'C:\Program Files\IVI Foundation\IVI\Bin\AqMD2_64.dll'
##netdllpath = r'C:\Program Files\IVI Foundation\IVI\Microsoft.NET\Framework64\v4.0.30319\Acqiris.AqMD3 3.4.2029\Acqiris.AqMD3.Fx40.dll'
#
##netdllfolder = r'C:\Program Files\IVI Foundation\IVI\Bin'
##netdllpath = os.path.join(netdllfolder, 'AqMD3_64')
#
#netdllfolder = r'C:\Program Files\IVI Foundation\IVI\Microsoft.NET\Framework64\v4.0.30319\Acqiris.AqMD3 3.4.2029'
#netdllpath = os.path.join(netdllfolder, 'Acqiris.AqMD3.Fx40')
##netdllpath = os.path.join(netdllfolder, 'Acqiris.AqMD3')
#
##netdllpath = r'C:\Program Files\IVI Foundation\IVI\Microsoft.NET\Framework64\v4.0.30319\Acqiris.AqMD3 3.4.2029\Acqiris.AqMD3.Fx40'
#
###clrpath = r'C:\Users\Kollarlab\Acqiris\pythonnet-2.4.0'
##clrpath = r'C:\Users\Kollarlab\anaconda3\lib\site-packages (2.4.0)'
##if not clrpath in sys.path:
##    sys.path.append(clrpath)
#
#if not netdllpath in sys.path:
#    sys.path.append(netdllpath)
#if not netdllfolder in sys.path:
#    sys.path.append(netdllfolder)
#    
#import clr
#
#temp = clr.AddReference(netdllpath)
##temp2 = clr.AddReference(cdllpath)
#
#
#hardwareAddress = 'PXI23::0::0::INSTR'
#idQuery = True
#reset   = True
##initOptions ='Simulate=False,  DriverSetup= model = SA220P'
##initOptions ='Simulate=True,  DriverSetup= model = SA220P'
#initOptions ='Simulate=True,  DriverSetup= model = U5309A'
##initOptions ='Simulate=True'
##
#
#
##temp.GetModule('Acqiris')
##temp.GetModule('AqMD3')
##temp.GetModule('Acqiris.AqMD3')
##temp.GetModule('Acqiris.AqMD3.Fx40')
#import Acqiris
#import Acqiris.AqMD3
##from Acqiris import *
##from Acqiris.AqMD3 import *
#import Ivi
#import Ivi.Driver
#
##from Acqiris import AqMD3
##from AqMD3 import AqMD3 as cardDriver
##
##card = cardDriver(hardwareAddress, idQuery, reset, initOptions)
#
##temp.GetModule('Acqiris')
##import Acqiris
##import Acqiris.AqMD3
##from Acqiris.AqMD3 import AqMD3
###import Acqiris.AqMD3.AqMD3
##
#
##card = AqMD3(hardwareAddress, idQuery, reset)
#card = Acqiris.AqMD3.AqMD3(hardwareAddress, idQuery, reset)
##card = Acqiris.AqMD3.AqMD3(hardwareAddress, idQuery, reset, initOptions)
##card = Acqiris.AqMD3.IAqMD3(hardwareAddress, idQuery, reset, initOptions)
#
#
#
#for chan in list(card.Channels):
#    print (chan.Name)
#
#
#
##driver = temp.CreateInstance()
##driver= temp.GetModule(netdllpath)
##driver = temp.CreateInstance('Acqiris.AqMD3.Fx40')
##driver = temp.GetModule('Acqiris.AqMD3.Fx40')
##driver = clr.clrModule(netdllpath)
##clr.Ivi(netdllpath)
##clr.Acqiris.AqMD3.init()
##driver = clr.Acqiris.AqMD3.(hardwareAddress, idQuery, reset, initOptions)
##driver = clr.Acqiris.AqMD3.Initiate(hardwareAddress, idQuery, reset, initOptions)
#
##lib = ctypes.CDLL(netdllpath)
###class _CMathematics(ctypes.Structure): 
###    pass
###
###CMathematics = ctypes.POINTER(_CMathematics)
##
##
##GetModule(dllpath)
#
#
