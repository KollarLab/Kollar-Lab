# -*- coding: utf-8 -*-
'''4/3/20
Modifying old version of the Acqiris object that ran off of the C driver to use 
the new python beta version of the driver from acqiris.

card = Acqiris(hardwareAddress)  '''

import ctypes


# import comtypes
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

#import pickle
#import datetime
#import itertools
import sys


import AqMD3 as driver #this is the Acqiris python driver



class Acqiris(object):
    def __init__(self, ResourceName, simulate = False):  
        #look to the right place.
        IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
        
        if not IVIbinPath in sys.path:
            sys.path.append(IVIbinPath)
        
#        self._prefix = "AqMD3"
#        self._lib = ctypes.cdll.LoadLibrary(library)
#        self.visession = ctypes.c_int()
#        #fill in and store the necessary driver-level attribute ID numbers
#        self._fillHardwareIDs()
            
        self._fillHardwareKeys()
        
        
        #set default values. Eventrually this should probable load from some default config file
        self._loadDefaultConfig()
        
        self.simulate = simulate
        self.hardwareAddress = ResourceName
        self.ReInitialize()
        
        
        #some flags for controlling the read
        self.armed = False
        self.acquisitionFinished = False
        
#    ##############################
#    #hardwae properties
#    
#    def get_samples(self):
#        return self.GetDriverAttribute('samples')
#    
#    def set_samples(self,val):
#        self.SetDriverAttribute('samples', val)
#    
#    samples = property()
#    
#    
#    
#    #############################
#        
            
    def SetParams(self):
        '''Autoinititalize function that will set up the card with all the driver
        parameters stored in this object. It will push all the settings and do cross checks.
        
        This function loads all the settings to the card, so it sets up the trigger and the
        data acquisition channels. Then it calls Configure Acquisition to set up things
        like averaging and sample/sample rate'''
        
        #configure the trigger
        self.ConfigureTrigger(Source = self.triggerSource, Level = self.triggerLevel, Slope = self.triggerSlope, Mode = self.triggerMode)
        
        #configure the channels:
        for chind in [1,2]:
            if chind in self.activeChannels:
                enabled = True
            else:
                enabled= False
            self.ConfigureChannel(channelNum = chind, Range = self.channelRange, offset = self.channelOffset, enabled = enabled)
            
        #configure the acquisition
        self.ConfigureAcquisition(self.samples, self.sampleRate, self.segments)
        
    def GetParams(self, driverInfo = False):
        '''Autoget function for params. Will querry the hardware and package all the settings, as well as 
        update the fields of the python software object.'''
        
        #going to need this
        for key in pickledict.keys():
                    setattr(self, key, pickledict[key])
        
        if self.verbose:
            print('    ')
        
        hardwareSettings = {}
        for key in self.hardwareIDs.keys():
            ID = self.hardwareIDs[key][0]
            driverInfoFlag =  (ID > self.driverBases['IVI_INHERENT_ATTR_BASE'] + 500) and  (ID < self.driverBases['IVI_INHERENT_ATTR_BASE'] + 600)
            serialNumberFlag = (ID == self.driverBases['IVI_SPECIFIC_ATTR_BASE'] + 8)
            if (driverInfoFlag or serialNumberFlag) and not driverInfo:
                #this should be driver and firmware info, not acquisition settings for the card.
                pass
            else:
                print(key)
                #channel properties need special handling
                if key[0:7] == 'channel':
                    #this is a channel specific variable and needs to be access with a repeat capability identifier
                    val = self.GetDriverAttribute(key, recap = 'Channel1')
                    val2 = self.GetDriverAttribute(key, recap = 'Channel2')
                    if not val == val2:
                        print('Warning: ' + key+ ' set differently on the two channels')
                        
                    if  key == 'channelEnabled':
                        flags = numpy.asarray([val,val2])
                        inds = numpy.where(flags == True)[0]
                        full = numpy.asarray([1,2])
                        actives = full[inds]
                        hardwareSettings['activeChannels'] = actives
                        if self.verbose:
                            print('activeChannels : ' + str(actives))
                    else:
                        hardwareSettings[key] = val
                        if self.verbose:
                            print(key + ' : ' + str(val) + ' , ' + str(val2))
                            
                #trigger properties need special handling        
                elif key[0:7] == 'trigger':
                    #trigger settings are mixed convnetional and recap identifier settings
                    if (key == 'triggerSource') or (key == 'triggerDelay'):
                        val = self.GetDriverAttribute(key)
                    else:
                        trigSource = self.GetDriverAttribute('triggerSource')
                        val = self.GetDriverAttribute(key, recap = trigSource)
                    hardwareSettings[key] = val
                else:
                    val = self.GetDriverAttribute(key)
                    hardwareSettings[key] = val
                
                #print all the results    
                if self.verbose:
                    print(key + ' : ' + str(val))
                
   
        if self.verbose:
            print('    ')
            
        return hardwareSettings
            

#    def SetDriverAttribute(self, keyword, val, recap = None):
#        '''Set function for atributes of the Cdriver.
#        takes in a python string keyword which needs to be a key of 
#        self.hardwareIDs. Then it uses thes stored information in there
#        to call the lower level _Set function. 
#        Note: if the variable is associated with a repeat capability identifier
#        then that needs to be handled separately. 
#        Also: this function will not check that the value is of the right type. User beware.'''
#        if keyword in self.hardwareIDs.keys():
#            #known driver attribute
#            [attrNum, attrType] = self.hardwareIDs[keyword]
#            self._SetAttribute(recap, attrNum, val, attrType)
#        else:
#            raise ValueError('Invalid keyword. Either a typo, or it needs to be looked up in AqMD3 and added.')

#    def _SetAttribute(self, RepCapIdentifier, Attribute, Value, AttributeType) :
#        '''Raw set attribute frunction from Acqiris '''
#
#        # Encode the RepCapIdentifier if needed
#        if (RepCapIdentifier != None):
#            RepCapIdentifier = ctypes.c_char_p(RepCapIdentifier.encode('utf-8'))
#                   
#        if (AttributeType == 'ViBoolean') :            
#            self.call('SetAttributeViBoolean',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute,
#                  ctypes.c_bool(Value))
#            
#        if (AttributeType == 'ViInt32') :
#            self.call('SetAttributeViInt32',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute,
#                 ctypes.c_int32(Value))
#        
#        if (AttributeType == 'ViInt64') :            
#            self.call('SetAttributeViInt64',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute,
#                  ctypes.c_int64(Value))
#
#        if (AttributeType == 'ViReal64') :
#            self.call('SetAttributeViReal64',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute,
#                  ctypes.c_double(Value))
#            
#        if (AttributeType == 'ViString') :
#            buffer = ctypes.create_string_buffer(256)
#            self.call('SetAttributeViString',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute,
#                  ctypes.c_char_p(Value.encode('utf-8')))
            
#    def GetDriverAttribute(self, keyword, recap = None):
#        '''Get function for atributes of the Cdriver.
#        takes in a python string keyword which needs to be a key of 
#        self.hardwareIDs. Then it uses thes stored information in there
#        to call the lower level _Get function. 
#        Note: if the variable is associated with a repeat capability identifier
#        then that needs to be handled separately. '''
#        if keyword in self.hardwareIDs.keys():
#            #known driver attribute
#            [attrNum, attrType] = self.hardwareIDs[keyword]
#            return self._GetAttribute(recap, attrNum, attrType)
#        else:
#            raise ValueError('Invalid keyword. Either a typo, or it needs to be looked up in AqMD3 and added.')
#
#    def _GetAttribute(self, RepCapIdentifier, Attribute, AttributeType) :
#        '''Raw get attribute function from Acqiris '''
#        
#        # Encode the RepCapIdentifier if needed
#        if (RepCapIdentifier != None):
#            RepCapIdentifier = ctypes.c_char_p(RepCapIdentifier.encode('utf-8'))
#                   
#        if (AttributeType == 'ViBoolean') :
#            buffer = ctypes.c_bool();
#            self.call('GetAttributeViBoolean',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute, 
#                  ctypes.byref(buffer))
#            return buffer.value      
#            
#        if (AttributeType == 'ViInt32') :
#            buffer = ctypes.c_int32()
#            self.call('GetAttributeViInt32',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute, 
#                  ctypes.byref(buffer))
#            return buffer.value
#        
#        if (AttributeType == 'ViInt64') :
#            buffer = ctypes.c_int64()
#            self.call('GetAttributeViInt64',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute, 
#                  ctypes.byref(buffer))
#            return buffer.value
#
#        if (AttributeType == 'ViReal64') :
#            buffer = ctypes.c_double()
#            self.call('GetAttributeViReal64',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute, 
#                  ctypes.byref(buffer))
#            return buffer.value
#            
#        if (AttributeType == 'ViString') :
#            buffer = ctypes.create_string_buffer(256)
#            self.call('GetAttributeViString',
#                  self.visession,
#                  RepCapIdentifier,
#                  Attribute, 
#                  256,
#                  ctypes.byref(buffer))
#            return buffer.value.decode('ascii')        
#           
#        return ''

#    def call(self, funcname, *args):
#        method_to_call = getattr(self._lib, self._prefix + '_' + funcname)
#        return_value = method_to_call(*args)        
#        if return_value:
#            error_message = ctypes.create_string_buffer(256)
#            self.call('error_message',
#                      self.visession,
#                      return_value, 
#                      ctypes.byref(error_message))
#            # raise IviCError(error_message.value)   
#            raise ValueError(error_message.value)      
            
            
    def close(self):    
        self.call('close', self.visession) 
        
    def SelfCalibrate(self):
        print('Calibrating...')
        self.driver.Calibration.SelfCalibrate()
#        self.call('SelfCalibrate', self.visession) 

    ##################################################
    
    def Arm(self):
        '''Initialize the acquisition and get ready to take and read data.
        Tells the card to get ready and wait for triggers.
        
        Unlike it big brother ArmAndWait, this function will allow python to 
        do stuff while the card is waiting for a trigger. E.g. start up the
        AWG so that it actually delivers said triggers.
        
        Note: If this is used instead of ArmAndWait, then WaitForAcquisition will
        need to be called at some point so that the card realizes that it is done.
        ReadData has been modified so that it should take care of this automatically,
        but FYI.
        
        '''
        try:
            self.InitiateAcquisition()
        except:
            if self.verbose:
                print('Initiate Failed. Trying to fix with a recalibrate.')
            self.SelfCalibrate()
            self.InitiateAcquisition()
    
    def ArmAndWait(self):
        '''Initialize the acquisition and get ready to take and read data.
        Tells the card to get ready and wait for triggers, and then waits 
        for it to finish.
        
        This version works fine if the AWG is already setup and running, but it will
        cause python to hang until the acqusition is finished. So, you will not be able
        to do anything, and that includes starting up the AWG to have it give the 
        card triggers. So, this function should only be used if the card is triggered
        off an external generator, not the AWG.
        '''
        self.Arm()
        self.WaitForAcquisition()
     
        
    def ConfigureAcquisition(self, samples, sampleRate, segments = 1):
        ''''Configure Acquistion function that gets everything ready for data taking.
#        Does more than the equivelent C function.
        
        Test version of configure acquisition. Trying to use the built in function of the C
        driver and not manually setting all the variables.'''
        
        self.sampleRate = sampleRate
        sampleRateC = ctypes.c_double(self.sampleRate)

        self.segments = segments
        numRecordsC = ctypes.c_int64(self.segments)
        

#        #software needs to know at this point so that it can adjust the number of samples
#        #I can't write to the driver here because old number of samples could still be
#        #large from a non-averaged run.
        
        #turn off averaging in hardware so that it doesn't muck with setting the number of samples.
        #This is necessary when switching back from averaging mode.
        self.SetDriverAttribute('averageMode', 0)
        self.SetDriverAttribute('averages', 1)
        #using this low-level hack isntead of ConfigureAveraging, because at this point, python needs
        #to know about averaging, but the hardware needs to always not be in averaging mode, so that 
        #it can take a large number of samples if the upcoming acqusition is going to be in regular mode.
        #whereas the software needs to adjust the limit on the number of samples if the upcoming acquisition
        # is going to be in averaging mode.
        
        
        if self.verbose:
            print('Setting acquisition parameters manually.')
            print('Sample will autoround to a multiple of 1024 and shorten if tto long for avg mode.')
        ###!!!!!!!!!!!!!
        #!!!!!!!!!empircally, record size must be 1024*n, and the largest that it can be is 513*1024
        #roughly 500 kS. Which is what the mnual says, but this sucks. Did the old Acqiris do this?
        #or is it just that the ME is only for non-averaging mode?
        multiple = int(numpy.ceil(self.samples/1024))
        if self.averages > 1:
            #check the number of channels
            numChans = len(self.activeChannels)
            if numChans == 1:
                maxSamples = 1024*1024
            else:
                maxSamples = 1024*512
                
            if multiple*1024 > maxSamples:
                print('Data is too long for averaging mode. Setting to max length: 512*1024')
                self.samples = int(maxSamples)
            else:
                #auto round to next multiple of 1024 up, regardless
                self.samples = int(multiple*1024) 
                #it is very important that this winds up an integer, or it least it was at some point
                
        numPointsPerRecordC = ctypes.c_int64(self.samples)
        #it looks like now that we have the order of operations right between configuring
        #the acquisition and turning averaging on and off, so that we can use the
        #built in C function for configuring, that we might not need this multiple of 1024 check.
        #But we do need the auto check that turns down the number of samples if its averaging mode
        
#        print('self.averageMode = ' + str(self.averageMode))
#        print('self.averages = ' + str(self.averages))
#        print('self.samples = ' + str(self.samples))
        
        self.call('ConfigureAcquisition', self.visession, numRecordsC, numPointsPerRecordC, sampleRateC)
        
        #configure the the actual averaging mode and restore number of averages
        self.ConfigureAveraging() #trip the flags to handle the averaging
        #(both in python and sending down to the hardware driver.)
        #Hopefully doing it after everything else is set will only try to flip to averager 
        #after the number of samples has been reduced.


    def ConfigureAveraging(self):
        if self.averages >1:
            self.averageMode = 1 #make sure this variable is conistent. I probably want to do away with it eventually
        else:
            self.averageMode = 0
#        self.SetDriverAttribute('averageMode', self.averageMode)
        self.driver.Acquisition.Mode = self.averageMode
#        self.SetDriverAttribute('averages', self.averages)
        self.driver.Acquisition.NumberOfAverages = self.averages

    def ConfigureChannel(self, channelNum = 1, Range = 2.5, offset = 0, enabled = True):
        if channelNum ==1:
            chanName = b'Channel1'
        elif channelNum ==2:
            chanName = b'Channel2'
        else:
            raise ValueError('Invalid channel. Should be 1 or 2')
        chanNameC = ctypes.c_char_p(chanName)
        
        if Range in [0.5, 2.5]:
            rangeC = ctypes.c_double(Range)
        else:
            raise ValueError('Range must be 0.5 or 2.5')
            
        self.channelRange = Range
        self.channelOffset = offset

        offsetC = ctypes.c_double(offset)

        couplingC = ctypes.c_int32(1) # 1 = DC coupling. Always. Card can't do AC.

        enabledC = ctypes.c_bool(enabled)

        self.call('ConfigureChannel', self.visession, chanNameC, rangeC, offsetC, couplingC, enabledC)
 
    
    def ConfigureTrigger(self, Source = 'External1', Level = 0, Slope = 'Falling', Mode = 'Edge'):
        if Mode != 'Edge':
            raise ValueError('This trigger mode not yet supported. Edge trigger only.')
        else:
            self.triggerMode = Mode
            self.triggerSource = Source
            self.triggerLevel = Level
            
            if Slope == 'Falling':
                self.triggerSlope = Slope
                triggerSlopeC = 0
            elif Slope == 'Rising':
                self.triggerSlope = Slope
                triggerSlopeC = 1
            else:
                raise ValueError("Edge trigger slope must be either 'Rising' or 'Falling'")
            
#            triggerSourceC = ctypes.create_string_buffer(self.triggerSource.encode('utf-8'))
#            triggerLevelC = ctypes.c_longdouble(self.triggerLevel)
            
#            self.SetDriverAttribute('triggerDelay', self.triggerDelay)
            self.driver.Trigger.Delay = self.triggerDelay
            



            #manually set trigger source because the configure function doesn't actually do it.
#            self.SetDriverAttribute('triggerSource', self.triggerSource)
            self.driver.Trigger.ActiveSource = self.triggerSource
            
#            self.SetDriverAttribute('triggerCoupling', self.triggerCoupling, recap = self.triggerSource)
            #!!!!!!!!haven't found this one
            #self..driver.Trigger.Sources[self.triggerSource].Coupling, but I think this always has to be 1
            self.driver.Trigger.Sources[self.triggerSource].Coupling = self.triggerCoupling
            
            #set trigger source before and after to make sure configure hits the right hardware channel
#            self.call('ConfigureEdgeTriggerSource', self.visession, triggerSourceC ,triggerLevelC, triggerSlopeC)
            activeTrigger = self.driver.Trigger.Sources[self.triggerSource]
            activeTrigger.Level = self.triggerLevel
            activeTrigger.Edge.Slope = triggerSlopeC 
            
            
            
    def _fillHardwareKeys(self):
        ''' Storing the names of the harware level vairables. '''
        
        self.driverKeys = []
        self.driverKeys.append('driverDescription')
        self.driverKeys.append('driverRevision')
        self.driverKeys.append('firmwareRevision')
        self.driverKeys.append('manufacturer')
        self.driverKeys.append('instrumentModel')
        self.driverKeys.append('serialNumber')
        
        
        
        
#        # Declaration of some attributes (values from AqMd3.h)
#        self.driverBases['IVI_ATTR_BASE']         = 1000000
#        self.driverBases['IVI_INHERENT_ATTR_BASE'] = self.driverBases['IVI_ATTR_BASE'] +  50000
#        self.driverBases['IVI_CLASS_ATTR_BASE']    = self.driverBases['IVI_ATTR_BASE'] +  250000
#        self.driverBases['IVI_LXISYNC_ATTR_BASE']  = self.driverBases['IVI_ATTR_BASE'] +  950000
#        self.driverBases['IVI_SPECIFIC_ATTR_BASE'] = self.driverBases['IVI_ATTR_BASE'] +  150000
#        
#        baseID = self.driverBases['IVI_ATTR_BASE']  
#        inherentID = self.driverBases['IVI_INHERENT_ATTR_BASE']
#        classID = self.driverBases['IVI_CLASS_ATTR_BASE'] 
#        lxisyncID = self.driverBases['IVI_LXISYNC_ATTR_BASE']
#        specificID = self.driverBases['IVI_SPECIFIC_ATTR_BASE']
            
#        #firmware
##        AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION          = IVI_INHERENT_ATTR_BASE + 514  # ViString, read-only 
#        self.hardwareIDs['driverDescription'] = [inherentID + 514, 'ViString']
##        AQMD3_ATTR_SPECIFIC_DRIVER_REVISION             = IVI_INHERENT_ATTR_BASE + 551  # ViString, read-only  
#        self.hardwareIDs['driverRevision'] = [inherentID + 551, 'ViString']
##        AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION         = IVI_INHERENT_ATTR_BASE + 510  # ViString, read-only
#        self.hardwareIDs['firmwareRevision'] = [inherentID + 510, 'ViString']
##        AQMD3_ATTR_INSTRUMENT_MANUFACTURER              = IVI_INHERENT_ATTR_BASE + 511  # ViString, read-only
#        self.hardwareIDs['manufacturer'] = [inherentID+511, 'ViString']
##        AQMD3_ATTR_INSTRUMENT_MODEL                     = IVI_INHERENT_ATTR_BASE + 512  # ViString, read-only
#        self.hardwareIDs['instrumentModel'] = [inherentID + 512, 'ViString']
##        AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING = IVI_SPECIFIC_ATTR_BASE + 8    # ViString, read-only 
#        self.hardwareIDs['serialNumber'] = [specificID +8, 'ViString']
        
        
        self.hardwareKeys = []
        self.hardwareKeys.append('channelOffset')
        self.hardwareKeys.append('channelRange')
        self.hardwareKeys.append('channelEnabled')
        
        self.hardwareKeys.append('triggerSource')
        self.hardwareKeys.append('triggerLevel')
        self.hardwareKeys.append('triggerMode')
        self.hardwareKeys.append('triggerDelay')
        self.hardwareKeys.append('triggerCoupling')
        
        self.hardwareKeys.append('averageMode')
        self.hardwareKeys.append('averages')
        
        self.hardwareKeys.append('samples')
        self.hardwareKeys.append('segments')
        
        self.hardwareKeys.append('sampleRate')
        
        self.hardwareKeys.append('simulateMode')
        
        
#        #channel
##        AQMD3_ATTR_VERTICAL_OFFSET                      = IVI_CLASS_ATTR_BASE + 25      # ViReal64, read-write 
##        self.hardwareIDs['verticalOffset'] = [classID + 25, 'ViReal64']
#        self.hardwareIDs['channelOffset'] = [classID + 25, 'ViReal64']
#        self.hardwareIDs['channelRange'] = [classID + 26, 'ViReal64']
#        self.hardwareIDs['channelEnabled'] = [classID + 2, 'ViBoolean']
#        
#        #trigger
##        AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE                = IVI_CLASS_ATTR_BASE + 1       # ViString, read-write
#        self.hardwareIDs['triggerSource'] = [classID+1, 'ViString']
##        self.hardwareIDs['triggerLevel'] = [classID+19, 'ViReal64']
#        self.hardwareIDs['triggerLevel'] = [classID+21, 'ViInt32']
#        self.hardwareIDs['triggerMode'] = [classID+23, 'ViInt32']
#        self.hardwareIDs['triggerDelay'] = [classID+17, 'ViReal64']
#        self.hardwareIDs['triggerCoupling'] = [classID+16, 'ViInt32']
#    
#        #averaging
##        AQMD3_ATTR_ACQUISITION_MODE                     = IVI_SPECIFIC_ATTR_BASE + 11   # ViInt32, read-write 
#        self.hardwareIDs['averageMode'] = [specificID + 11, 'ViInt32']
##        AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES       = IVI_SPECIFIC_ATTR_BASE + 69   # ViInt32, read-write 
#        self.hardwareIDs['averages'] = [specificID+69, 'ViInt32']
#    
#        #acquistion parameters
##        AQMD3_ATTR_RECORD_SIZE                          = IVI_CLASS_ATTR_BASE + 14      # ViTnt64, read-write
#        self.hardwareIDs['samples'] = [classID+14, 'ViInt64']
##        AQMD3_ATTR_NUM_RECORDS_TO_ACQUIRE               = IVI_CLASS_ATTR_BASE + 13      # ViInt64, read-write
#        self.hardwareIDs['segments'] = [classID+13, 'ViInt64']
#    
#        #sampling and clocking
##        AQMD3_ATTR_SAMPLE_RATE                          = IVI_CLASS_ATTR_BASE + 15       # ViReal64, read-write  
#        self.hardwareIDs['sampleRate'] = [classID+15, 'ViReal64']
#       
#        #simulation mode
##        AQMD3_ATTR_SIMULATE                             = IVI_INHERENT_ATTR_BASE + 5       # ViBoolean, read-write  
#        self.hardwareIDs['simulateMode'] = [inherentID + 5, 'ViBoolean']
        

    def InitiateAcquisition(self):
        self.call('InitiateAcquisition', self.visession)
        self.armed = True #throw flag so that python knows acquisition has been itiated
        #and system is waiting for trigger.
        
    def _loadDefaultConfig(self):
        #sampling
        self.samples = 1024
        self.sampleRate = 1*10**9
        
        #averaging
        self.averageMode = 0
        self.averages = 1
        
        #multisegment acquisition
        self.segments = 1
        
        #channel settings
        self.activeChannels = [1,2]
        self.channelRange = 2.5
        self.channelOffset = 0
        
        #trigger settings
        self.triggerSource = 'External1'
        self.triggerMode = 'Edge'
        self.triggerLevel = 0
        self.triggerSlope = 'Falling'
        self.triggerDelay = 0
        self.triggerCoupling = 1 #1 for DC, 0 for AC
        

        #diagnostic prints
        self.verbose = False
        
        #timeout of failed acquisitions
        self.timeout = 5  #seconds
        
    def ReadAllData(self, returnRaw = False):
        ''' Highest level read function. Will automatically read all active
        channels.
        
        Will trip flags when it is done that the data acquisition and read are complete
        '''
        #check if the card is done. 
        #NOTE: Card needs to be asked if its done, or it will get mad, even if it is actually done
        if self.armed:
            #there was an active acquisition
            if not self.acquisitionFinished:
                #need to confirm with the card that it is done.
                #(this will already have happened if you used ArmAndWait, but not if you use Arm)
                self.WaitForAcquisition()
        else:
            raise ValueError("Cannot read data. Card acquisition hasn't happened.")
        
        
        data1 = []
        data2 = []
        if len(self.activeChannels) > 1:
            data1 = self.ReadData(1, returnRaw = returnRaw)
            data2 = self.ReadData(2, returnRaw = returnRaw)
            #notify the python that this acquisition is done
            self.armed = False
            self.acquisitionFinished = False
            return data1, data2
        else:
            data =  self.ReadData(self.activeChannels[0], returnRaw = returnRaw)
            #notify the python that this acquisition and read are done
            self.armed = False
            self.acquisitionFinished = False
            return data

    def ReadData(self, chanNum, returnRaw = False):
        '''Single channel data read function.
        Will automatically call lowerlevel read functions of the driver,
        depending on how the card is set up.

        Somewhat intended as a hidden internal function of ReadAllData, so
        it may not do the safety checks as well as it's parent because they are 
        harder to manage when each channel is called separately.
        
        This single channel read function is not protected against reading the same dat twice,
        as least not at the Kollar-lab python level.
        
        '''
        #check if the card is done. (Ideally this happens in ReadAllData, but it channels are
        #being called individually, then it may need to happen here.)
        #NOTE: Card needs to be asked if its done, or it will get mad, even if it is actually done
        #WaitForAcquisition MUST happen before the hardware-level reads.
        if self.armed:
            #there was an active acquisition
            if not self.acquisitionFinished:
                #need to confirm with the card that it is done.
                #(this will already have happened if you used ArmAndWait, but not if you use Arm)
                self.WaitForAcquisition()
        else:
            raise ValueError("Cannot read data. Card acquisition hasn't happened.")
        
        
        if chanNum ==1:
            chanName = b'Channel1'
        elif chanNum ==2:
            chanName = b'Channel2'
        else:
            raise ValueError('Invalid channel. Should be 1 or 2')

        #memout = driverdll.AqMD3_QueryMinWaveformMemory(init_status.ViSession, c_int32(64), numRecordsC, c_int64(0), numPointsPerRecordC, byref(numSamples))
        numPointsPerRecordC = ctypes.c_int64(self.samples) 
        numRecordsC = ctypes.c_int64(self.segments)
        self.offsetWithinRecord = 0
#        self.offsetWithinRecord = 1
        offsetWithinRecordC = ctypes.c_int64(self.offsetWithinRecord) 
        numSamplesC = ctypes.c_int64()
        self.call('QueryMinWaveformMemory', self.visession, ctypes.c_int32(64), numRecordsC,\
             offsetWithinRecordC, numPointsPerRecordC, ctypes.byref(numSamplesC))
        if self.verbose:
            print('Memeory Allocation Determined')
            
#        print(chanNum)
#        print(offsetWithinRecordC.value)
            
        numSamples = int(numSamplesC.value)
        self.totalSamples = numSamples

        if self.verbose:
            print('Trying to Read')
        if self.averageMode:
            out = self._ReadAveragerData(chanName, returnRaw = returnRaw)
        else:
#            if self.multiseg:
            if self.segments > 1:
                out = self._ReadMultiSegData(chanName, returnRaw = returnRaw)
            else:
                out = self._ReadSingleSegData(chanName, returnRaw = returnRaw)
                

        return out

    def _ReadAveragerData(self,chanName, returnRaw = False):
        if self.verbose:
            print('reading averaged data')
        arraySize_int = int(self.totalSamples)
        segments_int = int(self.segments)

        firstRecordC = ctypes.c_int64(0)
        numRecordsC = ctypes.c_int64(self.segments)
        offsetWithinRecordC = ctypes.c_int64(self.offsetWithinRecord)
        numPointsPerRecordC = ctypes.c_int64(self.samples)

        WavefromArraySizeC = ctypes.c_int64(self.totalSamples)

        ActualAveragesC = ctypes.c_int32()
        ActualRecordsC = ctypes.c_int64()
        InitialXOffsetC = ctypes.c_longdouble()
        XIncrementC = ctypes.c_longdouble()

        # FlagsC = ctypes.c_int32*segments_int
        # ActualPointsC = ctypes.c_int64*segments_int
        # FirstValidPointsC = ctypes.c_int64*segments_int
        # InitialXTimeSecondsC = ctypes.c_longdouble*segments_int
        # InitialXTimeFractionC = ctypes.c_longdouble*segments_int

        # WaveformArrayC = ctypes.c_double*arraySize_int #this doesn't work
        class WAVEFORMHOLDER(ctypes.Structure):   
            _fields_ = [("WaveformArrayC", ctypes.c_longdouble*arraySize_int),
                        ("ActualPointsC", ctypes.c_int64*segments_int),
                        ("FirstValidPointsC", ctypes.c_int64*segments_int),
                        ("InitialXTimeSecondsC", ctypes.c_longdouble*segments_int),
                        ("InitialXTimeFractionC", ctypes.c_longdouble*segments_int),
                        ("FlagsC", ctypes.c_int64*segments_int)]
#        class WAVEFORMHOLDER(ctypes.Structure):   
#            _fields_ = [("WaveformArrayC", ctypes.c_int32*arraySize_int),
#                        ("ActualPointsC", ctypes.c_int64*segments_int),
#                        ("FirstValidPointsC", ctypes.c_int64*segments_int),
#                        ("InitialXTimeSecondsC", ctypes.c_longdouble*segments_int),
#                        ("InitialXTimeFractionC", ctypes.c_longdouble*segments_int),
#                        ("FlagsC", ctypes.c_int64*segments_int)]
        WaveHolder = WAVEFORMHOLDER()   # This DOES WORK! It can be sent in.
        #Successfully returns data to "WaveHolder.WaveformArrayC"

        

        chanNameC = ctypes.c_char_p(chanName)

        self.call('FetchAccumulatedWaveformReal64', self.visession, chanNameC, \
                    firstRecordC,\
                    numRecordsC,\
                    offsetWithinRecordC,\
                    numPointsPerRecordC,\
                    WavefromArraySizeC,
                    WaveHolder.WaveformArrayC,\
                    ctypes.byref(ActualAveragesC),\
                    ctypes.byref(ActualRecordsC),\
                    WaveHolder.ActualPointsC,\
                    WaveHolder.FirstValidPointsC,\
                    ctypes.byref(InitialXOffsetC),\
                    WaveHolder.InitialXTimeSecondsC,\
                    WaveHolder.InitialXTimeFractionC,\
                    ctypes.byref(XIncrementC),\
                    WaveHolder.FlagsC)
#        self.call('FetchAccumulatedWaveformInt32', self.visession, chanNameC, \
#                    firstRecordC,\
#                    numRecordsC,\
#                    offsetWithinRecordC,\
#                    numPointsPerRecordC,\
#                    WavefromArraySizeC,
#                    WaveHolder.WaveformArrayC,\
#                    ctypes.byref(ActualAveragesC),\
#                    ctypes.byref(ActualRecordsC),\
#                    WaveHolder.ActualPointsC,\
#                    WaveHolder.FirstValidPointsC,\
#                    ctypes.byref(InitialXOffsetC),\
#                    WaveHolder.InitialXTimeSecondsC,\
#                    WaveHolder.InitialXTimeFractionC,\
#                    ctypes.byref(XIncrementC),\
#                    WaveHolder.FlagsC)
        
        if self.verbose:
            print('Fetch Complete. Processing Data.')
        rawData = numpy.asarray(WaveHolder.WaveformArrayC)
#        rawData = numpy.asarray(WaveHolder.WaveformArrayC)
        dataRawSize = rawData.size
        dataActualSegments = int(ActualRecordsC.value)
        dataActualPoints_full = numpy.asarray(WaveHolder.ActualPointsC)
        dataActualPoints = dataActualPoints_full[0]
        dataFirstValidPoints = numpy.asarray(WaveHolder.FirstValidPointsC).astype('int64')
#        print(dataFirstValidPoints)
        if returnRaw:
            out = [rawData, dataActualPoints, dataFirstValidPoints, dataActualSegments]
            return out
        else:
            if dataActualPoints != self.samples:
                print("Warning. Data size doesn't match the number of samples. Something wierd happened.")

            if dataActualSegments == 1:
                startInd = dataFirstValidPoints[0]
                data = rawData[startInd:(startInd+dataActualPoints)]
            else:
                data = numpy.zeros((dataActualSegments,dataActualPoints))
                for segind in range(0,dataActualSegments ):
                    startInd = dataFirstValidPoints[segind]
                    
                    data[segind,:] = rawData[startInd:(startInd+dataActualPoints)]
            if self.verbose:
                print('Data Processed.')
            return data

    def _ReadMultiSegData(self,chanName, returnRaw = False):
        if self.verbose:
            print('reading non-averaged, multisegment data')
        arraySize_int = int(self.totalSamples)
        segments_int = int(self.segments)

        ActualRecordsC = ctypes.c_int64()
        ActualPointsC = ctypes.c_int64*segments_int
        FirstValidPointsC = ctypes.c_int64*segments_int
        InitialXOffsetC = ctypes.c_longdouble*segments_int
        InitialXTimeSecondsC = ctypes.c_longdouble*segments_int
        InitialXTimeFractionC = ctypes.c_longdouble*segments_int
        XIncrementC = ctypes.c_longdouble()

        # WaveformArrayC = ctypes.c_double*arraySize_int #this doesn't work
        class WAVEFORMHOLDER(ctypes.Structure):   
            _fields_ = [("WaveformArrayC", ctypes.c_longdouble*arraySize_int),
                        ("ActualPointsC", ctypes.c_int64*segments_int),
                        ("FirstValidPointsC", ctypes.c_int64*segments_int),
                        ("InitialXOffsetC", ctypes.c_longdouble*segments_int),
                        ("InitialXTimeSecondsC", ctypes.c_longdouble*segments_int),
                        ("InitialXTimeFractionC", ctypes.c_longdouble*segments_int)]
        WaveHolder = WAVEFORMHOLDER()   # This DOES WORK! It can be sent in.
        #Successfully returns data to "WaveHolder.WaveformArrayC"
        
        WavefromArraySizeC = ctypes.c_int64(self.totalSamples)

        chanNameC = ctypes.c_char_p(chanName)

        firstRecordC = ctypes.c_int64(0)
        offsetWithinRecordC = ctypes.c_int64(self.offsetWithinRecord)
        numRecordsC = ctypes.c_int64(self.segments)
        numPointsPerRecordC = ctypes.c_int64(self.samples)

        self.call('FetchMultiRecordWaveformReal64', self.visession, chanNameC,\
            firstRecordC,\
                numRecordsC,\
                    offsetWithinRecordC,\
                        numPointsPerRecordC,\
                            WavefromArraySizeC,\
                                 WaveHolder.WaveformArrayC,\
                                 ctypes.byref(ActualRecordsC),\
                                 WaveHolder.ActualPointsC,\
                                 WaveHolder.FirstValidPointsC,\
                                 WaveHolder.InitialXOffsetC,\
                                 WaveHolder.InitialXTimeSecondsC,\
                                 WaveHolder.InitialXTimeFractionC,\
                                 ctypes.byref(XIncrementC))
        if self.verbose:
            print('Fetch Complete. Processing Data.')
        rawData = numpy.asarray(WaveHolder.WaveformArrayC)
        dataRawSize = rawData.size
        dataActualSegments = int(ActualRecordsC.value)
        dataActualPoints_full = numpy.asarray(WaveHolder.ActualPointsC)
        dataActualPoints = dataActualPoints_full[0]
        dataFirstValidPoints = numpy.asarray(WaveHolder.FirstValidPointsC).astype('int64')
        if returnRaw:
            out = [rawData, dataActualPoints, dataFirstValidPoints, dataActualSegments]
            return out
        else:
            if dataActualPoints != self.samples:
                print("Warning. Data size doesn't match the number of samples. Something wierd happened.")

            data = numpy.zeros((dataActualSegments,dataActualPoints))
            for segind in range(0,dataActualSegments ):
                startInd = dataFirstValidPoints[segind]
                
                data[segind,:] = rawData[startInd:(startInd+dataActualPoints)]
            return data

    def _ReadSingleSegData(self,chanName, returnRaw = False):
        if self.verbose:
            print('reading non-averaged, single-segment data')
        arraySize_int = int(self.totalSamples)

        ActualPointsC = ctypes.c_int64()
        FirstValidPointC = ctypes.c_int64()
        InitialXOffsetC = ctypes.c_longdouble()
        InitialXTimeSecondsC = ctypes.c_longdouble()
        InitialXTimeFractionC = ctypes.c_longdouble()
        XIncrementC = ctypes.c_longdouble()

        # WaveformArrayC = ctypes.c_double*arraySize_int #this doesn't work
        class WAVEFORMHOLDER(ctypes.Structure):   
            _fields_ = [("WaveformArrayC", ctypes.c_longdouble*arraySize_int)]
        WaveHolder = WAVEFORMHOLDER()   # This DOES WORK! It can be sent in.
        #Successfully returns data to "WaveHolder.WaveformArrayC"

        WavefromArraySizeC = ctypes.c_int64(self.totalSamples)

        chanNameC = ctypes.c_char_p(chanName)

        self.call('FetchWaveformReal64', self.visession, chanNameC, WavefromArraySizeC,\
             WaveHolder.WaveformArrayC,\
            ctypes.byref(ActualPointsC),\
                 ctypes.byref(FirstValidPointC),\
                      ctypes.byref(InitialXOffsetC),\
                 ctypes.byref(InitialXTimeSecondsC),\
                      ctypes.byref(InitialXTimeFractionC),\
                       ctypes.byref(XIncrementC))
        
        if self.verbose:
            print('Fetch Complete. Processing Data.')
        rawData = numpy.asarray(WaveHolder.WaveformArrayC)
        dataRawSize = rawData.size
        dataActualPoints = int(ActualPointsC.value)
        dataFirstValidPoint = int(FirstValidPointC.value)
        if returnRaw:
            out = [rawData, dataActualPoints, dataFirstValidPoint]
            return out
        else:
            if dataActualPoints != self.samples:
                print("Warning. Data size doesn't match the number of samples. Something wierd happened.")

            startInd = dataFirstValidPoint
            data = rawData[startInd: (startInd + dataActualPoints)]
            return data
        
    def ReInitialize(self):
        '''Basic init function. Can also be used to wipe he settings if the card is very confused. 
        Will push the currently stored settings to the card and calibrate'''
        print('Initializing...')
        if self.simulate == True:
            strInitOptions =  'Simulate=True,  DriverSetup= model = SA220P'
        else:
            strInitOptions = 'Simulate=False,  DriverSetup= model = SA220P'
        
        self.driver =  AqMD3( self.hardwareAddress , False, False, strInitOptions)
        
##        with AqMD3( resourceString, False, False, initOptions ) as driver:
#        self.call('InitWithOptions',
#            ctypes.create_string_buffer(self.hardwareAddress.encode('utf-8')),
#            False,
#            False, 
#            strInitOptionsC,
#            ctypes.byref(self.visession))
        
       
        #push the currently stored settings to the card and calibrate. Hopefully ths will actually save hassle.
        #if this is running at _init_ then it will push the defaults
#        self.SetParams()  #!!!!!!!!!!!!!!!
        self.SelfCalibrate()

    def SetAverageMode(self, averages):
        ''' Intermediate function to switch the card over to averaging mode.
        It manually writes the necessary fields. 
        I expect that it will fall out of use ang instead be done in the general 
        SetParams'''
#        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 1, 'ViInt32') 
        self.SetDriverAttribute('averageMode', 1)
#        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, averages, 'ViInt32') 
        self.SetDriverAttribute('averages', averages)
        self.averageMode = True
#        self.multiseg = False

    def SetSingleMode(self):
        ''' Intermediate function to switch the card over to single acquisition mode.
        It manually writes the necessary fields. 
        I expect that it will fall out of use ang instead be done in the general 
        SetParams'''
#        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 0, 'ViInt32') 
        self.SetDriverAttribute('averageMode', 0)
#        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, 1, 'ViInt32') 
        self.SetDriverAttribute('averages', 1)
        self.averageMode = False
#        self.multiseg = False
        self.segments = 1

    def SetMultiMode(self):
        ''' Intermediate function to switch the card over to multisegment, non-avergaing mode.
        It manually writes the necessary fields. 
        I expect that it will fall out of use ang instead be done in the general 
        SetParams'''
#        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 0, 'ViInt32') 
        self.SetDriverAttribute('averageMode', 0)
#        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, 1, 'ViInt32') 
        self.SetDriverAttribute('averages', 1)
        self.averageMode = False
#        self.multiseg = True   #segment number will be configured elsewhere
        self.segments = 2

    def WaitForAcquisition(self):
        '''Timeout in seconds '''
        if self.verbose:
            print('Waiting until the acquisition is done (or checking that it is) and getting confirmation')
#        self.timeout = timeout
        timeoutC = ctypes.c_int32(self.timeout*1000)
        self.call('WaitForAcquisitionComplete', self.visession, timeoutC)
        self.acquisitionFinished = True #throw the internal flag to show that the card is done.



if __name__ == '__main__':

    hardwareAddress = "PXI23::0::0::INSTR"

    IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
    sys.path.append(IVIbinPath)
    
#    alwaysInit = False
    alwaysInit = True
    if alwaysInit:
        print("Initialize instrument @", hardwareAddress)
        card = Acqiris(hardwareAddress)
    else:
        try:
            #random driver call to see if things are set up
#            card.GetAttribute(None, AQMD3_ATTR_SIMULATE, 'ViBoolean')
            card.GetDriverAttribute('simulateMode')
            print('Card already initialized')
        except:
            print("Initialize instrument @", hardwareAddress)
            card = Acqiris(hardwareAddress)
    

#    #####################
#    #schose acquisition type
#    #####################
#
#    averageMode = True
##    averageMode = False
#
##    multisegMode = True
#    multisegMode = False
#    
#    
#    
#    ##########
#    #set settings to the card object
#    ### samples = 1*10**7
#    card.samples = 513*1024
#    card.sampleRate = 1*10**9
#    
#    
#    #swtich between test cases
#    if averageMode:
#        card.averageMode = 1
#        if multisegMode:
#            card.averages = 100
#            card.segments = 2
#        else:
#            card.averages = 100
#            card.segments = 1
#    else:
#        card.averageMode = 0
#        if multisegMode:
#            card.averages = 1
#            card.segments = 2
#        else:
#            card.averages = 1
#            card.segments = 1
#            
#    card.triggerSource = 'External1'
#    card.triggerLevel = 0
#    card.triggerSlope = 'Falling'
#    
#    card.activeChannels = [1,2]
##    card.activeChannels = [1] #it looks like single channel has to be channel 1
#    
#    card.timeOut = 2
#    
#    card.verbose = True #get lots of diagnostic print statements
#    
#    ##########
#    #push settings to hardware
#    card.SetParams()
##    card.ConfigureChannel(channelNum = 1, Range = 2.5, offset = 0, enabled = True)
##    card.ConfigureChannel(channelNum = 2, Range = 2.5, offset = 0, enabled = True)
##    card.ConfigureTrigger('External1', 0, 'Falling')
##    card.ConfigureAcquisition(samples, sampleRate, segments)
#
#    card.GetParams()
#    
#    ##########
#    #initiate acquisition and wait for it to finish
#    card.ArmAndWait()
##    card.Arm()
##    card.WaitForAcquisition()
##    card.WaitForAcquisition(timeOut)
#
#    print('Data Acquired (In Theory)')
#
#    
##    data = card.ReadData(1, returnRaw = False) #read channel 1
#    data = card.ReadData(2, returnRaw = False) #read channel 1
##    data, data2 = card.ReadAllData()
#    
#
#    ####
#    #plot data
#
#    numSamples = card.samples
##    sampleRate = card.GetDriverAttribute('sampleRate')
#    sampleRate = card.sampleRate
#    dt = 1/sampleRate
#    xaxis = scipy.arange(0, numSamples,1)*dt
#
#    if multisegMode:
#        dataSegments = data.shape[0]
#
#    print('Plotting.')
#
#    fig1 = pylab.figure(1)
#    pylab.clf()
#    ax = pylab.subplot(1,1,1)
#    if multisegMode:
#        pylab.title('Multiseg: Active Channel')
#        for segind in  range(0,dataSegments ):
#            labelstr = 'seg: ' + str(segind)
#            pylab.plot(xaxis*1000, data[segind,:] + segind*0.125, label = labelstr)
#    else:
#        if averageMode:
#            pylab.title('Averager: Active Channel')
#        else:
#            pylab.title('Singleseg: Active Channel')
#        labelstr = 'seg: 0' 
#        pylab.plot(xaxis*1000, data[:], label = labelstr)
#    ax.legend(loc = 'upper right')
#    pylab.ylabel('Voltage (waterfall)')
#    pylab.xlabel('Time (ms)')
#
#    pylab.show()
#
#
#    print('Done plotting.')
#    
##    # Close the instrument
##    print("\nClose the instrument")    
##    card.close()
#    
    

    
    
    
    
    
    
    
    
    
    ###################
    #older, cruder codes
    ################
    #    #    averageMode = True
#    averageMode = False
#
#    # multisegMode = True
#    multisegMode = False
#    
#    # samples = 1*10**7
#    samples = 513*1024
#    sampleRate = 1*10**9
    

#    if averageMode:
#        numAverages = 100
#        segments = 1
#        card.SetAverageMode(numAverages)
#    else:
#        if multisegMode:
#            numAverages = 1
#            segments = 2
#            card.SetMultiMode()
#        else:
#            numAverages = 1
#            segments = 1
#            card.SetSingleMode()
#
#
#    #configure channels
#    card.activeChannels = [1,2]
#    card.ConfigureChannel(channelNum = 1, Range = 2.5, offset = 0, enabled = True)
#    card.ConfigureChannel(channelNum = 2, Range = 2.5, offset = 0, enabled = True)
#
#    # Set active trigger source to External1
#    print('\nSet active trigger source to External1')
#    card.ConfigureTrigger('External1', 0, 'Falling')
#
#    # aconfigout = driverdll.AqMD3_ConfigureAcquisition(init_status.ViSession, numRecordsC, numPointsPerRecordC, sampleRateC)
#    card.ConfigureAcquisition(samples, sampleRate, segments)
    
#    # acqout = driverdll.AqMD3_InitiateAcquisition(init_status.ViSession)
#    # waitout = driverdll.AqMD3_WaitForAcquisitionComplete(init_status.ViSession,c_int32(timeouttime) )
#    try:
#        card.InitiateAcquisition()
#    except:
#        print('Initiate Failed. Trying a Calibration and Reinnitiate')
#        card.SelfCalibrate()
#        card.InitiateAcquisition()
#
#    timeOut = 10
#    card.WaitForAcquisition(timeOut)
#
#    print('Data Acquired (In Theory)')
#
#
#    # out = card.ReadData(1, returnRaw = True)
#    data = card.ReadData(1, returnRaw = False)
#
#
#    ####
#    #plot data
#
#    numSamples = card.samples
#    #sampleRate = card.GetAttribute('Channel1', AQMD3_ATTR_SAMPLE_RATE, 'ViReal64')
#    sampleRate = card.GetDriverAttribute('sampleRate')
#    dt = 1/sampleRate
#    xaxis = scipy.arange(0, numSamples,1)*dt
#
#    if multisegMode:
#        dataSegments = data.shape[0]
#
#    print('Plotting.')
#
#    fig1 = pylab.figure(1)
#    pylab.clf()
#    ax = pylab.subplot(1,1,1)
#    if multisegMode:
#        pylab.title('Multiseg: Active Channel')
#        for segind in  range(0,dataSegments ):
#            labelstr = 'seg: ' + str(segind)
#            pylab.plot(xaxis*1000, data[segind,:] + segind*0.125, label = labelstr)
#    else:
#        if averageMode:
#            pylab.title('Averager: Active Channel')
#        else:
#            pylab.title('Singleseg: Active Channel')
#        labelstr = 'seg: 0' 
#        pylab.plot(xaxis*1000, data[:], label = labelstr)
#    ax.legend(loc = 'upper right')
#    pylab.ylabel('Voltage (waterfall)')
#    pylab.xlabel('Time (ms)')
#
#    pylab.show()
#
#
#    print('Done plotting.')
#    
##    # Close the instrument
##    print("\nClose the instrument")    
##    card.close()    
    
    
    
    
    