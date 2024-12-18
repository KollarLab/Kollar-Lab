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
import matplotlib.pyplot as plt
#import tarfile
#import struct
#import glob
import numpy
import time

#import pickle
#import datetime
#import itertools
import sys
import warnings


#import AqMD3 as driver #this is the Acqiris python driver
import AqMD3

InstrumentFolderPath = r'C:\Users\Kollarlab\Desktop\Kollar-Lab\Control'
if not InstrumentFolderPath in sys.path:
    sys.path.append(InstrumentFolderPath)


from userfuncs import freeze

@freeze
class Acqiris(object):
    
    def __init__(self, ResourceName, simulate = False):  
        self.instrument_type = 'Acqiris'
        
        #look to the right place.
        IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
        
        if not IVIbinPath in sys.path:
            sys.path.append(IVIbinPath)
        
        self.settingsCurrent = False
            
        
        self._fillHardwareKeys()
        
        
        self.simulate = simulate #python for whether we want to simulate or not
        #perperty simulateMode will come from the driver.
        
        self.hardwareAddress = ResourceName
#        self.ReInitialize()
        
        #in the python version, there seems to be an issue with reinitializing an already
        #existing instnace of the driver. It causes an error that the sample clocks are unstable
        #or something about mising TDC curves. Running close before reinitializing it seems to allow
        #things to work. So, putting this try in here in case the object card already exists.
        try:
            self.Close()
        except:
            pass
        self.InitializeDriver()
        
        #set default values. Eventually this should probably load from some default config file
        self._loadDefaultConfig()
        
        #push the currently stored settings to the card and calibrate. Hopefully ths will actually save hassle.
        #if this is running at _init_ then it will push the defaults
        self.SetParams()
        try:
            self.SelfCalibrate()
        except:
            #most likely cause of failure at this point is either that funny
            #pll synching issue, or the absence of an external clock.
            if self._clockSource == 'External':
                warnings.warn('External Clock Failure. Booting in Internal Mode.', RuntimeWarning )
                self.clockSource = 'Internal'
                self.SelfCalibrate()
            else:
                warnings.warn('Self Calibrate Failed: Cause Unkown.', RuntimeWarning )
        
        
        #some flags for controlling the read
        self.armed = False
        self.acquisitionFinished = False
        
#        self._fill_slots()
        #declare some attributes that will be needed later
        #because class is prozen they need to exist from the get go.
        self.offsetWithinRecord = 0
        self.totalSamples = self.samples
        self._simulateMode = self.simulate 
        temp  = self.GetParams()
        temp['timeout'] = self.timeout
        temp['verbose'] = self.verbose
        self.settings = temp
        
#    ##############################
    #hardware properties
    
    
    #master hardware settings parameter, to be used for
    #saving, and printing
    #includes what yout would need to reset to a particular configuration
    #but not changeable things like run time flags or hardware address.
    @property
    def settings(self):
        ''' getts all the settings needed to save or initialize a python card
        object. So, it calls GetParams to get the hardware level settings. 
        It adds a couple more like timeout and verbose that don't exist at the lowest level.
        
        For looking at amtching function names, this makes something more like a configuration
        getten from LoadConfig or _laodDefaultConfig
        
        '''
        paramsDict = self.GetParams()
        paramsDict['timeout'] = self.timeout #timeout is an argument 
        #of a driver level function, not a driver level parameter, so 
        #it is not gotten in get params
        paramsDict['verbose'] = self.verbose #added this so that this
        #has exactly the same fields as the default config.
        #Python object can therfore boot from this drictionary.
        return paramsDict 
    @settings.setter
    def settings(self,val):
        self.LoadConfig(val) #load config checks for invalid fields
    
    #icky properties that have to be batch set
    @property
    def sampleRate(self): #this one might not need to be here, but I did it already
        val = self.driver.Acquisition.SampleRate
        self._sampleRate = val
        return val
    @sampleRate.setter
    def sampleRate(self,val):
        self.settingsCurrent = False
        self._sampleRate = val 
    
    @property
    def samples(self):
        val = self.driver.Acquisition.RecordSize
        self._samples = val
        return val
    @samples.setter
    def samples(self,val):
        self.settingsCurrent = False
        self._samples = int(val) 

    @property
    def segments(self):
        val = self.driver.Acquisition.NumberOfRecordsToAcquire
        self._segments = val
        return val
    @segments.setter
    def segments(self,val):
        self.settingsCurrent = False
        self._segments = int(val) 
    
    @property
    def averageMode(self):
        val = self.driver.Acquisition.Mode
        self._averageMode = val
        return val
    @averageMode.setter
    def averageMode(self,val):
        self.settingsCurrent = False
        self._averageMode = val 

    @property
    def averages(self):
        val = self.driver.Acquisition.NumberOfAverages
        self._averages = val
        return val
    @averages.setter
    def averages(self,val):
        self.settingsCurrent = False
        self._averages = int(val) 
        if val > 1:
            self.averageMode = 1
        else:
            self.averageMode = 0
            
            
    #I hope that this one can be a regular property, but I'm not sure
    #might run into trouble with limits in avergaing mode
    @property
    def activeChannels(self):
        #have to check what it enabled and convert
        val = self.driver.Channels[0].Enabled
        val2 = self.driver.Channels[1].Enabled
        flags = numpy.asarray([val,val2])
        inds = numpy.where(flags == True)[0]
        full = numpy.asarray([1,2])
        actives = full[inds]
        self._activeChannels = actives
        return actives
    @activeChannels.setter
    def activeChannels(self,val):
        if (type(val) == int) or (type(val) == float):
            val = [val]
        if not val[0] == 1:
            raise ValueError('Active channels needs to be either 1, [1], or [1,2].')
        if len(val) > 1:
            if  not val[1] == 2:
                raise ValueError('Active channels needs to be either 1, [1], or [1,2].')
            if len(val)>2:
                raise ValueError('Only two channels')
        self.settingsCurrent = False
        self._activeChannels = val


    
    
    #normal properties that I think can be set reasonably.
    @property
    def triggerSource(self):
        val = self.driver.Trigger.ActiveSource
        self._triggerSource = val
        return val
    @triggerSource.setter
    def triggerSource(self,val):
        self._triggerSource = val
        self.driver.Trigger.ActiveSource = val 

    @property
    def triggerLevel(self):
        activeTrigger = self.driver.Trigger.Sources[self.triggerSource]
        val = activeTrigger.Level
        self._triggerLevel = val
        return val
    @triggerLevel.setter
    def triggerLevel(self,val):
        activeTrigger = self.driver.Trigger.Sources[self.triggerSource]
        self._triggerLevel = val
        activeTrigger.Level = val 
    
    @property
    def triggerSlope(self):
        activeTrigger = self.driver.Trigger.Sources[self.triggerSource]
        temp = activeTrigger.Edge.Slope
#        self._triggerSlope = temp
        if temp == 0:
            val = 'Falling'
        elif temp == 1:
            val = 'Rising'
        else:
            raise ValueError("Unkown trigger slope returned")
        self._triggerSlope = val
        return val
    @triggerSlope.setter
    def triggerSlope(self,val):
        if val == 'Falling':
            triggerSlopeC = 0
        elif val == 'Rising':
            triggerSlopeC = 1
        else:
            raise ValueError("Edge trigger slope must be either 'Rising' or 'Falling'")
#        self._triggerSlope = triggerSlopeC
        self._triggerSlope = val
        activeTrigger = self.driver.Trigger.Sources[self.triggerSource]
        activeTrigger.Edge.Slope = triggerSlopeC 
        
    @property
    def triggerDelay(self):
        val = self.driver.Trigger.Delay
        self._triggerDelay = val
        return val
    @triggerDelay.setter
    def triggerDelay(self,val):
        self._triggerDelay = val
        self.driver.Trigger.Delay = val 
        
    @property
    def triggerCoupling(self):
        activeTrigger = self.driver.Trigger.Sources[self.triggerSource]
        val = activeTrigger.Coupling
        self._triggerCoupling = val
        return val
    @triggerCoupling.setter
    def triggerCoupling(self,val):
         activeTrigger = self.driver.Trigger.Sources[self.triggerSource]
         self._triggerCoupling = val
         activeTrigger.Coupling = val
         
         
    
#    @property
#    def channelRange(self):
#        #have to check what it enabled and convert
#        val = self.driver.Channels[0].Range
#        val2 = self.driver.Channels[1].Range
#        if not val == val2:
#            print('Warning: Range set differently on the two channels')
#        return val
#    @channelRange.setter
#    def channelRange(self,val):
#        if val in [0.5, 2.5]:
#            self.driver.Channels[0].Range = val
#            self.driver.Channels[1].Range = val
#        else:
#            raise ValueError('Range must be 0.5 or 2.5')
    
#    @property
#    def channelOffset(self):
#        #have to check what it enabled and convert
#        val = self.driver.Channels[0].Offset
#        val2 = self.driver.Channels[1].Offset
#        if not val == val2:
#            print('Warning: Offset set differently on the two channels')
#        return val
#    @channelOffset.setter
#    def channelOffset(self,val):
#        self.driver.Channels[0].Offset = val
#        self.driver.Channels[1].Offset = val
         
    @property
    def channelRange(self):
        #have to check what it enabled and convert
        val = self.driver.Channels[0].Range
        val2 = self.driver.Channels[1].Range
        self._channelRange = numpy.asarray([val, val2]) #store the basic range
        #for both channels under the hood
        if val == val2:
            #two channels set the same, use one only
            return val
        else:
            return numpy.asarray([val, val2])
    @channelRange.setter
    def channelRange(self,val):
        if type(val) == float:
            #only single value specified
            self._channelRange = numpy.asarray([val, val])
            if val in [0.5, 2.5]:
                self.driver.Channels[0].Range = val
                self.driver.Channels[1].Range = val
            else:
                raise ValueError('Range must be 0.5 or 2.5')   
        elif (type(val) == float) or (type(val) == numpy.ndarray) :
            #given an array or list of values
            if len(val) == 2:
                if (val[0] in [0.5, 2.5]) and (val[1] in [0.5, 2.5]):
                    self.driver.Channels[0].Range = val[0]
                    self.driver.Channels[1].Range = val[1]
                else:
                    raise ValueError('Range must be 0.5 or 2.5')  
            else:
                raise ValueError('Array too long. Only two channels.') 

    @property
    def channelOffset(self):
        #have to check what it enabled and convert
        val = self.driver.Channels[0].Offset
        val2 = self.driver.Channels[1].Offset
        self._channelOffset = numpy.asarray([val,val2])
        if not val == val2:
            print('Warning: Offset set differently on the two channels')
        return val
    @channelOffset.setter
    def channelOffset(self,val):
        self.driver.Channels[0].Offset = val
        self.driver.Channels[1].Offset = val
        self._channelOffset = numpy.asarray([val,val])
        
        
        
    @property
    def simulateMode(self):
        val = self.driver.DriverOperation.Simulate
        self._simulateMode = val
        return val
    @simulateMode.setter
    def simulateMode(self,val):  #be very wary just trrying to set this blindly.
        self.driver.DriverOperation.Simulate = val 
        self._simulateMode = val
        

    @property
    def clockSource(self):
        temp = self.driver.ReferenceOscillator.Source
        if temp == 0:
            val = 'Internal'
        elif temp == 1:
            val = 'External'
        else:
            raise ValueError("Unkown clock source returned")
        self._clockSource = val
        return val
    @clockSource.setter
    def clockSource(self,val):
        if val == 'Internal':
            sourceC = 0
        elif val == 'External':
            sourceC = 1
        else:
            raise ValueError("Edge clock source. Must be either 'Internal' or 'External'")
        self._clockSource = val
        self.driver.ReferenceOscillator.Source = sourceC
        
    @property
    def clockFrequency(self):
        val = self.driver.ReferenceOscillator.ExternalFrequency
        if self.verbose:
            if self.clockSource == 'Internal':
                print('Internal Clock. Frequency Setting Ignored (I think).')
            if not val == 10**7:
                print('Warning: Clock frequency is not set to 10 MHz.')
        self._clockFrequency = val
        return val
    @clockFrequency.setter
    def clockFrequency(self,val):  #be very wary just trrying to set this blindly.
        if self.verbose:
            if self.clockSource == 'Internal':
                print('Internal Clock. Frequency Setting Ignored (I think).')
        if not val == 10**7:
            raise ValueError('Clock frequency must be 10 MHz?')
        self._clockFrequency = val
        self.driver.ReferenceOscillator.ExternalFrequency = val 


    
#    
#    #############################
#        
            
    def SetParams(self, params = {}):
        '''Autoinititalize function that will set up the card with all the driver
        parameters stored in this object. It will push all the settings and do cross checks.
        
        This function loads all the settings to the card, so it sets up the trigger and the
        data acquisition channels. Then it calls Configure Acquisition to set up things
        like averaging and sample/sample rate
        
        params is an optional settings dictionary that can be used set everything at the python level'''
        
        fields = params.keys()
        if not len(fields) == 0:
            self.LoadConfig(params)
        
        #configure the trigger
#        self.ConfigureTrigger(Source = self.triggerSource, Level = self.triggerLevel, Slope = self.triggerSlope, Mode = self.triggerMode)
        self.ConfigureTrigger(Source = self._triggerSource, Level = self._triggerLevel, Slope = self._triggerSlope, Mode = self.triggerMode)
        
        #configure the channels:
        for chind in [1,2]:
            if chind in self._activeChannels: #need the underscore here so it takes the python value and then sets up hardware to match
                enabled = True
            else:
                enabled= False
            self.ConfigureChannel(channelNum = chind, Range = self._channelRange, offset = self.channelOffset, enabled = enabled)
            
        #configure the acquisition
#        self.ConfigureAcquisition(self.samples, self.sampleRate, self.segments)
        self.ConfigureAcquisition(self._samples, self._sampleRate, self._segments)
        
        self.settingsCurrent = True
        
    def GetParams(self):
        '''Autoget function for params. Will querry the hardware and package all the settings, as well as 
        update the fields of the python software object.
        
        In order to fully initialize or save a python Acqiris object, you will need the @settings
        property. This incorporates the hardware-level settings and a couple more. (but still no
        run time flags or things like that.)'''
        
        hardwareSettings= {}
        
        if self.verbose:
            print('    ')
        
#        hardwareSettings = {}
        for key in self.hardwareKeys:
            
            #If these are all properties, then I can do this like this
            val = getattr(self,key)
            
            #store everything away
            setattr(self, key, val)
            hardwareSettings[key] = val
            
            #print all the results    
            if self.verbose:
                print(key + ' : ' + str(val))
   
        if self.verbose:
            print('    ')
            
#        return
        return hardwareSettings
            
            
    def Close(self):    
        self.driver.Close()
        
    def SelfCalibrate(self):
        print('Calibrating...')
        self.driver.Calibration.SelfCalibrate()

    ##################################################
    
    def Abort(self):
        '''Abort a hung acquisition '''
        self.driver.Acquisition.Abort()
    
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
        #check if the settings are up to date
        if not self.settingsCurrent: 
            print('Warning: Python contains stale settings.')
            print('Pushing current settings.')
            self.SetParams()
        
        if self.driver.Calibration.IsRequired:
             if self.verbose:
                print('Calibration needed before acquisition. Doing it.')
             self.SelfCalibrate()
        self.InitiateAcquisition()
        
#        time.sleep(0.1) #putting in a pause to see if that stops the clocks from being unhappy. Nope.
#        try:
#            self.InitiateAcquisition()
#        except:
#            time.sleep(5)
#            print('Waiting for clocks to stabilize?') #hopefully this fixes stuff
#            self.InitiateAcquisition()
            
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

        self.segments = segments
        

#        #software needs to know at this point so that it can adjust the number of samples
#        #I can't write to the driver here because old number of samples could still be
#        #large from a non-averaged run.
        
        #turn off averaging in hardware so that it doesn't muck with setting the number of samples.
        #This is necessary when switching back from averaging mode.
        self.driver.Acquisition.Mode = 0
        self.driver.Acquisition.NumberOfAverages = 1
        #using this low-level hack isntead of ConfigureAveraging, because at this point, python needs
        #to know about averaging, but the hardware needs to always not be in averaging mode, so that 
        #it can take a large number of samples if the upcoming acqusition is going to be in regular mode.
        #whereas the software needs to adjust the limit on the number of samples if the upcoming acquisition
        # is going to be in averaging mode.
        
        
        if self.verbose:
            print('Setting acquisition parameters manually.')
            print('Samples will shorten if too long for avg mode.')
#            print('Sample will autoround to a multiple of 1024 and shorten if to long for avg mode.')
        ###!!!!!!!!!!!!!
        #!!!!!!!!!empircally, record size must be 1024*n, and the largest that it can be is 513*1024
        #roughly 500 kS. Which is what the mnual says, but this sucks. Did the old Acqiris do this?
        #or is it just that the ME is only for non-averaging mode?
        multiple = int(numpy.ceil(self._samples/1024))
        if self._averages > 1:
            #check the number of channels
            numChans = len(self._activeChannels)
            if numChans == 1:
                maxSamples = 1024*1024
            else:
                maxSamples = 1024*512
                
            if self._samples > maxSamples:
                print('Data is too long for averaging mode. Setting to max length: 512*1024')
                self.samples = int(maxSamples)
                
#            if multiple*1024 > maxSamples:
#                print('Data is too long for averaging mode. Setting to max length: 512*1024')
#                self.samples = int(maxSamples)
#            else:
#                #auto round to next multiple of 1024 up, regardless
#                self.samples = int(multiple*1024) 
#                #it is very important that this winds up an integer, or it least it was at some point
                
        #it looks like now that we have the order of operations right between configuring
        #the acquisition and turning averaging on and off, so that we can use the
        #built in C function for configuring, that we might not need this multiple of 1024 check.
        #But we do need the auto check that turns down the number of samples if its averaging mode
        
#        self.driver.Acquisition.ConfigureAcquisition(int(self.segments), int(self.samples), self.sampleRate)
        self.driver.Acquisition.ConfigureAcquisition(int(self._segments), int(self._samples), self._sampleRate)
        
        #configure the the actual averaging mode and restore number of averages
        self.ConfigureAveraging() #trip the flags to handle the averaging
        #(both in python and sending down to the hardware driver.)
        #Hopefully doing it after everything else is set will only try to flip to averager 
        #after the number of samples has been reduced.


    def ConfigureAveraging(self):
        if self._averages >1:
            self.averageMode = 1 #make sure this variable is conistent. I probably want to do away with it eventually
        else:
            self.averageMode = 0
#        self.driver.Acquisition.Mode = self.averageMode 
#        self.driver.Acquisition.NumberOfAverages = self.averages
        self.driver.Acquisition.Mode = self._averageMode 
        self.driver.Acquisition.NumberOfAverages = self._averages

    def ConfigureChannel(self, channelNum = 1, Range = 2.5, offset = 0, enabled = True):
        '''Configure channel function.
        
        Meant to becalled by configure acquisiton, but can be called manually.
        If it is called mnaually, it will write it's inputs to the driver settings.
        '''

        if channelNum in [1,2]:
            pass
        else:
            raise ValueError('Invalid channel. Should be 1 or 2')
        
#        if Range in [0.5, 2.5]: #changed to a property, so this is checked elsewhere
#            self.channelRange = Range
#        else:
#            raise ValueError('Range must be 0.5 or 2.5')
            
        self.channelRange = Range
        self.channelOffset = offset

        couplingC = 1 # 1 = DC coupling. Always. Card can't do AC.

        chan = self.driver.Channels[int(channelNum-1)] #python channels are index from 0
        chan.Configure(self._channelRange[int(channelNum-1)], self.channelOffset, couplingC, enabled)
    
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
            
            self.driver.Trigger.Delay = self.triggerDelay


            #manually set trigger source because the configure function doesn't actually do it.
            self.driver.Trigger.ActiveSource = self.triggerSource
            

            #self.driver.Trigger.Sources[self.triggerSource].Coupling, but I think this always has to be 1
#            self.driver.Trigger.Sources[self.triggerSource].Coupling = self.triggerCoupling
            
            #couldn't find a configure trigger function in the python. Just set all the fields instead
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
        
        
        

        self.hardwareKeys = []
        self.hardwareKeys.append('channelOffset')
        self.hardwareKeys.append('channelRange')
#        self.hardwareKeys.append('channelEnabled') #removed for better syntax
        self.hardwareKeys.append('activeChannels')
        
        self.hardwareKeys.append('triggerSource')
        self.hardwareKeys.append('triggerLevel')
#        self.hardwareKeys.append('triggerMode') #I don't think we can change this one
        self.hardwareKeys.append('triggerSlope')
        self.hardwareKeys.append('triggerDelay')
        self.hardwareKeys.append('triggerCoupling')
        
        self.hardwareKeys.append('averageMode')
        self.hardwareKeys.append('averages')
        
        self.hardwareKeys.append('samples')
        self.hardwareKeys.append('segments')
        
        self.hardwareKeys.append('sampleRate')
        
        self.hardwareKeys.append('simulateMode')
        
        self.hardwareKeys.append('clockFrequency')
        self.hardwareKeys.append('clockSource')
        
#        self.hardwareKeys.append('timeout') #i don't know how to get this from the hardware, so it's a software setting
    
    def _fill_slots(self):
        firstRound = list(self.__dict__.keys()) #this gets all of the normal attributes, but it doesn't handle the properties
        propertyKeys = self.hardwareKeys
        
        self.__slots__ = firstRound + propertyKeys
    
    def _generateConfig(self):
        '''Make a dictionary of the hardware settings and a couple of useful flags.
        
        As it currently stands it will try to get the settings from the hardware
        and ignore the python shadow copies.
        
        This will hopefully be the basis for saving and loading.
        '''
        params = {}
        for key in self.hardwareKeys:
            val = getattr(self,key)
            params[key] = val
        
        #a couple of the python software settings need to go too.
        params['verbose'] = self.verbose
        params['timeout'] = self.timeout
        params['simulate'] = self.simulate
        return params
        
    
    def InitializeDriver(self):
        '''Basic init function. Creates the driver.'''
        print('Initializing...')
        if self.simulate == True:
            strInitOptions =  'Simulate=True,  DriverSetup= model = SA220P'
        else:
            strInitOptions = 'Simulate=False,  DriverSetup= model = SA220P'
        
        self.driver =  AqMD3.AqMD3( self.hardwareAddress , False, False, strInitOptions)
        

    def InitiateAcquisition(self):
        ''' Tells the hardware to start looking for data.
        This is a relatively low-level hardware function.
        User-friendly version with more checks and balances
        is the wrappers Arm and ArmAndWait.'''
        self.driver.Acquisition.Initiate()
        self.armed = True #throw flag so that python knows acquisition has been itiated
        #and system is waiting for trigger.
        
    def LoadConfig(self, params):
        ''' Takes in a dictionary of settings and loads them into this python
        class. Those that are normal properties will get set to hardware.
        The icky properties will not get set until a call to SetParams is made,
        either by the user, or when the flag is tripped in Arm
        '''
        for key in params.keys():
            bool1 = key in self.__dict__.keys()
            bool2 = ('_' + key) in self.__dict__.keys()
            if not (bool1 or bool2):
#                raise ValueError('invalid setting key : ' + str(key))
                warnings.warn('Warning:  Invalid setting key : ' + key, RuntimeWarning)
#                print('Warning:  setting key : ' + key)
            else:
                setattr(self, key, params[key])
        
    def _loadDefaultConfig(self):
        #diagnostic prints
        self.verbose = False   #this needs to be first because it effects how properties below it act
        
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
        self.triggerMode = 'Edge' #this thing should maybe go away.
        self.triggerLevel = 0
        self.triggerSlope = 'Rising'
        self.triggerDelay = 0
        self.triggerCoupling = 1 #1 for DC, 0 for AC. I think it must always be DC, at least on external
        
        #clock settings
        self.clockSource = 'External'
#        self.clockSource = 'Internal'
        self.clockFrequency = 10**7
        
        #timeout of failed acquisitions
        self.timeout = 5  #seconds
        
    def Print(self):
        '''Print the most important settings. 
        Will not show the run time flags and stuff like this.'''
        params = self._generateConfig()
        print()
        for key in params.keys():
            print(key + ' : ' + str(params[key]))
        print()
        
    def ReadAllData(self, returnRaw = False, returnTwoArrays = True):
        ''' Highest level read function. Will automatically read all active
        channels.
        
        Will trip flags when it is done that the data acquisition and read are complete
        
        if returnTwoArrays is true, then it will alwyas return two data arrays, even if channel 
        two is turned off. But it will return an empty array for the inactive channel
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
        
        
#        data1 = []
#        data2 = []
        if len(self._activeChannels) > 1:
            data1 = self.ReadData(1, returnRaw = returnRaw)
            data2 = self.ReadData(2, returnRaw = returnRaw)
            #notify the python that this acquisition is done
            self.armed = False
            self.acquisitionFinished = False
            return data1, data2
        else:
            data =  self.ReadData(self._activeChannels[0], returnRaw = returnRaw)
            data_blank = numpy.zeros((0,0))
            #notify the python that this acquisition and read are done
            self.armed = False
            self.acquisitionFinished = False
            if returnTwoArrays:
                return data, data_blank
            else:
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
        
        if chanNum in [1,2]:
            pass
        else:
            raise ValueError('Invalid channel. Should be 1 or 2')
        
        self.offsetWithinRecord = 0
        numSamples = self.driver.Acquisition.QueryMinWaveformMemory(64, int(self.segments), int(self.offsetWithinRecord), \
                                                       int(self.samples))
        if self.verbose:
            print('Memeory Allocation Determined')
            
        self.totalSamples = numSamples

        if self.verbose:
            print('Trying to Read')
        if self.averageMode:
            out = self._ReadAveragerData(chanNum, returnRaw = returnRaw)
        else:
            if self.segments > 1:
                out = self._ReadMultiSegData(chanNum, returnRaw = returnRaw)
            else:
                out = self._ReadSingleSegData(chanNum, returnRaw = returnRaw)
                

        return out
    
    def _ReadAveragerData(self,chanNum, returnRaw = False):
        if self.verbose:
            print('reading averaged data')
            
        channel = self.driver.Channels[int(chanNum - 1)]
        
        
        waveformObj = channel.Measurement.FetchAccumulatedWaveform(firstRecord = 0,\
                                                                      numberOfRecords = int(self.segments),\
                                                                      offsetWithinRecord = int(self.offsetWithinRecord), \
                                                                      numberOfPointsPerRecord = int(self.samples))   

        if self.verbose:
            print('Fetch Complete. Processing Data.')
            
#        if self.segments == 1:
#            rawData = numpy.zeros( self.samples)
#            waveform  = waveformObj[0]
#            rawData[:] = waveform.Samples* waveform.ScaleFactor + waveform.ScaleOffset
#        else:
#            rawData = numpy.zeros((self.segments, self.samples))
#            for segind in range(0, self.segments):
#                waveform  = waveformObj[segind]
#                rawData[segind,:] = waveform.Samples* waveform.ScaleFactor + waveform.ScaleOffset
            
        #always get matrix shaped data
        rawData = numpy.zeros((self.segments, self.samples))
        for segind in range(0, self.segments):
            waveform  = waveformObj[segind]
            rawData[segind,:] = waveform.Samples* waveform.ScaleFactor + waveform.ScaleOffset
        
        
        if waveformObj[0].ActualSamples != self.samples:
            print("Warning. Data size doesn't match the number of samples. Something wierd happened.")
        
        if returnRaw:
            out = [rawData, waveformObj]
            return out
        else:
            data = rawData
            if self.verbose:
                print('Data Processed.')
            return data
        
    def _ReadMultiSegData(self,chanNum, returnRaw = False):
        '''Python function for reading non-averaged multisegment data. '''
        if self.verbose:
            print('reading non-averaged, multisegment data')
            
        if not (chanNum in[1,2]):
            raise ValueError('Channel number must be 1 or 2.')
            
        channel = self.driver.Channels[int(chanNum - 1)]
        
        
        waveformObj = channel.MultiRecordMeasurement.FetchMultiRecordWaveform(firstRecord = 0,\
                                                                      numberOfRecords = int(self.segments),\
                                                                      offsetWithinRecord = int(self.offsetWithinRecord), \
                                                                      numberOfPointsPerRecord = int(self.samples))
        if self.verbose:
            print('Fetch Complete. Processing Data.')
            
        rawData = numpy.zeros((self.segments, self.samples))
        for segind in range(0, self.segments):
            waveform  = waveformObj[segind]
            rawData[segind,:] = waveform.Samples* waveform.ScaleFactor + waveform.ScaleOffset
        
        if waveformObj[0].ActualSamples != self.samples:
            print("Warning. Data size doesn't match the number of samples. Something wierd happened.")
        
        if returnRaw:
            out = [rawData, waveformObj]
            return out
        else:
            data = rawData
            if self.verbose:
                print('Data Processed.')
            return data

    def _ReadSingleSegData(self,chanNum, returnRaw = False):
        '''Python read function for single segment data that is not averaged '''
        if self.verbose:
            print('reading non-averaged, single-segment data')
             
        if not (chanNum in[1,2]):
            raise ValueError('Channel number must be 1 or 2.')
            
        channel = self.driver.Channels[int(chanNum - 1)]
        
        waveformObj = channel.Measurement.FetchWaveform()
        
        if self.verbose:
            print('Fetch Complete. Processing Data.')
            
            
        rawData = numpy.asarray(waveformObj.Samples*waveformObj.ScaleFactor+ waveformObj.ScaleOffset)
        
        if waveformObj.ActualSamples != self.samples:
            print("Warning. Data size doesn't match the number of samples. Something wierd happened.")
            
        if returnRaw:
            out = [rawData, waveformObj]
            return out
        else:
            #always make matrix shaped data to be compatible with the other read functions
            data = numpy.zeros((1, len(rawData)))
            data[0,:] = rawData
            if self.verbose:
                print('Data Processed.')
            return data
  
    def ReInitialize(self):
        '''Basic init function. Can also be used to wipe he settings if the card is very confused. 
        Will push the currently stored settings to the card and calibrate'''
        try:
            #try to pull the current values of all the settings
            params = self._generateConfig()
            gotSettings = True
        except:
            gotSettings = False
            pass


        #in the python version, there seems to be an issue with reinitializing an already
        #existing instnace of the driver. It causes an error that the sample clocks are unstable
        #or something about mising TDC curves. Running close before reinitializing it seems to allow
        #things to work. So, putting this try in here to try to catch things
        try:
            self.Close()
        except:
            pass
        
        self.InitializeDriver() #this puts the driver back in a default
        #and then some of the properties act funny, because they are only stored in hardware
        #so they revert to the driver default.
        #trying to patch this by pulling the current settings, but if it's vry confuzzled,
        #this may not work. So, if that fails, I will just load the default.
        if gotSettings:
            #load rescued settings
            self.LoadConfig(params)
        else:
            #load the default if you can't get the current ones.
            self._loadDefaultConfig()

       
        #push the currently stored settings to the card and calibrate. Hopefully ths will actually save hassle.
        #if this is running at _init_ then it will push the defaults
        self.SetParams()
        self.SelfCalibrate()

    def WaitForAcquisition(self):
        '''Timeout in seconds '''
        if self.verbose:
            print('Waiting until the acquisition is done (or checking that it is) and getting confirmation')
        timeoutC = self.timeout*1000 #convert to ms
#        self.driver.Acquisition.WaitForAcquisitionComplete(timeoutC)
#        self.acquisitionFinished = True #throw the internal flag to show that the card is done.
        try:
            self.driver.Acquisition.WaitForAcquisitionComplete(timeoutC)
            self.acquisitionFinished = True #throw the internal flag to show that the card is done.
        except (RuntimeError):
            print('Acquisition probably timed out.')
            if self.verbose:
                print('Aborting to clear hung acquisition.')
            self.Abort()
            self.acquisitionFinished = False #throw the internal flag to indicate failure
        



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
    
    
    ######
    #testing raw functions
#    card.ConfigureTrigger()
#    card.ConfigureAveraging()
#    card.ConfigureChannel(channelNum = 1, Range = 2.5, offset = 0, enabled = True)
#    card.SetParams()
    
    

    #####################
    #chose acquisition type
    #####################

    averageMode = True
#    averageMode = False

    multisegMode = True
#    multisegMode = False
    
    
    
    ##########
    #set settings to the card object
    ### samples = 1*10**7
    card.samples = 513*1024
    card.sampleRate = 1*10**9
    
    
    #swtich between test cases
    if averageMode:
        card.averageMode = 1
        if multisegMode:
            card.averages = 100
            card.segments = 2
        else:
            card.averages = 100
            card.segments = 1
    else:
        card.averageMode = 0
        if multisegMode:
            card.averages = 1
            card.segments = 2
        else:
            card.averages = 1
            card.segments = 1
            
    card.triggerSource = 'External1'
    card.triggerLevel = 0
    card.triggerSlope = 'Falling'
    
    card.activeChannels = [1,2]
#    card.activeChannels = [1] #it looks like single channel has to be channel 1
    
    card.timeOut = 2
    
    card.verbose = True #get lots of diagnostic print statements
    
#    card.channelRange = 0.5
#    card.channelRange = [0.5, 2.5]
    card.channelRange = numpy.asarray([0.5, 2.5])
    
    card.channelOffset = 0.1
    
    ##########
    #push settings to hardware
    card.SetParams()
#    card.ConfigureChannel(channelNum = 1, Range = 2.5, offset = 0, enabled = True)
#    card.ConfigureChannel(channelNum = 2, Range = 2.5, offset = 0, enabled = True)
#    card.ConfigureTrigger('External1', 0, 'Falling')
#    card.ConfigureAcquisition(samples, sampleRate, segments)

    #card.GetParams()
    
    ##########
    #initiate acquisition and wait for it to finish
    card.ArmAndWait()
#    card.Arm()
#    card.WaitForAcquisition()
#    card.WaitForAcquisition(timeOut)

    print('Data Acquired (In Theory)')


    
    
#    data = card.ReadData(1, returnRaw = False) #read channel 1
#    data = card.ReadData(2, returnRaw = False) #read channel 1
    data, data2 = card.ReadAllData()
    

    ####
    #plot data

    numSamples = card.samples
    sampleRate = card.sampleRate
    dt = 1/sampleRate
    xaxis = scipy.arange(0, numSamples,1)*dt

    if multisegMode:
        dataSegments = data.shape[0]

    print('Plotting.')

    fig1 = plt.figure(1)
    plt.clf()
    ax = plt.subplot(1,1,1)
    if multisegMode:
        plt.title('Multiseg: Active Channel')
        for segind in  range(0,dataSegments ):
            labelstr = 'seg: ' + str(segind)
            plt.plot(xaxis*1000, data[segind,:] + segind*0.125, label = labelstr)
    else:
        if averageMode:
            plt.title('Averager: Active Channel')
        else:
            plt.title('Singleseg: Active Channel')
        labelstr = 'seg: 0' 
        plt.plot(xaxis*1000, data[:], label = labelstr)
    ax.legend(loc = 'upper right')
    plt.ylabel('Voltage (waterfall)')
    plt.xlabel('Time (ms)')

    plt.show()

#
#    print('Done plotting.')
#    
#    # Close the instrument
#    print("\nClose the instrument")    
#    card.close()
#    
    

    
    
  #    card.offsetWithinRecord = 0
#    channel = card.driver.Channels[0]
#    if averageMode:
#        waveformObj = channel.Measurement.FetchAccumulatedWaveform(firstRecord = 0,\
#                                                       numberOfRecords = int(card.segments),\
#                                                       offsetWithinRecord = int(card.offsetWithinRecord),\
#                                                       numberOfPointsPerRecord = int(card.samples),\
#                                                       dtype=numpy.int32)
#    else:
#        if multisegMode:
#            waveformObj = channel.MultiRecordMeasurement.FetchMultiRecordWaveform(firstRecord = 0,\
#                                                                          numberOfRecords = int(card.segments),\
#                                                                          offsetWithinRecord = int(card.offsetWithinRecord), \
#                                                                          numberOfPointsPerRecord = int(card.samples))
#        else:
#            waveformObj = channel.Measurement.FetchWaveform()  
    
    
    
    
    
    
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
#    fig1 = plt.figure(1)
#    plt.clf()
#    ax = plt.subplot(1,1,1)
#    if multisegMode:
#        plt.title('Multiseg: Active Channel')
#        for segind in  range(0,dataSegments ):
#            labelstr = 'seg: ' + str(segind)
#            plt.plot(xaxis*1000, data[segind,:] + segind*0.125, label = labelstr)
#    else:
#        if averageMode:
#            plt.title('Averager: Active Channel')
#        else:
#            plt.title('Singleseg: Active Channel')
#        labelstr = 'seg: 0' 
#        plt.plot(xaxis*1000, data[:], label = labelstr)
#    ax.legend(loc = 'upper right')
#    plt.ylabel('Voltage (waterfall)')
#    plt.xlabel('Time (ms)')
#
#    plt.show()
#
#
#    print('Done plotting.')
#    
##    # Close the instrument
##    print("\nClose the instrument")    
##    card.close()    
    
    
    
    
    