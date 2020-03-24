# -*- coding: utf-8 -*-
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





class AqMd3(object):
    def __init__(self, library, ResourceName, simulate = False):        
        self._prefix = "AqMD3"
        self._lib = ctypes.cdll.LoadLibrary(library)
              
        self.visession = ctypes.c_int()
        
        self.simulate = simulate
        self.hardwareAddress = ResourceName
        self.ReInitialize()
            
#        self.call('InitWithOptions',
#            ctypes.create_string_buffer(ResourceName.encode('utf-8')),
#            False,
#            False, 
#            '',
#            ctypes.byref(self.visession))

    def SetAttribute(self, RepCapIdentifier, Attribute, Value, AttributeType) :

        # Encode the RepCapIdentifier if needed
        if (RepCapIdentifier != None):
            RepCapIdentifier = ctypes.c_char_p(RepCapIdentifier.encode('utf-8'))
                   
        if (AttributeType == 'ViBoolean') :            
            self.call('SetAttributeViBoolean',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                  ctypes.c_bool(Value))
            
        if (AttributeType == 'ViInt32') :
            self.call('SetAttributeViInt32',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                 ctypes.c_int32(Value))
        
        if (AttributeType == 'ViInt64') :            
            self.call('SetAttributeViInt64',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                  ctypes.c_int64(Value))

        if (AttributeType == 'ViReal64') :
            self.call('SetAttributeViReal64',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                  ctypes.c_double(Value))
            
        if (AttributeType == 'ViString') :
            buffer = ctypes.create_string_buffer(256)
            self.call('SetAttributeViString',
                  self.visession,
                  RepCapIdentifier,
                  Attribute,
                  ctypes.c_char_p(Value.encode('utf-8')))

    def GetAttribute(self, RepCapIdentifier, Attribute, AttributeType) :
        
        # Encode the RepCapIdentifier if needed
        if (RepCapIdentifier != None):
            RepCapIdentifier = ctypes.c_char_p(RepCapIdentifier.encode('utf-8'))
                   
        if (AttributeType == 'ViBoolean') :
            buffer = ctypes.c_bool();
            self.call('GetAttributeViBoolean',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  ctypes.byref(buffer))
            return buffer.value      
            
        if (AttributeType == 'ViInt32') :
            buffer = ctypes.c_int32()
            self.call('GetAttributeViInt32',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  ctypes.byref(buffer))
            return buffer.value
        
        if (AttributeType == 'ViInt64') :
            buffer = ctypes.c_int64()
            self.call('GetAttributeViInt64',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  ctypes.byref(buffer))
            return buffer.value

        if (AttributeType == 'ViReal64') :
            buffer = ctypes.c_double()
            self.call('GetAttributeViReal64',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  ctypes.byref(buffer))
            return buffer.value
            
        if (AttributeType == 'ViString') :
            buffer = ctypes.create_string_buffer(256)
            self.call('GetAttributeViString',
                  self.visession,
                  RepCapIdentifier,
                  Attribute, 
                  256,
                  ctypes.byref(buffer))
            return buffer.value.decode('ascii')        
           
        return ''

    def call(self, funcname, *args):
        method_to_call = getattr(self._lib, self._prefix + '_' + funcname)
        return_value = method_to_call(*args)        
        if return_value:
            error_message = ctypes.create_string_buffer(256)
            self.call('error_message',
                      self.visession,
                      return_value, 
                      ctypes.byref(error_message))
            # raise IviCError(error_message.value)   
            raise ValueError(error_message.value)      
            
            
    def close(self):    
        self.call('close', self.visession) 
        
    def SelfCalibrate(self):    
        self.call('SelfCalibrate', self.visession) 

    def ConfigureAcquisition(self, samples, sampleRate, segments = 1):
        
        self.samples = samples
        numPointsPerRecordC = ctypes.c_int64(self.samples)
        self.sampleRate = sampleRate
        sampleRateC = ctypes.c_double(self.sampleRate)

        #trying to handle the vairous acquistion types. Not quite sure how this logic will 
        #eventually go.
        if self.averager == True:
            print('Warning: Manual says multiseg should work with this, not sure I have the right read function though.')
            self.segments = segments
            numRecordsC = ctypes.c_int64(self.segments)

            print('Trying to set acquisition parameters manually.')
            ###!!!!!!!!!!!!!
            #!!!!!!!!!empircally, record size must be 1024*n, and the largest that it can be is 513*1024
            #roughly 500 kS. Which is what the mnual says, but this sucks. Did the old Acqiris do this?
            #or is it just that the ME is only for non-averaging mode?
            multiple = numpy.ceil(self.samples/1024)
            if multiple > 512:
                print('Data is too long for averaging mode. Setting to max length: 512*1024')
                multiple = 512
                self.samples = multiple*1024
                numPointsPerRecordC = ctypes.c_int64(self.samples)

            self.SetAttribute(None, AQMD3_ATTR_RECORD_SIZE, int(1024*multiple), 'ViInt64')
            self.SetAttribute(None, AQMD3_ATTR_NUM_RECORDS_TO_ACQUIRE, int(self.segments), 'ViInt64')
            self.SetAttribute(None, AQMD3_ATTR_SAMPLE_RATE, self.sampleRate, 'ViReal64')

            # self.call('ConfigureAcquisition', self.visession, numRecordsC, numPointsPerRecordC, sampleRateC)

        else:
            if self.multiseg == False:
                self.segments = 1
                if segments > 1:
                    print('Single Segment Mode: Ignoring requested extra segments')
            else:
                self.segments = segments

            numRecordsC = ctypes.c_int64(self.segments)
            print('Using auto config function for the acquisition.')
            self.call('ConfigureAcquisition', self.visession, numRecordsC, numPointsPerRecordC, sampleRateC)
            #I think this configure funtion only works in non-averaging mode. !!!?????

        #stupid testing
        # card.SetSingleMode()
        # self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 0, 'ViInt32') 
        # self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, 1, 'ViInt32') 

        # print('Configuring')
        # # self.call('ConfigureAcquisition', self.visession, numRecordsC, numPointsPerRecordC, sampleRateC)
        # self.call('ConfigureAcquisition', self.visession, ctypes.c_int64(numAverages), numPointsPerRecordC, sampleRateC)
        # print('Setting averaging mode')
        # self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 1, 'ViInt32') 
        # print('Setting averages')
        # self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, 4, 'ViInt32') 

        # print('going to averager mode')
        # Avs = 4
        # self.SetAverageMode(Avs)

        # print('going to single mode')
        # self.SetSingleMode()
        # print('configuring')
        # self.call('ConfigureAcquisition', self.visession, ctypes.c_int64(1), numPointsPerRecordC, sampleRateC)

        # print('going to averager mode')
        # Avs = 4
        # self.SetAverageMode(Avs)
        # print('configuring')
        # self.call('ConfigureAcquisition', self.visession, ctypes.c_int64(Avs), numPointsPerRecordC, sampleRateC)

        # print('Configuring')
        # self.call('ConfigureAcquisition', self.visession, numRecordsC, numPointsPerRecordC, sampleRateC)
        #ConfigureAcquisition(init_status.ViSession, numRecordsC, numPointsPerRecordC, sampleRateC)

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

        offsetC = ctypes.c_double(offset)

        couplingC = ctypes.c_int32(1) # 1 = DC coupling. Always. Card can't do AC.

        enabledC = ctypes.c_bool(enabled)

        self.call('ConfigureChannel', self.visession, chanNameC, rangeC, offsetC, couplingC, enabledC)

    # def ConfigureEdgeTriggerSource:
    #     #call configure edge trigger

    #     #actually set the trigger source
    #     pass
    
    def ConfigureTrigger(self, Source = 'External1', Level = 0, Slope = 'Falling', Mode = 'Edge'):
        if Mode != 'Edge':
            raise ValueError('This trigger mode not yet supported. Edge trigger only.')
        else:
            self.triggerMode = Mode
            self.triggerSource = Source
            self.triggerLevel = Level
            
            print(Slope)
            if Slope == 'Falling':
                self.triggerSlope = Slope
                triggerSlopeC = ctypes.c_int32(0)
            elif Slope == 'Rising':
                self.triggerSlope = Slope
                triggerSlopeC = ctypes.c_int32(1)
            else:
                raise ValueError("Edge trigger slope must be either 'Rising' or 'Falling'")
            
            triggerSourceC = ctypes.create_string_buffer(self.triggerSource.encode('utf-8'))
            triggerLevelC = ctypes.c_longdouble(self.triggerLevel)

            self.call('ConfigureEdgeTriggerSource', self.visession, triggerSourceC ,triggerLevelC, triggerSlopeC)
            #manually set trigger source because the configure function doesn't actually do it.
            self.SetAttribute(None,AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, self.triggerSource, 'ViString')

    def InitiateAcquisition(self):
        self.call('InitiateAcquisition', self.visession)

    def ReadData(self, chanNum, returnRaw = False):
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
        offsetWithinRecordC = ctypes.c_int64(self.offsetWithinRecord) 
        numSamplesC = ctypes.c_int64()
        self.call('QueryMinWaveformMemory', self.visession, ctypes.c_int32(64), numRecordsC,\
             offsetWithinRecordC, numPointsPerRecordC, ctypes.byref(numSamplesC))
        print('Memeory Allocation Determined')
        numSamples = int(numSamplesC.value)
        self.totalSamples = numSamples

        print('Trying to Read')
        if self.averager:
            out = self._ReadAveragerData(chanName, returnRaw = returnRaw)
        else:
            if self.multiseg:
                out = self._ReadMultiSegData(chanName, returnRaw = returnRaw)
            else:
                out = self._ReadSingleSegData(chanName, returnRaw = returnRaw)
        return out

    def _ReadAveragerData(self,chanName, returnRaw = False):
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

            if dataActualSegments == 1:
                startInd = dataFirstValidPoints[0]
                data = rawData[startInd:(startInd+dataActualPoints)]
            else:
                data = numpy.zeros((dataActualSegments,dataActualPoints))
                for segind in range(0,dataActualSegments ):
                    startInd = dataFirstValidPoints[segind]
                    
                    data[segind,:] = rawData[startInd:(startInd+dataActualPoints)]
            print('Data Processed.')
            return data

    def _ReadMultiSegData(self,chanName, returnRaw = False):
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
            _fields_ = [("WaveformArrayC", ctypes.c_longdouble*arraySize_int)]
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
                                 ActualPointsC,\
                                 FirstValidPointsC,\
                                 InitialXOffsetC,\
                                 InitialXTimeSecondsC,\
                                 InitialXTimeFractionC,\
                                 ctypes.byref(XIncrementC))

        print('Fetch Complete. Processing Data.')
        rawData = numpy.asarray(WaveHolder.WaveformArrayC)
        dataRawSize = rawData.size
        dataActualSegments = numpy.asarray(ActualRecordsC)[0]
        dataActualPoints_full = numpy.asarray(ActualPointsC)
        dataActualPoints = dataActualPoints_full[0]
        dataFirstValidPoints = numpy.asarray(FirstValidPointsC).astype('int64')
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
        '''Basic init function. Can also be used to wipe he settings if the card is very confused. '''
        if self.simulate == True:
            strInitOptionsC = ctypes.c_char_p(b'Simulate=True,  DriverSetup= model = SA220P')
        else:
            strInitOptionsC = ctypes.c_char_p(b'Simulate=False,  DriverSetup= model = SA220P')
            
        self.call('InitWithOptions',
            ctypes.create_string_buffer(self.hardwareAddress.encode('utf-8')),
            False,
            False, 
            strInitOptionsC,
            ctypes.byref(self.visession))

    def SetAverageMode(self, averages):
        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 1, 'ViInt32') 
        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, numAverages, 'ViInt32') 
        self.averager = True
        self.multiseg = False

    def SetSingleMode(self):
        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 0, 'ViInt32') 
        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, 1, 'ViInt32') 
        self.averager = False
        self.multiseg = False

    def SetMultiMode(self):
        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 0, 'ViInt32') 
        self.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, 1, 'ViInt32') 
        self.averager = False
        self.multiseg = True   #segment number will be configured elsewhere

    def WaitForAcquisition(self, timeOut = 5):
        '''Timeout in seconds '''
        timeOutC = ctypes.c_int32(timeOut*1000)
        self.call('WaitForAcquisitionComplete', self.visession, timeOutC)



if __name__ == '__main__':
    
    # Declaration of some attributes (values from AqMd3.h)
    IVI_ATTR_BASE          = 1000000
    IVI_INHERENT_ATTR_BASE = IVI_ATTR_BASE +  50000
    IVI_CLASS_ATTR_BASE    = IVI_ATTR_BASE +  250000
    IVI_LXISYNC_ATTR_BASE  = IVI_ATTR_BASE +  950000
    IVI_SPECIFIC_ATTR_BASE = IVI_ATTR_BASE +  150000
        
    #firmware
    AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION          = IVI_INHERENT_ATTR_BASE + 514  # ViString, read-only 
    AQMD3_ATTR_SPECIFIC_DRIVER_REVISION             = IVI_INHERENT_ATTR_BASE + 551  # ViString, read-only   
    AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION         = IVI_INHERENT_ATTR_BASE + 510  # ViString, read-only
    AQMD3_ATTR_INSTRUMENT_MANUFACTURER              = IVI_INHERENT_ATTR_BASE + 511  # ViString, read-only
    AQMD3_ATTR_INSTRUMENT_MODEL                     = IVI_INHERENT_ATTR_BASE + 512  # ViString, read-only
    AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING = IVI_SPECIFIC_ATTR_BASE + 8    # ViString, read-only 
    
    #channel
    AQMD3_ATTR_VERTICAL_OFFSET                      = IVI_CLASS_ATTR_BASE + 25      # ViReal64, read-write 
    
    #trigger
    AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE                = IVI_CLASS_ATTR_BASE + 1       # ViString, read-write

    #averaging
    AQMD3_ATTR_ACQUISITION_MODE                     = IVI_SPECIFIC_ATTR_BASE + 11   # ViInt32, read-write 
    AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES       = IVI_SPECIFIC_ATTR_BASE + 69   # ViInt32, read-write 

    #acquistion parameters
    AQMD3_ATTR_RECORD_SIZE                          = IVI_CLASS_ATTR_BASE + 14      # ViTnt64, read-write
    AQMD3_ATTR_NUM_RECORDS_TO_ACQUIRE               = IVI_CLASS_ATTR_BASE + 13      # ViInt64, read-write

    #sampling and clocking
    AQMD3_ATTR_SAMPLE_RATE                          = IVI_CLASS_ATTR_BASE + 15       # ViReal64, read-write  
   
    #simulation mode
    AQMD3_ATTR_SIMULATE                             = IVI_INHERENT_ATTR_BASE + 5       # ViBoolean, read-write  
    
    
    
    
    
    ResourceName = "PXI23::0::0::INSTR"

    IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
    sys.path.append(IVIbinPath)
    
#    alwaysInit = False
    alwaysInit = True
    if alwaysInit:
        print("Initialize instrument @", ResourceName)
        card = AqMd3("AqMd3_64.dll", ResourceName)
    else:
        try:
            #random driver call to see if things are set up
            card.GetAttribute(None, AQMD3_ATTR_SIMULATE, 'ViBoolean')
            print('Card already initialized')
        except:
            print("Initialize instrument @", ResourceName)
            card = AqMd3("AqMd3_64.dll", ResourceName)
    
    
    
    # Read some attributes
    print('Manufacturer  : {}'.format(card.GetAttribute(None, AQMD3_ATTR_INSTRUMENT_MANUFACTURER, 'ViString')))
    print('Description   : {}'.format(card.GetAttribute(None, AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION, 'ViString')))
    print('Revision      : {}'.format(card.GetAttribute(None, AQMD3_ATTR_SPECIFIC_DRIVER_REVISION, 'ViString')))
    print('Model         : {}'.format(card.GetAttribute(None, AQMD3_ATTR_INSTRUMENT_MODEL, 'ViString')))
    print('Firmware Rev. : {}'.format(card.GetAttribute(None, AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION, 'ViString')))
    print('Serial #      : {}'.format(card.GetAttribute(None, AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING, 'ViString')))
    
    # # Read Channel1 offset
    # print('\nRead Channel1 Offset: {}'.format(drv.GetAttribute('Channel1', AQMD3_ATTR_VERTICAL_OFFSET, 'ViReal64')))
    
    # # Set Channel1 offset to 0.5V
    # print('\nChange Channel1 offset to 0.05 V')
    # drv.SetAttribute('Channel1', AQMD3_ATTR_VERTICAL_OFFSET, 0.05, 'ViReal64')

    # # Read Channel1 offset
    # print('\nRead Channel1 Offset: {}'.format(drv.GetAttribute('Channel1', AQMD3_ATTR_VERTICAL_OFFSET, 'ViReal64')))    
    
    # # Read Channel1 offset
    # print('\nRead Channel1 Offset: {}'.format(drv.GetAttribute('Channel1', AQMD3_ATTR_VERTICAL_OFFSET, 'ViReal64')))    
    
    # # Read active trigger source
    # print('\nRead active trigger source : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, 'ViString')))    

    # # Set active trigger source to External1
    # print('\nSet active trigger source to External1')
    # drv.SetAttribute(None, AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, 'External1', 'ViString')
    
    # # Read active trigger source
    # print('\nRead active trigger source : {}'.format(drv.GetAttribute(None, AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, 'ViString')))  

    # # Set acquisition mode to averager
    # print('\nSet acquisition mode to averager')
    # drv.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 1, 'ViInt32')   
 
    # # Read acquisition mode
    # print('\nRead acquisition mode: {}'.format(drv.GetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 'ViInt32')))
    
    # # Perform a self-calibration
    # print("\nPerform internal calibrations")    
    # card.SelfCalibrate()
    



    #####################
    #set up acquisition and try
    #####################

#    averageMode = True
    averageMode = False

    # multisegMode = True
    multisegMode = False

    if averageMode:
        numAverages = 100
        segments = 1
        card.SetAverageMode(numAverages)
        # card.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 1, 'ViInt32') 
        # card.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, numAverages, 'ViInt32') 
    else:
        if multisegMode:
            numAverages = 1
            segments = 2
            card.SetMultiMode()
        else:
            numAverages = 1
            segments = 1
            card.SetSingleMode()

        # card.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 0, 'ViInt32') 
        # card.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, numAverages, 'ViInt32') 

    #configure channels
    card.ConfigureChannel(channelNum = 1, Range = 2.5, offset = 0, enabled = True)
    card.ConfigureChannel(channelNum = 2, Range = 2.5, offset = 0, enabled = True)

    # Set active trigger source to External1
    print('\nSet active trigger source to External1')
#    card.SetAttribute(None, AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, 'External1', 'ViString')
    card.ConfigureTrigger('External1', 0, 'Falling')


    # samples = 1*10**7
    samples = 513*1024
    sampleRate = 1*10**9
    # aconfigout = driverdll.AqMD3_ConfigureAcquisition(init_status.ViSession, numRecordsC, numPointsPerRecordC, sampleRateC)
    card.ConfigureAcquisition(samples, sampleRate, segments)

    # acqout = driverdll.AqMD3_InitiateAcquisition(init_status.ViSession)
    # waitout = driverdll.AqMD3_WaitForAcquisitionComplete(init_status.ViSession,c_int32(timeouttime) )
    try:
        card.InitiateAcquisition()
    except:
        print('Initiate Failed. Trying a Calibration and Reinnitiate')
        card.SelfCalibrate()
        card.InitiateAcquisition()

    timeOut = 10
    card.WaitForAcquisition(timeOut)

    print('Data Acquired (In Theory)')


    # out = card.ReadData(1, returnRaw = True)
    data = card.ReadData(1, returnRaw = False)


    ####
    #plot data

    numSamples = card.samples
    #sampleRate = card.GetAttribute('Channel1', AQMD3_ATTR_SAMPLE_RATE, 'ViReal64')
    dt = 1/sampleRate
    xaxis = scipy.arange(0, numSamples,1)*dt

    if multisegMode:
        dataSegments = data.shape[0]

    print('Plotting.')

    fig1 = pylab.figure(1)
    pylab.clf()
    ax = pylab.subplot(1,1,1)
    if multisegMode:
        pylab.title('Multiseg: Active Channel')
        for segind in  range(0,dataSegments ):
            labelstr = 'seg: ' + str(segind)
            pylab.plot(xaxis*1000, data[segind,:] + segind*0.125, label = labelstr)
    else:
        if averageMode:
            pylab.title('Averager: Active Channel')
        else:
            pylab.title('Singleseg: Active Channel')
        labelstr = 'seg: 0' 
        pylab.plot(xaxis*1000, data[:], label = labelstr)
    ax.legend(loc = 'upper right')
    pylab.ylabel('Voltage (waterfall)')
    pylab.xlabel('Time (ms)')

    pylab.show()


    print('Done plotting.')
    
#    # Close the instrument
#    print("\nClose the instrument")    
#    card.close()