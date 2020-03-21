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
    def __init__(self, library, ResourceName):        
        self._prefix = "AqMD3"
        self._lib = ctypes.cdll.LoadLibrary(library)
              
        self.visession = ctypes.c_int()
        
        self.call('InitWithOptions',
            ctypes.create_string_buffer(ResourceName.encode('utf-8')),
            False,
            False, 
            '',
            ctypes.byref(self.visession))

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
            raise IviCError(error_message.value)        
            
            
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
            print('Warning: I suspect mutisegment may not work with averager on.')
            self.segments = segments

        if self.multiseg == False:
            self.segments = 1
            if segments > 1:
                print('Single Segment Mode: Ignoring requested extra segments')
        else:
            self.segments = segments

        numRecordsC = ctypes.c_int64(self.segments)
        self.call('ConfigureAcquisition', self.visession, numRecordsC, numPointsPerRecordC, sampleRateC)
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

        
        WaveformArrayC = ctypes.c_longdouble*arraySize_int
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
                                 WaveformArrayC,\
                                 ctypes.byref(ActualRecordsC),\
                                 ActualPointsC,\
                                 FirstValidPointsC,\
                                 InitialXOffsetC,\
                                 InitialXTimeSecondsC,\
                                 InitialXTimeFractionC,\
                                 ctypes.byref(XIncrementC))

        print('Fetch Complete. Processing Data.')
        rawData = numpy.asarray(WaveformArrayC)
        dataRawSize = rawData.size
        dataActualSegments = numpy.asarray(ActualRecordsC)[0]
        dataActualPoints_full = numpy.asarray(ActualPointsC)
        dataActualPoints = dataActualPoints_full[0]
        dataFirstValidPoints = numpy.asarray(FirstValidPointsC).astype('int64')
        if returnRaw:
            out = [rawData, dataActualPoints, dataFirstValidPoints. dataActualSegments]
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

        WaveformArrayC = ctypes.c_double*arraySize_int
        # WaveformArrayC = ctypes.c_longdouble*arraySize_int
        # tempArray = numpy.zeros(arraySize_int)
        # WaveformArrayC.value = tempArray
        WavefromArraySizeC = ctypes.c_int64(self.totalSamples)

        chanNameC = ctypes.c_char_p(chanName)

        self.call('FetchWaveformReal64', self.visession, chanNameC, WavefromArraySizeC, WaveformArrayC,\
            ctypes.byref(ActualPointsC),\
                 ctypes.byref(FirstValidPointC),\
                      ctypes.byref(InitialXOffsetC),\
                 ctypes.byref(InitialXTimeSecondsC),\
                      ctypes.byref(InitialXTimeFractionC),\
                       ctypes.byref(XIncrementC))

        # self.call('FetchWaveformInt32', self.visession, chanNameC, WavefromArrayActualSizeC, WaveformArrayC,\
        #     ctypes.byref(ActualPointsC),\
        #          ctypes.byref(FirstValidPointC),\
        #               ctypes.byref(InitialXOffsetC),\
        #          ctypes.byref(InitialXTimeSecondsC),\
        #               ctypes.byref(InitialXTimeFractionC),\
        #                ctypes.byref(XIncrementC))

        print('Fetch Complete. Processing Data.')
        rawData = numpy.asarray(WaveformArrayC)
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
    
    ResourceName = "PXI23::0::0::INSTR"

    IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
    sys.path.append(IVIbinPath)
    
    print("Initialize instrument @", ResourceName)
    card = AqMd3("AqMd3_64.dll", ResourceName)
    
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

    # averageMode = True
    averageMode = False
    if averageMode:
        numAverages = 40
        segments = 1
        card.SetAverageMode(numAverages)
        # card.SetAttribute(None, AQMD3_ATTR_ACQUISITION_MODE, 1, 'ViInt32') 
        # card.SetAttribute(None, AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, numAverages, 'ViInt32') 
    else:
        # numAverages = 1
        # segments = 2
        # card.SetMultiMode()

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
    card.SetAttribute(None, AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, 'External1', 'ViString')


    samples = 1*10**7
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


    out = card.ReadData(1, returnRaw = True)





    
#    # Close the instrument
#    print("\nClose the instrument")    
#    card.close()