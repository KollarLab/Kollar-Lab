from ctypes import *
import os
import time
#import subprocess

#import re
#import scipy
import pylab
#import tarfile
#import struct
#import glob
import numpy
import time

#import pickle
#import datetime
#import itertools
#import sys

#os.chdir(r'C:\Users\Kollarlab\Acqiris\testing')

##############################
###### error message check
class ERRORMSG(Structure):
#    _fields_ = [("ViStatus", c_ulong),
#                ("ErrorMessage", c_char*256)]
    _fields_ = [("ErrorMessage", c_char*256)]
    
def get_error_message(outC):
    err_status = ERRORMSG()
    errorout = driverdll.AqMD3_error_message(init_status.ViSession, outC, byref(err_status))
    print (errorout)
    #print err_status.ViStatus
    print (err_status.ErrorMessage)
    return
        

class INTVALUE(Structure):
    _fields_ = [("AttributeValue", c_int64)]       

class REAL64VALUE(Structure):
    _fields_ = [("AttributeValue", c_longdouble)]    

class BOOLVALUE(Structure):
    _fields_ = [("AttributeValue", c_bool)]  
    
class STRINGVALUE(Structure):
    _fields_ = [("AttributeValue",c_char*256)]   
    
attr_intval = INTVALUE() 
attr_realval = REAL64VALUE()
attr_boolval = BOOLVALUE()
attr_stringval = STRINGVALUE()



#class INITSTATUS(Structure):
#    #    _fields_ = [("ViStatus", c_ulong),
#    #                ("ViSession", c_ulong)]    
#    #    _fields_ = [("ViSession", c_ulong),
#    #                ("ViStatus", c_ulong)]  
#        _fields_ = [("ViSession", c_ulong)]  
#init_status = INITSTATUS()  
    
#class CHANELSTATUS(Structure):
#    _fields_ = [("ViStatus", c_ulong),
#                ]    
#chan_status = CHANELSTATUS()  


#class TRIGSTATUS(Structure):
#    _fields_ = [("ViStatus", c_ulong)]    
#trig_status = TRIGSTATUS() 


def get_int32Val(IDnum):
    getout = driverdll.AqMD3_GetAttributeViInt32(init_status.ViSession, None, c_int32(IDnum), byref(attr_intval) )
    return attr_intval.AttributeValue

def get_realVal(IDnum):
    getout = driverdll.AqMD3_GetAttributeViReal64(init_status.ViSession, None, c_int32(IDnum), byref(attr_realval) )
    return attr_realval.AttributeValue

def get_stringVal(IDnum, bufferSize = 256):
    getout = driverdll.AqMD3_GetAttributeViString(init_status.ViSession, None, c_int32(IDnum),c_int32(bufferSize), byref(attr_stringval) )
    return attr_stringval.AttributeValue

def set_stringVal(IDnum, strr):
    setout = driverdll.AqMD3_SetAttributeViString(init_status.ViSession, None, c_int32(IDnum), strr )
    return setout

##############################


driverdll = WinDLL(r'AqMD3_64.dll')


hardwareAddressC = c_char_p(b'PXI23::0::0::INSTR') #needs to e a byte string for python 3. Used to be regular string would work.


# strInitOptionsC = c_char_p(b'Simulate=False, DriverSetup= CAL=0, Trace=false, model = SA220P')
#strInitOptionsC = c_char_p(b'Simulate=True, DriverSetup= CAL=0, Trace=false, model = SA220P')
strInitOptionsC = c_char_p(b'Simulate=False')

idQueryC = c_bool(True)
resetC   = c_bool(True)

###############################
#init = True
init = False
if init:
    print ('initializing ...')
    ######initialize card  
    class INITSTATUS(Structure):
        #    _fields_ = [("ViStatus", c_ulong),
        #                ("ViSession", c_ulong)]    
        #    _fields_ = [("ViSession", c_ulong),
        #                ("ViStatus", c_ulong)]  
            _fields_ = [("ViSession", c_ulong)]  
    init_status = INITSTATUS()  
    initout = driverdll.AqMD3_InitWithOptions(hardwareAddressC,idQueryC, resetC , strInitOptionsC,  byref(init_status))            
    if initout == 0:
        print('Driver Initialized')
    else:
        print("Drivier initialization error = " + str(initout))
        get_error_message(initout)
    
    print ('calibrating ...')    
    calout = driverdll.AqMD3_SelfCalibrate(init_status.ViSession)            
    if calout == 0:
        print('Calibrated')
    else:
        print("Calibration error = " + str(calout))
        get_error_message(calout)
        
    
print ('setting channels ...')
##### check the channels
class CHANELSTATUS(Structure):
    _fields_ = [("ViStatus", c_ulong),
                ]    
chan_status = CHANELSTATUS()  

rangeVal =2.5
if rangeVal in [0.5, 2.5]:
    rangeC = c_double(rangeVal)
else:
    raise ValueError('Range must be 0.5 or 2.5')

offset = 0
offsetC = c_double(offset)

#coupling always DC? Yes.
couplingC = c_int32(1) # 1 = DC coupling


chanStr = b'Channel1'
chanC = c_char_p(chanStr)
chanout1 = driverdll.AqMD3_ConfigureChannel(init_status.ViSession, chanC,rangeC,offsetC, couplingC , c_bool(True))  


chanStr = b'Channel2'
chanC = c_char_p(chanStr)
chanout2 = driverdll.AqMD3_ConfigureChannel(init_status.ViSession, chanC,rangeC,offsetC, couplingC , c_bool(True))  
if chanout1 ==0 and chanout2 == 0:
    print ('Channels Initialized')
else:
    print ("Channel error = " + str(chanout1) + ' or ' + str(chanout2))
    get_error_message(chanout1)
    get_error_message(chanout2)
    


#is it a real boy?
getout = driverdll.AqMD3_GetAttributeViBoolean(init_status.ViSession, None, c_int32(1050005), byref(attr_boolval) )
print ('simulated = ' + str(attr_boolval.AttributeValue))




############## calibration seems to be needed before taking data?
#print ('calibrating ...')    
#calout = driverdll.AqMD3_SelfCalibrate(init_status.ViSession)            
#if calout == 0:
#    print('Calibrated')
#else:
#    print("Calibration error = " + str(calout))
#    get_error_message(calout)


#############set the acquisition
#numPointsPerRecord = c_int64(10000000)
#numRecords = c_int64(1)
numPointsPerRecordC = c_int64(1*10**7)
numRecordsC = c_int64(2)
sampleRateC = c_longdouble(1*10**9) #not currently being used correctly


aconfigout = driverdll.AqMD3_ConfigureAcquisition(init_status.ViSession, numRecordsC, numPointsPerRecordC, sampleRateC)
if aconfigout == 0:
    print ('acquisition set up')
else:
    print ("acquisition setup error = " + str(aconfigout))
    get_error_message(aconfigout)
        
    
    
    
########### set up the trigger
class TRIGSTATUS(Structure):
    _fields_ = [("ViStatus", c_ulong)]    
trig_status = TRIGSTATUS() 

#softwareTrigger = True
softwareTrigger = False
if softwareTrigger:
    trigout = driverdll.AqMD3_ConfigureEdgeTriggerSource(init_status.ViSession, c_char_p(b'Software') , c_double(0), c_int32(0))
    set_stringVal(1250001, b'Software')    
else:
    trigout = driverdll.AqMD3_ConfigureEdgeTriggerSource(init_status.ViSession, c_char_p(b'External1') , c_double(0), c_int32(0))
    set_stringVal(1250001, b'External1')    
#    trigout = driverdll.AqMD3_ConfigureEdgeTriggerSource(init_status.ViSession, c_char_p(b'Internal1') , c_double(0), c_int32(0))

if trigout == 0:
    print ('trigger set up')
else:
    print ("trigger error = " + str(trigout))
    get_error_message(trigout)

    
    
 
##############get ready to allocte memeory
class MINMEMQUERRY(Structure):
#    _fields_ = [("NumSamples", c_ulong)]  
#    _fields_ = [("NumSamples", c_int64)]  
#    _fields_ = [("NumSamples", c_int64),
#                ("ViStatus", c_ulong)]  
    _fields_ = [("NumSamples", c_int64)]  
numSamples =  MINMEMQUERRY() 
memout = driverdll.AqMD3_QueryMinWaveformMemory(init_status.ViSession, c_int32(64), numRecordsC, c_int64(0), numPointsPerRecordC, byref(numSamples))
#driverdll.AgMD1_QueryMinWaveformMemory(byref(numSamples), c_int64(64), numRecords, c_int64(64), numPointsPerRecord)
if memout == 0:
    print ('Tried to get the min memory')
else:
    print ("mem error = " + str(memout))
    get_error_message(memout)

print (numSamples.NumSamples)




############## calibration seems to be needed before taking data?
#print ('calibrating ...')    
#calout = driverdll.AqMD3_SelfCalibrate(init_status.ViSession)            
#if calout == 0:
#    print('Calibrated')
#else:
#    print("Calibration error = " + str(calout))
#    get_error_message(calout)

##############prep for acquisition
class MEASSTATUS(Structure):
    _fields_ = [("Measuring" , c_int32)]    
meas_status = MEASSTATUS() 


measout = driverdll.AqMD3_IsMeasuring(init_status.ViSession,byref(meas_status))

if measout == 0:
    print ('querried measureent status')
else:
    print ("error querrying measurment status = " + str(measout))
    get_error_message(measout)
   



    
#################initiate acquisition
print ('starting acquisition')
acqout = driverdll.AqMD3_InitiateAcquisition(init_status.ViSession)
if acqout  == 0:
    print ('Initialized acquisition')
else:
    print ("arming error = " + str(acqout))
    get_error_message(acqout)
    
    #trying to recover from its sometimes need for recalibration
    abortout = driverdll.AqMD3_Abort(init_status.ViSession)
    
    print ('trying to recalbirate')
    calout = driverdll.AqMD3_SelfCalibrate(init_status.ViSession)  
    if calout == 0:
        print('Calibrated')
    else:
        print("Calibration error = " + str(calout))
        get_error_message(calout)
        
    print ('trying to restart acquisition')
    acqout = driverdll.AqMD3_InitiateAcquisition(init_status.ViSession)
    if acqout  == 0:
        print ('Initialized acquisition')
    else:
        print ("arming error = " + str(acqout))
        get_error_message(acqout)
    


if softwareTrigger:
    for steps in range(0, numRecordsC.value+1):
        
        print ('sending software trigger')
        softtrigout = driverdll.AqMD3_SendSoftwareTrigger(init_status.ViSession)
        if softtrigout == 0:
            print ('sent software trigger')
        else:
            print ("software trigger error = " + str(softtrigout))
            get_error_message(softtrigout)
        
        time.sleep(0.1)
else:
    #wait for acquisition to be done
    dt = 0.01
    maxT = 5
    for ind in range(0,int(numpy.ceil(maxT/dt))):
        print (str(ind))
        measout = driverdll.AqMD3_IsMeasuring(init_status.ViSession,byref(meas_status))
        flag = meas_status.Measuring
        if flag == 1:
            print ('Still Measureing')
            time.sleep(0.01)
        else:
            print ('Done Measuring')
            break
    
    
#arraySize = int(numSamples.NumSamples)
arraySize = int(numRecordsC.value*numPointsPerRecordC.value) #forced it to take the right amount of data for now.
temp2 = int(numRecordsC.value)*1+ 1
#temp2 = int(numRecords.value)*1

samplesToTake = arraySize *temp2
#samplesToTake = numSamples.NumSamples


class FETCHPARAMS(Structure):
    _fields_ = [("WaveformArray", c_longdouble*arraySize),
                ("WaveformArrayActualSize", c_int64*1),
                ("ActualRecords", c_int64*1),
                ("ActualPoints", c_int64*temp2),
                ("FirstValidPoint", c_int64*temp2),
                ("InitialXOffset", c_longdouble*temp2),
                ("InitialXTimeSeconds", c_longdouble*temp2),
                ("InitialXTimeFraction", c_longdouble*temp2),
                ("XIncrement", c_longdouble*1),
                ("ScaleFactor", c_longdouble*1),
                ("ScaleOffset", c_longdouble*1)]  
#    _fields_ = [("WaveformArray", c_int16*arraySize),
#                ("WaveformArrayActualSize", c_int64*1),
#                ("ActualRecords", c_int64*1),
#                ("ActualPoints", c_int64*temp2),
#                ("FirstValidPoint", c_int64*temp2),
#                ("InitialXOffset", c_longdouble*temp2),
#                ("InitialXTimeSeconds", c_longdouble*temp2),
#                ("InitialXTimeFraction", c_longdouble*temp2),
#                ("XIncrement", c_longdouble*1),
#                ("ScaleFactor", c_longdouble*1),
#                ("ScaleOffset", c_longdouble*1),
#                ("ViStatus", c_ulong)]  
ta = time.time()
fetch_params = FETCHPARAMS()
tb = time.time()
print ('preallocation time = ' + str(tb-ta))


t0 = time.time()
fetchout = driverdll.AqMD3_FetchMultiRecordWaveformReal64(init_status.ViSession,\
                                                         c_char_p(b'Channel2'),\
                                                         c_int64(0),\
                                                         numRecordsC,\
                                                         c_int64(0),\
                                                         numPointsPerRecordC,\
                                                         c_int64(samplesToTake),\
                                                         fetch_params.WaveformArray,\
                                                         byref(fetch_params.WaveformArrayActualSize),\
                                                         byref(fetch_params.ActualRecords),\
                                                         fetch_params.ActualPoints,\
                                                         fetch_params.FirstValidPoint,\
                                                         fetch_params.InitialXOffset,\
                                                         fetch_params.InitialXTimeSeconds,\
                                                         fetch_params.InitialXTimeFraction,\
                                                         byref(fetch_params.XIncrement),\
                                                         byref(fetch_params.ScaleFactor),\
                                                         byref(fetch_params.ScaleOffset))
t1 = time.time()
print ('time to pull data as ctypes = ' + str(t1-t0))

#fetchout = driverdll.AqMD3_FetchMultiRecordWaveformReal64(init_status.ViSession,\
#                                                         c_char_p(unicode('Channel2')),\
#                                                         c_int64(0),\
#                                                         numRecords,\
#                                                         c_int64(0),\
#                                                         numPointsPerRecordC,\
#                                                         c_int64(numSamples.NumSamples),\
#                                                         fetch_params.WaveformArray,\
#                                                         byref(fetch_params.WaveformArrayActualSize),\
#                                                         byref(fetch_params.ActualRecords),\
#                                                         fetch_params.ActualPoints,\
#                                                         fetch_params.FirstValidPoint,\
#                                                         fetch_params.InitialXOffset,\
#                                                         fetch_params.InitialXTimeSeconds,\
#                                                         fetch_params.InitialXTimeFraction,\
#                                                         byref(fetch_params.XIncrement),\
#                                                         byref(fetch_params.ScaleFactor),\
#                                                         byref(fetch_params.ScaleOffset))
#class SINGLEFETCHPARAMS(Structure):
#    _fields_ = [("WaveformArray", c_int16*arraySize),
#                ("WaveformArrayActualSize", c_int64),
#                ("ActualPoints", c_int64),
#                ("FirstValidPoint", c_int64),
#                ("InitialXOffset", c_longdouble),
#                ("InitialXTimeSeconds", c_longdouble),
#                ("InitialXTimeFraction", c_longdouble),
#                ("XIncrement", c_longdouble*1),
#                ("ScaleFactor", c_longdouble*1),
#                ("ScaleOffset", c_longdouble*1)]  
#fetch_params_single = SINGLEFETCHPARAMS()
#t0 = time.time()
#fetchout = driverdll.AgMD1_FetchWaveformInt16(init_status.ViSession,\
#                                                         c_char_p(unicode('Channel2')),\
#                                                         c_int64(numSamples.NumSamples),\
#                                                         fetch_params.WaveformArray,\
#                                                         fetch_params.ActualPoints,\
#                                                         fetch_params.FirstValidPoint,\
#                                                         fetch_params.InitialXOffset,\
#                                                         fetch_params.InitialXTimeSeconds,\
#                                                         fetch_params.InitialXTimeFraction,\
#                                                         byref(fetch_params.XIncrement),\
#                                                         byref(fetch_params.ScaleFactor),\
#                                                         byref(fetch_params.ScaleOffset))
#t1 = time.time()
#print 'time to pull data as ctypes = ' + str(t1-t0)








#print (fetchout)
if fetchout  == 0:
    print ('Data Fetched')
else:
    print ("Fetching error = " + str(fetchout))
    get_error_message(fetchout)


#dataHolder = numpy.zeros(arraySize)
t2 = time.time()
temp = numpy.asarray(fetch_params.WaveformArray)
#temp = numpy.asarray(fetch_params.WaveformArray)*numpy.asarray(fetch_params.ScaleFactor) + numpy.asarray(fetch_params.ScaleOffset)
#temp = numpy.asarray(fetch_params.WaveformArray).astype('double')
#temp = numpy.asarray(fetch_params.WaveformArray).astype('float16')
#dataHolder = numpy.asarray(fetch_params.WaveformArray).astype('double')
t3 = time.time()
print ('time to cast to numpy array = ' + str(t3-t2))
print ('total time = ' + str(t3-t2 + t1-t0 + tb-ta))
print ('number of samples  = ' + str(int(numPointsPerRecordC.value)*int(numRecordsC.value)))
pylab.figure(1)
pylab.clf()
pylab.plot(temp)
pylab.show()
    
    

#def get_realVal(IDnum):
#    getout = driverdll.AqMD3_GetAttributeViReal64(init_status.ViSession, None, c_int32(IDnum), byref(attr_realval) )
#    return attr_realval.AttributeValue

sampRate = get_realVal(1250015)
print ('sample rate = ' + str(sampRate))    
    

#getout = driverdll.AqMD3_GetAttributeViString(init_status.ViSession, None, c_int32(1250001),c_int32(256), byref(attr_stringval) )
##getout = driverdll.AqMD3_GetAttributeViString(init_status.ViSession, c_char_p(b'TriggerSource'), c_int32(1250001), byref(attr_stringval) )
#print ('trigger source = ' + str(attr_stringval.AttributeValue))

#def get_stringVal(IDnum, bufferSize = 256):
#    getout = driverdll.AqMD3_GetAttributeViString(init_status.ViSession, None, c_int32(IDnum),c_int32(bufferSize), byref(attr_stringval) )
#    return attr_stringval.AttributeValue
trigSour = get_stringVal(1250001)
print ('trigger source = ' + str(trigSour))


#def set_stringVal(IDnum, strr):
#    setout = driverdll.AqMD3_SetAttributeViString(init_status.ViSession, None, c_int32(IDnum), strr )
#    return setout

#set_stringVal(1250001, b'Software')    
#trigSour = get_stringVal(1250001)
#print ('trigger source = ' + str(trigSour))  
    



#############abort if it isn't done yet
print ('aborting acquisition')
abortout = driverdll.AqMD3_Abort(init_status.ViSession)
if abortout == 0:
    print ('Aborting acquisition')
else:
    print ("Aborting error = " + str(abortout))
    get_error_message(abortout)





