from ctypes import *
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

#import pickle
#import datetime
#import itertools
#import sys

#os.chdir(r'C:\Users\Kollarlab\Acqiris\testing')

#%% helper functions
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



def get_boolVal(IDnum, recapID = None):
    getout = driverdll.AqMD3_GetAttributeViBoolean(init_status.ViSession, recapID, c_int32(IDnum), byref(attr_boolval) )
    return attr_boolval.AttributeValue

def set_boolVal(IDnum, val, recapID = None):
    setout = driverdll.AqMD3_SetAttributeViBoolean(init_status.ViSession, recapID, c_int32(IDnum), val )
    return setout





def get_int32Val(IDnum):
    getout = driverdll.AqMD3_GetAttributeViInt32(init_status.ViSession, None, c_int32(IDnum), byref(attr_intval) )
    return attr_intval.AttributeValue

def set_int32Val(IDnum, val):
    setout = driverdll.AqMD3_SetAttributeViInt32(init_status.ViSession, None, c_int32(IDnum), val )
    return setout

def get_int64Val(IDnum):
    getout = driverdll.AqMD3_GetAttributeViInt64(init_status.ViSession, None, c_int32(IDnum), byref(attr_intval) )
    return attr_intval.AttributeValue

def set_int64Val(IDnum, val):
    setout = driverdll.AqMD3_SetAttributeViInt64(init_status.ViSession, None, c_int32(IDnum), val )
    return setout



def get_realVal(IDnum):
    getout = driverdll.AqMD3_GetAttributeViReal64(init_status.ViSession, None, c_int32(IDnum), byref(attr_realval) )
    return attr_realval.AttributeValue

def set_realVal(IDnum, recapID = None):
    setout = driverdll.AqMD3_SetAttributeViReal64(init_status.ViSession, recapID, c_int32(IDnum), val )
    return setout




def get_stringVal(IDnum, bufferSize = 256):
    getout = driverdll.AqMD3_GetAttributeViString(init_status.ViSession, None, c_int32(IDnum),c_int32(bufferSize), byref(attr_stringval) )
    return attr_stringval.AttributeValue

def set_stringVal(IDnum, strr):
    setout = driverdll.AqMD3_SetAttributeViString(init_status.ViSession, None, c_int32(IDnum), strr )
    return setout


class INITSTATUS(Structure):
        #    _fields_ = [("ViStatus", c_ulong),
        #                ("ViSession", c_ulong)]    
        #    _fields_ = [("ViSession", c_ulong),
        #                ("ViStatus", c_ulong)]  
            _fields_ = [("ViSession", c_ulong)]  
#init_status = INITSTATUS()  #essentially a global that everyone needs. Don't reset.
    
def initializeCard(hardwareAddressC, strInitOptionsC = c_char_p(b'Simulate=False') ):
    print ('initializing ...')
    ######initialize card  
#    class INITSTATUS(Structure):
#        #    _fields_ = [("ViStatus", c_ulong),
#        #                ("ViSession", c_ulong)]    
#        #    _fields_ = [("ViSession", c_ulong),
#        #                ("ViStatus", c_ulong)]  
#            _fields_ = [("ViSession", c_ulong)]  
    global init_status 
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

    return

def check_channel_status(chanNum):
    if chanNum == 1:
        chanNameC = b'Channel1'
    elif chanNum == 2:
        chanNameC = b'Channel2'
    else:
        raise ValueError('Invalid channel number. Should be 1 or 2.')
    
    IDnum = 1250002
    
    val = get_boolVal(IDnum, recapID = chanNameC)
    
    return val

def get_card_options():
    IDnum = 1150023 #AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS
    val = get_stringVal(IDnum, bufferSize = 532)
    return val

##############################



#%% Initialization of card
print ('##############################')
print ('##############################')

driverdll = WinDLL(r'AqMD3_64.dll')
#driverdll = WinDLL(r'AqMD2_64.dll')


#hardwareAddress = u'PXI23::0::0::INSTR'
#hardwareAddressC = c_char_p(bytes(hardwareAddress))
hardwareAddressC = c_char_p(b'PXI23::0::0::INSTR') #needs to e a byte string for python 3. Used to be regular string would work.
#hardwareAddressC = c_char_p(b'PXI21::0::0::INSTR') #junk address


# strInitOptionsC = c_char_p(b'Simulate=False, DriverSetup= CAL=0, Trace=false, model = SA220P')
#strInitOptionsC = c_char_p(b'Simulate=True, DriverSetup= CAL=0, Trace=false, model = SA220P')

#strInitOptions = 'Simulate=False'
#strInitOptionsC = c_char_p(bytes(strInitOptions))


#strInitOptionsC = c_char_p(b'Simulate=False')
strInitOptionsC = c_char_p(b'Simulate=False,  DriverSetup= model = SA220P')
#strInitOptionsC = c_char_p(b'Simulate=True,  DriverSetup= model = SA220P')
#strInitOptionsC = c_char_p(b'Simulate=True, DriverSetup= model = SA220P')

idQueryC = c_bool(True)
resetC   = c_bool(True)
#resetC   = c_bool(False)

###############################


###try to auto do the initialize.
###fake call to see if the card is happy
try:
    getout = driverdll.AqMD3_GetAttributeViBoolean(init_status.ViSession, None, c_int32(1050005), byref(attr_boolval) )
    
    if getout == -1074134951:  #error = -1074134951 'Unexpected responce from the instrument' - needs full reinit.
        get_error_message(getout)
        print ('card is not happy. Reinitializing.')
        initializeCard(hardwareAddressC, strInitOptionsC)
    
except:
    print ('couldnt call. Initializing')
    initializeCard(hardwareAddressC, strInitOptionsC)



###mem error = -1074134951 'Unexpected responce from the instrument' - needs full reinit.
#init = True
##init = False
#if init:
#    initializeCard(hardwareAddressC, strInitOptionsC)
#    
#    
##    print ('initializing ...')
##    ######initialize card  
##    class INITSTATUS(Structure):
##        #    _fields_ = [("ViStatus", c_ulong),
##        #                ("ViSession", c_ulong)]    
##        #    _fields_ = [("ViSession", c_ulong),
##        #                ("ViStatus", c_ulong)]  
##            _fields_ = [("ViSession", c_ulong)]  
##    init_status = INITSTATUS()  
##    initout = driverdll.AqMD3_InitWithOptions(hardwareAddressC,idQueryC, resetC , strInitOptionsC,  byref(init_status))            
##    if initout == 0:
##        print('Driver Initialized')
##    else:
##        print("Drivier initialization error = " + str(initout))
##        get_error_message(initout)
##    
##    print ('calibrating ...')    
##    calout = driverdll.AqMD3_SelfCalibrate(init_status.ViSession)            
##    if calout == 0:
##        print('Calibrated')
##    else:
##        print("Calibration error = " + str(calout))
##        get_error_message(calout)
#    
    

 
a = 1
print ('##############################')
print ('   ')
print ('   ')       

       
    
#%% configuring channels and calibrating
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


activeChannels = numpy.asarray([1,2])
#activeChannels = numpy.asarray([1]) #lloks like 1 channel mode has to be channel1?

chanOuts = numpy.zeros(2)
for ind in [1,2]:

    if ind in activeChannels:
        enabledC = c_bool(True) #override default off state
    else:
        enabledC = c_bool(False)
    
    if ind == 1:
        chanStr = b'Channel1'
        
    elif ind == 2:
        chanStr = b'Channel2'
        
    chanC = c_char_p(chanStr)
    chanout = driverdll.AqMD3_ConfigureChannel(init_status.ViSession, chanC,rangeC,offsetC, couplingC , enabledC) 
    
    chanOuts[ind-1] = chanout
if chanOuts[0] ==0 and chanOuts[1] == 0:
    print ('Channels Initialized')
else:
    print ("Channel error = " + str(chanOuts[0] ) + ' or ' + str(chanOuts[1]))
    get_error_message(int(chanOuts[0]))
    get_error_message(int(chanOuts[1]))

#chanStr = b'Channel1'
#chanC = c_char_p(chanStr)
#chanout1 = driverdll.AqMD3_ConfigureChannel(init_status.ViSession, chanC,rangeC,offsetC, couplingC , c_bool(True))  
#
#
#chanStr = b'Channel2'
#chanC = c_char_p(chanStr)
#chanout2 = driverdll.AqMD3_ConfigureChannel(init_status.ViSession, chanC,rangeC,offsetC, couplingC , c_bool(True))  
#if chanout1 ==0 and chanout2 == 0:
#    print ('Channels Initialized')
#else:
#    print ("Channel error = " + str(chanout1) + ' or ' + str(chanout2))
#    get_error_message(chanout1)
#    get_error_message(chanout2)
#    
print ('##############################')
print ('   ')
print ('   ')    



##is it a real boy?
#getout = driverdll.AqMD3_GetAttributeViBoolean(init_status.ViSession, None, c_int32(1050005), byref(attr_boolval) )
#print ('simulated = ' + str(attr_boolval.AttributeValue))




############## calibration seems to be needed before taking data?
#print ('calibrating ...')    
#calout = driverdll.AqMD3_SelfCalibrate(init_status.ViSession)            
#if calout == 0:
#    print('Calibrated')
#else:
#    print("Calibration error = " + str(calout))
#    get_error_message(calout)


#%% setting the acquisition and the trigger
#############set the acquisition
#numPointsPerRecord = c_int64(10000000)
#numRecords = c_int64(1)
numPointsPerRecordC = c_int64(1*10**7)
#numRecordsC = c_int64(2)
sampleRateC = c_longdouble(1*10**9) #not currently being used correctly

#averageMode = True
averageMode = False
if averageMode:
    averageModeC = c_int32(1)
    averagesC = c_int32(80)
    numRecordsC = c_int64(1)
    
#    numRecordsC = c_int64(40)
#    testout = set_int64Val(1250013, numRecordsC)
else:
    averageModeC = c_int32(0)
    averagesC = c_int32(0)
    numRecordsC = c_int64(2)


#this doesn't seem to do anything
applyout = driverdll.AqMD3_ApplySetup(init_status.ViSession)
if applyout == 0:
    print ('setup applied')
else:
    print ("setup application error = " + str(applyout))
    get_error_message(applyout)
    


#st whether it should be averaging. 
##setting ACQUISITION_MODE
#avgmodesetout = set_int32Val(1150011, averageModeC) #1 for average. 0 for not.  4 might be interesting
#avgmodesetout =0
if averageMode:
    avgnumsetout = set_int32Val(1150069, averagesC)
else:
    avgnumsetout = set_int32Val(1150069, c_int32(4)) #4 is the minimum
    
numAvgs_gotten = get_int32Val(1150069)

#setting ACQUISITION_MODE
avgmodesetout = set_int32Val(1150011, averageModeC) #1 for average. 0 for not.  4 might be interesting    
    
    
if (avgmodesetout ==0) & (avgnumsetout == 0):
    print ('averaging succesfully configured or turned off')
else:
    print ('error configuring averaging = ' + str(avgmodesetout) + ' and ' + str(avgnumsetout))
    get_error_message(avgmodesetout)
    get_error_message(avgnumsetout)



aconfigout = driverdll.AqMD3_ConfigureAcquisition(init_status.ViSession, numRecordsC, numPointsPerRecordC, sampleRateC)
if aconfigout == 0:
    print ('acquisition set up')
else:
    print ("acquisition setup error = " + str(aconfigout))
    get_error_message(aconfigout)
        


##st whether it should be averaging. 
##setting ACQUISITION_MODE
##avgmodesetout = set_int32Val(1150011, averageModeC) #1 for average. 0 for not.  4 might be interesting
#avgmodesetout =0
#if averageMode:
#    avgnumsetout = set_int32Val(1150069, averagesC)
#else:
#    avgnumsetout = set_int32Val(1150069, c_int32(4)) #4 is the minimum
#if (avgmodesetout ==0) & (avgnumsetout == 0):
#    print ('averaging succesfully configured or turned off')
#else:
#    print ('error configuring averaging = ' + str(avgmodesetout) + ' and ' + str(avgnumsetout))
#    get_error_message(avgmodesetout)
#    get_error_message(avgnumsetout)
    
    
    
    
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




print ('##############################')
print ('   ')
print ('   ')    


    
 
#%% setting up memory  
 
###############get ready to allocte memeory
#class MINMEMQUERRY(Structure):
##    _fields_ = [("NumSamples", c_ulong)]  
##    _fields_ = [("NumSamples", c_int64)]  
##    _fields_ = [("NumSamples", c_int64),
##                ("ViStatus", c_ulong)]  
#    _fields_ = [("NumSamples", c_int64)]  
#numSamples =  MINMEMQUERRY() 
#memout = driverdll.AqMD3_QueryMinWaveformMemory(init_status.ViSession, c_int32(64), numRecordsC, c_int64(0), numPointsPerRecordC, byref(numSamples))
##driverdll.AgMD1_QueryMinWaveformMemory(byref(numSamples), c_int64(64), numRecords, c_int64(64), numPointsPerRecord)
#if memout == 0:
#    print ('Ran min memory function')
#else:
#    print ("mem error = " + str(memout))
#    get_error_message(memout)
#
#print (numSamples.NumSamples)




############## calibration seems to be needed before taking data?
#print ('calibrating ...')    
#calout = driverdll.AqMD3_SelfCalibrate(init_status.ViSession)            
#if calout == 0:
#    print('Calibrated')
#else:
#    print("Calibration error = " + str(calout))
#    get_error_message(calout)

###############prep for acquisition
#class MEASSTATUS(Structure):
#    _fields_ = [("Measuring" , c_int32)]    
#meas_status = MEASSTATUS() 


#measout = driverdll.AqMD3_IsMeasuring(init_status.ViSession,byref(meas_status))
#if measout == 0:
#    print ('querried measureent status')
#else:
#    print ("error querrying measurment status = " + str(measout))
#    get_error_message(measout)
   



#%% acquire data, hopefully    
#################initiate acquisition
class MEASSTATUS(Structure):
    _fields_ = [("Measuring" , c_int32)]    
meas_status = MEASSTATUS() 



print ('starting acquisition')
acqout = driverdll.AqMD3_InitiateAcquisition(init_status.ViSession)
if acqout  == 0:
    print ('Initialized acquisition')
else:
    print ("arming error = " + str(acqout))
    get_error_message(acqout)
    
    #trying to recover from its sometimes need for recalibration
    abortout = driverdll.AqMD3_Abort(init_status.ViSession)
    
    #error code -1074118652 = needs calibration
    if acqout == -1074118647:
        print ('Cabilbration needed. Trying to recalbirate')
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




print ('##############################')
print ('   ')
print ('   ')    




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
    timeouttime = 10*1000 #ms
    print ('waiting for measurement to finish')
    waitout = driverdll.AqMD3_WaitForAcquisitionComplete(init_status.ViSession,c_int32(timeouttime) )
    if waitout == 0:
        print ('successfully done waiting')
    else:
        print ("waiting error = " + str(waitout))
        get_error_message(waitout)
#    dt = 0.01
#    maxT = 5
#    for ind in range(0,int(numpy.ceil(maxT/dt))):
#        print (str(ind))
#        measout = driverdll.AqMD3_IsMeasuring(init_status.ViSession,byref(meas_status))
#        flag = meas_status.Measuring
#        if flag == 1:
#            print ('Still Measureing')
#            time.sleep(0.01)
#        else:
#            print ('Done Measuring')
#            break
    

print ('##############################')
print ('   ')
print ('   ')    





#%% querry the memory
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
    print ('Ran min memory function')
else:
    print ("mem error = " + str(memout))
    get_error_message(memout)

print (numSamples.NumSamples)



print ('##############################')
print ('   ')
print ('   ')   

#%% fetch data    

 
#temp2 = int(numRecordsC.value)*1+ 1
##temp2 = int(numRecords.value)*1
#samplesToTake = arraySize *temp2 # ignore the correct memory function
#arraySize = int(numRecordsC.value*numPointsPerRecordC.value) #forced it to take the right amount of data for now.

arraySize = int(numSamples.NumSamples)
samplesToTake = numSamples.NumSamples
segmentsNum = int(numRecordsC.value)*1


#class FETCHPARAMS(Structure):
#    _fields_ = [("WaveformArray", c_longdouble*arraySize),
#                ("WaveformArrayActualSize", c_int64*1),
#                ("ActualAverages", c_int32*1),
#                ("ActualRecords", c_int64*1),
#                ("ActualPoints", c_int64*segmentsNum),
#                ("FirstValidPoint", c_int64*segmentsNum),
#                ("InitialXOffset", c_longdouble*segmentsNum),
#                ("InitialXTimeSeconds", c_longdouble*segmentsNum),
#                ("InitialXTimeFraction", c_longdouble*segmentsNum),
#                ("XIncrement", c_longdouble*1),
#                ("ScaleFactor", c_longdouble*1),
#                ("ScaleOffset", c_longdouble*1)]  
    
class FETCHPARAMS(Structure):
    _fields_ = [("WaveformArray", c_longdouble*arraySize),
                ("ActualRecords", c_int64*1),
                ("ActualPoints", c_int64*segmentsNum),
                ("FirstValidPoint", c_int64*segmentsNum),
                ("InitialXOffset", c_longdouble*segmentsNum),
                ("InitialXTimeSeconds", c_longdouble*segmentsNum),
                ("InitialXTimeFraction", c_longdouble*segmentsNum),
                ("XIncrement", c_longdouble*1),
                ("ScaleFactor", c_longdouble*1),
                ("ScaleOffset", c_longdouble*1)]  
    
class ACCUMULATEDFETCHPARAMS(Structure):
    _fields_ = [("WaveformArray", c_longdouble*arraySize),
                ("ActualAverages", c_int32*1),
                ("ActualRecords", c_int64*1),
                ("ActualPoints", c_int64*segmentsNum),
                ("FirstValidPoint", c_int64*segmentsNum),
                ("InitialXOffset", c_longdouble*segmentsNum),
                ("InitialXTimeSeconds", c_longdouble*segmentsNum),
                ("InitialXTimeFraction", c_longdouble*segmentsNum),
                ("XIncrement", c_longdouble*1),
                ("Flags", c_int32*segmentsNum)]  
    
ta = time.time()
if averageMode:
    fetch_params = ACCUMULATEDFETCHPARAMS()  
else:
    fetch_params = FETCHPARAMS()
tb = time.time()
print ('preallocation time = ' + str(tb-ta))


t0 = time.time()


chanNameC  = c_char_p(b'Channel1')
firstRecordC = c_int64(0)
#numRecordsC
OffsetWithinRecordC = c_int64(0)
#numPointsPerRecordC
WaveformArrayPlannedSizeC = c_int64(samplesToTake)
#byref waveformarray

# class WAVEFORMHOLDER(Structure):   
#     _fields_ = [("WaveformArrayC", c_longdouble*arraySize)]
# WaveHolder = WAVEFORMHOLDER()   # This DOES WORK! It can be sent in.
# #Successfully returns data to "WaveHolder.WaveformArrayC"


# WaveformArrayC= c_longdouble*arraySize ## this syntax doesn't work. It can't be cast.

if averageMode:
    fetchout = driverdll.AqMD3_FetchAccumulatedWaveformReal64(init_status.ViSession,\
                                                             chanNameC,\
                                                             firstRecordC,\
                                                             numRecordsC,\
                                                             OffsetWithinRecordC,\
                                                             numPointsPerRecordC,\
                                                             WaveformArrayPlannedSizeC,\
                                                             fetch_params.WaveformArray,\
                                                             byref(fetch_params.ActualAverages),\
                                                             byref(fetch_params.ActualRecords),\
                                                             fetch_params.ActualPoints,\
                                                             fetch_params.FirstValidPoint,\
                                                             byref(fetch_params.InitialXOffset),\
                                                             fetch_params.InitialXTimeSeconds,\
                                                             fetch_params.InitialXTimeFraction,\
                                                             byref(fetch_params.XIncrement),\
                                                             fetch_params.Flags)
else:
    # fetchout = driverdll.AqMD3_FetchMultiRecordWaveformReal64(init_status.ViSession,\
    #                                                          chanNameC,\
    #                                                          firstRecordC,\
    #                                                          numRecordsC,\
    #                                                          OffsetWithinRecordC,\
    #                                                          numPointsPerRecordC,\
    #                                                          WaveformArrayPlannedSizeC,\
    #                                                          WaveHolder.WaveformArrayC,\
    #                                                          byref(fetch_params.ActualRecords),\
    #                                                          fetch_params.ActualPoints,\
    #                                                          fetch_params.FirstValidPoint,\
    #                                                          fetch_params.InitialXOffset,\
    #                                                          fetch_params.InitialXTimeSeconds,\
    #                                                          fetch_params.InitialXTimeFraction,\
    #                                                          byref(fetch_params.XIncrement))
    fetchout = driverdll.AqMD3_FetchMultiRecordWaveformReal64(init_status.ViSession,\
                                                             chanNameC,\
                                                             firstRecordC,\
                                                             numRecordsC,\
                                                             OffsetWithinRecordC,\
                                                             numPointsPerRecordC,\
                                                             WaveformArrayPlannedSizeC,\
                                                             fetch_params.WaveformArray,\
                                                             byref(fetch_params.ActualRecords),\
                                                             fetch_params.ActualPoints,\
                                                             fetch_params.FirstValidPoint,\
                                                             fetch_params.InitialXOffset,\
                                                             fetch_params.InitialXTimeSeconds,\
                                                             fetch_params.InitialXTimeFraction,\
                                                             byref(fetch_params.XIncrement))



t1 = time.time()
print ('time to pull data as ctypes = ' + str(t1-t0))


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



#%% diagnostics

#dataHolder = numpy.zeros(arraySize)
t2 = time.time()
rawData = numpy.asarray(fetch_params.WaveformArray)
# rawData = numpy.asarray(WaveHolder.WaveformArrayC)

#temp = numpy.asarray(fetch_params.WaveformArray)*numpy.asarray(fetch_params.ScaleFactor) + numpy.asarray(fetch_params.ScaleOffset)
#temp = numpy.asarray(fetch_params.WaveformArray).astype('double')
#temp = numpy.asarray(fetch_params.WaveformArray).astype('float16')
#dataHolder = numpy.asarray(fetch_params.WaveformArray).astype('double')
t3 = time.time()
print ('time to cast to numpy array = ' + str(t3-t2))
print ('total time = ' + str(t3-t2 + t1-t0 + tb-ta))
print ('number of samples  = ' + str(int(numPointsPerRecordC.value)*int(numRecordsC.value)))


sampRate = get_realVal(1250015)
print ('sample rate = ' + str(sampRate))    


#
#dataActualSizeIsh = numpy.asarray(fetch_params.WaveformArrayActualSize) #set to zero and never modified. Nevermind.
dataRawSize = rawData.size
dataActualSegments = numpy.asarray(fetch_params.ActualRecords)[0]
dataActualPoints_full = numpy.asarray(fetch_params.ActualPoints)
dataActualPoints = dataActualPoints_full[0]
dataFirstValidPoints = numpy.asarray(fetch_params.FirstValidPoint).astype('int64')


data = numpy.zeros((dataActualSegments,dataActualPoints))
for segind in range(0,dataActualSegments ):
    startInd = dataFirstValidPoints[segind]
    
    data[segind,:] = rawData[startInd:(startInd+dataActualPoints)]


    
numSamples = dataActualPoints
dt = 1/sampRate
xaxis = scipy.arange(0, numSamples,1)*dt

pylab.figure(1)
pylab.clf()
ax = pylab.subplot(1,1,1)
for segind in  range(0,dataActualSegments ):
    labelstr = 'seg: ' + str(segind)
    pylab.plot(xaxis*1000, data[segind,:] + segind*0.125, label = labelstr)
ax.legend(loc = 'upper right')
pylab.ylabel('Voltage (waterfall)')
pylab.xlabel('Time (ms)')
pylab.title('Active Channel')
pylab.show()
#    
    

print ('##############################')
print ('   ')
print ('   ')    






#def get_realVal(IDnum):
#    getout = driverdll.AqMD3_GetAttributeViReal64(init_status.ViSession, None, c_int32(IDnum), byref(attr_realval) )
#    return attr_realval.AttributeValue

#sampRate = get_realVal(1250015)
#print ('sample rate = ' + str(sampRate))    
#    

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
    


#%% abort
#############abort if it isn't done yet

#abort = True
abort = False

if abort:
    print ('aborting acquisition')
    abortout = driverdll.AqMD3_Abort(init_status.ViSession)
    if abortout == 0:
        print ('Aborting acquisition')
    else:
        print ("Aborting error = " + str(abortout))
        get_error_message(abortout)
    
    
    print ('##############################')
    print ('##############################')      
    print ('   ')
    print ('   ')    



#print ('aborting acquisition')
#abortout = driverdll.AqMD3_Abort(init_status.ViSession)
#if abortout == 0:
#    print ('Aborting acquisition')
#else:
#    print ("Aborting error = " + str(abortout))
#    get_error_message(abortout)
#
#
#print ('##############################')
#print ('##############################')      
#print ('   ')
#print ('   ')

