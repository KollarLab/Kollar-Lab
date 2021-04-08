#
# Acqiris IVI-Python AqMD3 Driver Example Program
#
# Creates a driver object, reads a few Identity interface properties, and performs a simple
# acquisition.
#
# See http://www.ivifoundation.org/resources/ for additional programming information.
#
# Runs in simulation mode without an instrument.
#
# Requires a working installation of the IVI-C driver.
#
import time
import scipy
import sys
hardwareAddress = "PXI23::0::0::INSTR"
IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
if not IVIbinPath in sys.path:
    sys.path.append(IVIbinPath)


import numpy as np
from AqMD3 import *

# Edit resource and options as needed. Resource is ignored if option Simulate=True.
resourceString = hardwareAddress
initOptions = "Simulate=false, DriverSetup= Model=SA220P"
idquery = False
reset = False

rinds = scipy.arange(0,10,1.)

try:
    print("IVI-Python SimpleAcquisition")
    print()

    with AqMD3( resourceString, False, False, initOptions ) as driver:
        print("Driver initialized")
        for rind in rinds:
            print()
            print('Round : ' + str(rind))
    
#            # Print a few IIviDriverIdentity properties.
#            #print("Driver identifier:  ", driver.Identity.Identifier)
#            #print("Driver revision:    ", driver.Identity.Revision)
#            #print("Driver vendor:      ", driver.Identity.Vendor)
#            #print("Driver description: ", driver.Identity.Description)
#            print("Instrument manufact:", driver.Identity.InstrumentManufacturer)
#            print("Instrument model:   ", driver.Identity.InstrumentModel)
#            print("Firmware revision:  ", driver.Identity.InstrumentFirmwareRevision)
#            print("Serial number:      ", driver.InstrumentInfo.SerialNumberString)
#            print("Options:            ", driver.InstrumentInfo.Options)
#            print("Simulate:           ", driver.DriverOperation.Simulate)
#            print()
    
            # Configure channel properties.
            range = 2.5
            offset = 0.0
            coupling = VerticalCoupling.DC
            print("Configuring channel properties")
#            print("Range:              ", range)
#            print("Offset:             ", offset)
#            print("Coupling:           ", coupling)
            for channel in driver.Channels:
#                print("Applying on ", channel.Name)
                channel.Configure(range, offset, coupling, True)
    
            # Configure the acquisition.
#            print()
            print("Configuring acquisition")
            numRecords = 1
            numPointsPerRecord = 1601
            numAverages = 8
            driver.Acquisition.ConfigureAcquisition(numRecords, numPointsPerRecord, 2.0e9)
#            print("Number of records:  ", numRecords)
#            print("Number of Averages: ", numAverages)
#            print("Record size:        ", numPointsPerRecord)
            driver.Acquisition.NumberOfRecordsToAcquire = numRecords;
            driver.Acquisition.RecordSize = numPointsPerRecord;
            driver.Acquisition.Mode = AcquisitionMode.Averager;
            driver.Acquisition.NumberOfAverages = numAverages;
            driver.Acquisition.Mode = AcquisitionMode.Averager;
          
            # Configure the trigger.
            sourceName = "External1"
            level = 0.0
            slope = TriggerSlope.Positive
    
#            print()
            print("Configuring trigger")
#            print("Active source:      ", sourceName)
            driver.Trigger.ActiveSource = sourceName
            activeTrigger = driver.Trigger.Sources[sourceName]
#            print("Level:              ", level)
            activeTrigger.Level = level
#            print("Slope:              ", slope)
            activeTrigger.Edge.Slope = slope
    
            # Calibrate the instrument.
#            print()
            if rind == 0:
                print("Performing self-calibration")
                t0 = time.time()
                driver.Calibration.SelfCalibrate()
                t1 = time.time()
                print('Calibration Time : ' + str(t1-t0))
#            print("Performing self-calibration")
#            t0 = time.time()
#            driver.Calibration.SelfCalibrate()
#            t1 = time.time()
#            print('Calibration Time : ' + str(t1-t0))
    
            # Perform the acquisition.
#            print()
            print("Performing acquisition")
            driver.Acquisition.Initiate()
            timeoutInMs = 1000
            driver.Acquisition.WaitForAcquisitionComplete(timeoutInMs)
            print("Acquisition completed")
    
            firstRecord = 0
            offsetWithinRecord = 0
            waveforms = None
            for channel in driver.Channels:
#                print("Fetching data from ", channel.Name)
                
                # Fetch the acquired data
                waveforms = channel.Measurement.FetchAccumulatedWaveform(firstRecord, numRecords, offsetWithinRecord, numPointsPerRecord, dtype=np.int32)
#                print("Actual averages: ", waveforms.ActualAverages)
#                print("Processing data fetched from ", channel.Name)
                
                # Convert data to Volts.
                for waveform in waveforms:        
                    for sampleIndex, sampleRaw in enumerate(waveform):
                        sampleInVolts = sampleRaw * waveform.ScaleFactor + waveform.ScaleOffset
                        if sampleIndex<10:
#                            print(sampleRaw, sampleInVolts)
                            pass
    
    
            print("Processing completed. ")

    # Automaticaly close the driver at end of with clause.
    print("Driver closed")

except RuntimeError as e:
    print( e )

print("\nDone - Press enter to exit")
print()

