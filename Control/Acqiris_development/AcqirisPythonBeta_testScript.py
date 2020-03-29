# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 16:44:14 2020

@author: Kollarlab
"""

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


#hardwareAddress = "PXI23::0::0::INSTR"
#
IVIbinPath = "C:\\Program Files\\IVI Foundation\\IVI\\Bin\\"
sys.path.append(IVIbinPath)

###############
#modified from:
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

from AqMD3 import *

# Edit resource and options as needed. Resource is ignored if option Simulate=True.
resourceString = "PXI23::0::0::INSTR"
initOptions = "Simulate=True, DriverSetup= Model=U5303A"
idquery = False
reset = False

try:
    print("IVI-Python SimpleAcquisition")
    print()

    with AqMD3( resourceString, False, False, initOptions ) as driver:

        print("Driver initialized")

        # Print a few IIviDriverIdentity properties.
        #print("Driver identifier:  ", driver.Identity.Identifier)
        #print("Driver revision:    ", driver.Identity.Revision)
        #print("Driver vendor:      ", driver.Identity.Vendor)
        #print("Driver description: ", driver.Identity.Description)
        print("Instrument manufact:", driver.Identity.InstrumentManufacturer)
        print("Instrument model:   ", driver.Identity.InstrumentModel)
        print("Firmware revision:  ", driver.Identity.InstrumentFirmwareRevision)
        print("Serial number:      ", driver.InstrumentInfo.SerialNumberString)
        print("Options:            ", driver.InstrumentInfo.Options)
        print("Simulate:           ", driver.DriverOperation.Simulate)
        print()

        # Configure channel properties.
        range = 1.0
        offset = 0.0
        coupling = VerticalCoupling.DC
        print("Configuring channel properties")
        print("Range:              ", range)
        print("Offset:             ", offset)
        print("Coupling:           ", coupling)
        for channel in driver.Channels:
            print("Applying on ", channel.Name)
            channel.Configure(range, offset, coupling, True)

        # Configure the acquisition.
        numPointsPerRecord = 1000000

        print()
        print("Configuring acquisition")
        print("Record size:        ", numPointsPerRecord)
        driver.Acquisition.RecordSize = numPointsPerRecord

        # Configure the trigger.
        sourceName = "Internal1"
        level = 0.0
        slope = TriggerSlope.Positive

        print()
        print("Configuring trigger")
        print("Active source:      ", sourceName)
        driver.Trigger.ActiveSource = sourceName
        activeTrigger = driver.Trigger.Sources[sourceName]
        print("Level:              ", level)
        activeTrigger.Level = level
        print("Slope:              ", slope)
        activeTrigger.Edge.Slope = slope

        # Calibrate the instrument.
        print()
        print("Performing self-calibration")
        driver.Calibration.SelfCalibrate()

        # Perform the acquisition.
        print()
        print("Performing acquisition")
        driver.Acquisition.Initiate()
        timeoutInMs = 1000
        driver.Acquisition.WaitForAcquisitionComplete(timeoutInMs)
        print("Acquisition completed")

        waveform = None
        for channel in driver.Channels:
            print()
            print("Fetching data from ", channel.Name)

            # Fetch the acquired data
            waveform = channel.Measurement.FetchWaveform(waveform)
            print("First samples:", *waveform[:10])

            print("Processing data fetched from ", channel.Name)
            # Convert data to Volts.
            for sampleRaw in waveform:
                # the same sample in volts
                sampleInVolts = sampleRaw * waveform.ScaleFactor + waveform.ScaleOffset

        print("Processing completed. ")

    # Automaticaly close the driver at end of with clause.
    print("Driver closed")

except RuntimeError as e:
    print( e )

print("\nDone - Press enter to exit")
print()


