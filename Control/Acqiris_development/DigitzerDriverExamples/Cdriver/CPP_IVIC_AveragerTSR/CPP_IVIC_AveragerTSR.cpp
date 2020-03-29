///
/// Acqiris IVI-C Driver Example Program
///
/// Initializes the driver, reads a few Identity interface properties, and performs
/// accumulated acquisitions using TSR.
///
/// For additional information on programming with IVI drivers in various IDEs, please see
/// http://www.ivifoundation.org/resources/
///
/// WARNING:
/// The Averager TSR features are not supported in simulation mode. You will have to update
/// the resource string (resource[]) to match your configuration and disable the
/// simulation mode (options[] - set Simulate=false) to be able to run this example
/// successfully.
///

#include "AqMD3.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::hex;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <stdexcept>
using std::runtime_error;


// Edit resource and options as needed. Resource is ignored if option has Simulate=true.
// An input signal is necessary if the example is run in non simulated mode, otherwise
// the acquisition will time out.
ViChar resource[] = "PXI40::0::0::INSTR";
ViChar options[] = "Simulate=true, DriverSetup= Model=U5303A";


ViInt64 const numRecords = 1;
ViInt32 const numAverages = 80;


#define checkApiCall( f ) do { ViStatus s = f; testApiCall( s, #f ); } while( false )


// Utility function to check status error during driver API call.
void testApiCall(ViStatus status, char const * functionName)
{
    ViInt32 ErrorCode;
    ViChar ErrorMessage[256];

    if (status>0) // Warning occurred.
    {
        AqMD3_GetError(VI_NULL, &ErrorCode, sizeof(ErrorMessage), ErrorMessage);
        cerr << "** Warning during " << functionName << ": 0x" << hex << ErrorCode << ", " << ErrorMessage << "\n";

    }
    else if (status<0) // Error occurred.
    {
        AqMD3_GetError(VI_NULL, &ErrorCode, sizeof(ErrorMessage), ErrorMessage);
        cerr << "** ERROR during " << functionName << ": 0x" << hex << ErrorCode << ", " << ErrorMessage << "\n";
        throw runtime_error(ErrorMessage);
    }
}

int main()
{
    cout << "\nStarting Averager+TSR\n";

    // Initialize the driver. See driver help topic "Initializing the IVI-C Driver" for additional information.
    ViSession session;
    ViBoolean idQuery = VI_FALSE;
    ViBoolean reset   = VI_FALSE;
    checkApiCall( AqMD3_InitWithOptions( resource, idQuery, reset, options, &session ) );

    cout << "\nDriver initialized\n";

    // Abort execution if instrument is still in simulated mode.
    ViBoolean simulate;
    checkApiCall(AqMD3_GetAttributeViBoolean(session, "", AQMD3_ATTR_SIMULATE, &simulate));
    if (simulate == VI_TRUE)
    {
        cout << "\nThe Averager features are not supported in simulated mode.\n";
        cout << "Please update the resource string (resource[]) to match your configuration,";
        cout << " and update the init options string (options[]) to disable simulation.\n";

        AqMD3_close(session);

        return 1;
    }
    cout << "\nSimulate:           " << (simulate == VI_TRUE ? "true" : "false") << "\n";

    // Check the instrument contains the required AVG module option.
    ViChar str[128] = { '\0' };
    checkApiCall(AqMD3_GetAttributeViString(session, "", AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS, sizeof(str), str));
    if (string(str).find("AVG") == string::npos || string(str).find("TSR") == string::npos)
    {
        cout << "Both AVG and TSR module options are required on the instrument.\n";

        AqMD3_close(session);

        return 1;
    }


    // Read and output a few attributes.
    checkApiCall(AqMD3_GetAttributeViString(session, "", AQMD3_ATTR_SPECIFIC_DRIVER_PREFIX, sizeof(str), str));
    cout << "Driver prefix:      " << str << "\n";
    checkApiCall(AqMD3_GetAttributeViString(session, "", AQMD3_ATTR_SPECIFIC_DRIVER_REVISION, sizeof(str), str));
    cout << "Driver revision:    " << str << "\n";
    checkApiCall(AqMD3_GetAttributeViString(session, "", AQMD3_ATTR_SPECIFIC_DRIVER_VENDOR, sizeof(str), str));
    cout << "Driver vendor:      " << str << "\n";
    checkApiCall(AqMD3_GetAttributeViString(session, "", AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION, sizeof(str), str));
    cout << "Driver description: " << str << "\n";
    checkApiCall(AqMD3_GetAttributeViString(session, "", AQMD3_ATTR_INSTRUMENT_MODEL, sizeof(str), str));
    cout << "Instrument model:   " << str << "\n";
    string const instrModel = str;
    checkApiCall(AqMD3_GetAttributeViString(session, "", AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS, sizeof(str), str));
    cout << "Instrument options: " << str << '\n';
    checkApiCall(AqMD3_GetAttributeViString(session, "", AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION, sizeof(str), str));
    cout << "Firmware revision:  " << str << "\n";
    checkApiCall(AqMD3_GetAttributeViString(session, "", AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING, sizeof(str), str));
    cout << "Serial number:      " << str << "\n";

    // Configure the acquisition.
    cout << "\nConfiguring acquisition\n";
	ViInt64 const recordSize = 1600;
    ViInt64 const numRecords = 8;
    ViInt32 const numAverages = 8;
    ViReal64 const range = 1.0;
    ViReal64 const offset = 0.0;
    ViInt32 const coupling = AQMD3_VAL_VERTICAL_COUPLING_DC;
    cout << "Range:              " << range << "\n";
    cout << "Offset:             " << offset << "\n";
    cout << "Coupling:           " << (coupling ? "DC" : "AC") << "\n";
    checkApiCall(AqMD3_ConfigureChannel(session, "Channel1", range, offset, coupling, VI_TRUE));
    cout << "Record size:        " << recordSize << "\n";
    checkApiCall(AqMD3_SetAttributeViInt64(session, "", AQMD3_ATTR_RECORD_SIZE, recordSize));
    cout << "Number of averages: " << numAverages << "\n";
    checkApiCall(AqMD3_SetAttributeViInt32(session, "", AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, numAverages));
    // Have to enable time-interleaving in order to enable TSR with AVG
    if(instrModel == "U5303A")
        checkApiCall(AqMD3_SetAttributeViString(session, "Channel1", AQMD3_ATTR_TIME_INTERLEAVED_CHANNEL_LIST, "Channel2"));
    cout << "Mode:               Averager\n";
    checkApiCall(AqMD3_SetAttributeViInt32(session, "", AQMD3_ATTR_ACQUISITION_MODE, AQMD3_VAL_ACQUISITION_MODE_AVERAGER));
    cout << "TSR:                Enabled\n";
    checkApiCall(AqMD3_SetAttributeViBoolean(session, "", AQMD3_ATTR_TSR_ENABLED, VI_TRUE));
    cout << "Number of records:  " << numRecords << "\n";
    checkApiCall(AqMD3_SetAttributeViInt64(session, "", AQMD3_ATTR_NUM_RECORDS_TO_ACQUIRE, numRecords));
    checkApiCall(AqMD3_ApplySetup(session));

    // Configure the trigger.
    cout << "\nConfiguring trigger (self-trigger)\n";
	double squareWaveFrequency = 10e3; // 10 KHz
    double squareWaveDutyCycle = 10.2; // 10.2% of the period - maximum autorized with frequency of 10 KHz
    ViInt32 squareWaveSlope = AQMD3_VAL_TRIGGER_SLOPE_POSITIVE;
    checkApiCall(AqMD3_SetAttributeViString(session, "", AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, "SelfTrigger"));
    checkApiCall(AqMD3_SetAttributeViString(session, "ControlIO3", AQMD3_ATTR_CONTROL_IO_SIGNAL, "Out-AveragerAwg")); // Shunt self-trigger signal to the control IO 3 output.
    checkApiCall(AqMD3_SetAttributeViInt32(session, "SelfTrigger", AQMD3_ATTR_SELF_TRIGGER_MODE, AQMD3_VAL_SELF_TRIGGER_MODE_SQUARE_WAVE));
    cout << "Frequency:          " << squareWaveFrequency << "\n";
    cout << "Duty cycle:         " << squareWaveDutyCycle << "\n";
    cout << "Slope:              " << (squareWaveSlope == AQMD3_VAL_TRIGGER_SLOPE_POSITIVE ? "Positive" : "Negative") << "\n\n";
    checkApiCall(AqMD3_SelfTriggerSquareWaveConfigure(session, "SelfTrigger", squareWaveFrequency, squareWaveDutyCycle, squareWaveSlope));

    // Calibrate the instrument.
    cout << "\nPerforming self-calibration\n";
    checkApiCall(AqMD3_SelfCalibrate(session));

    // Prepare the buffer to fetch data.
    ViInt64 arraySize = 0;
    checkApiCall(AqMD3_QueryMinWaveformMemory(session, 32, numRecords, 0, recordSize, &arraySize));
    vector<ViInt32> dataArray(arraySize);
    ViInt32 actualAverages = 0;
    ViInt64 actualRecords = 0;
    ViInt64 actualPoints[numRecords], firstValidPoint[numRecords];
    ViReal64 initialXTimeSeconds[numRecords], initialXTimeFraction[numRecords];
    ViReal64 initialXOffset = 0.0, xIncrement = 0.0, scaleFactor = 0.0, scaleOffset = 0.0;
    ViInt32 flags[numRecords];
    ViBoolean abort = VI_FALSE;

    // Start the acquisition loop.
    cout << "\nPerforming acquisition\n";
    checkApiCall(AqMD3_SelfTriggerInitiateGeneration(session, "SelfTrigger"));
    checkApiCall(AqMD3_InitiateAcquisition(session));
    while (abort == VI_FALSE)
    {
        ViBoolean isComplete = VI_FALSE;
        ViBoolean overflow = VI_FALSE;

        // Check whether acquisition is complete.
        while (!isComplete)
            checkApiCall(AqMD3_GetAttributeViBoolean(session, "", AQMD3_ATTR_TSR_IS_ACQUISITION_COMPLETE, &isComplete));
        printf("Acquisition completed\n");

        // Fetch the acquired data in array.
        checkApiCall(AqMD3_FetchAccumulatedWaveformInt32(session, "Channel1", 0, numRecords, 0, recordSize, arraySize, &dataArray[0],
            &actualAverages, &actualRecords, actualPoints, firstValidPoint,
            &initialXOffset, initialXTimeSeconds, initialXTimeFraction,
            &xIncrement, &scaleFactor, &scaleOffset, flags));

        // Release the memory to the instrument.
        checkApiCall(AqMD3_TSRContinue(session));

        // Check for trigger overflow.
        checkApiCall(AqMD3_GetAttributeViBoolean(session, "", AQMD3_ATTR_TSR_MEMORY_OVERFLOW_OCCURRED, &overflow));
        if (overflow)
        {
            cerr << "Memory Overflow occurred\n";
            abort = VI_TRUE;
        }

        // Convert data to Volts.
        cout << "Processing data\n";
        for (ViInt64 currentRecord = 0; currentRecord < numRecords; ++currentRecord)
        {
            for (ViInt64 currentPoint = 0; currentPoint < actualPoints[currentRecord]; ++currentPoint)
            {
                // The data returned by the driver is unsigned. Since the API is limited to signed data type, a conversion is required
                // to convert from raw value to Volts.
                ViReal64 valueInVolts = static_cast<ViUInt32>(dataArray[firstValidPoint[currentRecord] + currentPoint]) * scaleFactor + scaleOffset;
            }
        }

        // Display trigger-time delta (in seconds)
        for (ViInt64 i = 1; i < numRecords; ++i)
        {
            double triggerTimePrevious = initialXTimeSeconds[i-1] + initialXTimeFraction[i-1] + initialXOffset;
            double triggerTimeCurrent = initialXTimeSeconds[i] + initialXTimeFraction[i] + initialXOffset;
            std::cout << "    - Trigger time difference ( record#"<<i<<" - record#" <<i-1<<") = "<< triggerTimeCurrent- triggerTimePrevious<<" s"<<std::endl;
        }

        cout << "Processing completed\n";
    }

    // Stop the SelfTrigger signal
    checkApiCall(AqMD3_SelfTriggerAbortGeneration(session, "SelfTrigger"));
    // Stop the Acquisition
    checkApiCall(AqMD3_Abort(session));

    // Close the driver.
    checkApiCall(AqMD3_close(session));
    printf("Driver closed \n");

    return 0;
}

