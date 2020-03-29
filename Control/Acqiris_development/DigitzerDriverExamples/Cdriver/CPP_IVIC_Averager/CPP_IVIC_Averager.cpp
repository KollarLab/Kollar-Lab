///
///  IVI-C Driver Example Program
///
/// Initializes the driver, reads a few Identity interface properties, and performs an
/// accumulated acquisition.
///
/// For additional information on programming with IVI drivers in various IDEs, please see
/// http://www.ivifoundation.org/resources/
///
/// WARNING:
/// The Averager features are not supported in simulation mode. You will have to update
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
ViChar options[]  = "Simulate=true, DriverSetup= Model=U5303A";


ViInt64 const numRecords = 1;
ViInt32 const numAverages = 80;


#define checkApiCall( f ) do { ViStatus s = f; testApiCall( s, #f ); } while( false )


// Utility function to check status error during driver API call.
void testApiCall( ViStatus status, char const * functionName )
{
    ViInt32 ErrorCode;
    ViChar ErrorMessage[256];

    if( status>0 ) // Warning occurred.
    {
        AqMD3_GetError( VI_NULL, &ErrorCode, sizeof(ErrorMessage), ErrorMessage );
        cerr << "** Warning during " << functionName << ": 0x" << hex << ErrorCode << ", " << ErrorMessage << "\n";

    }
    else if( status<0 ) // Error occurred.
    {
        AqMD3_GetError( VI_NULL, &ErrorCode, sizeof(ErrorMessage), ErrorMessage );
        cerr << "** ERROR during " << functionName << ": 0x" << hex << ErrorCode << ", " << ErrorMessage << "\n";
        throw runtime_error( ErrorMessage );
    }
}


int main()
{
    cout << "\nStarting Averager\n";

    // Initialize the driver. See driver help topic "Initializing the IVI-C Driver" for additional information.
    ViSession session;
    ViBoolean idQuery = VI_FALSE;
    ViBoolean reset   = VI_FALSE;
    checkApiCall( AqMD3_InitWithOptions( resource, idQuery, reset, options, &session ) );

    cout << "\nDriver initialized\n";

    // Abort execution if instrument is still in simulated mode.
    ViBoolean simulate;
    checkApiCall( AqMD3_GetAttributeViBoolean( session, "", AQMD3_ATTR_SIMULATE, &simulate ) );
    if( simulate==VI_TRUE )
    {
        cout << "\nThe Averager features are not supported in simulated mode.\n";
        cout << "Please update the resource string (resource[]) to match your configuration,";
        cout << " and update the init options string (options[]) to disable simulation.\n";

        AqMD3_close( session );

        return 1;
    }
    cout << "\nSimulate:           " << ( simulate==VI_TRUE ? "true" : "false" ) << "\n";

    // Check the instrument contains the required AVG module option.
    ViChar str[128] = { '\0' };
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS, sizeof(str), str ) );
    if( string( str ).find( "AVG" )==string::npos )
    {
        cout << "The required AVG module option is missing from the instrument.\n";

        AqMD3_close( session );

        return 1;
    }

    // Read and output a few attributes.
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_PREFIX, sizeof(str), str ) );
    cout << "Driver prefix:      " << str << "\n";
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_REVISION, sizeof(str), str ) );
    cout << "Driver revision:    " << str << "\n";
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_VENDOR, sizeof(str), str ) );
    cout << "Driver vendor:      " << str << "\n";
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION, sizeof(str), str ) );
    cout << "Driver description: " << str << "\n";
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_MODEL, sizeof(str), str ) );
    cout << "Instrument model:   " << str << "\n";
    string const instrModel = str;
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS, sizeof(str), str) );
    cout << "Instrument options: " << str << '\n';
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION, sizeof(str), str ) );
    cout << "Firmware revision:  " << str << "\n";
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING, sizeof(str), str ) );
    cout << "Serial number:      " << str << "\n";

    // Configure the acquisition.
    cout << "\nConfiguring acquisition\n";
    ViInt64 const recordSize = 1600;
    ViReal64 const range = 1.0;
    ViReal64 const offset = 0.0;
    ViInt32 const coupling = AQMD3_VAL_VERTICAL_COUPLING_DC;
    cout << "Range:              " << range << "\n";
    cout << "Offset:             " << offset << "\n";
    cout << "Coupling:           " << ( coupling ? "DC" : "AC" ) << "\n";
    checkApiCall( AqMD3_ConfigureChannel( session, "Channel1", range, offset, coupling, VI_TRUE ) );
    cout << "Record size:        " << recordSize << "\n";
    checkApiCall( AqMD3_SetAttributeViInt64( session, "", AQMD3_ATTR_RECORD_SIZE, recordSize ) );
    cout << "Number of averages: " << numAverages << "\n";
    checkApiCall( AqMD3_SetAttributeViInt32( session, "", AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, numAverages ) );
    // Have to disable "Channel2" in order to set Averager mode for U5309A
    if( instrModel=="U5309A" )
        checkApiCall( AqMD3_SetAttributeViBoolean( session, "Channel2", AQMD3_ATTR_CHANNEL_ENABLED, VI_FALSE ) );
    checkApiCall( AqMD3_SetAttributeViInt32( session, "", AQMD3_ATTR_ACQUISITION_MODE, AQMD3_VAL_ACQUISITION_MODE_AVERAGER ) );
    checkApiCall( AqMD3_ApplySetup( session ) );

    // Configure the trigger.
    cout << "\nConfiguring trigger (self-trigger)\n";
    double squareWaveFrequency = 10e3; // 10 KHz
    double squareWaveDutyCycle = 10.2; // 10.2% of the period - maximum autorized with frequency of 10 KHz
    ViInt32 squareWaveSlope = AQMD3_VAL_TRIGGER_SLOPE_POSITIVE;
    checkApiCall( AqMD3_SetAttributeViString( session, "", AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, "SelfTrigger" ) );
    checkApiCall( AqMD3_SetAttributeViString( session, "ControlIO3", AQMD3_ATTR_CONTROL_IO_SIGNAL, "Out-AveragerAwg" ) ); // Shunt self-trigger signal to the control IO 3 output.
    checkApiCall( AqMD3_SetAttributeViInt32( session, "SelfTrigger", AQMD3_ATTR_SELF_TRIGGER_MODE, AQMD3_VAL_SELF_TRIGGER_MODE_SQUARE_WAVE ) );
    cout << "Frequency:          " << squareWaveFrequency << "\n";
    cout << "Duty cycle:         " << squareWaveDutyCycle << "\n";
    cout << "Slope:              " << ( squareWaveSlope==AQMD3_VAL_TRIGGER_SLOPE_POSITIVE ? "Positive" : "Negative" ) << "\n\n";
    checkApiCall( AqMD3_SelfTriggerSquareWaveConfigure( session, "SelfTrigger", squareWaveFrequency, squareWaveDutyCycle, squareWaveSlope ) );
    checkApiCall( AqMD3_SelfTriggerInitiateGeneration( session, "SelfTrigger" ) );

    // Calibrate the instrument.
    cout << "\nPerforming self-calibration\n";
    checkApiCall( AqMD3_SelfCalibrate( session ) );

    // Perform the acquisition.
    cout << "\nPerforming acquisition\n";
    ViInt32 timeoutInMs = 1000;
    checkApiCall( AqMD3_InitiateAcquisition( session ) );
    checkApiCall( AqMD3_WaitForAcquisitionComplete( session, timeoutInMs ) );
    cout << "Acquisition completed\n";

    // Fetch the acquired data in array.
    ViInt64 arraySize = 0;
    checkApiCall( AqMD3_QueryMinWaveformMemory( session, 32, 1, 0, recordSize, &arraySize ) );

    vector<ViInt32> dataArray( arraySize );
    ViInt32 actualAverages = 0;
    ViInt64 actualRecords = 0;
    ViInt64 actualPoints[numRecords], firstValidPoint[numRecords];
    ViReal64 initialXTimeSeconds[numRecords], initialXTimeFraction[numRecords];
    ViReal64 initialXOffset = 0.0, XIncrement = 0.0, scaleFactor = 0.0, scaleOffset = 0.0;
    ViInt32 flags[numRecords];
    checkApiCall( AqMD3_FetchAccumulatedWaveformInt32( session, "Channel1", 0, 1, 0, recordSize, arraySize, &dataArray[0],
        &actualAverages, &actualRecords, actualPoints, firstValidPoint,
        &initialXOffset, initialXTimeSeconds, initialXTimeFraction,
        &XIncrement, &scaleFactor, &scaleOffset, flags ) );

    // Convert data to Volts.
    cout << "\nProcessing data\n";
    for( ViInt64 currentRecord=0 ; currentRecord<actualRecords ; ++currentRecord )
    {
        for( ViInt64 currentPoint=0 ; currentPoint<actualPoints[currentRecord] ; ++currentPoint )
        {
            ViInt32 const valueRaw = dataArray[firstValidPoint[currentRecord]+currentPoint];
            ViReal64 const valueInVolts = ViReal64( valueRaw )*scaleFactor + scaleOffset;
            (void)valueInVolts; // Use it!

        }
    }
    cout << "\nProcessing completed\n";

    // Close the driver.
    checkApiCall( AqMD3_close( session ) );
    cout << "Driver closed \n";

    return 0;
}


