///
/// Acqiris IVI-C Driver Example Program
///
/// Initializes the driver, reads a few Identity interface properties, and performs an
/// acquisition with peak detection.
///
/// For additional information on programming with IVI drivers in various IDEs, please see
/// http://www.ivifoundation.org/resources/
///
/// WARNING:
/// The PeakDetection features are not supported in simulation mode. You will have to update
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
#include <stdexcept>
using std::runtime_error;
#include <string>
using std::string;

#define checkApiCall( f ) do { ViStatus s = f; testApiCall( s, #f ); } while( false )

// Edit resource and options as needed. Resource is ignored if option Simulate=true.
// An input signal is necessary if the example is run in non simulated mode, otherwise
// the acquisition will time out. Set the Simulate option to false to disable simulation.
char resource[] = "PXI40::0::0::INSTR";
char options[] = "Simulate=true, DriverSetup= Model=U5303A";

ViInt64 const numRecords = 1;
ViInt32 const numAverages = 80;

// Utility function to check status error during driver API call.
void testApiCall( ViStatus status, char const * functionName )
{
    ViInt32 ErrorCode;
    ViChar ErrorMessage[256];

    if( status>0 ) // Warning occurred.
    {
        AqMD3_GetError( VI_NULL, &ErrorCode, sizeof( ErrorMessage ), ErrorMessage );
        cerr << "** Warning during " << functionName << ": 0x" << hex << ErrorCode << ", " << ErrorMessage << '\n';

    }
    else if( status<0 ) // Error occurred.
    {
        AqMD3_GetError( VI_NULL, &ErrorCode, sizeof( ErrorMessage ), ErrorMessage );
        cerr << "** ERROR during " << functionName << ": 0x" << hex << ErrorCode << ", " << ErrorMessage << '\n';
        throw runtime_error( ErrorMessage );
    }
}

int main()
{
    cout << "PeakDetection\n\n";

    ViSession session;
    ViBoolean const idQuery = VI_FALSE;
    ViBoolean const reset   = VI_FALSE;

    // Initialize the driver. See driver help topic "Initializing the IVI-C Driver" for additional information.
    checkApiCall( AqMD3_InitWithOptions( resource, idQuery, reset, options, &session ) );

    cout << "Driver initialized \n";

    // Abort execution if instrument is still in simulated mode.
    ViBoolean simulate;
    checkApiCall( AqMD3_GetAttributeViBoolean( session, "", AQMD3_ATTR_SIMULATE, &simulate ) );

    if ( simulate == VI_TRUE )
    {
        cout << "The PeakDetection features are not supported in simulated mode.";
        cout << "\nPlease update the resource string (resource[]) to match your configuration, and update the init options string (options[]) to disable simulation.";

        return -1;
    }

    // Check the instrument contains the required PKD module option.
    ViChar str[128] = {'\0'};
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS, sizeof( str ), str ) );
    if ( string( str ).find( "PKD" ) == string::npos )
    {
        cout << "The required PKD module option is missing from the instrument.";
        return -1;
    }

    // Read and output a few attributes.
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_PREFIX, sizeof( str ), str ) );
    cout << "Driver prefix:      " << str << '\n';
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_REVISION, sizeof( str), str ) );
    cout << "Driver revision:    " << str << '\n';
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_VENDOR, sizeof( str ), str ) );
    cout << "Driver vendor:      " << str << '\n';
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION, sizeof( str ), str ) );
    cout << "Driver description: " << str << '\n';
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_MODEL, sizeof( str ), str) );
    cout << "Instrument model:   " << str << '\n';
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS, sizeof( str ), str ) );
    cout << "Instrument options: " << str << '\n';
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION, sizeof( str ), str ) );
    cout << "Firmware revision:  " << str << '\n';
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING, sizeof( str ), str ) );
    cout << "Serial number:      " << str << '\n';

    cout << "\nSimulate:           " << ( simulate?"True":"False" ) << '\n';

    // Configure the acquisition.
    ViInt64 const recordSize = 1600;
    ViReal64 const range = 1.0;
    ViReal64 const offset = 0.0;
    ViInt32 const coupling = AQMD3_VAL_VERTICAL_COUPLING_DC;
    cout << "Configuring acquisition\n";
    cout << "Range:              " << range << '\n';
    cout << "Offset:             " << offset << '\n';
    cout << "Coupling:           " << ( coupling?"DC":"AC" ) << '\n';
    checkApiCall( AqMD3_ConfigureChannel( session, "Channel1", range, offset, coupling, VI_TRUE ) );
    cout << "Number of records:  " << numRecords << '\n';
    cout << "Record size:        " << recordSize << '\n';
    checkApiCall( AqMD3_SetAttributeViInt64( session, "", AQMD3_ATTR_NUM_RECORDS_TO_ACQUIRE, numRecords ) );
    checkApiCall( AqMD3_SetAttributeViInt64( session, "", AQMD3_ATTR_RECORD_SIZE, recordSize ) );
    cout << "Number of averages: " << numAverages << "\n\n";
    checkApiCall( AqMD3_SetAttributeViInt32( session, "", AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, numAverages ) );
    checkApiCall( AqMD3_SetAttributeViInt32( session, "", AQMD3_ATTR_ACQUISITION_MODE, AQMD3_VAL_ACQUISITION_MODE_PEAK_DETECTION ) );

    // Configure the trigger.
    cout << "Configuring trigger\n";
    checkApiCall( AqMD3_SetAttributeViString( session, "", AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, "Internal1" ) );

    // Calibrate the instrument.
    cout << "Performing self-calibration\n";
    checkApiCall( AqMD3_SelfCalibrate( session ) );

    // Perform the acquisition.
    ViInt32 const timeoutInMs = 1000;
    cout << "Performing acquisition\n";
    checkApiCall( AqMD3_InitiateAcquisition( session ) );
    checkApiCall( AqMD3_WaitForAcquisitionComplete( session, timeoutInMs ) );
    cout << "Acquisition completed\n";

    // Fetch the acquired data in array.
    ViInt64 arraySize = 0;
    checkApiCall( AqMD3_QueryMinWaveformMemory( session, 32, numRecords, 0, recordSize, &arraySize ) );

    vector<ViInt32> dataArray( arraySize );
    ViInt32 actualAverages = 0;
    ViInt64 actualRecords = 0, waveformArrayActualSize = 0;
    ViInt64 actualPoints[numRecords], firstValidPoint[numRecords];
    ViReal64 initialXOffset, initialXTimeSeconds[numRecords], initialXTimeFraction[numRecords];
    ViReal64 xIncrement = 0.0, scaleFactor = 0.0, scaleOffset = 0.0;
    ViInt32 flags[numRecords];

    checkApiCall( AqMD3_FetchAccumulatedWaveformInt32( session, "Channel1",
        0, numRecords, 0, recordSize, arraySize, &dataArray[0],
        &actualAverages, &actualRecords, actualPoints, firstValidPoint,
        &initialXOffset, initialXTimeSeconds, initialXTimeFraction,
        &xIncrement, &scaleFactor, &scaleOffset, flags ) );

    // Convert data to Volts.
    cout << "Processing data\n";
    for ( int currentRecord = 0; currentRecord < actualRecords; ++currentRecord )
    {
        for ( int currentPoint = 0; currentPoint < actualPoints[currentRecord]; ++currentPoint )
        {
            ViReal64 valueInVolts = dataArray[firstValidPoint[currentRecord] + currentPoint] * scaleFactor + scaleOffset;
        }
    }
    cout << "Processing completed\n";

    // Close the driver.
    checkApiCall( AqMD3_close( session ) );
    cout << "Driver closed\n";

    return 0;
}

