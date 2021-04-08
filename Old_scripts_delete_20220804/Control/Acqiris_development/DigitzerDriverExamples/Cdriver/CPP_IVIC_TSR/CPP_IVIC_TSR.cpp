///
/// Acqiris IVI-C Driver Example Program
///
/// Initializes the driver, reads a few Identity interface properties, and performs an
/// acquisition with the TSR multi-record mode.
///
/// For additional information on programming with IVI drivers in various IDEs, please see
/// http://www.ivifoundation.org/resources/
///
/// WARNING:
/// The TSR features are not supported in simulation mode. You will have to update
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
#include <chrono>
using std::chrono::milliseconds;
#include <thread>
using std::this_thread::sleep_for;


// Edit resource and options as needed. Resource is ignored if option has Simulate=true.
// An input signal is necessary if the example is run in non simulated mode, otherwise
// the acquisition will time out.
ViChar resource[] = "PXI40::0::0::INSTR";
ViChar options[]  = "Simulate=true, DriverSetup= Model=U5303A";

#define checkApiCall( f ) do { ViStatus s = f; testApiCall( s, #f ); } while( false )

ViInt64 const numRecords = 20;
ViInt64 const recordSize = 1600;

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
    cout << "TsrModeAcquisition\n\n";

    // Initialize the driver. See driver help topic "Initializing the IVI-C Driver" for additional information.
    ViSession session;
    ViBoolean const idQuery = VI_FALSE;
    ViBoolean const reset   = VI_FALSE;
    checkApiCall( AqMD3_InitWithOptions( resource, idQuery, reset, options, &session ) );
    cout << "Driver initialized \n";

    // Abort execution if instrument is still in simulated mode.
    ViBoolean simulate;
    checkApiCall( AqMD3_GetAttributeViBoolean( session, "", AQMD3_ATTR_SIMULATE, &simulate ) );
    if ( simulate == VI_TRUE )
    {
        cout << "The TSR features are not supported in simulated mode.";
        cout << "\nPlease update the resource string (resource[]) to match your configuration, and update the init options string (options[]) to disable simulation.";

        return -1;
    }

    // Check the instrument contains the required TSR module option.
    ViChar str[128] = { '\0' };
    checkApiCall(AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS, sizeof( str ), str ) );
    if (string( str ).find( "TSR" ) == string::npos)
    {
        cout << "The required TSR module option is missing from the instrument.";
        return -1;
    }

    // Read and output a few attributes.
    checkApiCall(AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_PREFIX, sizeof( str ), str ) );
    cout << "Driver prefix:      " << str << '\n';
    checkApiCall(AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_REVISION, sizeof( str ), str ) );
    cout << "Driver revision:    " << str << '\n';
    checkApiCall(AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_VENDOR, sizeof( str ), str ) );
    cout << "Driver vendor:      " << str << '\n';
    checkApiCall(AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION, sizeof( str ), str ) );
    cout << "Driver description: " << str << '\n';
    checkApiCall(AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_MODEL, sizeof( str ), str ) );
    cout << "Instrument model:   " << str << '\n';
    checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS, sizeof( str ), str ) );
    cout << "Instrument options: " << str << '\n';
    checkApiCall(AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION, sizeof( str ), str ) );
    cout << "Firmware revision:  " << str << '\n';
    checkApiCall(AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING, sizeof( str ), str ) );
    cout << "Serial number:      " << str << '\n';

    cout << "\nSimulate:           " << ( simulate?"True":"False") << '\n';

    // Configure the acquisition.
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
    checkApiCall( AqMD3_SetAttributeViInt64( session, "", AQMD3_ATTR_RECORD_SIZE, recordSize ) );
    checkApiCall( AqMD3_SetAttributeViInt64( session, "", AQMD3_ATTR_NUM_RECORDS_TO_ACQUIRE, numRecords ) );
    checkApiCall( AqMD3_SetAttributeViBoolean( session, "", AQMD3_ATTR_TSR_ENABLED, true ) );


    // Configure the trigger.
    cout << "\nConfiguring trigger\n";
    checkApiCall(AqMD3_SetAttributeViString( session, "", AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, "Internal1" ) );

    // Calibrate the instrument.
    cout << "\nPerforming self-calibration\n";
    checkApiCall( AqMD3_SelfCalibrate( session ) );

    // Initiate the acquisition.
    cout << "\nInitiating acquisition\n";
    checkApiCall( AqMD3_InitiateAcquisition( session ) );

    // Recover acquired data.
    ViConstString channelName = "Channel1";
    ViInt64 firstRecord = 0;
    ViInt64 offsetWithinRecord = 0;

    ViInt64 arraySize = 0;

    // Retrieve the array size needed to fetch the waveforms
    checkApiCall( AqMD3_QueryMinWaveformMemory( session, 16, numRecords, 0, recordSize, &arraySize ) );
    vector<ViInt16> dataArray( arraySize );

    ViInt64 actualRecords = 0;
    ViInt64 actualPoints[numRecords];
    ViInt64 firstValidPoint[numRecords];

    ViReal64 initialXOffset[numRecords];
    ViReal64 initialXTimeSeconds[numRecords];
    ViReal64 initialXTimeFraction[numRecords];

    ViReal64 xIncrement = 0.0;
    ViReal64 scaleFactor = 0.0;
    ViReal64 scaleOffset = 0.0;

    int const numberOfLoops = 50;
    for ( int loop = 0; loop < numberOfLoops; ++loop )
    {
        // Check for memory overflow.
        ViBoolean tsrOverflow;
        checkApiCall(AqMD3_GetAttributeViBoolean( session, "", AQMD3_ATTR_TSR_MEMORY_OVERFLOW_OCCURRED, &tsrOverflow ) );
        if ( tsrOverflow )
        {
            cout << "A memory overflow occurred during TSR acquisition.\n";
            return -1;
        }

        // Poll for trigger events.
        int const timeoutInMs = 1000;
        int timeCounter = 0;

        while ( timeCounter < timeoutInMs )
        {
            ViBoolean IsTsrAcquisitionComplete;
            checkApiCall( AqMD3_GetAttributeViBoolean( session, "", AQMD3_ATTR_TSR_IS_ACQUISITION_COMPLETE, &IsTsrAcquisitionComplete ) );
            if ( IsTsrAcquisitionComplete )
                break;

            sleep_for(milliseconds(1000) );
            ++timeCounter;
        }

        // Check for timeout.
        if ( timeCounter >= timeoutInMs )
        {
            cout << "A timeout occurred while waiting for trigger during TSR acquisition.\n";

            return -1;
        }

        // Fetch acquired data.
        checkApiCall( AqMD3_FetchMultiRecordWaveformInt16(
            session,
            channelName,
            firstRecord,
            numRecords,
            offsetWithinRecord,
            recordSize,
            arraySize,
            &dataArray[0],
            &actualRecords,
            actualPoints,
            firstValidPoint,
            initialXOffset,
            initialXTimeSeconds,
            initialXTimeFraction,
            &xIncrement,
            &scaleFactor,
            &scaleOffset
            )
        );

        // Release the corresponding memory and mark it as available for new acquisitions.
        checkApiCall( AqMD3_TSRContinue( session ) );

        // Convert data to Volts for each record.
        for( ViInt64 currentRecord = 0; currentRecord<actualRecords; ++currentRecord )
        {
            for( ViInt64 currentPoint = 0; currentPoint<actualPoints[currentRecord]; ++currentPoint )
            {
                ViReal64 valueInVolts = dataArray[firstValidPoint[currentRecord]+currentPoint]*scaleFactor + scaleOffset;
            }
        }

    }
    cout << "Processing completed\n";

    // Stop the acquisition.
    cout << "Stopping acquisition\n";
    checkApiCall( AqMD3_Abort( session ) );

    // Close the driver.
    checkApiCall( AqMD3_close( session ) );
    cout << "Driver closed\n";

    return 0;
}

