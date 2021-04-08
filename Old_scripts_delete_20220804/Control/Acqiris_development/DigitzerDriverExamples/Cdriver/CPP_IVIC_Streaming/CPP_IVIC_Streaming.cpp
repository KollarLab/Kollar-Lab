///
/// Acqiris IVI-C Driver Example Program
///
/// Initializes the driver, reads a few Identity interface properties, and performs a
/// streaming acquisition.
///
/// For additional information on programming with IVI drivers in various IDEs, please see
/// http://www.ivifoundation.org/resources/
///
/// The Example requires a real instrument having CST and input signal on "Channel1". It
/// also requires an AVG option to enable Averager acquisition mode.
///

#include "LibTool.h"
#include "AqMD3.h"

#include <iomanip>
#include <iostream>
using std::cout;
using std::cerr;
using std::hex;
#include <vector>
using std::vector;
#include <stdexcept>
using std::runtime_error;
#include <chrono>
using std::chrono::minutes;
using std::chrono::seconds;
using std::chrono::milliseconds;
using std::chrono::system_clock;
using LibTool::ToString;
#include <thread>
using std::this_thread::sleep_for;
#include <fstream>
#include <algorithm>
#include <fstream>
#include <queue>


#define checkApiCall( f ) do { ViStatus s = f; testApiCall( s, #f ); } while( false )

typedef std::vector<int32_t> FetchBuffer;

//! Validate success status of the given functionName.
void testApiCall( ViStatus status, char const * functionName );

//! Perform a fetch of exact 'nbrElementsToFetch' elements into the given fetch "buffer".
/*! The function throws std::runtime_error exception if it cannot fetch the requested number of elements within the given time delay 'timeout'.
    It might also perform two separate fetch operations to handle the case of memory end boundary (i.e. not enough elements are available at
    the end of memory). In such situation, a second fetch is required to read the remaining elements from the beginning of the memory.*/
LibTool::ArraySegment<int32_t> FetchElements(ViSession session, ViConstString streamName, ViInt64 nbrElementsToFetch, FetchBuffer& buffer, milliseconds timeout=milliseconds(2000));

//! Save record information in output steam
void SaveRecord(LibTool::TriggerMarker const& triggerMarker, ViInt64 recordElements, LibTool::ArraySegment<int32_t> const& sampleBuffer, double sampleInterval, std::ostream& output);

// name-space gathering all user-configurable parameters
namespace
{
    // Edit resource and options as needed. Resource is ignored if option has Simulate=true.
    // An input signal is necessary if the example is run in non simulated mode, otherwise
    // the acquisition will time out.
    ViChar resource[] = "PXI40::0::0::INSTR";
    ViChar options[]  = "Simulate=true, DriverSetup= Model=SA220P";

    // Acquisition configuration parameters
    ViReal64 const sampleRate = 2.0e9;
    ViReal64 const sampleInterval = 1.0 / sampleRate;
    ViInt64 const recordSize = 1024;
    ViInt32 const streamingMode = AQMD3_VAL_STREAMING_MODE_TRIGGERED;
    ViInt32 const acquisitionMode = AQMD3_VAL_ACQUISITION_MODE_NORMAL;
    ViInt32 const nbrAverages = 16;

    // Channel configuration parameters
    ViReal64 const range = 2.5;
    ViReal64 const offset = 0.0;
    ViInt32 const coupling = AQMD3_VAL_VERTICAL_COUPLING_DC;

    // Trigger configuration parameters
    ViConstString triggerSource = "Internal1";
    ViReal64 const triggerLevel = 0.0;
    ViInt32 const triggerSlope =  AQMD3_VAL_TRIGGER_SLOPE_POSITIVE;

    // Fetch parameters
    ViConstString sampleStreamName = "StreamCh1";
    ViConstString markerStreamName = "MarkersCh1";
    ViInt64 const nbrOfRecordsToFetchAtOnce = 8096 ;

    // duration of the streaming session
    auto const streamingDuration = seconds(60);

    // Output file
    std::string const outputFileName("Streaming.log");
}

int main()
{
    cout << "Triggered Streaming \n\n";

    // Initialize the driver. See driver help topic "Initializing the IVI-C Driver" for additional information.
    ViSession session = VI_NULL;
    ViBoolean const idQuery = VI_FALSE;
    ViBoolean const reset   = VI_FALSE;

    try
    {
        checkApiCall( AqMD3_InitWithOptions( resource, idQuery, reset, options, &session ) );

        cout << "Driver initialized \n";

        // Read and output a few attributes.
        ViChar str[128];
        checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_PREFIX,               sizeof( str ), str ) );
        cout << "Driver prefix:      " << str << '\n';
        checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_REVISION,             sizeof( str ), str ) );
        cout << "Driver revision:    " << str << '\n';
        checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_VENDOR,               sizeof( str ), str ) );
        cout << "Driver vendor:      " << str << '\n';
        checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_SPECIFIC_DRIVER_DESCRIPTION,          sizeof( str ), str ) );
        cout << "Driver description: " << str << '\n';
        checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_MODEL,                     sizeof( str ), str ) );
        cout << "Instrument model:   " << str << '\n';
        checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_OPTIONS,              sizeof( str ), str ) );
        cout << "Instrument options: " << str << '\n';
        checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION,         sizeof( str ), str ) );
        cout << "Firmware revision:  " << str << '\n';
        checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_INFO_SERIAL_NUMBER_STRING, sizeof( str ), str ) );
        cout << "Serial number:      " << str << '\n';
        cout << '\n';

        // Abort execution if instrument is still in simulated mode.
        ViBoolean simulate;
        checkApiCall( AqMD3_GetAttributeViBoolean( session, "", AQMD3_ATTR_SIMULATE, &simulate ) );
        if( simulate==VI_TRUE )
        {
            cout << "\nThe Streaming features are not supported in simulated mode.\n";
            cout << "Please update the resource string (resource[]) to match your configuration,";
            cout << " and update the init options string (options[]) to disable simulation.\n";

            AqMD3_close( session );

            return 1;
        }

        // Get full sampling rate.
        ViReal64 fullSampleRate = 0.0;
        checkApiCall( AqMD3_GetAttributeViReal64( session, "", AQMD3_ATTR_SAMPLE_RATE, &fullSampleRate ) );
        ViReal64 const fullRateSampleInterval = 1.0 / fullSampleRate;

        // Configure the acquisition in triggered mode with ZeroSuppress enabled.
        cout << "Configuring Acquisition\n";
        cout << "  Record size :        " << recordSize << '\n';
        cout << "  Streaming mode :     " << streamingMode << '\n';
        cout << "  SampleRate:          " << sampleRate << '\n';
        cout << "  Acquisition mode:    " << acquisitionMode << '\n';
        cout << "  Number of averages:  " << nbrAverages << '\n';
        checkApiCall( AqMD3_SetAttributeViInt32( session, "", AQMD3_ATTR_ACQUISITION_NUMBER_OF_AVERAGES, nbrAverages) );
        checkApiCall( AqMD3_SetAttributeViInt32( session, "", AQMD3_ATTR_STREAMING_MODE, streamingMode) );
        checkApiCall( AqMD3_SetAttributeViReal64( session, "", AQMD3_ATTR_SAMPLE_RATE, sampleRate ) );
        checkApiCall( AqMD3_SetAttributeViInt32( session, "", AQMD3_ATTR_ACQUISITION_MODE, acquisitionMode) );
        checkApiCall( AqMD3_SetAttributeViInt64( session, "", AQMD3_ATTR_RECORD_SIZE, recordSize) );

        // Configure the channels.
        cout << "Configuring Channel1\n";
        cout << "  Range:              " << range << '\n';
        cout << "  Offset:             " << offset << '\n';
        cout << "  Coupling:           " << ( coupling?"DC":"AC" ) << '\n';
        checkApiCall( AqMD3_ConfigureChannel( session, "Channel1", range, offset, coupling, VI_TRUE ) );

        // Configure the trigger.
        cout << "Configuring Trigger\n";
        cout << "  ActiveSource:       " << triggerSource << '\n';
        cout << "  Level:              " << triggerLevel << "\n";
        cout << "  Slope:              " << (triggerSlope ? "Positive" : "Negative") << "\n";
        checkApiCall( AqMD3_SetAttributeViString( session, "", AQMD3_ATTR_ACTIVE_TRIGGER_SOURCE, triggerSource ) );
        checkApiCall( AqMD3_SetAttributeViReal64( session, triggerSource, AQMD3_ATTR_TRIGGER_LEVEL, triggerLevel ) );
        checkApiCall( AqMD3_SetAttributeViInt32( session, triggerSource, AQMD3_ATTR_TRIGGER_SLOPE, triggerSlope ) );

        // Calibrate the instrument.
        cout << "\nApply setup and run self-calibration\n";
        checkApiCall( AqMD3_ApplySetup( session ) );
        checkApiCall( AqMD3_SelfCalibrate( session ) );

        // Prepare readout buffer
        ViInt64 sampleStreamGrain = 0;
        ViInt64 markerStreamGrain = 0;
        checkApiCall( AqMD3_GetAttributeViInt64( session, sampleStreamName , AQMD3_ATTR_STREAM_GRANULARITY_IN_BYTES, &sampleStreamGrain) );
        checkApiCall( AqMD3_GetAttributeViInt64( session, markerStreamName , AQMD3_ATTR_STREAM_GRANULARITY_IN_BYTES, &markerStreamGrain) );
        ViInt64 const sampleStreamGrainElements = sampleStreamGrain / sizeof(int32_t);
        ViInt64 const markerStreamGrainElements = markerStreamGrain / sizeof(int32_t);

        ViInt64 const nbrSamplesPerElement = (acquisitionMode == AQMD3_VAL_ACQUISITION_MODE_AVERAGER) ? 1
                                           : (acquisitionMode == AQMD3_VAL_ACQUISITION_MODE_NORMAL)   ? 2
                                           : throw std::logic_error("Unexpected acquisition mode");

        ViInt64 const recordElements = recordSize / nbrSamplesPerElement;
        ViInt64 const acquisitionElements = recordElements * nbrOfRecordsToFetchAtOnce;
        ViInt64 const markerElements = LibTool::MarkerStreamDecoder::NbrTriggerMarkerElements * nbrOfRecordsToFetchAtOnce;

        ViInt64 const sampleStreamBufferSize = acquisitionElements           // required elements
                                             + acquisitionElements/2         // unwrapping overhead (only in single channel mode)
                                             + sampleStreamGrainElements - 1;// alignment overhead

        ViInt64 const markerStreamBufferSize = markerElements                // required elements
                                             + markerStreamGrainElements - 1;// alignment overhead

        FetchBuffer sampleStreamBuffer(sampleStreamBufferSize);
        FetchBuffer markerStreamBuffer(markerStreamBufferSize);

        // Expected values and statistics
        double minXtime = 0.0; // InitialXTime is the time of the very first sample in the record.

        LibTool::MarkerTag const expectedTag = (acquisitionMode == AQMD3_VAL_ACQUISITION_MODE_AVERAGER) ? LibTool::MarkerTag::TriggerAverager
                                             : (acquisitionMode == AQMD3_VAL_ACQUISITION_MODE_NORMAL)   ? LibTool::MarkerTag::TriggerNormal
                                             : throw std::logic_error("Unexpected acquisition mode");

        ViInt64 expectedRecordIndex = 0;

        // Count the total volume of fetched markers and elements.
        ViInt64 totalSampleElements = 0;
        ViInt64 totalMarkerElements = 0;

        // Start the acquisition.
        cout << "\nInitiating acquisition\n";
        checkApiCall( AqMD3_InitiateAcquisition( session ) );
        cout << "Acquisition is running\n\n";

        std::ofstream outputFile(outputFileName);

        auto const endTime = system_clock::now() + streamingDuration;
        while( system_clock::now() < endTime )
        {
            // Fetch markers of requested records
            LibTool::ArraySegment<int32_t> markerArraySegment = FetchElements(session, markerStreamName, markerElements, markerStreamBuffer);
            totalMarkerElements += markerArraySegment.Size();

            // Fetch all samples of requested records
            LibTool::ArraySegment<int32_t> sampleArraySegment = FetchElements(session, sampleStreamName, acquisitionElements, sampleStreamBuffer);
            totalSampleElements += sampleArraySegment.Size();

            // Process acquired records
            for(ViInt64 i = 0; i < nbrOfRecordsToFetchAtOnce; ++i)
            {
                // 1. decode trigger marker from marker stream
                LibTool::TriggerMarker const nextTriggerMarker = LibTool::MarkerStreamDecoder::DecodeTriggerMarker(markerArraySegment);

                // 2. Validate marker consistency: tag, incrementing record index, increasing xtime.
                if(expectedTag !=  nextTriggerMarker.tag)
                    throw std::runtime_error("Unexpected trigger marker tag: got "+ToString(int(nextTriggerMarker.tag))+", expected "+ToString(int(expectedTag)));

                // 2.1 Check that the record descriptor holds the expected record index.
                if ((expectedRecordIndex & LibTool::TriggerMarker::RecordIndexMask) != nextTriggerMarker.recordIndex)
                    throw std::runtime_error("Unexpected record index: expected="+ToString(expectedRecordIndex)+", got " + ToString(nextTriggerMarker.recordIndex));

                // 2.2 initialXTime (time of first sample in record) must increase.
                ViReal64 const xtime = nextTriggerMarker.GetInitialXTime(fullRateSampleInterval);
                if(xtime <= minXtime)
                    throw std::runtime_error("InitialXTime not increasing: minimum expected="+ToString(minXtime)+", got " + ToString(xtime));

                // 3. Process the record "acquisitionElements" with "nextTriggerMarker"
                SaveRecord(nextTriggerMarker, recordElements, sampleArraySegment, fullRateSampleInterval, outputFile);

                // 3.1 remove record elements from the segment and advance to elements of the next record
                sampleArraySegment.PopFront(recordElements);

                ++expectedRecordIndex;
                minXtime = xtime;
            }
        }
        outputFile.close();

        ViInt64 const totalSampleData = totalSampleElements * sizeof(ViInt32);
        ViInt64 const totalMarkerData = totalMarkerElements * sizeof(ViInt32);
        cout << "\nTotal sample data read: " << (totalSampleData/(1024*1024)) << " MBytes.\n";
        cout << "Total marker data read: " << (totalMarkerData/(1024*1024)) << " MBytes.\n";
        cout << "Duration: " << (streamingDuration/seconds(1)) << " seconds.\n";
        ViInt64 const totalData = totalSampleData + totalMarkerData;
        cout << "Data rate: " << (totalData)/(1024*1024)/(streamingDuration/seconds(1)) << " MB/s.\n";

        // Stop the acquisition.
        cout << "\nStopping acquisition\n";
        checkApiCall( AqMD3_Abort( session ) );

        // Close the driver.
        checkApiCall( AqMD3_close( session ) );
        cout << "\nDriver closed\n";
        return 0;
    }
    catch(std::exception const& exc)
    {
        std::cerr << "Example Code Error: " << exc.what() << std::endl;

        if(session != VI_NULL)
            checkApiCall( AqMD3_close( session ) );
        cout << "\nDriver closed\n";

        return(1);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Definition of local functions
//

// Utility function to check status error during driver API call.
void testApiCall( ViStatus status, char const * functionName )
{
    ViInt32 ErrorCode;
    ViChar ErrorMessage[512];

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


// Helper function to wait untill the data available on the module reaches the requested volume within a specific time-frame.
void FetchWithTimeout(ViSession session, ViConstString streamName, ViInt64 nbrElementsToFetch, ViInt64 dataSize, int32_t* data, ViInt64& firstValidElement, ViInt64& actualElements, milliseconds timeout)
{
    auto endTime = system_clock::now() + timeout;
    do
    {
        ViInt64 availableElements;

        // Try to fetch acquired data.
        checkApiCall( AqMD3_StreamFetchDataInt32( session, streamName, nbrElementsToFetch, dataSize, (ViInt32*)data, &availableElements, &actualElements, &firstValidElement ) );

        if (actualElements == 0)
        {
            std::cout << ".";
            sleep_for( milliseconds(5) );
            continue;
        }

        cout << streamName << " fetch: " << actualElements << " actual elements, remaining elements: " << availableElements << "\n";
        return; // successful fetch.

    }while(system_clock::now() <= endTime);

    throw std::runtime_error("Maximum time exceeded while waiting for fetch of " + ToString(nbrElementsToFetch) + " elements from stream " + streamName +
                              ". Please verify the source of trigger and its rate.");
}

LibTool::ArraySegment<int32_t> FetchElements(ViSession session, ViConstString streamName, ViInt64 nbrElementsToFetch, FetchBuffer& buffer, milliseconds timeout)
{
    int32_t* bufferData = buffer.data();
    ViInt64 const bufferSize = buffer.size();

    if(bufferSize < nbrElementsToFetch)
        throw std::invalid_argument("Buffer size is smaller than the requested elements to fetch");

    ViInt64 firstElement = 0;
    ViInt64 actualElements = 0;
    FetchWithTimeout(session, streamName, nbrElementsToFetch, bufferSize, bufferData, firstElement, actualElements, timeout);

    if(nbrElementsToFetch <= actualElements)
    {
        return LibTool::ArraySegment<int32_t>(buffer, firstElement, actualElements);
    }
    else
    {
        /* Handle the case of memory-end boundary where the fetch might return less elements than requested. Perform a second fetch to read the remaining elements. */
        ViInt64 nextFetchFirstElement = 0;
        ViInt64 nextFetchActualElements = 0;
        int32_t* nexFetchtBufferData = bufferData + firstElement + actualElements;
        ViInt64 const nextFetchBufferSize = bufferSize - firstElement - actualElements;

        ViInt64 const remainingElementsToFetch = nbrElementsToFetch - actualElements;
        if(nextFetchBufferSize < remainingElementsToFetch)
            throw std::invalid_argument("Buffer size is smaller than the remaining elements to fetch");

        FetchWithTimeout(session, streamName, remainingElementsToFetch, nextFetchBufferSize, nexFetchtBufferData, nextFetchFirstElement, nextFetchActualElements, timeout);

        // given buffer is expected to be aligned (thanks to the first fetch). firstElement must be equal to 0.
        if(nextFetchFirstElement != 0)
            throw std::runtime_error("First valid point index is different than zero " +ToString(nextFetchFirstElement));

        return LibTool::ArraySegment<int32_t>(buffer, firstElement, actualElements+nextFetchActualElements);
    }
}


void SaveRecord(LibTool::TriggerMarker const& triggerMarker, ViInt64 recordElements, LibTool::ArraySegment<int32_t> const& sampleBuffer, double fullRateSampleInterval, std::ostream& output)
{
    double const xTime = triggerMarker.GetInitialXTime(fullRateSampleInterval);
    double const xOffset = triggerMarker.GetInitialXOffset(fullRateSampleInterval);
    output << "# record index                 : " << std::dec << triggerMarker.recordIndex << '\n';
    output << "# Absolute Time of First Sample: " << std::setprecision(12) << xTime << '\n';
    output << "# Absolute Time of Trigger     : " << std::setprecision(12) << xTime+xOffset << '\n';

    output << "Elements(" << std::dec << recordElements << ") = [ "; // An element= 2 samples in Normal mode, and 1 sample in Averager mode.

    // Print all elements of small records.
    if (recordElements <= 10)
    {
        for (int i = 0; i < recordElements; ++i)
            output << sampleBuffer[i] <<" ";
    }
    else
    {
        // print first 5 elements and last 2 elements
        output << sampleBuffer[0] << " "
               << sampleBuffer[1] << " "
               << sampleBuffer[2] << " "
               << sampleBuffer[3] << " "
               << sampleBuffer[4] << " "
               << "... "
               << sampleBuffer[recordElements-2] << " "
               << sampleBuffer[recordElements-1] << " ";
    }

    output << "]\n\n";
}
