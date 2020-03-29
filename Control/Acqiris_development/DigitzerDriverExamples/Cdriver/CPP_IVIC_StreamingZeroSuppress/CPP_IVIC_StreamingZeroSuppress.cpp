///
/// Acqiris IVI-C Driver Example Program
///
/// Initializes the driver, reads a few Identity interface properties, and performs a
/// streaming acquisition.
///
/// For additional information on programming with IVI drivers in various IDEs, please see
/// http://www.ivifoundation.org/resources/
///
/// Requires a real instrument having CST and SZ1 options and an input signal on "Channel1".
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
#include<fstream>
#include <algorithm>


#define checkApiCall( f ) do { ViStatus s = f; testApiCall( s, #f ); } while( false )
void testApiCall( ViStatus status, char const * functionName);

// Get the sample with 'index' as 16-bit aligned integer from the stream of elements (32-bit signed integers).
ViInt16 GetSample(LibTool::AlignedCircularBuffer<int32_t> const& sampleBuffer, size_t index);

//! Unpack the record samples into 'output' stream from 'sampleBuffer' according to the record descriptor 'recordDesc'
void UnpackRecord(LibTool::RecordDescriptor const& recordDesc, LibTool::AlignedCircularBuffer<int32_t> const& sampleBuffer, LibTool::ProcessingParameters const& processingParams, double sampleInterval, std::ostream& output);

//! Perform a stream fetch of at maximum 'nbrElementsToFetch' elements into the aligned circular buffer "buffer"
ViInt64 FetchAligned(ViSession session, ViConstString streamName, ViInt64 nbrElementsToFetch, LibTool::AlignedCircularBuffer<int32_t>& buffer, milliseconds timeout=milliseconds(2000));

// name-space gathering all user-configurable parameters
namespace
{
    // Edit resource and options as needed. Resource is ignored if option has Simulate=true.
    // An input signal is necessary if the example is run in non simulated mode, otherwise
    // the acquisition will time out.
    ViChar resource[] = "PXI40::0::0::INSTR";
    ViChar options[]  = "Simulate=false, DriverSetup= Model=SA220P";

    // Acquisition configuration parameters
    ViReal64 const sampleRate = 2.0e9;
    ViReal64 const sampleInterval = 1.0 / sampleRate;
    ViInt64 const recordSize = 1024;
    ViInt32 const streamingMode = AQMD3_VAL_STREAMING_MODE_TRIGGERED;
    ViInt32 const dataReductionMode = AQMD3_VAL_ACQUISITION_DATA_REDUCTION_MODE_ZERO_SUPPRESS;

    // Channel configuration parameters
    ViConstString channel = "Channel1";
    ViReal64 const range = 2.5;
    ViReal64 const offset = 0.0;
    ViInt32 const coupling = AQMD3_VAL_VERTICAL_COUPLING_DC;

    // Channel ZeroSuppress configuration parameters
    ViInt32 const zsThreshold = 0;
    ViInt32 const zsHysteresis = 300;
    ViInt32 const zsPreGateSamples = 0;
    ViInt32 const zsPostGateSamples = 0;

    // Trigger configuration
    ViConstString triggerSource = "Internal1";
    ViReal64 const triggerLevel = 0.0;
    ViInt32 const triggerSlope =  AQMD3_VAL_TRIGGER_SLOPE_POSITIVE;

    // Readout parameters
    ViInt64 const processingBlockSamples = 8;
    ViConstString sampleStreamName = "StreamCh1";
    ViConstString markerStreamName = "MarkersCh1";
    ViInt64 const nbrSampleElementsToFetch = 4 * 1024 * 1024;
    ViInt64 const nbrMarkerElementsToFetch = 64 * 1024;
    ViInt64 const nbrSampleBufferElements = nbrSampleElementsToFetch * 4;
    ViInt64 const nbrMarkerBufferElements = nbrMarkerElementsToFetch * 4;

    // Streaming session parameters
    seconds streamingDuration = seconds(5);

    // Output file
    std::string const outputFileName("StreamingZeroSuppress.log");
}

int main()
{
    cout << "Triggered Streaming with ZeroSuppress\n\n";

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

        // Configure the acquisition in triggered mode with ZeroSuppress enabled.
        cout << "Configuring Acquisition\n";
        cout << "  Record size :        " << recordSize << '\n';
        cout << "  Streaming mode :     " << streamingMode << '\n';
        cout << "  SampleRate:          " << sampleRate << '\n';
        cout << "  Data Reduction mode: " << dataReductionMode << '\n';
        checkApiCall( AqMD3_SetAttributeViInt32( session, "", AQMD3_ATTR_STREAMING_MODE, streamingMode) );
        checkApiCall( AqMD3_SetAttributeViReal64( session, "", AQMD3_ATTR_SAMPLE_RATE, sampleRate ) );
        checkApiCall( AqMD3_SetAttributeViInt32( session, "", AQMD3_ATTR_ACQUISITION_DATA_REDUCTION_MODE, dataReductionMode) );
        checkApiCall( AqMD3_SetAttributeViInt64( session, "", AQMD3_ATTR_RECORD_SIZE, recordSize) );

        // Configure the channels.
        cout << "Configuring " << channel << "\n";
        cout << "  Range:              " << range << '\n';
        cout << "  Offset:             " << offset << '\n';
        cout << "  Coupling:           " << ( coupling?"DC":"AC" ) << '\n';
        checkApiCall( AqMD3_ConfigureChannel( session, channel, range, offset, coupling, VI_TRUE ) );

        // Configure ZeroSuppress
        cout << "Configuring ZeroSuppress\n";
        cout << "  Threshold:          " << zsThreshold << '\n';
        cout << "  Hysteresis:         " << zsHysteresis << '\n';
        cout << "  PreGate Samples:    " << zsPreGateSamples << '\n';
        cout << "  PostGate Samples:   " << zsPostGateSamples << '\n';
        checkApiCall( AqMD3_SetAttributeViInt32( session, channel, AQMD3_ATTR_CHANNEL_ZERO_SUPPRESS_HYSTERESIS, zsHysteresis) );
        checkApiCall( AqMD3_SetAttributeViInt32( session, channel, AQMD3_ATTR_CHANNEL_ZERO_SUPPRESS_THRESHOLD, zsThreshold) );
        checkApiCall( AqMD3_SetAttributeViInt32( session, channel, AQMD3_ATTR_CHANNEL_ZERO_SUPPRESS_PRE_GATE_SAMPLES, zsPreGateSamples) );
        checkApiCall( AqMD3_SetAttributeViInt32( session, channel, AQMD3_ATTR_CHANNEL_ZERO_SUPPRESS_POST_GATE_SAMPLES, zsPostGateSamples) );

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
        checkApiCall( AqMD3_GetAttributeViString( session, "", AQMD3_ATTR_INSTRUMENT_FIRMWARE_REVISION,         sizeof( str ), str ) );
        cout << "Firmware revision:  " << str << '\n';
        checkApiCall( AqMD3_SelfCalibrate( session ) );

        // Prepare readout buffers
        ViInt64 sampleStreamGrain = 0;
        ViInt64 markerStreamGrain = 0;
        checkApiCall( AqMD3_GetAttributeViInt64( session, sampleStreamName , AQMD3_ATTR_STREAM_GRANULARITY_IN_BYTES, &sampleStreamGrain) );
        checkApiCall( AqMD3_GetAttributeViInt64( session, markerStreamName , AQMD3_ATTR_STREAM_GRANULARITY_IN_BYTES, &markerStreamGrain) );

        LibTool::AlignedCircularBuffer<int32_t> sampleStreamBuffer(nbrSampleBufferElements, sampleStreamGrain);
        LibTool::AlignedCircularBuffer<int32_t> markerStreamBuffer(nbrMarkerBufferElements, markerStreamGrain);
        LibTool::MarkerStreamDecoder markerDecoder(markerStreamBuffer);

        // Start the acquisition.
        cout << "\nInitiating acquisition\n";
        checkApiCall( AqMD3_InitiateAcquisition( session ) );
        cout << "Acquisition is running\n\n";

        ViInt64 expectedRecordIndex = 0;
        ViInt64 totalSampleElements = 0;
        ViInt64 totalMarkerElements = 0;

        LibTool::ProcessingParameters const processingParams(sampleStreamGrain / sizeof(ViInt16), processingBlockSamples);

        std::ofstream outputFile(outputFileName);
        auto endTime = system_clock::now() + streamingDuration;
        while( system_clock::now() < endTime )
        {
            LibTool::RecordDescriptor recordDesc;

            // Decode the next record descriptor from marker stream (do fetch if necessary)
            while(!markerDecoder.DecodeNextRecordDescriptor(recordDesc))
                totalMarkerElements += FetchAligned(session, markerStreamName, nbrMarkerElementsToFetch, markerStreamBuffer);

            // Check that the record descriptor holds the expected record index and the expected tag
            if ((expectedRecordIndex & LibTool::TriggerMarker::RecordIndexMask) != recordDesc.GetTriggerMarker().recordIndex)
                throw std::runtime_error("Unexpected record index: expected="+ToString(expectedRecordIndex)+", got " + ToString(recordDesc.GetTriggerMarker().recordIndex));

            if(LibTool::MarkerTag::TriggerNormal !=  recordDesc.GetTriggerMarker().tag)
                throw std::runtime_error("Expected normal trigger tag, got "+ToString(uint8_t(recordDesc.GetTriggerMarker().tag)));

            // Compute the volume of packed data and make sure the circular buffer has the required storage capacity.
            size_t const packedRecordSamples = recordDesc.GetRequiredMemorySamples(processingParams);
            size_t const packedRecordElements = packedRecordSamples / 2; // One 32-bit element encodes two 16-bit samples.

            if (sampleStreamBuffer.Capacity() < packedRecordElements)
                throw std::runtime_error("Packed record elements "+ToString(packedRecordElements)+" exceeds the capacity of the circular buffer " + ToString(sampleStreamBuffer.Capacity()));

            // Fetch all the volume of samples associated with the current packed record.
            while (sampleStreamBuffer.Size() < packedRecordElements)
                totalSampleElements += FetchAligned(session, sampleStreamName, nbrSampleElementsToFetch, sampleStreamBuffer);

            // Rebuild the record and do some processing into output file
            UnpackRecord(recordDesc, sampleStreamBuffer, processingParams, 1.0/fullSampleRate, outputFile);

            // Consume all the volume of elements associated with the current record.
            sampleStreamBuffer.Pop(packedRecordElements);

            ++expectedRecordIndex;
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
        std::cerr << "Unexpected error: " << exc.what() << std::endl;

        if(session != VI_NULL)
            checkApiCall( AqMD3_close( session ) );
        cout << "\nDriver closed\n";

        return(1);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//
// definition of Helper functions
//

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

ViInt16 GetSample(LibTool::AlignedCircularBuffer<int32_t> const& sampleBuffer, size_t index)
{
    size_t const elementIndex = index / 2;
    int32_t const element = sampleBuffer[elementIndex];
    if((index % 2) == 0)
        return ViInt16(element & 0xffff);
    else
        return ViInt16((element >> 16) & 0xffff);
}

void UnpackRecord(
      LibTool::RecordDescriptor const& recordDesc
    , LibTool::AlignedCircularBuffer<int32_t> const& sampleBuffer
    , LibTool::ProcessingParameters const& processingParams
    , double fullRateSampleInterval
    , std::ostream& output)
{
    LibTool::TriggerMarker const& trig = recordDesc.GetTriggerMarker();
    double const xTime = trig.GetInitialXTime(fullRateSampleInterval);
    double const xOffset = trig.GetInitialXOffset(fullRateSampleInterval);
    output << "# record index                 : " << std::dec << trig.recordIndex << '\n';
    output << "# Absolute Time of First Sample: " << std::setprecision(12) << xTime << '\n';
    output << "# Absolute Time of Trigger     : " << std::setprecision(12) << xTime+xOffset << '\n';

    output << std::dec;
    size_t nextGateOffset = 0;
    for (LibTool::GateMarker const& gate : recordDesc.GetGateList())
    {
        output << " - Gate: start=" << gate.start.GetPositionSample(processingParams)
               << ", stop=" << gate.stop.GetPositionSample(processingParams)
               << ", data=[";

        size_t const startPosition = gate.start.GetPositionSample(processingParams); // sample position in record
        size_t const stopPosition = gate.stop.GetPositionSample(processingParams);
        size_t const nbrGateSamples = stopPosition-startPosition;

        if(nbrGateSamples <= 32) // print all samples of small gates
        {
            size_t sampleIndex=nextGateOffset; // sample index in the stream of samples
            for (size_t samplePosition = startPosition; samplePosition < stopPosition; ++samplePosition, ++sampleIndex)
            {
                output << GetSample(sampleBuffer, sampleIndex);
                if(samplePosition < gate.stop.GetPositionSample(processingParams) - 1)
                    output << ", ";
            }
        }
        else
        {
            // print first 3 samples and last 2 samples of the gate.
            output << GetSample(sampleBuffer, nextGateOffset+0) << ", "
                   << GetSample(sampleBuffer, nextGateOffset+1) << ", "
                   << GetSample(sampleBuffer, nextGateOffset+2) << ", "
                   << "... , "
                   << GetSample(sampleBuffer, nextGateOffset+nbrGateSamples-2) << ", "
                   << GetSample(sampleBuffer, nextGateOffset+nbrGateSamples-1);
        }

        output << "]\n";
        nextGateOffset += gate.GetRequiredMemorySamples(processingParams);
    }

    output << '\n';
}

ViInt64 FetchAligned(ViSession session, ViConstString streamName, ViInt64 nbrElementsToFetch, LibTool::AlignedCircularBuffer<int32_t>& buffer, milliseconds timeout)
{
    // Request the next chunk of contiguous space in the circular buffer, and check it is not empty.
    LibTool::ArraySegment<int32_t> nextSegment = buffer.GetNextContiguousSpace();

    if(nextSegment.size() == 0)
        throw std::invalid_argument("No contiguous space left in circular buffer");

    // compute the effective number of elements to fetch according to the capacity of 'nextSegment'.
    ViInt64 elementsToFetch = std::min(ViInt64(nextSegment.size()), nbrElementsToFetch);

    auto endTime = system_clock::now() + timeout;
    do
    {
        ViInt64 availableElements;
        ViInt64 actualElements;
        ViInt64 firstValidElement;

        // Try to fetch acquired data.
        checkApiCall( AqMD3_StreamFetchDataInt32( session, streamName, nbrElementsToFetch, nextSegment.size(), (ViInt32*)nextSegment.GetData(),
                  &availableElements, &actualElements, &firstValidElement ) );

        if (actualElements == 0)
        {
            std::cout << ".";
            sleep_for( milliseconds(5) );
            continue;
        }

        // the circular buffer is assumed to be aligned according to stream granularity. So, 'firstValidElement' must be 0
        if(firstValidElement != 0)
            throw std::runtime_error("First valid point index is different than zero " +ToString(firstValidElement));

        // Append the circular buffer with the number of fetched elements and return.
        buffer.Append(size_t(actualElements));
        cout << streamName << " fetch: " << actualElements << " actual elements, remaining elements: " << availableElements << "\n";
        return actualElements;

    }while(system_clock::now() <= endTime);

    throw std::runtime_error("Maximum time exceeded while waiting for fetch of " + ToString(elementsToFetch) + " elements from stream " + streamName +
                              ". Please verify the source of trigger and its rate.");
}


