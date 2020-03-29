////////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (C) Acqiris SA 2019
//--------------------------------------------------------------------------------------------------
// LibTool: header only library for AqMD3 IVI-C examples.
//
////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef LIBTOOL_H
#define LIBTOOL_H

#include <sstream>
#include <vector>
#include <numeric>
#include <limits>

namespace LibTool
{
    //! Convert a value to a string.
    /*! Uses the stream insertion operator as defined for the given type to do the conversion. */
    template <typename T>
    inline std::string ToString( T const & value )
    {
        std::ostringstream strm;
        strm << value;
        return strm.str();
    }

    //! Divide `value` by `divider`, always rounding up.
    /*! `value` and `divider` must be integers, and `divider` must be positive.
        \return the smallest integer which is not smaller than `value` divided by `divider`. */
    template <class T>
    inline T CeilDiv(T value, T divider)
    {
        static_assert(std::numeric_limits<T>::is_integer, "Requires integer type");
        if (divider <= 0)
            throw std::invalid_argument("Divider must be positive; got " + ToString(divider));

        if (value > 0 && value > (std::numeric_limits<T>::max() - divider + 1))
            throw std::overflow_error("Integer overflow in CeilDiv");

        return (value + ( value > 0 ? divider - 1 : 0)) / divider;
    }

    //! Align a value to the next higher integer multiple of `grain`.
    template <class T>
    inline T AlignUp(T value, T grain)
    {
        static_assert(std::numeric_limits<T>::is_integer, "Requires integer type");

        if (value > std::numeric_limits<T>::max() - (grain-1))
            throw std::logic_error("Aligned up value exceeding the numeric limit.");

        return CeilDiv(value, grain) * grain;
    }

    //! Represent processing and storage parameters
    struct ProcessingParameters
    {
        size_t storageBlockSamples;     //!< number of samples in a memory block
        size_t processingBlockSamples;  //!< number of samples in a processing block

        ProcessingParameters(size_t storageSamples, size_t processingSamples)
            : storageBlockSamples(storageSamples)
            , processingBlockSamples(processingSamples)
        {}
    };

    //! Represent available tags for markers.
    enum class MarkerTag : uint8_t
    {
        None            = 0x00,
        TriggerNormal   = 0x01, // 512-bit: Trigger marker standard Normal acquisition mode
        TriggerAverager = 0x02, // 512-bit: Trigger marker Averager acquisition mode.
        TriggerExtended = 0x03, // reserved.
        GateStartCst    = 0x04, //  64-bit: Gate start marker in CST mode.
        GateStopCst     = 0x05, //  64-bit: Gate stop marker in CST mode.
        GateStartCsr    = 0x06, // reserved.
        GateStopCsr     = 0x07, // reserved
        DummyGate       = 0x08, //  64-bit: Dummy gate marker.
    };

    //! Represent a marker of a trigger
    struct TriggerMarker
    {
        static uint32_t const RecordIndexMask = 0x00ffffff;

        MarkerTag tag;                //!< marker tag.
        double triggerTimeSamples;    /*!< represent the time difference (in sample interval) between trigger and next sampling time. Values are in [0,1[.
                                           This field does not include trigger delay.*/
        uint64_t absoluteSampleIndex; //!< the absolute index (since module init/reset) of the very first sample of acquisition.
        uint32_t recordIndex;         //!< index of the record.

        //! Return the absolute time of the very first sample of record.
        double GetInitialXTime(double fullRateSampleInterval ) const
        { return double(absoluteSampleIndex)*fullRateSampleInterval ; }

        //! Return the time difference between the very first sample of record and trigger event.
        double GetInitialXOffset(double fullRateSampleInterval , double triggerDelay=0.0) const
        { return triggerTimeSamples*fullRateSampleInterval + triggerDelay; }
    };

    //! Represent a marker of a gate (start & stop)
    struct GateMarker
    {
        //! Represent a marker of a gate stop
        struct GateStopMarker
        {
            uint64_t position;      //!< gate stop position (in processing blocks)
            int32_t startPosition;  //!< reserved
            int32_t thresholdIndex; //!< reserved

            //! Return the stop position in samples
            size_t GetPositionSample(ProcessingParameters const& params) const { return size_t(position-1) * params.processingBlockSamples; }
        };

        //! Represent a marker of a gate start
        struct GateStartMarker
        {
            uint64_t position;  //!< gate start position (in processing blocks)
            int32_t counter;    //!< reserved

            //! Return the start position in samples
            size_t GetPositionSample(ProcessingParameters const& params) const { return size_t(position-1) * params.processingBlockSamples; }
        };

        GateStartMarker start;  //!< gate start marker
        GateStopMarker stop;    //!< gate stop marker

        // Return the number of samples required to store data associated with the gate in memory (including samples added for padding).
        size_t GetRequiredMemorySamples(ProcessingParameters const& params) const
        {
            size_t const gateBlocks = size_t(stop.position - start.position);
            return AlignUp<size_t>(gateBlocks * params.processingBlockSamples, params.storageBlockSamples);
        }

    };

    //! Represent a record descriptor with trigger marker and gate marker list.
    class RecordDescriptor
    {
    public:
        typedef std::vector<GateMarker> GateList;

        RecordDescriptor(RecordDescriptor const&) = default;

        explicit RecordDescriptor() : m_trigger(), m_gateList() {}

        //! Setter for the trigger marker
        void SetTriggerMarker(TriggerMarker const& t) {m_trigger=t;}

        //! Getter for the trigger marker.
        TriggerMarker const& GetTriggerMarker() const { return m_trigger;}

        //! Add a gate marker to gate marker list.
        void AddGate(GateMarker const& gate) {m_gateList.push_back(gate);}

        //! Return a const reference to the list of gate markers.
        GateList const& GetGateList() const { return m_gateList;}

        // Return the number of samples required to store all gate data in memory (including samples added for padding).
        size_t GetRequiredMemorySamples(ProcessingParameters const& params) const
        {
            auto const& fct = [&params](size_t size, GateMarker const& gate) { return size + gate.GetRequiredMemorySamples(params); };
            return std::accumulate(m_gateList.begin(), m_gateList.end(), size_t(0), fct);
        }


    private:
        TriggerMarker m_trigger;    //!< trigger marker.
        GateList m_gateList;        //!< list of gate markers.
    };

    //! Represents a subsegment of a read-only array.
    /*! The class receives a const reference to a read-only instance of 'std::vector'. The class does not own this given instance, and the
        later must not be resized and/or destroyed until the destruction of all associated 'ArraySegment' instances.*/
    template <typename T> class ArraySegment
    {
    private:
        using value_type = T;
        using pointer = T*;
        using container_type = std::vector<T>;

    public:
        //! Build a segment of 'count' elements starting at 'offset' from 'data' container.
        /*! 'data' must not be resized and/or destroyed during the life of all ArraySegment objects referencing it.*/
        explicit ArraySegment(container_type const& data, size_t offset, size_t count);

        //! Return the size of the segment.
        size_t size() const { return m_size; }

        //! Return a const reference to the item associated with the given index in the segment.
        value_type const& operator[](size_t index) const
        { return m_data[m_offset+index];}

        //! Return a pointer of the first element in the segment.
        pointer GetData() {return pointer(&m_data[m_offset]);}

    private:
        container_type const& m_data;   //!< Const reference to the underlaying data container
        size_t m_offset;                //!< offset in 'm_data' of the very first item of the segment.
        size_t m_size;                  //!< size of the segment
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // ArraySegment class member definitions
    //

    template <typename T>
    inline ArraySegment<T>::ArraySegment(container_type const& data, size_t offset, size_t count)
        : m_data(data)
        , m_offset(offset)
        , m_size(count)
    {
        if (data.size() < offset + count)
        {
            throw std::logic_error("Array segment definition exceeds array size: "
                                    "offset=" + ToString(offset) +
                                    ", count=" + ToString(count) +
                                    ", array size=" + ToString(m_data.size())
                                    );
        }
    }

    //! Represent a circular buffer with both start address and size are aligned the given "grainBytes" constraints.
    /*! The class creates an internal buffer with a slightly larger capacity in order the satisfy the alignment constraint.
        It creates an array segment that starts at an address aligned to 'grainBytes' with the user-requested capacity.

        The circular behavior is achieved by defining two pointers 'start' which points to the first valid element, and 'end'
        pointing to one-past last valid element. Both 'start' and 'end' wrap-around when they reach the end of the linear
        buffer. This mechanism requires a sentinel item to distinguish between two situations: the buffer is empty
        (start == end) and when it is full (end+1=start).*/
    template <class T>
    class AlignedCircularBuffer
    {
    private:
        using container_type=std::vector<T>;

    public:
        //! Creates an aligned circular buffer with maximum capacity 'capacity' aligned to 'grainBytes'.
        /*! Produces an empty circular buffer aligned to 'grainBytes' with a maximum capacity of 'capacity'.
            \param[capacity] requested buffer capacity. It must satisfy the alignment constraint
                             (i.e. capacity*sizeof(T) must be multiple of 'grainBytes')
            \param[grainBytes] alignment constraint in Bytes.*/
        explicit AlignedCircularBuffer(size_t capacity, size_t grainBytes);

        //! Return the capacity of the buffer
        size_t Capacity() const { return m_capacity; }

        //! Return the number of valid elements in the buffer.
        /*! Elements are validated with 'Append()' and invalidated with 'Pop()'*/
        size_t Size() const { return ((m_end - m_start) + Capacity()) % Capacity(); }

        //! Return 'true' when the buffer is empty, return 'false' otherwise.
        bool IsEmpty() const { return m_end == m_start; }

        //! Return 'true' when the buffer is full, return 'false' otherwise.
        bool IsFull() const { return SuccessorIndex(m_end) == m_start; }

        //! Return a const reference to the element associated with the given index.
        T const& operator[](size_t index) const;

        //! Invalidate the 'count' first elements in the buffer.
        void Pop(size_t count);

        //! Validate the 'count' elements in the back of the buffer.
        void Append(size_t count);

        //! Return a segment to the next free contiguous space in the circular buffer.
        ArraySegment<T> GetNextContiguousSpace() const;

        ~AlignedCircularBuffer(){ if(m_alignedData) delete(m_alignedData);}

    private:
        //! Return the index of the one-past the latest element of range [offset, offset+size-1].
        size_t RangeEnd(size_t offset, size_t size) const
        { return (offset + size) % Capacity(); }

        //! Return the index of successor of 'currentIndex'
        size_t SuccessorIndex(size_t currentIndex) const
        { return RangeEnd(currentIndex, 1); }

    private:
        container_type m_data;          //!< internal raw storage with extra items for alignments and sentinal.
        size_t m_capacity;              //!< capacity of the circular buffer.
        size_t m_firstAlignedIndex;     //!< index of the very first item of "m_data" satisfying the alignment constraint.

        ArraySegment<T>* m_alignedData; //!< aligned internal storage.
        size_t m_start;                 //!< start index (in m_alignedData) of valid items.
        size_t m_end;                   //!< index of the one-past the latest item (in m_alignedData) of valid data.
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // member function definitions
    //

    template <class T>
    inline AlignedCircularBuffer<T>::AlignedCircularBuffer(size_t capacity, size_t grainBytes)
        : m_data(capacity + (grainBytes / sizeof(T) - 1) + 1) // user requested capacity + alignment constraint + sentinel item
        , m_capacity(capacity)
        , m_firstAlignedIndex(0)
        , m_alignedData(nullptr)
        , m_start(0)
        , m_end(0)
    {
        if(grainBytes == 0)
            throw std::invalid_argument("Grain must be strict positive power of 2");

        if((grainBytes & (grainBytes-1)) != 0)
            throw std::invalid_argument("Grain must be a power of 2");

        if((grainBytes % sizeof(T)) != 0)
            throw std::invalid_argument("Grain must be multiple of item bytes " + ToString(sizeof(T)));

        size_t const grainElements = grainBytes / sizeof(T);

        if((capacity % grainElements) != 0)
            throw std::invalid_argument("Capacity must be multiple of grain elements " + ToString(grainElements));

        // find the very first item in 'm_data' satisfying alignment constraint.
        uint64_t const alignmentMask = uint64_t(grainBytes - 1);
        for (size_t i = 0; i < grainElements; ++i)
        {
            uint64_t const ptr = reinterpret_cast<uint64_t>(&m_data[i]);
            if ((ptr & alignmentMask) == 0)
            {
                m_firstAlignedIndex = i;
                m_alignedData = new ArraySegment<T>(m_data, m_firstAlignedIndex, m_capacity);
                return;
            }
        }

        throw std::runtime_error("Cannot find any item satisfying the alignment constraint");
    }

    template <class T>
    inline T const& AlignedCircularBuffer<T>::operator[](size_t index) const
    {
        if( Size() <= index )
            throw std::invalid_argument("Index is out of range: " + ToString(index) + ", size=" + ToString(Size()));

        size_t const itemIndex = (m_start + index) % Capacity();
        return (*m_alignedData)[itemIndex];
    }

    template <class T>
    inline void AlignedCircularBuffer<T>::Pop(size_t count)
    {
        if(Size() < count || count == 0)
            throw std::invalid_argument("Cannot pop "+ ToString(count) + " item(s) out from circular buffer of " + ToString(Size()));

        m_start = (m_start + count) % Capacity();
    }

    template <class T>
    inline void AlignedCircularBuffer<T>::Append(size_t count)
    {
        if(count == 0)
            throw std::logic_error("Cannot append buffer with 0 elements");

        if(Capacity() - Size() < count)
            throw std::invalid_argument("Cannot append buffer with the requested number of elements");

        m_end = RangeEnd(m_end, count);
    }

    template <class T>
    inline ArraySegment<T> AlignedCircularBuffer<T>::GetNextContiguousSpace() const
    {
        size_t const remainingSpace = (m_start <= m_end) ? Capacity() - m_end
                                                         : m_start - m_end;

        return ArraySegment<T>(m_data, m_firstAlignedIndex + m_end, remainingSpace);
    }

    //! Represent a decode of a marker stream. The class is designed to work in perpetual streaming context.
    /*! The marker decoder works along with a circular input buffer which reference is taken upon construction. It has one public method
        `DecodeNextRecordDescriptor` which decodes the next record descriptor out from the input stream. The method returns 'true' on success
        and 'false' when the input stream doe not contain enough elements to build a complete record description. In such case, the user
        owns the responsibility of filling the input marker stream and calling the method again. The class is designed to work in the
        following fashion

            AlignedCircularBuffer<ViInt32> inputStream;
            MarkerStreamDecoder decoder(inputStream);

            for(;;)
            {
                RecordDescriptor nextRecord;

                while(!decoder.DecodeNextRecordDescriptor(nextRecord))
                    FillStreamWithElements(inputStream);

                Use( nextRecord );
            }*/
    class MarkerStreamDecoder
    {
    private:
        //! States of the decoder
        enum class State
        {
            ExpectTrigger,              //!< Only expect a trigger marker from input stream (new record descriptor)
            ExpectGateOrDummyOrTrigger  //!< Expect gate marker, dummy gate marker or a trigger marker (currently building a record descriptor).
        };

    public:
        using MarkerStream = AlignedCircularBuffer<int32_t>;

        explicit MarkerStreamDecoder(MarkerStream& markerStream)
            : m_state(State::ExpectTrigger)
            , m_markerStream(markerStream)
            , m_currentRecord()
        {}

        //! Decode the next record descriptor from input stream.
        /*! It processes the input stream which reference is kept in `m_markerStream` until a complete record
            descriptor is built in member `m_currentRecord`, or when processing reaches the end of the marker stream.
            If the input stream is completely consumed without building a complete record descriptor, the user has to
            fill the input buffer with data and call the method again.
            \param[record] output parameter with the resulting record descriptor. The output parameter is updated
                           only when the method returns 'true'. Otherwise, the parameter is left as-is.
            \return 'true' when the input stream has enough data to build a complete record descriptor.*/
        bool DecodeNextRecordDescriptor(RecordDescriptor& record);

    private:
        static RecordDescriptor DecodeRecordMarker(MarkerStream&);
        //! Expect a trigger marker from the input marker stream, decode it and return it as a result.
        static TriggerMarker DecodeTriggerMarker(MarkerStream&);
        //! Expect a gate start marker from the input marker stream, decode it and return it as a result.
        static GateMarker::GateStartMarker DecodeGateStartMarker(MarkerStream&);
        //! Expect a gate stop marker from the input marker stream, decode it and return it as a result.
        static GateMarker::GateStopMarker DecodeGateStopMarker(MarkerStream&);
        //! Expect a gate marker (start+stop) from the input marker stream, decode it and return it as a result.
        static GateMarker DecodeGateMarker(MarkerStream&);
        //! Expect a dummy gate marker and consume it from the input marker stream.
        static void WalkthroughDummyMarker(MarkerStream&);

    private: // helper static methods
        static bool IsTriggerMarkerTag(MarkerTag tag)
        { return tag == MarkerTag::TriggerNormal || tag == MarkerTag::TriggerAverager || tag == MarkerTag::TriggerExtended; }

        static MarkerTag ExtractTag(uint32_t element)
        { return MarkerTag(element & 0xff); }

        // Extract gate position from two 32-bit integers (lsb & msb)
        static uint64_t ExtractPosition(uint32_t lsb, uint32_t msb);

    private:
        State m_state;                      //!< the current state of the decoder.
        MarkerStream& m_markerStream;       //!< a reference to the input marker stream.
        RecordDescriptor m_currentRecord;   //!< the current record descriptor.
    };


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // MarkerStreamDecoder member definitions
    //

    inline bool MarkerStreamDecoder::DecodeNextRecordDescriptor(RecordDescriptor& record)
    {
        while(!m_markerStream.IsEmpty())
        {
            if (m_state == State::ExpectTrigger)
            {
                m_currentRecord.SetTriggerMarker(DecodeTriggerMarker(m_markerStream));
                m_state = State::ExpectGateOrDummyOrTrigger;
            }

            if(m_markerStream.IsEmpty())
                return false;

            if (m_state == State::ExpectGateOrDummyOrTrigger)
            {
                uint32_t const header = m_markerStream[0];
                MarkerTag const tag = ExtractTag(header);

                if(IsTriggerMarkerTag(tag)) // trigger tag is associated with the next record descriptor
                {
                    record = m_currentRecord;
                    m_currentRecord = RecordDescriptor();
                    m_state = State::ExpectTrigger;
                    return true;
                }

                switch(tag)
                {
                    case MarkerTag::DummyGate:
                        WalkthroughDummyMarker(m_markerStream);
                        break;
                    case MarkerTag::GateStartCst:
                        m_currentRecord.AddGate(DecodeGateMarker(m_markerStream));
                        break;
                    default:
                        throw std::runtime_error("Unexpected tag: " + ToString(int(tag)));
                }
            }
        }
        return false;
    }

    inline uint64_t MarkerStreamDecoder::ExtractPosition(uint32_t lsb, uint32_t msb)
    {
        uint64_t const low = (lsb >> 24) & 0xff;
        uint64_t const high = uint64_t(msb & 0xffffff) << 8;
        return high | low;
    }

    inline void MarkerStreamDecoder::WalkthroughDummyMarker(MarkerStream& stream)
    {
        MarkerTag const tag = ExtractTag(stream[0]);
        if( tag != MarkerTag::DummyGate)
            throw std::logic_error("Expected Dummy gate marker, got "+ToString(int(tag)));

        stream.Pop(2);
    }

    inline GateMarker MarkerStreamDecoder::DecodeGateMarker(MarkerStream& stream)
    {
        GateMarker result;
        result.start = DecodeGateStartMarker(stream);
        result.stop = DecodeGateStopMarker(stream);
        return result;
    }

    inline GateMarker::GateStartMarker MarkerStreamDecoder::DecodeGateStartMarker(MarkerStream& stream)
    {
        uint32_t const header = stream[0];
        uint32_t const next = stream[1];
        MarkerTag const tag = ExtractTag(header);
        if(tag != MarkerTag::GateStartCst)
            throw std::logic_error("Expected gate start marker, got "+ToString(int(tag)));

        stream.Pop(2);

        GateMarker::GateStartMarker result;
        result.counter = (header >> 8) & 0xffff;
        result.position = ExtractPosition(header, next);

        return result;
    }

    inline GateMarker::GateStopMarker MarkerStreamDecoder::DecodeGateStopMarker(MarkerStream& stream)
    {
        uint32_t const header = stream[0];
        uint32_t const next = stream[1];
        MarkerTag const tag = ExtractTag(header);
        if(tag != MarkerTag::GateStopCst)
            throw std::logic_error("Expected gate stop marker from header, got "+ToString(int(tag)));

        stream.Pop(2);

        GateMarker::GateStopMarker result;
        result.startPosition = (header >> 8) & 0xff;
        result.thresholdIndex = (header >> 16) & 0xff;
        result.position = ExtractPosition(header, next);

        return result;
    }

    inline TriggerMarker MarkerStreamDecoder::DecodeTriggerMarker(MarkerStream& stream)
    {
        uint32_t const header = stream[0];
        MarkerTag const tag = ExtractTag(header);

        if (!IsTriggerMarkerTag(tag))
            throw std::runtime_error( "Expected trigger marker, got "+ToString(int(tag)) );

        uint32_t const low = stream[1];
        uint32_t const high = stream[2];

        TriggerMarker result;

        result.tag = tag;
        result.recordIndex = (header >> 8) & TriggerMarker::RecordIndexMask;

        result.triggerTimeSamples = -(double(low & 0x000000ff) / 256.0);
        uint64_t const timestampLow = (low >> 8) & 0x0000000000ffffffL;
        uint64_t const timestampHigh = uint64_t(high) << 24;
        result.absoluteSampleIndex = timestampHigh | timestampLow;

        stream.Pop(16);

        return result;
    }

}

#endif
