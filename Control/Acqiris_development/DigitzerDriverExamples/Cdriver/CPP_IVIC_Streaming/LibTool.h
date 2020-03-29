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

    //! Represent available tags for markers.
    enum class MarkerTag : uint8_t
    {
        None            = 0x00,
        TriggerNormal   = 0x01, // 512-bit: Trigger marker standard Normal acquisition mode
        TriggerAverager = 0x02, // 512-bit: Trigger marker Averager acquisition mode.
        TriggerExtended = 0x03, // reserved.
        GateStartCst    = 0x04, // reserved
        GateStopCst     = 0x05, // reserved
        GateStartCsr    = 0x06, // reserved.
        GateStopCsr     = 0x07, // reserved
        DummyGate       = 0x08, // reserved.
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
        double GetInitialXOffset(double fullRateSampleInterval, double triggerDelay=0.0) const
        { return triggerTimeSamples*fullRateSampleInterval + triggerDelay; }
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
        size_t Size() const { return m_size; }

        //! Return a const reference to the item associated with the given index in the segment.
        value_type const& operator[](size_t index) const
        { return m_data[m_offset+index];}

        //! Return a pointer of the first element in the segment.
        pointer GetData() {return pointer(&m_data[m_offset]);}

        //! Skip the first 'nbrElements' elements from the array segment. Size is reduced accordingly, elements are not destroyed.
        void PopFront(size_t nbrElements);

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

    template <typename T>
    inline void ArraySegment<T>::PopFront(size_t nbrElements)
    {
        if(Size() < nbrElements)
            throw std::logic_error("Cannot pop "+ToString(nbrElements)+" elements out from a segment of "+ToString(Size())+".");

        m_offset += nbrElements;
        m_size -= nbrElements;
    }

    //! Special namespace gathering helper function associated with marker stream decoding.
    namespace MarkerStreamDecoder
    {
        using MarkerStream = ArraySegment<int32_t>;

        //! Expect a trigger marker from the input marker stream, decode it and return it as a result.
        TriggerMarker DecodeTriggerMarker(MarkerStream&);

        //! Tell whether the given tag corresponds to a trigger marker.
        bool IsTriggerMarkerTag(MarkerTag tag)
        { return tag == MarkerTag::TriggerNormal || tag == MarkerTag::TriggerAverager || tag == MarkerTag::TriggerExtended; }

        //! Extract the marker tag from the given header.
        MarkerTag ExtractTag(uint32_t element)
        { return MarkerTag(element & 0xff); }

        static constexpr size_t NbrTriggerMarkerElements = 16; // 512-bit (16 elements of 32-bit)
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // MarkerStreamDecoder member definitions
    //

    inline TriggerMarker MarkerStreamDecoder::DecodeTriggerMarker(MarkerStream& stream)
    {
        uint32_t const header = stream[0];
        MarkerTag const tag = ExtractTag(header);

        if (!IsTriggerMarkerTag(tag))
            throw std::runtime_error( "Unexpected trigger marker: "+ToString(int(tag)) );

        uint32_t const low = stream[1];
        uint32_t const high = stream[2];

        TriggerMarker result;

        result.tag = tag;
        result.recordIndex = (header >> 8) & TriggerMarker::RecordIndexMask;

        result.triggerTimeSamples = -(double(low & 0x000000ff) / 256.0);
        uint64_t const timestampLow = (low >> 8) & 0x0000000000ffffffL;
        uint64_t const timestampHigh = uint64_t(high) << 24;
        result.absoluteSampleIndex = timestampHigh | timestampLow;

        stream.PopFront(NbrTriggerMarkerElements);

        return result;
    }

}

#endif
