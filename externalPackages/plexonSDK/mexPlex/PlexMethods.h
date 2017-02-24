//
// PlexMethods.h
//
// Copyright (c) 1998-2013 Plexon, Inc, Dallas, Texas 75206
//     www.plexon.com, support@plexon.com
//
// This code is provided as-is, in order to assist Plexon users who wish to 
// use non-Windows platforms such as Linux. Plexon cannot provide support 
// specific to any non-Windows platform, system, compiler, etc., nor guarantee 
// correct operation or performance in such environments.
//
// Plexon users are free to modify this code as necessary for their personal 
// use, but this code is not open-source. Only code downloaded directly from 
// plexon.com is guaranteed to be the latest version.
//

#ifndef _PLEXMETHODS_H_INCLUDED
#define _PLEXMETHODS_H_INCLUDED

#include "PlexonFiles.h"
#include <mex.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

// the code assumes little-endian byte ordering
// the code also assumes the following is true regarding int-related data types:
// int is 32 bits
// long long is 64 bits
// both MSVC and GCC interpret int and long long this way in 32-bit and 64-bit builds

// size_t corresponds to the integral data type returned by the language operator sizeof and 
// it usually has 32 bits in 32-bit builds and 64 bits in 64-bit builds

// we use 64-bit integers for timestamps
typedef unsigned long long TimeStampType;

// convert timestamp parts to a full 64-bit timestamp
#define toTimeStampType( upper8, lower32 ) \
    ( (static_cast<TimeStampType>(upper8)<<32) + static_cast<TimeStampType>(lower32) )

// max number of unit entries in tscounts and wfcounts matrices returned by plx_info
#define MAX_NUM_UNITS               (26)

#define MAX_NUMBER_OF_WORDS_FOLLOWING_DATA_BLOCK_HEADER (512)

// with huge .plx files, we may have a situation when matlab array requires too much memory
// for example, if you try this command in Matlab: zeros(1024*1024*1024,1) on a 64-bit W7 with 6 GB RAM,
// the system will practically freeze for several minutes.
// I found that on a 64-bit W7 with 6 GB RAM zeros(1024*1024*512,1) works OK. 
// you may want to adjust this number according to your system
#define MAX_NUMBER_OF_VALUES_IN_ONE_MATLAB_MATRIX (1024*1024*512)

// the trick to figure out maximum value of size_t:
// we cast minus one to size_t:
// size_t is guaranteed to be an unsigned type, 
// and a signed integer is converted to unsigned by taking it
// modulo 2^n, where n is the number of bits in the unsigned type.
#define MAX_SIZE_T_VALUE (( size_t )- 1)

#define MEXPLEX_VERSION_NUMBER (182)

// supported mexPlex functions
enum SupportedFunctions {
    FIRST_FUNCTION          = 0,

    DDT                     = 1,
    PLX_AD                  = 2,
    PLX_EVENT_TS            = 3,
    PLX_INFO                = 4,
    PLX_TS                  = 5,
    PLX_WAVES               = 6,
    PLX_AD_SPAN             = 7,
    PLX_CHAN_GAINS          = 8,
    PLX_CHAN_THRESHOLD      = 9,
    PLX_CHAN_FILTERS        = 10,
    PLX_ADCHAN_GAINS        = 11,
    PLX_ADCHAN_FREQS        = 12,
    PLX_INFORMATION         = 13,
    PLX_CHAN_NAMES          = 14,
    PLX_ADCHAN_NAMES        = 15,
    PLX_EVENT_NAMES         = 16,
    PLX_AD_V                = 17,
    PLX_AD_SPAN_V           = 18,
    PLX_WAVES_V             = 19,
    DDT_V                   = 20,
    PLX_VT_INTERPRET        = 21,
    PLX_CLOSE               = 22,
    PLX_ADCHAN_SAMPLECOUNTS = 23,
    PLX_AD_GAP_INFO         = 24,
    DDT_WRITE               = 25,
    PLX_CHANMAP             = 26,
    PLX_AD_CHANMAP          = 27,
    PLX_EV_CHANMAP          = 28,
    GET_VERSION             = 29,

    LAST_FUNCTION               // dummy
};

// helper template function to return the value from a map
template< typename MapType, typename KeyArgType, typename ValueArgType >
bool GetMapValue( const MapType& theMap, const KeyArgType& theKey, ValueArgType& theValue )
{
    typename MapType::const_iterator it = theMap.find( theKey );
    if ( it != theMap.end() ) {
        theValue = it->second;
        return true;
    } else {
        return false;
    }
}

// file reader interface
class IReader
{
public:
    IReader() {}
    virtual ~IReader() {}
    virtual bool Open( const char* fileName ) = 0;
    virtual unsigned long long GetFileSize() = 0;
    virtual bool Seek( unsigned long long seekPositionFromStartOfFile ) = 0;
    virtual void Close() = 0;
    virtual unsigned int Read( void* buffer, unsigned int numberOfBytesToRead ) = 0;
    virtual bool IsAtEndOfFile() = 0;
};

// simple file reader implementation using FILE*
class Reader : public IReader
{
public:
    Reader();
    virtual ~Reader();
    virtual bool Open( const char* fileName );
    virtual unsigned long long GetFileSize();
    virtual bool Seek( unsigned long long seekPositionFromStartOfFile );
    virtual void Close();
    virtual unsigned int Read( void* buffer, unsigned int numberOfBytesToRead );
    virtual bool IsAtEndOfFile();
protected:
    FILE* m_pFile;
};

// this class contains information about a single .plx file
// we have a single global instance of this class in mexPlex library
class CPlexMethods
{
private:
    // note that .ddt file read and write methods use local file variables, not m_Reader class member
    bool DoOpenFile( const char* fileName );

    // internal methods to read data from .plx files
    bool SeekToDataStart();
    bool ReadNextBlock( PL_DataBlockHeader& m_db, std::vector<short>& buf );
    bool ValidateDataBlockHeader( const PL_DataBlockHeader &db ) const;
    double GetTimestamp( const PL_DataBlockHeader &db ) const;

    // current file and its properties
    IReader* m_Reader;
    PL_FileHeader m_FileHeader;
    std::string m_CurrentFileName;
    bool m_bTallied;
    unsigned long long m_PlxFileSize;
    unsigned int m_DataStartFilePosition;

    // file channel headers
    std::vector<PL_ChanHeader>  m_SpikeChannelHeaders;
    std::vector<PL_EventHeader> m_EventChannelHeaders;
    std::vector<PL_SlowChannelHeader> m_AnalogChannelHeaders;

    // mapping from channel number to header index
    std::map<int, int> m_SpikeChannelNumberToSpikeChannelHeaderIndex;
    std::map<int, int> m_EventChannelNumberToEventChannelHeaderIndex;
    std::map<int, int> m_AnalogChannelNumberToAnalogChannelHeaderIndex;

    // here we want to report that there is no header for the specified channel
    // therefore, we return true, if there is a header and false, if there is no header
    bool GetSpikeHeaderIndexFromSpikeChannelNumber( int channelNumber, int& index ) const  {
        return GetMapValue<std::map<int, int>, int, int>( m_SpikeChannelNumberToSpikeChannelHeaderIndex, channelNumber, index );
    }

    bool GetEventHeaderIndexFromEventChannelNumber( int channelNumber, int& index ) const {
        return GetMapValue<std::map<int, int>, int, int>( m_EventChannelNumberToEventChannelHeaderIndex, channelNumber, index );
    }

    bool GetAnalogHeaderIndexFromAnalogChannelNumber( int channelNumber, int& index ) const {
        return GetMapValue<std::map<int, int>, int, int>( m_AnalogChannelNumberToAnalogChannelHeaderIndex, channelNumber, index );
    }

    // channels are 0-based in these maps, unit 0 = unsorted
    // channel is the first member of the pair, unit is the second
    std::map<std::pair<int, int>, long long> m_WaveformCounts;
    std::map<std::pair<int, int>, long long> m_TimestampCounts;

    // these maps are all indexed by the raw channel number (the channel number in the data blocks)
    std::map<int, long long > m_AnalogSampleCounts;
    std::map<int, long long > m_EventCounts;

    // if there is no count for the specified [channel,unit] pair, we simply return zero
    long long GetWaveformCount( int channel, int unit ) const {
        long long theCount = 0;
        if ( GetMapValue<std::map<std::pair<int, int>, long long>, std::pair<int, int>, long long>( m_WaveformCounts, std::make_pair( channel, unit ), theCount ) ) {
            return theCount;
        } else {
            return 0;
        }
    }

    long long GetTimestampCount( int channel, int unit ) const {
        long long theCount = 0;
        if ( GetMapValue<std::map<std::pair<int, int>, long long>, std::pair<int, int>, long long>( m_TimestampCounts, std::make_pair( channel, unit ), theCount ) ) {
            return theCount;
        } else {
            return 0;
        }
    }

    long long GetAnalogChannelSampleCount( int channelNumber ) const {
        long long theCount = 0;
        if ( GetMapValue<std::map<int, long long >, int, long long>( m_AnalogSampleCounts, channelNumber, theCount ) ) {
            return theCount;
        } else {
            return 0;
        }
    }

    long long  GetEventChannelSampleCount( int channelNumber ) const {
        long long theCount = 0;
        if ( GetMapValue<std::map<int, long long >, int, long long>( m_EventCounts, channelNumber, theCount ) ) {
            return theCount;
        } else {
            return 0;
        }
    }

    // conversion coefficient from ADC samples to mV
    double  calculateRawToMilliVoltsCoeffForAnalogChannel( int analogChannelHeaderIndex ) const;
    // conversion coefficient from spike samples to mV
    double  calculateRawToMilliVoltsCoeffForSpikeChannel( int spikeChannelHeaderIndex ) const;

public:
    CPlexMethods();
    ~CPlexMethods();

    enum Unit { NONE, RAW, MV };

    int GetNumberOfPointsInWaveform() { return m_FileHeader.NumPointsWave;};
    const char* GetCurrentFileName() { return m_CurrentFileName.c_str(); }

    // public methods that are called from mexFunction 

    // .ddt methods 
    void ddt_read( const char*, mxArray **ppmxDdt, bool mV, int& numberOfChannels, unsigned int& dataCount, int& ADfrequency );
    int  ddt_write( const char* filename, int nch, int npoints, double freq, const mxArray** d );

    // .plx methods
    bool plx_open( const char* fileName );
    void plx_close();
    bool tallyPlxFile( );
    void plx_info( mxArray **ppmxTs, mxArray **ppmxWf, mxArray **ppmxEv, mxArray **ppmxCont, int fullRead );

    // methods that retrieve actual data from the file return the number of 
    // data items retrieved as size_t. On 64-bit systems, size_t is 64-bit wide 
    size_t plx_ts( int channel_to_extract, int unit_to_extract, mxArray **ppmxTs );
    size_t plx_waves( int channel_to_extract, int unit_to_extract, mxArray **ppmxWf, mxArray **ppmxTs, bool mV = false );
    size_t plx_ad( int channel_to_extract, mxArray **ppmxTs, mxArray **ppmxFn, mxArray **ppmxWf, Unit mV, int& adFrequency );
    size_t plx_ad_span( int t, long long, long long, mxArray **ppmxWf, bool mV, int& ADFrequency );
    size_t plx_event_ts( int event_to_extract, mxArray **ppmxTs, mxArray **ppmxEv );

    void plx_information( int*, int*, char**, int*, int*, int*, int*, int*, int*, int*, double*, char*, size_t DateTimeLen );
    int  plx_chan_gains( mxArray **ppmxGains );
    int  plx_chan_thresholds( mxArray **ppmxTreshs );
    int  plx_chan_filters( mxArray **ppmxFilters );
    int  plx_adchan_gains( mxArray **ppmxGains );
    int  plx_adchan_freqs( mxArray **ppmxFreqs );
    int  plx_adchan_samplecounts( mxArray **ppmxSamples );
    int  plx_chan_names( mxArray **ppmxNames );
    int  plx_adchan_names( mxArray **ppmxNames );
    int  plx_event_names( mxArray **ppmxNames );
    int  plx_chanmap( mxArray **ppmxChanmap );
    int  plx_ad_chanmap( mxArray **ppmxChanmap );
    int  plx_ev_chanmap( mxArray **ppmxChanmap );

    void vt_interpret(
        const mxArray* ts,
        const mxArray* sv,
        mxArray** nCoords,
        mxArray** nDim,
        mxArray** vtMode,
        mxArray** coords
        );
};

#endif /* _PLEXMETHODS_H_INCLUDED */