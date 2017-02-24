#ifndef _PLEXON_H_INCLUDED
#define _PLEXON_H_INCLUDED


///////////////////////////////////////////////////////////////////////////////
// Plexon Client API Definitions
///////////////////////////////////////////////////////////////////////////////


#define PL_SingleWFType         (1)
#define PL_StereotrodeWFType    (2)     // reserved
#define PL_TetrodeWFType        (3)     // reserved
#define PL_ExtEventType         (4)
#define PL_ADDataType           (5)
#define PL_StrobedExtChannel    (257)
#define PL_StartExtChannel      (258)   // delineates frames, sent for resume also
#define PL_StopExtChannel       (259)   // delineates frames, sent for pause also
#define PL_Pause                (260)   // not used
#define PL_Resume               (261)   // not used

#define MAX_WF_LENGTH           (56)
#define MAX_WF_LENGTH_LONG      (120)



// If the server closes the connection, dll sends WM_CONNECTION_CLOSED message to hWndMain
#define WM_CONNECTION_CLOSED    (WM_USER + 401)



//
// PL_Event is used in PL_GetTimestampStructures(...)
//
struct PL_Event
{
    char    Type;                       // PL_SingleWFType, PL_ExtEventType or PL_ADDataType
    char    NumberOfBlocksInRecord;     // reserved   
    char    BlockNumberInRecord;        // reserved 
    unsigned char    UpperTS;           // Upper 8 bits of the 40-bit timestamp
    unsigned long    TimeStamp;         // Lower 32 bits of the 40-bit timestamp
    short   Channel;                    // Channel that this came from, or Event number
    short   Unit;                       // Unit classification, or Event strobe value
    char    DataType;                   // reserved
    char    NumberOfBlocksPerWaveform;  // reserved
    char    BlockNumberForWaveform;     // reserved
    char    NumberOfDataWords;          // number of shorts (2-byte integers) that follow this header 
}; // 16 bytes


//
// The same as PL_Event above, but with Waveform added
//
struct PL_Wave 
{
    char    Type;                       // PL_SingleWFType, PL_ExtEventType or PL_ADDataType
    char    NumberOfBlocksInRecord;     // reserved   
    char    BlockNumberInRecord;        // reserved 
    unsigned char    UpperTS;           // Upper 8 bits of the 40-bit timestamp
    unsigned long    TimeStamp;         // Lower 32 bits of the 40-bit timestamp
    short   Channel;                    // Channel that this came from, or Event number
    short   Unit;                       // Unit classification, or Event strobe value
    char    DataType;                   // reserved
    char    NumberOfBlocksPerWaveform;  // reserved
    char    BlockNumberForWaveform;     // reserved
    char    NumberOfDataWords;          // number of shorts (2-byte integers) that follow this header 
    short   WaveForm[MAX_WF_LENGTH];    // The actual waveform data
}; // size should be 128

//
// An extended version of PL_Wave for longer waveforms
//
struct PL_WaveLong 
{
    char    Type;                       // PL_SingleWFType, PL_ExtEventType or PL_ADDataType
    char    NumberOfBlocksInRecord;     // reserved   
    char    BlockNumberInRecord;        // reserved 
    unsigned char    UpperTS;           // Upper 8 bits of the 40-bit timestamp
    unsigned long    TimeStamp;         // Lower 32 bits of the 40-bit timestamp
    short   Channel;                    // Channel that this came from, or Event number
    short   Unit;                       // Unit classification, or Event strobe value
    char    DataType;                   // reserved
    char    NumberOfBlocksPerWaveform;  // reserved
    char    BlockNumberForWaveform;     // reserved
    char    NumberOfDataWords;          // number of shorts (2-byte integers) that follow this header 
    short   WaveForm[MAX_WF_LENGTH_LONG];   // The actual long waveform data
}; // size should be 256





// PL_InitClient - initialize client
// In: 
//      type  -- client type. SC registers with type = 256, electrode client with type = 1 
//      hWndList -- handle to the listbox type window.
//                  if hWndList is not null,  PlexClient.dll will send 
//                  LB_ADDSTRING messages with error or debug strings to this window
// Effect:
//      Initializes PlexClient.dll for a client. Opens MMF's and registers 
//                  the client with the server.
//
extern "C" int      WINAPI PL_InitClient(int type, HWND hWndList);



// PL_InitClientEx2 - initialize client
// In: 
//      type  -- client type. SC registers with type = 256, electrode client with type = 1 
//      hWndMain -- handle to the main client window.
//                  if hWndMain is not null, the server sends 
//                  WM_COPYDATA broadcas messages to this window,  
// Effect:
//      Initializes PlexClient.dll for a client. Opens MMF's and registers 
//                  the client with the server.
//
extern "C" int      WINAPI PL_InitClientEx2(int type, HWND hWndMain);


// PL_InitClientEx3 - initialize client
// In: 
//      type  -- client type. SC registers with type = 256, electrode client with type = 1 
//      hWndList -- handle to the listbox type window.
//                  if hWndList is not null,  PlexClient.dll will send 
//                  LB_ADDSTRING messages with error or debug strings to this window
//      hWndMain -- handle to the main client window.
//                  if hWndMain is not null, the server sends 
//                  WM_COPYDATA broadcas messages to this window,  
// Effect:
//      Initializes PlexClient.dll for a client. Opens MMF's and registers 
//                  the client with the server.
extern "C" int      WINAPI PL_InitClientEx3(int type, HWND hWndList, HWND hWndMain);



// PL_CloseClient - closes client connection to the server
// Effect: 
//      Cleans up PlexClient.dll (deletes CClient object) and
//      Sends ClientDisconnected command to the server. 
//      The server decrements the counter for the number of connected clients 
extern "C" void     WINAPI PL_CloseClient();



// PL_IsLongWaveMode - is serevr using long wave mode? 
// Returns:
//      1, if the server uses long waves, 0 otherwise
// Effect: 
//      none
extern "C" int      WINAPI PL_IsLongWaveMode();




// PL_GetTimeStampArrays - get recent timestamps
// In: 
//      *pnmax  - maximum number of timestamps to transfer
// Out:
//      *pnmax - actual number of timestamps transfered
//      type - array of types (PL_SingleWFType, PL_ExtEventType or PL_ADDataType)
//      ch - array of channel numbers
//      cl - array of unit numbers
//      ts - array of timestamps
// Effect: 
//      Copies the timestamps that the server transferred to MMF since
//          any of the PL_GetTimeStamp* or PL_GetWave* was called last time
extern "C" void     WINAPI PL_GetTimeStampArrays(int* pnmax, short* type, short* ch,
                                              short* cl, int* ts);



// PL_GetTimeStampStructures - get recent timestamp structures
// In: 
//      *pnmax  - maximum number of timestamp structures to transfer
// Out:
//      *pnmax - actual number of timestamp structures transferred
//      events - array of PL_Event structures filled with new data
// Effect: 
//      Copies the timestamp structures that the server transferred to MMF since
//          any of the PL_GetTimeStamp* or PL_GetWave* was called last time
extern "C" void     WINAPI PL_GetTimeStampStructures(int* pnmax, 
                                                        PL_Event* events);




// PL_GetTimeStampStructuresEx - get recent timestamp structures
// In: 
//      *pnmax  - maximum number of timestamp structures to transfer
// Out:
//      *pnmax - actual number of timestamp structures transferred
//      events - array of PL_Event structures filled with new data
//      *pollhigh - high DWORD of the perf. counter at the time of HB poll
//      *pollhigh - low DWORD of the perf. counter at the time of HB poll
// Effect: 
//      Copies the timestamp structures that the server transferred to MMF since
//          any of the PL_GetTimeStamp* or PL_GetWave* was called last time
extern "C" void     WINAPI PL_GetTimeStampStructuresEx(int* pnmax, 
                                                    PL_Event* events,
                                                    int* pollhigh,
                                                    int* polllow);


// PL_GetTimeStampStructuresEx2 - get recent timestamp structures
// In: 
//      *pnmax  - maximum number of timestamp structures to transfer
//      includeContinuous - if zero, only spike and external event data is returned
// Out:
//      *pnmax - actual number of timestamp structures transferred
//      events - array of PL_Event structures filled with new data
// Effect: 
//      Copies the timestamp structures that the server transferred to MMF since
//          any of the PL_GetTimeStamp* or PL_GetWave* was called last time
extern "C" void     WINAPI PL_GetTimeStampStructuresEx2(int* pnmax, 
                                                        PL_Event* events,
                                                        int includeContinuous);
                                                        

// PL_GetWaveFormStructures - get recent waveform structures
// In: 
//      *pnmax  - maximum number of waveform structures to transfer
// Out:
//      *pnmax - actual number of waveform structures transferred
//      waves - array of PL_Wave structures filled with new data
// Effect: 
//      Copies the waveform structures that the server transferred to MMF since
//          any of the PL_GetTimeStamp* or PL_GetWave* was called last time
extern "C" void     WINAPI PL_GetWaveFormStructures(int* pnmax, 
                                                        PL_Wave* waves);



// PL_GetWaveFormStructuresEx - get recent waveform structures
// In: 
//      *pnmax  - maximum number of waveform structures to transfer
// Out:
//      *pnmax - actual number of waveform structures transferred
//      waves - array of PL_Wave structures filled with new data
//      *serverdropped - number of waveforms that were dropped in MXI transfer
//      *mmfdropped - number of waveforms that were dropped in MMF->client transfer
// Effect: 
//      Copies the waveform structures that the server transferred to MMF since
//          any of the PL_GetTimeStamp* or PL_GetWave* was called last time
extern "C" void     WINAPI PL_GetWaveFormStructuresEx(int* pnmax, 
                                        PL_Wave* waves, 
                                        int* serverdropped,
                                        int* mmfdropped);




// PL_GetLongWaveFormStructures - get recent long waveform structures
// In: 
//      *pnmax  - maximum number of waveform structures to transfer
// Out:
//      *pnmax - actual number of waveform structures transferred
//      waves - array of PL_WaveLong structures filled with new data
//      *serverdropped - number of waveforms that were dropped in MXI transfer
//      *mmfdropped - number of waveforms that were dropped in MMF->client transfer
// Effect: 
//      Copies the waveform structures that the server transferred to MMF since
//          any of the PL_GetTimeStamp* or PL_GetWave* was called last time
extern "C" void     WINAPI PL_GetLongWaveFormStructures(int* pnmax, 
                                        PL_WaveLong* waves, 
                                        int* serverdropped,
                                        int* mmfdropped);



// PL_GetWaveFormStructuresEx2 - get recent waveform structures
// In: 
//      *pnmax  - maximum number of waveform structures to transfer
// Out:
//      *pnmax - actual number of waveform structures transferred
//      waves - array of PL_Wave structures filled with new data
//      *serverdropped - number of waveforms that were dropped in MXI transfer
//      *mmfdropped - number of waveforms that were dropped in MMF->client transfer
//      *pollhigh - high DWORD of the perf. counter at the time of HB poll
//      *pollhigh - low DWORD of the perf. counter at the time of HB poll
// Effect: 
//      Copies the waveform structures that the server transferred to MMF since
//          any of the PL_GetTimeStamp* or PL_GetWave* was called last time
extern "C" void     WINAPI PL_GetWaveFormStructuresEx2(int* pnmax, PL_Wave* waves,
                                                  int* serverdropped, 
                                                  int* mmfdropped,
                                                  int* pollhigh,
                                                  int* polllow);



// PL_GetLongWaveFormStructuresEx2 - get recent long waveform structures
// In: 
//      *pnmax  - maximum number of waveform structures to transfer
// Out:
//      *pnmax - actual number of waveform structures transferred
//      waves - array of PL_WaveLong structures filled with new data
//      *serverdropped - number of waveforms that were dropped in MXI transfer
//      *mmfdropped - number of waveforms that were dropped in MMF->client transfer
//      *pollhigh - high DWORD of the perf. counter at the time of HB poll
//      *pollhigh - low DWORD of the perf. counter at the time of HB poll
// Effect: 
//      Copies the waveform structures that the server transferred to MMF since
//          any of the PL_GetTimeStamp* or PL_GetWave* was called last time
extern "C" void     WINAPI PL_GetLongWaveFormStructuresEx2(int* pnmax, PL_WaveLong* waves,
                                                  int* serverdropped, 
                                                  int* mmfdropped,
                                                  int* pollhigh,
                                                  int* polllow);


// PL_SendUserEvent - send an external user event to the server. 
// In:
//		channel -- event channel
// Effect: 
//		Sends an external user event to the server. Used by SC to send start/stop
//			commands and to send keyboard events. The server inserts the event into the
//			data stream, i.e. copies it to the MMF.
// Note: this routine will NOT work across PlexNet. A process calling this routine
// must be running on the same machine as Server.
extern "C" void		WINAPI PL_SendUserEvent(int channel);


// PL_SendUserEventWord - send an external strobed-word user event to the server. 
// In:
//		w -- strobed event word value
// Effect: 
//		Same as PL_SendUserEvent, except that a strobed event with a user-specified
//    strobed-word value is sent.
extern "C" void		WINAPI PL_SendUserEventWord(WORD w);


// "get" commands
extern "C" void     WINAPI PL_GetOUTInfo(int* out1, int* out2);
extern "C" void     WINAPI PL_GetOUTInfoEx(int* out1, int* out2, int* out3, int* out4);
extern "C" void     WINAPI PL_GetSlowInfo(int* freq, int* channels, int* gains);
extern "C" void     WINAPI PL_GetSlowInfo64(int* freq, int* channels, int* gains);
extern "C" void     WINAPI PL_GetNumNIDAQCards(int* numcards);
extern "C" void     WINAPI PL_GetSlowInfo256(int* freqs, int* channels, int* gains);
extern "C" void     WINAPI PL_GetNIDAQCardSlow4(int* IsSlow);
extern "C" int      WINAPI PL_GetTIMClockFreq(void);
extern "C" int      WINAPI PL_GetNIDAQBandwidth(void);
extern "C" int      WINAPI PL_GetActiveChannel();
extern "C" int      WINAPI PL_IsElClientRunning();
extern "C" int      WINAPI PL_IsSortClientRunning();
extern "C" int      WINAPI PL_IsNIDAQEnabled();
extern "C" int      WINAPI PL_IsDSPProgramLoaded();
extern "C" int      WINAPI PL_GetTimeStampTick();
extern "C" void     WINAPI PL_GetGlobalPars(int* numch, int* npw, int* npre, int* gainmult);
extern "C" void     WINAPI PL_GetGlobalParsEx(int* numch, int* npw, int* npre, int* gainmult, int* maxwflength);
extern "C" void     WINAPI PL_GetChannelInfo(int* nsig, int* ndsp, int* nout);
extern "C" void     WINAPI PL_GetSIG(int* sig);
extern "C" void     WINAPI PL_GetFilter(int* filter);
extern "C" void     WINAPI PL_GetGain(int* gain);
extern "C" void     WINAPI PL_GetMethod(int* method);
extern "C" void     WINAPI PL_GetThreshold(int* thr);
extern "C" void     WINAPI PL_GetNumUnits(int* numunits);
extern "C" void     WINAPI PL_GetTemplate(int ch, int unit, int* t);
extern "C" void     WINAPI PL_GetNPointsSort(int* npts);
extern "C" int      WINAPI PL_SWHStatus();      // 1 if SWH board is present, otherwise 0
extern "C" int      WINAPI PL_GetPollingInterval();
// this now returns the number of channels per NIDAQ card, not the total number of NIDAQ chans
extern "C" int      WINAPI PL_GetNIDAQNumChannels();
extern "C" void     WINAPI PL_EnableExtLevelStartStop(int enable);
extern "C" int      WINAPI PL_IsNidaqServer();
extern "C" int      WINAPI PL_GetNIDAQBitsPerSample();
extern "C" int      WINAPI PL_IsNIDAQmx();

extern "C" void     WINAPI PL_GetName(int ch1x, char* name);
extern "C" void     WINAPI PL_GetEventName(int ch1x, char* name);
extern "C" void     WINAPI PL_SetSlowChanName(int ch0x, char* name);
extern "C" void     WINAPI PL_GetSlowChanName(int ch0x, char* name);

// not implemented in the verison 09.98
extern "C" void     WINAPI PL_GetValidPCA(int* num);
extern "C" void     WINAPI PL_GetTemplateFit(int ch, int* fit);
extern "C" void     WINAPI PL_GetBoxes(int ch, int* b);
extern "C" void     WINAPI PL_GetPC(int ch, int unit, float* pc);
extern "C" void     WINAPI PL_GetMinMax(int ch,  float* mm);
extern "C" void     WINAPI PL_GetGlobalWFRate(int* t);
extern "C" void     WINAPI PL_GetWFRate(int* t);













///////////////////////////////////////////////////////////////////////////////
// Plexon .plx File Structure Definitions
///////////////////////////////////////////////////////////////////////////////


#define LATEST_PLX_FILE_VERSION 107

#define PLX_HDR_LAST_SPIKE_CHAN     128     // max spike channel number with counts in TSCounts and WFCounts arrays
#define PLX_HDR_LAST_UNIT           4       // max unit number supported by PL_FileHeader information

#define PLX_HDR_LAST_EVENT_CHAN     299     // max digital event number that will be counted in EVCounts

#define PLX_HDR_FIRST_CONT_CHAN_IDX 300     // index in EVCounts for analog channel 0
#define PLX_HDR_LAST_CONT_CHAN      211     // max (0-based) analog channel number that has counts in EVCounts, starting at [300]


// file header (is followed by the channel descriptors)
struct  PL_FileHeader 
{
    unsigned int MagicNumber;   // = 0x58454c50;

    int     Version;            // Version of the data format; determines which data items are valid
    char    Comment[128];       // User-supplied comment 
    int     ADFrequency;        // Timestamp frequency in hertz
    int     NumDSPChannels;     // Number of DSP channel headers in the file
    int     NumEventChannels;   // Number of Event channel headers in the file
    int     NumSlowChannels;    // Number of A/D channel headers in the file
    int     NumPointsWave;      // Number of data points in waveform
    int     NumPointsPreThr;    // Number of data points before crossing the threshold

    int     Year;               // Time/date when the data was acquired
    int     Month; 
    int     Day; 
    int     Hour; 
    int     Minute; 
    int     Second; 

    int     FastRead;           // reserved
    int     WaveformFreq;       // waveform sampling rate; ADFrequency above is timestamp freq 
    double  LastTimestamp;      // duration of the experimental session, in ticks
    
    // The following 6 items are only valid if Version >= 103
    char    Trodalness;                 // 1 for single, 2 for stereotrode, 4 for tetrode
    char    DataTrodalness;             // trodalness of the data representation
    char    BitsPerSpikeSample;         // ADC resolution for spike waveforms in bits (usually 12)
    char    BitsPerSlowSample;          // ADC resolution for slow-channel data in bits (usually 12)
    unsigned short SpikeMaxMagnitudeMV; // the zero-to-peak voltage in mV for spike waveform adc values (usually 3000)
    unsigned short SlowMaxMagnitudeMV;  // the zero-to-peak voltage in mV for slow-channel waveform adc values (usually 5000)
    
    // Only valid if Version >= 105
    unsigned short SpikePreAmpGain;     // usually either 1000 or 500

    // Only valid if Version >= 106
    char    AcquiringSoftware[18];      // name and version of the software that originally created/acquired this data file
    char    ProcessingSoftware[18];     // name and version of the software that last processed/saved this data file



    char    Padding[10];        // so that this part of the header is 256 bytes
    
    
    // Counters for the number of timestamps and waveforms in each channel and unit.
    // Note that even though there may be more than 4 (MAX_HDR_COUNTS_UNITS) units on any 
    // channel, these arrays only record the counts for the first 4 units in each channel.
    // Likewise, starting with .plx file format version 107, there may be more than 128 
    // (MAX_HDR_COUNTS_SPIKE_CHANS) spike channels, but these arrays only record the  
    // counts for the first 128 channels.
    // Channel and unit numbers are 1-based - channel entries at [0] and [129] are 
    // unused, and unit entries at [0] are unused.
    int     TSCounts[130][5]; // number of timestamps[channel][unit]
    int     WFCounts[130][5]; // number of waveforms[channel][unit]

    // Starting at index 300, this array also records the number of samples for the 
    // continuous channels.  Note that since EVCounts has only 512 entries, continuous 
    // channels above channel 211 do not have sample counts.
    int     EVCounts[512];    // number of timestamps[event_number]
};


struct PL_ChanHeader 
{
    char    Name[32];       // Name given to the DSP channel
    char    SIGName[32];    // Name given to the corresponding SIG channel
    int     Channel;        // DSP channel number, 1-based
    int     WFRate;         // When MAP is doing waveform rate limiting, this is limit w/f per sec divided by 10
    int     SIG;            // SIG channel associated with this DSP channel 1 - based
    int     Ref;            // SIG channel used as a Reference signal, 1- based
    int     Gain;           // actual gain divided by SpikePreAmpGain. For pre version 105, actual gain divided by 1000. 
    int     Filter;         // 0 or 1
    int     Threshold;      // Threshold for spike detection in a/d values
    int     Method;         // Method used for sorting units, 1 - boxes, 2 - templates
    int     NUnits;         // number of sorted units
    short   Template[5][64];// Templates used for template sorting, in a/d values
    int     Fit[5];         // Template fit 
    int     SortWidth;      // how many points to use in template sorting (template only)
    short   Boxes[5][2][4]; // the boxes used in boxes sorting
    int     SortBeg;        // beginning of the sorting window to use in template sorting (width defined by SortWidth)
    char    Comment[128];   // Version >=105
    unsigned char SrcId;    // Version >=106, Omniplex Source ID for this channel
    unsigned char reserved; // Version >=106
    unsigned short ChanId;  // Version >=106, Omniplex Channel ID within the Source for this channel
    int     Padding[10];
};

struct PL_EventHeader 
{
    char    Name[32];       // name given to this event
    int     Channel;        // event number, 1-based
    char    Comment[128];   // Version >=105
    unsigned char SrcId;    // Version >=106, Omniplex Source ID for this channel
    unsigned char reserved; // Version >=106
    unsigned short ChanId;  // Version >=106, Omniplex Channel ID within the Source for this channel
    int     Padding[32];
};

struct PL_SlowChannelHeader 
{
    char    Name[32];       // name given to this channel
    int     Channel;        // channel number, 0-based
    int     ADFreq;         // digitization frequency
    int     Gain;           // gain at the adc card
    int     Enabled;        // whether this channel is enabled for taking data, 0 or 1
    int     PreAmpGain;     // gain at the preamp

    // As of Version 104, this indicates the spike channel (PL_ChanHeader.Channel) of
    // a spike channel corresponding to this continuous data channel. 
    // <=0 means no associated spike channel.
    int     SpikeChannel;

    char    Comment[128];   // Version >=105
    unsigned char SrcId;    // Version >=106, Omniplex Source ID for this channel
    unsigned char reserved; // Version >=106
    unsigned short ChanId;  // Version >=106, Omniplex Channel ID within the Source for this channel
    int     Padding[27];
};

// The header for the data record used in the datafile (*.plx)
// This is followed by NumberOfWaveforms*NumberOfWordsInWaveform
// short integers that represent the waveform(s)

struct PL_DataBlockHeader
{
    short   Type;                       // Data type; 1=spike, 4=Event, 5=continuous
    unsigned short   UpperByteOf5ByteTimestamp; // Upper 8 bits of the 40 bit timestamp
    unsigned long    TimeStamp;                 // Lower 32 bits of the 40 bit timestamp
    short   Channel;                    // Channel number
    short   Unit;                       // Sorted unit number; 0=unsorted
    short   NumberOfWaveforms;          // Number of waveforms in the data to folow, usually 0 or 1
    short   NumberOfWordsInWaveform;    // Number of samples per waveform in the data to follow
}; // 16 bytes










///////////////////////////////////////////////////////////////////////////////
// Plexon continuous data file (.DDT) File Structure Definitions
///////////////////////////////////////////////////////////////////////////////

#define LATEST_DDT_FILE_VERSION 103

struct DigFileHeader 
{
    int     Version;        // Version of the data format; determines which data items are valid
    int     DataOffset;     // Offset into the file where the data starts
    double  Freq;           // Digitization frequency
    int     NChannels;      // Number of recorded channels; for version 100-101, this will always
                            // be the same as the highest channel number recorded; for versions >= 102,
                            // NChannels is the same as the number of enabled channels, i.e. channels
                            // whose entry in the ChannelGain array is not 255.

    int     Year;           // Time/date when the data was acquired
    int     Month;
    int     Day;
    int     Hour;
    int     Minute;
    int     Second;
    
    int     Gain;           // As of Version 102, this is the *preamp* gain, not ADC gain
    char    Comment[128];   // User-supplied comment 
    unsigned char BitsPerSample;    // ADC resolution, usually either 12 or 16. Added for ddt Version 101    
    unsigned char ChannelGain[64];  // Gains for each channel; 255 means channel was disabled (not recorded). 
									// The gain for Channel n is located at ChannelGain[n-1]
									// Added for ddt Version 102 
    unsigned char Unused;           // padding to restore alignment 
    short         MaxMagnitudeMV;   // ADC max input voltage in millivolts: 5000 for NI, 2500 for ADS64
                                    // Added for ddt version 103
    unsigned char Padding[188];
};





#endif
