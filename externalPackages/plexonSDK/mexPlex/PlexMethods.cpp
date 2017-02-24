//
// PlexMethods.cpp
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


/******************************************************************
*
* PlexMethods.cpp contains implementation of Plexon .plx and .ddt file
* reading functions in Matlab.
* You may also want to look at the corresponding m-files.
*
*******************************************************************/
#include "PlexMethods.h"

#include <time.h>
#include <memory.h>
#include <limits>
#include <list>

#ifdef _MSC_VER
  // Visual C++
  #define fopenEx fopen
  #define fseekEx _fseeki64
  #define ftellEx _ftelli64
#elif defined( __APPLE__ )
  // Apple
  #include "TargetConditionals.h"
  #if TARGET_OS_MAC
    #define _FILE_OFFSET_BITS 64
    #define fopenEx fopen
    #define fseekEx fseeko
    #define ftellEx ftello
  #else
    // other Apple platform: use default functions
    #define fopenEx fopen
    #define fseekEx fseek
    #define ftellEx ftell
  #endif
#elif defined( __GNUC__ )
  // gcc on Linux
  #define _FILE_OFFSET_BITS 64
  #define fopenEx fopen64
  #define fseekEx fseeko
  #define ftellEx ftello
#else 
  // unknown platform: use default functions
  #define fopenEx fopen
  #define fseekEx fseek
  #define ftellEx ftell
#endif

using namespace std;

// make sure we have min defined
#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

// expected structure sizes
//assert(sizeof(PL_FileHeader) == 256 + 2*130*5*4 + 512*4);
//assert(sizeof(PL_ChanHeader) == 1020);
//assert(sizeof(PL_EventHeader) == 296);
//assert(sizeof(PL_SlowChannelHeader) == 296);
//assert(sizeof(PL_DataBlockHeader) == 16);

// we have a single CPlexMethods object in the library
static CPlexMethods s_plxMethods;

// helper functions

// combines error message with file name
void ReportError( const char* whatFailed, const char* fileName )
{
    string s = whatFailed;
    s += " '";
    s += fileName;
    s += "'.";
    mexErrMsgTxt( s.c_str() );
}

bool FunctionRequiresFileName( int functionIndex )
{
    return !( functionIndex == PLX_VT_INTERPRET || functionIndex == PLX_CLOSE  || functionIndex == GET_VERSION );
}

bool IsDdtFunction( int functionIndex )
{
    return ( functionIndex == DDT || functionIndex == DDT_V || functionIndex == DDT_WRITE );
}

// expectedNumInputPars is the number of parameters in plx_*.m call. 
// the actual call to mexFunction will have one extra first parameter - function index
// for example, plx_info(filename, fullread) is translated into 
// mexPlex(4, filename, fullread)
bool ExpectedNumberOfArguments( int expectedNumOutputPars, int expectedNumInputPars, 
                               int actualNumberOfLeftHandSideArguments, int actualNumberOfRightHandSideArguments )
{
    if ( actualNumberOfRightHandSideArguments != expectedNumInputPars + 1 ) {
        char buf[256];
        sprintf( buf, "Expected %d input argument(s).", expectedNumInputPars );
        mexErrMsgTxt( buf );
        return false;
    }
    if ( actualNumberOfLeftHandSideArguments != expectedNumOutputPars ) {
        char buf[256];
        sprintf( buf, "Expected %d output argument(s).", expectedNumInputPars );
        mexErrMsgTxt( buf );
        return false;
    }
    return true;
}

/******************************************************************
*
* mexFunction is the exported method which gives access to all
* functions in the DLL.
*
*****************************************************************/
extern "C" void mexFunction( int nlhs, mxArray *plhs[],
                            int nrhs, const mxArray*prhs[] )
{
    int iFunc;
    std::string fileName;
    bool bDdtFile = false;

    // assign all left hand side arguments to minus one
    // to avoid Matlab error messages about unassigned outputs
    for ( int i = 0; i < nlhs; i++ ) {
        plhs[i] = mxCreateDoubleScalar( -1 );
    }

	// Check if the number of right hand side parameters is less than 1.
    // The first parameter is the function number
    if ( nrhs < 1 ) {
        mexErrMsgTxt( "mex function number is not specified" );
        return;
    }

    // Extract the function index from the first parameter
    iFunc = ( int )mxGetScalar( prhs[0] );
    if ( ( iFunc >= LAST_FUNCTION ) || ( iFunc <= FIRST_FUNCTION ) ) {
        mexErrMsgTxt( "Function not supported." );
        return;
    }

    if ( iFunc == GET_VERSION ) {
        if ( nlhs > 0 ) {
            plhs[0] = mxCreateDoubleScalar( MEXPLEX_VERSION_NUMBER );
        }
        return;
    }

    // Check if the number of right hand side parameters is less than 2.
    // The first parameter is function number, the second -- file name
    // plx_*.m scripts add the first argument - the function index
    if ( nrhs < 2 ) {
        mexErrMsgTxt( "At least one input argument required." );
        return;
    }

    if ( FunctionRequiresFileName( iFunc ) ) {
        // the first argument in plx_*.m call (we see it as a second arg)
        // should be a string with the file name
        if ( !( mxIsChar( prhs[1] ) ) ) {
            mexErrMsgTxt( "The first argument must be of type string." );
            return;
        }
        char* fileNameFromArgument = mxArrayToString( prhs[1] );
        if ( fileNameFromArgument == NULL ) {
            mexErrMsgTxt( "The first argument must be of type string." );
            return;
        }
        fileName = fileNameFromArgument;
    }

    // check that file extensions are .plx or .ddt
    bool bReadOk = true;
    if ( FunctionRequiresFileName( iFunc ) ) {
        if ( fileName.length() == 0 ) {
            mexErrMsgTxt( "Invalid (empty) file name." );
            return;
        }
        if ( fileName.length() < 4 ) {
            mexErrMsgTxt( "Invalid file extension." );
            return;
        }
        // get file extension
        std::string ext = fileName.substr( fileName.length() - 4, 4 );
        std::transform( ext.begin(), ext.end(), ext.begin(), ::tolower );
        if ( ext == ".plx" ) {
            if ( IsDdtFunction( iFunc ) ) {
                mexErrMsgTxt( "Invalid file extension." );
                return;
            }
            // This will read the file headers into memory if they are not already there.
            // Global objects are persisted in memory between calls to mexPlex.
            // Note that we call plx_open for every plx_funcion requiring file name.
            // plx_open will do nothing if currently opened file name is fileName
            // we then call s_plxMethods.plx_ad etc. without the file name
            bReadOk = s_plxMethods.plx_open( fileName.c_str() );
            // check if there is file reading error.
            if ( !bReadOk ) {
                ReportError( "Error while reading file", fileName.c_str() );
                return;
            }
        } else if ( ext == ".ddt" ) {
            if ( !IsDdtFunction( iFunc ) ) {
                mexErrMsgTxt( "Invalid file extension." );
                return;
            }
            bDdtFile = true;
        } else {
            mexErrMsgTxt( "Invalid file extension." );
            return;
        }
    }

    // This switch statement controls calling of all
    // plx_* and ddt_* methods implemented in this file.
    switch ( iFunc ) {
        case DDT:  {
            // read ddt file
            // [nch, npoints, freq, d] = ddt(filename)
            if ( !ExpectedNumberOfArguments( 4, 1, nlhs, nrhs ) ) {
                return;
            }
            if ( bDdtFile ) {
                int numberOfChannels = 0;
                unsigned int dataCount = 0;
                int ADfrequency = 0;
                s_plxMethods.ddt_read( fileName.c_str(), &plhs[3], false, numberOfChannels, dataCount, ADfrequency );
                plhs[0] = mxCreateDoubleScalar( numberOfChannels );
                plhs[1] = mxCreateDoubleScalar( ( double )dataCount );
                plhs[2] = mxCreateDoubleScalar( ADfrequency );
            } else {
                mexErrMsgTxt( "Invalid file extension." );
                return;
            }
                   }
                   break;
        case DDT_WRITE: {
            // write ddt file
            // [errCode] = ddt_write_v(filename, nch, npoints, freq, d)
            if ( !ExpectedNumberOfArguments( 1, 5, nlhs, nrhs ) ) {
                return;
            }
            int nch = static_cast<int>( mxGetScalar( prhs[2] ) );
            if ( nch < 1 ) {
                mexErrMsgTxt( "Number of channels should be positive." );
                return;
            }
            if ( nch > 64 ) {
                mexErrMsgTxt( "Number of channels cannot be more than 64." );
                return;
            }
            int npoints = static_cast<int>( mxGetScalar( prhs[3] ) );
            if ( npoints < 1 ) {
                mexErrMsgTxt( "Number of points should be positive." );
                return;
            }
            double freq = mxGetScalar( prhs[4] );
            if ( freq <= 0 ) {
                mexErrMsgTxt( "Frequency should be positive." );
                return;
            }
            if ( bDdtFile ) {
                plhs[0] = mxCreateDoubleScalar( s_plxMethods.ddt_write( fileName.c_str(), nch, npoints, freq, &prhs[5] ) );
            } else {
                mexErrMsgTxt( "Invalid file extension: .ddt was expected." );
                return;
            }
                        }
                        break;
        case PLX_AD: {
            // function plx_ad.m
            // [adfreq, n, ts, fn, ad] = plx_ad(filename, ch)
            if ( !ExpectedNumberOfArguments( 5, 2, nlhs, nrhs ) ) {
                return;
            }
            int channel = static_cast<int>( mxGetScalar( prhs[2] ) );
            int ADFrequency = 0;
            size_t nCount = s_plxMethods.plx_ad( channel, &plhs[2], &plhs[3], &plhs[4], CPlexMethods::RAW, ADFrequency );
            plhs[0] = mxCreateDoubleScalar( ADFrequency );
            plhs[1] = mxCreateDoubleScalar( ( double )nCount );
                     }
                     break;
        case PLX_AD_V: {
            // function plx_ad_v.m
            // [adfreq, n, ts, fn, ad] = plx_ad_v(filename, ch)
            if ( !ExpectedNumberOfArguments( 5, 2, nlhs, nrhs ) ) {
                return;
            }
            int channel = static_cast<int>( mxGetScalar( prhs[2] ) );
            int ADFrequency = 0;
            size_t nCount = s_plxMethods.plx_ad( channel, &plhs[2], &plhs[3], &plhs[4], CPlexMethods::MV, ADFrequency );
            plhs[0] = mxCreateDoubleScalar( ADFrequency );
            plhs[1] = mxCreateDoubleScalar( ( double )nCount );
                       }
                       break;
        case PLX_AD_GAP_INFO: {
            // function plx_ad_gap_info.m
            // [adfreq, n, ts, fn] = plx_ad_gap_info(filename, ch)
            if ( !ExpectedNumberOfArguments( 4, 2, nlhs, nrhs ) ) {
                return;
            }
            int channel = static_cast<int>( mxGetScalar( prhs[2] ) );
            int ADFrequency = 0;
            size_t nCount = s_plxMethods.plx_ad( channel, &plhs[2], &plhs[3], NULL, CPlexMethods::NONE, ADFrequency );
            plhs[0] = mxCreateDoubleScalar( ADFrequency );
            plhs[1] = mxCreateDoubleScalar( ( double )nCount );
                              }
                              break;
        case PLX_EVENT_TS: {
            // function plx_event_ts.m
            // [n, ts, sv] = plx_event_ts(filename, ch)
            if ( !ExpectedNumberOfArguments( 3, 2, nlhs, nrhs ) ) {
                return;
            }
            int channel = static_cast<int>( mxGetScalar( prhs[2] ) );
            size_t nCount = s_plxMethods.plx_event_ts( channel, &plhs[1], &plhs[2] );
            plhs[0] = mxCreateDoubleScalar( ( double )nCount );
                           }
                           break;
        case PLX_INFO: {
            // function plx_info.m
            // [tscounts, wfcounts, evcounts, contcounts] = plx_info(filename, fullread)
            if ( !ExpectedNumberOfArguments( 4, 2, nlhs, nrhs ) ) {
                return;
            }
            int fullRead = static_cast<int>( mxGetScalar( prhs[2] ) );
            s_plxMethods.plx_info( &plhs[0], &plhs[1], &plhs[2], &plhs[3], fullRead );
                       }
                       break;
        case PLX_CHANMAP: {
            // function plx_chanmap.m
            // [n, dspchans] = plx_chanmap(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_chanmap( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                          }
                          break;
        case PLX_AD_CHANMAP: {
            // function plx_ad_chanmap.m
            // [n, adchans] = plx_ad_chanmap(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_ad_chanmap( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                             }
                             break;
        case PLX_EV_CHANMAP: {
            // function plx_ev_chanmap.m
            // [n, evchans] = plx_event_chanmap(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_ev_chanmap( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                             }
                             break;
        case PLX_TS: {
            // function plx_ts.m
            // [n, ts] = plx_ts(filename, channel, unit)
            if ( !ExpectedNumberOfArguments( 2, 3, nlhs, nrhs ) ) {
                return;
            }
            int channel = static_cast<int>( mxGetScalar( prhs[2] ) );
            int unit = static_cast<int>( mxGetScalar( prhs[3] ) );
            size_t nCount = s_plxMethods.plx_ts( channel, unit, &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( ( double )nCount );
                     }
                     break;
        case PLX_WAVES: {
            // function plx_waves.m
            // [n, npw, ts, wave] = plx_waves(filename, channel, unit)
            if ( !ExpectedNumberOfArguments( 4, 3, nlhs, nrhs ) ) {
                return;
            }
            int channel = static_cast<int>( mxGetScalar( prhs[2] ) );
            int unit = static_cast<int>( mxGetScalar( prhs[3] ) );
            size_t nCount = s_plxMethods.plx_waves( channel, unit, &plhs[2], &plhs[3] );
            plhs[0] = mxCreateDoubleScalar( ( double )nCount );
            plhs[1] = mxCreateDoubleScalar( s_plxMethods.GetNumberOfPointsInWaveform() );
                        }
                        break;
        case PLX_WAVES_V: {
            // function plx_waves_v.m
            // [n, npw, ts, wave] = plx_waves_v(filename, channel, unit)
            if ( !ExpectedNumberOfArguments( 4, 3, nlhs, nrhs ) ) {
                return;
            }
            int channel = static_cast<int>( mxGetScalar( prhs[2] ) );
            int unit = static_cast<int>( mxGetScalar( prhs[3] ) );
            size_t nCount = s_plxMethods.plx_waves( channel, unit, &plhs[2], &plhs[3], true );
            plhs[0] = mxCreateDoubleScalar( ( double )nCount );
            plhs[1] = mxCreateDoubleScalar( s_plxMethods.GetNumberOfPointsInWaveform() );
                          }
                          break;
        case PLX_AD_SPAN: {
            // function plx_ad_span.m
            // [adfreq, n, ad] = plx_ad_span(filename, channel, startCount, endCount)
            if ( !ExpectedNumberOfArguments( 3, 4, nlhs, nrhs ) ) {
                return;
            }
            int channel = static_cast<int>( mxGetScalar( prhs[2] ) );
            long long startCount = static_cast<long long>( mxGetScalar( prhs[3] ) );
            long long endCount = static_cast<long long>( mxGetScalar( prhs[4] ) );
            int ADFrequency = 0;
            size_t nCount = s_plxMethods.plx_ad_span( channel, startCount, endCount, &plhs[2], false, ADFrequency );
            plhs[0] = mxCreateDoubleScalar( ADFrequency );
            plhs[1] = mxCreateDoubleScalar( ( double )nCount );
                          }
                          break;
        case PLX_AD_SPAN_V: {
            // function plx_ad_span_v.m
            // [adfreq, n, ad] = plx_ad_span_v(filename, channel, startCount,endCount)
            if ( !ExpectedNumberOfArguments( 3, 4, nlhs, nrhs ) ) {
                return;
            }
            int channel = static_cast<int>( mxGetScalar( prhs[2] ) );
            long long startCount = static_cast<long long>( mxGetScalar( prhs[3] ) );
            long long endCount = static_cast<long long>( mxGetScalar( prhs[4] ) );
            int ADFrequency = 0;
            size_t nCount = s_plxMethods.plx_ad_span( channel, startCount, endCount, &plhs[2], true, ADFrequency );
            plhs[0] = mxCreateDoubleScalar( ADFrequency );
            plhs[1] = mxCreateDoubleScalar( ( double )nCount );
                            }
                            break;
        case PLX_CHAN_GAINS: {
            // function plx_chan_gains.m
            // [n, gains] = plx_chan_gains(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_chan_gains( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                             }
                             break;
        case PLX_CHAN_THRESHOLD: {
            // function plx_chan_threshold.m
            // [n, thresholds] = plx_chan_thresholds(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_chan_thresholds( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                                 }
                                 break;
        case PLX_CHAN_FILTERS: {
            // function plx_chan_filters.m
            // [n, filters] = plx_chan_filters(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_chan_filters( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                               }
                               break;
        case PLX_ADCHAN_GAINS: {
            // function plx_adchan_gains.m
            // [n, gains] = plx_adchan_gains(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_adchan_gains( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                               }
                               break;
        case PLX_ADCHAN_FREQS: {
            // function plx_adchan_freqs.m
            // [n, freqs] = plx_adchan_freqs(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_adchan_freqs( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                               }
                               break;
        case PLX_ADCHAN_SAMPLECOUNTS: {
            // function plx_adchan_samplecounts.m
            // [n, samplecounts] = plx_adchan_samplecounts(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_adchan_samplecounts( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                                      }
                                      break;
        case PLX_INFORMATION: {
            // function plx_information.m
            // [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(filename)
            if ( !ExpectedNumberOfArguments( 13, 1, nlhs, nrhs ) ) {
                return;
            }
            char DateTime[50];
            int Version, Freq, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits;
            int SlowPeakV, SlowADResBits;
            double Duration;
            char* szComment;

            s_plxMethods.plx_information( &Version, &Freq, &szComment, &Trodalness,
                &NPW, &PreTresh, &SpikePeakV, &SpikeADResBits, &SlowPeakV,
                &SlowADResBits, &Duration, DateTime, sizeof( DateTime ) );

            plhs[0] = mxCreateString( s_plxMethods.GetCurrentFileName() );
            plhs[1] = mxCreateDoubleScalar( Version );
            plhs[2] = mxCreateDoubleScalar( Freq );
            plhs[3] = mxCreateString( szComment );
            plhs[4] = mxCreateDoubleScalar( Trodalness );
            plhs[5] = mxCreateDoubleScalar( NPW );
            plhs[6] = mxCreateDoubleScalar( PreTresh );
            plhs[7] = mxCreateDoubleScalar( SpikePeakV );
            plhs[8] = mxCreateDoubleScalar( SpikeADResBits );
            plhs[9] = mxCreateDoubleScalar( SlowPeakV );
            plhs[10] = mxCreateDoubleScalar( SlowADResBits );
            plhs[11] = mxCreateDoubleScalar( Duration );
            plhs[12] = mxCreateString( DateTime );

                              }
                              break;
        case PLX_CHAN_NAMES: {
            // function plx_chan_names.m
            // [n, names] = plx_chan_names(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_chan_names( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                             }
                             break;
        case PLX_ADCHAN_NAMES: {
            // function plx_adchan_names.m
            // [n, names] = plx_adchan_names(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_adchan_names( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                               }
                               break;
        case PLX_EVENT_NAMES: {
            // function plx_event_names.m
            // [n, names] = plx_event_names(filename)
            if ( !ExpectedNumberOfArguments( 2, 1, nlhs, nrhs ) ) {
                return;
            }
            int nCount = s_plxMethods.plx_event_names( &plhs[1] );
            plhs[0] = mxCreateDoubleScalar( nCount );
                              }
                              break;

        case DDT_V:  {
            // read ddt file
            // [nch, npoints, freq, d] = ddt_v(filename)
            if ( !ExpectedNumberOfArguments( 4, 1, nlhs, nrhs ) ) {
                return;
            }
            if ( bDdtFile ) {
                int numberOfChannels = 0;
                unsigned int dataCount = 0;
                int ADfrequency = 0;
                s_plxMethods.ddt_read( fileName.c_str(), &plhs[3], true, numberOfChannels, dataCount, ADfrequency );
                plhs[0] = mxCreateDoubleScalar( numberOfChannels );
                plhs[1] = mxCreateDoubleScalar( ( double )dataCount );
                plhs[2] = mxCreateDoubleScalar( ADfrequency );
            } else {
                mexErrMsgTxt( "Invalid file extension." );
                return;
            }
                     }
                     break;
        case PLX_VT_INTERPRET: {
            // [nCoords, nDim, nVTMode, c] = plx_vt_interpret(ts, sv)
            // translated into:
            // [nCoords, nDim, nVTMode, c] = mexPlex(21, '', ts, sv);
            if ( nrhs != 4 ) {
                mexErrMsgTxt( "Expected 2 input arguments." );
                break;
            }
            if ( nlhs != 4 ) {
                mexErrMsgTxt( "Expected 4 output arguments." );
                break;
            }
            s_plxMethods.vt_interpret( prhs[2], prhs[3], &plhs[0], &plhs[1], &plhs[2], &plhs[3] );
                               }
                               break;
        case PLX_CLOSE: {
            // [n] = plx_close(filename)
            if ( !ExpectedNumberOfArguments( 1, 1, nlhs, nrhs ) ) {
                return;
            }
            s_plxMethods.plx_close();
            plhs[0] = mxCreateDoubleScalar( 0.0 );
                        }
                        break;
    }
}

// CPlexMethods implementation

CPlexMethods::CPlexMethods() : m_bTallied( false ), m_Reader( NULL )
{
}

CPlexMethods::~CPlexMethods()
{
    if ( m_Reader != NULL ) {
        delete m_Reader;
        m_Reader = NULL;
    }
}

bool CPlexMethods::DoOpenFile( const char* fileName )
{
    if ( m_Reader == NULL ) {
        m_Reader = new Reader();
    }
    else {
        m_Reader->Close();
    }
    return m_Reader->Open(fileName);
}

// helper to transform data values from one unit to another
// (example: from ADC values to milliVolts)
template< typename Result = double >
struct  TransHelper {
    double  a;
    TransHelper( double v = 1 ) : a( v )  {}
    template< typename X >
    Result  operator()( X x ) const    { return static_cast<Result>( a * x ); }
};

bool CPlexMethods::plx_open( const char* fileName )
{
    // open file, if file is new or reader is null
    if ( ( m_CurrentFileName != fileName ) || m_Reader == NULL ) {
        m_CurrentFileName = "";
        if ( !DoOpenFile( fileName ) ) {
            ReportError( "Cannot open plx file", fileName );
            return false;
        } 

        m_PlxFileSize = m_Reader->GetFileSize();
        if ( m_PlxFileSize == 0 ) {
            ReportError( "plx file is empty or cannot get plx file size", fileName );
            return false;
        }

        // read file header
        if ( m_Reader->Read( &m_FileHeader, sizeof( m_FileHeader ) ) != sizeof( m_FileHeader ) ) {
            ReportError( "Cannot read plx file header", fileName );
            return false;
        }

        // check file header
        // we use short to represent channel number, so there should not be more than 32k channels of each type
        int maxNumberOfChannelHeaders = 32 * 1024;
        if ( m_FileHeader.MagicNumber != 0x58454c50
            || m_FileHeader.NumDSPChannels < 0
            || m_FileHeader.NumDSPChannels > maxNumberOfChannelHeaders
            || m_FileHeader.NumEventChannels < 0
            || m_FileHeader.NumEventChannels > maxNumberOfChannelHeaders
            || m_FileHeader.NumSlowChannels < 0
            || m_FileHeader.NumSlowChannels > maxNumberOfChannelHeaders
            || m_FileHeader.ADFrequency <= 0
            || m_FileHeader.NumPointsWave < 0
            ) {
                ReportError( "plx file header is invalid", fileName );
                return false;
        }

        // read channel headers
        if ( m_FileHeader.NumDSPChannels > 0 ) {
            m_SpikeChannelHeaders.resize( m_FileHeader.NumDSPChannels );
            if ( m_Reader->Read( &m_SpikeChannelHeaders[0], m_FileHeader.NumDSPChannels * sizeof( PL_ChanHeader ) ) != m_FileHeader.NumDSPChannels * sizeof( PL_ChanHeader ) ) {
                ReportError( "Cannot read plx spike channel headers", fileName );
                return false;
            }
        }
        if ( m_FileHeader.NumEventChannels > 0 ) {
            m_EventChannelHeaders.resize( m_FileHeader.NumEventChannels );
            if ( m_Reader->Read( &m_EventChannelHeaders[0], m_FileHeader.NumEventChannels * sizeof( PL_EventHeader ) ) != m_FileHeader.NumEventChannels * sizeof( PL_EventHeader ) ) {
                ReportError( "Cannot read plx event channel headers", fileName );
                return false;
            }
        }
        if ( m_FileHeader.NumSlowChannels > 0 ) {
            m_AnalogChannelHeaders.resize( m_FileHeader.NumSlowChannels );
            if ( m_Reader->Read( &m_AnalogChannelHeaders[0], m_FileHeader.NumSlowChannels * sizeof( PL_SlowChannelHeader ) ) !=  m_FileHeader.NumSlowChannels * sizeof( PL_SlowChannelHeader ) ) {
                ReportError( "Cannot read plx analog channel headers", fileName );
                return false;
            }
        }

        // TODO: verify channel headers

        // calculate channel number to header index maps
        m_SpikeChannelNumberToSpikeChannelHeaderIndex.clear();
        for ( int i = 0; i < ( int )m_SpikeChannelHeaders.size(); ++i ) {
            if ( m_SpikeChannelHeaders[i].Channel < 0 ) continue;
            m_SpikeChannelNumberToSpikeChannelHeaderIndex[ m_SpikeChannelHeaders[i].Channel ] = i;
        }

        m_EventChannelNumberToEventChannelHeaderIndex.clear();
        for ( int i = 0; i < ( int )m_EventChannelHeaders.size(); ++i ) {
            if ( m_EventChannelHeaders[i].Channel < 0 ) continue;
            m_EventChannelNumberToEventChannelHeaderIndex[ m_EventChannelHeaders[i].Channel ] = i;
        }

        m_AnalogChannelNumberToAnalogChannelHeaderIndex.clear();
        for ( int i = 0; i < ( int )m_AnalogChannelHeaders.size(); ++i ) {
            if ( m_AnalogChannelHeaders[i].Channel < 0 ) continue;
            m_AnalogChannelNumberToAnalogChannelHeaderIndex[ m_AnalogChannelHeaders[i].Channel ] = i;
        }

        // transfer into m_TSCounts, m_WFCounts so that things work if fullread==0
        for ( int i = 0; i < PLX_HDR_LAST_SPIKE_CHAN; i++ ) {
            for ( int j = 0; j <= PLX_HDR_LAST_UNIT; j++ ) {
                // note the difference in channel basing between m_XXCounts and m_FileHeader.XXCounts
                m_WaveformCounts[make_pair( i, j )] = m_FileHeader.WFCounts[i + 1][j];
                m_TimestampCounts[make_pair( i, j )] = m_FileHeader.TSCounts[i + 1][j];
            }
        }

        // calculate the data offset
        m_DataStartFilePosition = sizeof( m_FileHeader ) + m_FileHeader.NumDSPChannels * sizeof( PL_ChanHeader )
            + m_FileHeader.NumEventChannels * sizeof( PL_EventHeader )
            + m_FileHeader.NumSlowChannels * sizeof( PL_SlowChannelHeader );

        m_bTallied = false;

        // if open file succeed, save the name as current file name
        m_CurrentFileName = fileName;
    }
    return true;
}

void CPlexMethods::plx_close()
{
    if( m_Reader != NULL ) {
        m_Reader->Close();
        delete m_Reader;
        m_Reader = NULL;
    }
    m_CurrentFileName = "";
}

bool CPlexMethods::ValidateDataBlockHeader( const PL_DataBlockHeader &db ) const
{
    if ( !( db.Type == PL_SingleWFType || db.Type == PL_ExtEventType  || db.Type == PL_ADDataType ) ) {
        mexPrintf( "\ninvalid data block header type\n" );
        return false;
    }
    if ( db.Channel < 0 || (db.Type == PL_SingleWFType && db.Unit < 0)
        || db.NumberOfWaveforms < 0 
        || db.NumberOfWordsInWaveform < 0 
        || db.NumberOfWaveforms * db.NumberOfWordsInWaveform > MAX_NUMBER_OF_WORDS_FOLLOWING_DATA_BLOCK_HEADER ) {
            mexPrintf( "\ninvalid data block header\n" );
            return false;
    }
    return true;
}

bool CPlexMethods::ReadNextBlock( PL_DataBlockHeader& db, std::vector<short>& buf )
{
    buf.clear();
    if ( m_Reader == NULL ) {
        return false;
    }
    if ( m_Reader->IsAtEndOfFile() ) {
        return false;
    }
    // read data block header
    if ( m_Reader->Read( &db, sizeof( db ) ) != sizeof( db ) ) { 
        return false;
    }

    if ( !ValidateDataBlockHeader( db ) ) {
        return false;
    }

    // read the waveform after the block header
    if ( db.NumberOfWaveforms > 0 ) { 
        if ( m_Reader->IsAtEndOfFile() ) {
            return false;
        }
        unsigned int nbuf = db.NumberOfWaveforms * db.NumberOfWordsInWaveform;
        buf.resize( nbuf );
        if ( m_Reader->Read( &buf[0], nbuf * 2 ) != nbuf * 2 ) {
            return false;
        }
    }
    return true;
}

bool CPlexMethods::SeekToDataStart()
{
    if ( m_Reader == NULL ) {
        return false;
    }
    return m_Reader->Seek( m_DataStartFilePosition );
}

bool CPlexMethods::tallyPlxFile( )
{
    m_WaveformCounts.clear();
    m_TimestampCounts.clear();
    m_AnalogSampleCounts.clear();
    m_EventCounts.clear();

    if ( !SeekToDataStart() ) {
        return false;
    }

    // Start reading the data
    PL_DataBlockHeader db;
    std::vector<short> buf;
    buf.reserve( 512 );

    while ( ReadNextBlock( db, buf ) ) { 
        int type = db.Type;
        int chan = db.Channel;  // 1 based
        int unit = db.Unit;     // 0 means unsorted

        if ( type == PL_SingleWFType ) {
            m_TimestampCounts[make_pair( chan - 1, unit )]++;
            if ( buf.size() > 0 ) {
                m_WaveformCounts[make_pair( chan - 1, unit )]++;
            }
        }     
        else if ( type == PL_ExtEventType ) {
            // these arrays are indexed by the channel number in the data block
            m_EventCounts[chan]++;
        } else if ( type == PL_ADDataType ) {
            // likewise, these arrays are indexed by the channel number in the data block
            m_AnalogSampleCounts[chan] += buf.size();
        }
    }

    m_bTallied = true;
    return true;
}

void CPlexMethods::plx_info( mxArray **ppmxTs, mxArray **ppmxWf, mxArray **ppmxEv, mxArray **ppmxCont, int fullRead )
{
    double *pdTsData, *pdWfData, *pdEvData, *pdContData;
    int nUnits;

    nUnits = MAX_NUM_UNITS + 1;
    if ( fullRead == 0 ) nUnits = 5;

    if ( ( fullRead == 1 ) && !m_bTallied ) {
        tallyPlxFile();
    }

    // put timestamp data into the mxArray
    *ppmxTs = mxCreateDoubleMatrix( nUnits, m_SpikeChannelHeaders.size() + 1, mxREAL );
    pdTsData = mxGetPr( *ppmxTs );

    // Note that the tscounts and wfcounts arrays are returned indexed by header offset
    int k = 0;
    for ( int i = 0; i < ( int )m_SpikeChannelHeaders.size(); ++i ) {
        int dspchan = m_SpikeChannelHeaders[i].Channel - 1;
        for ( int j = 0; j < nUnits; ++j, ++k ) {
            // the +nUnits to is offset the channel numbers by one, for backwards compatibility
            *( pdTsData + k + nUnits ) = ( double )GetTimestampCount( dspchan, j );
        }
    }

    // put waveforms data into the mxArray
    *ppmxWf = mxCreateDoubleMatrix( nUnits, m_SpikeChannelHeaders.size() + 1, mxREAL );
    pdWfData = mxGetPr( *ppmxWf );

    k = 0;
    for ( int i = 0; i < ( int )m_SpikeChannelHeaders.size(); ++i ) {
        for ( int j = 0; j < nUnits; ++j, ++k ) {
            // the +nUnits to is offset the channel numbers by one, for backwards compatibility
            // GetWaveformCount expects zero-based channel number. spike channels are one-based, so we subtract 1
            *( pdWfData + k + nUnits ) = ( double )GetWaveformCount( m_SpikeChannelHeaders[i].Channel - 1, j );
        }
    }

    // put event data into the mxArray
    *ppmxEv = mxCreateDoubleMatrix( 1, m_EventChannelHeaders.size(), mxREAL );
    pdEvData = mxGetPr( *ppmxEv );
    for ( int i = 0; i < ( int )m_EventChannelHeaders.size(); ++i ) {
        // find the event number (from the data blocks) that corresponds to this index
        *( pdEvData + i ) = ( double )GetEventChannelSampleCount( m_EventChannelHeaders[i].Channel );
    }

    // and analog data
    *ppmxCont = mxCreateDoubleMatrix( 1, m_AnalogChannelHeaders.size(), mxREAL );
    pdContData = mxGetPr( *ppmxCont );

    for ( int i = 0; i < ( int )m_AnalogChannelHeaders.size(); ++i ) {
        // find the channel number (from the data blocks) that corresponds to this index
        int analogChannel = m_AnalogChannelHeaders[i].Channel;
        *( pdContData + i ) = ( double )GetAnalogChannelSampleCount( analogChannel );
    }
}

int CPlexMethods::plx_chanmap( mxArray **ppmxChanmap )
{
    int nCount = ( int )m_SpikeChannelHeaders.size();
    if ( nCount == 0 ) {
        mexPrintf( "\nplx_chanmap: file contains no spike channels\n" );
        return 0;
    }

    // one entry in the array for each channel header
    *ppmxChanmap = mxCreateDoubleMatrix( 1, nCount, mxREAL );
    double* pdChanmapData = mxGetPr( *ppmxChanmap );

    for ( int i = 0; i < nCount; ++i ) {
        *( pdChanmapData + i ) = m_SpikeChannelHeaders[i].Channel;
    }
    return nCount;
}

int CPlexMethods::plx_ad_chanmap( mxArray **ppmxChanmap )
{
    int nCount = ( int ) m_AnalogChannelHeaders.size();
    if ( nCount == 0 ) {
        mexPrintf( "\nplx_ad_chanmap: file contains no A/D channels\n" );
        return 0;
    }

    // one entry in the array for each channel header
    *ppmxChanmap = mxCreateDoubleMatrix( 1, nCount, mxREAL );
    double* pdChanmapData = mxGetPr( *ppmxChanmap );

    for ( int i = 0; i < nCount; ++i ) {
        *( pdChanmapData + i ) = m_AnalogChannelHeaders[i].Channel;
    }
    return nCount;
}

int CPlexMethods::plx_ev_chanmap( mxArray **ppmxChanmap )
{
    int nCount = ( int )m_EventChannelHeaders.size();
    if ( nCount == 0 ) {
        mexPrintf( "\nplx_ev_chanmap: file contains no event channels\n" );
        return 0;
    }

    // one entry in the array for each channel header
    *ppmxChanmap = mxCreateDoubleMatrix( 1, nCount, mxREAL );
    double* pdChanmapData = mxGetPr( *ppmxChanmap );

    for ( int i = 0; i < nCount; ++i ) {
        *( pdChanmapData + i ) = m_EventChannelHeaders[i].Channel;
    }
    return nCount;
}

double CPlexMethods::GetTimestamp( const PL_DataBlockHeader &db ) const
{
    return ( ( double )toTimeStampType( db.UpperByteOf5ByteTimestamp, db.TimeStamp ) ) / m_FileHeader.ADFrequency;
}

size_t CPlexMethods::plx_ts( int channel_to_extract, int unit_to_extract, mxArray **ppmxTs )
{
    // this routine relies on accurate tallies
    if ( !m_bTallied ) {
        if ( !tallyPlxFile()){
            return 0;
        }
    }

    long long longCount = GetTimestampCount( channel_to_extract - 1, unit_to_extract );

    if ( longCount == 0 ) {
        mexPrintf( "\nplx_ts: no timestamps for channel %d and unit %d.\n", channel_to_extract, unit_to_extract );
        return 0;
    }

    if ( longCount > MAX_NUMBER_OF_VALUES_IN_ONE_MATLAB_MATRIX || longCount > MAX_SIZE_T_VALUE ) {
        mexPrintf( "\nplx_ts: too many timestamp values for channel %d, unit %d.\n", channel_to_extract, unit_to_extract );
        return 0;
    }

    if ( !SeekToDataStart() ) {
        return 0;
    }

    size_t nCount = ( size_t )longCount;

    *ppmxTs = mxCreateDoubleMatrix( nCount, 1, mxREAL );
    double* pdTsData = mxGetPr( *ppmxTs );

    PL_DataBlockHeader db;
    std::vector<short> buf;
    buf.reserve( 512 );
    size_t i = 0;

    while ( ReadNextBlock( db, buf ) ) { 
        if ( db.Type == PL_SingleWFType ) {
            if ( ( db.Channel == channel_to_extract ) && ( db.Unit == unit_to_extract ) ) {
                *( pdTsData + i ) = GetTimestamp( db );
                i++;
            }
        }
    }
    return nCount;
}

size_t CPlexMethods::plx_waves( int channel_to_extract, int unit_to_extract, mxArray **ppmxTs, mxArray **ppmxWf, bool mV )
{
    double* pdTsData;
    double* pdWfData;
    bool createdMatrix = false;

    // this routine relies on accurate tallies
    if ( !m_bTallied ) {
        if ( !tallyPlxFile() ){
            return 0;
        }
    }

    int index = -1;
    if ( !GetSpikeHeaderIndexFromSpikeChannelNumber( channel_to_extract, index ) ) {
        mexPrintf( "\nplx_waves: no spike channel header for channel %d.\n", channel_to_extract );
        return 0;
    }

    long long longCount = GetWaveformCount( channel_to_extract - 1, unit_to_extract );

    if ( longCount == 0 ) {
        mexPrintf( "\nplx_waves: no waveforms for channel %d and unit %d.\n", channel_to_extract, unit_to_extract );
        return 0;
    }

    if ( longCount * m_FileHeader.NumPointsWave > MAX_NUMBER_OF_VALUES_IN_ONE_MATLAB_MATRIX || longCount * m_FileHeader.NumPointsWave > MAX_SIZE_T_VALUE ) {
        mexPrintf( "\nplx_waves: too many waveform values for channel %d and unit %d.\n", channel_to_extract, unit_to_extract );
        return 0;
    }

    if ( !SeekToDataStart() ) {
        return 0;
    }

    size_t nCount = ( size_t )longCount;

    *ppmxTs = mxCreateDoubleMatrix( nCount, 1, mxREAL );
    pdTsData = mxGetPr( *ppmxTs );

    double  to_mV = mV ? calculateRawToMilliVoltsCoeffForSpikeChannel( index ) : 1.0;

    PL_DataBlockHeader db;
    std::vector<short> buf;
    buf.reserve( 512 );
    size_t j = 0;

    while ( ReadNextBlock( db, buf ) ) { 
        if ( db.Type == PL_SingleWFType && buf.size() > 0 ) {
            if ( db.Channel == channel_to_extract && db.Unit == unit_to_extract ) {
                if ( !createdMatrix ) {
                    *ppmxWf = mxCreateDoubleMatrix( nCount, m_FileHeader.NumPointsWave , mxREAL );
                    pdWfData = mxGetPr( *ppmxWf );
                    createdMatrix = true;
                }

                *( pdTsData + j ) = GetTimestamp( db );

                for ( int i = 0; i < m_FileHeader.NumPointsWave ; i++ ) {
                    *( pdWfData + j + ( i * nCount ) ) = to_mV * buf[i];
                }
                j++;
            }
        }
    }
    return nCount;
}

size_t CPlexMethods::plx_event_ts( int event_to_extract, mxArray **ppmxTs, mxArray **ppmxEv )
{
    // this routine relies on accurate tallies
    if ( !m_bTallied ) {
        if ( !tallyPlxFile() ){
            return 0;
        }
    }

    long long longCount = GetEventChannelSampleCount( event_to_extract );

    if ( longCount == 0 ) {
        mexPrintf( "\nplx_event_ts: event data for channel %d not found.\n", event_to_extract );
        return 0;
    }

    if ( longCount > MAX_NUMBER_OF_VALUES_IN_ONE_MATLAB_MATRIX || longCount > MAX_SIZE_T_VALUE ) {
        mexPrintf( "\nplx_event_ts: too many data values for event channel %d.\n", event_to_extract );
        return 0;
    }

    if ( !SeekToDataStart() ) {
        return 0;
    }

    size_t nCount = ( size_t )longCount;

    *ppmxEv = mxCreateDoubleMatrix( nCount, 1, mxREAL );
    double* pdEvData = mxGetPr( *ppmxEv );
    *ppmxTs = mxCreateDoubleMatrix( nCount, 1, mxREAL );
    double* pdTsData = mxGetPr( *ppmxTs );

    PL_DataBlockHeader db;
    std::vector<short> buf;
    buf.reserve( 512 );
    size_t dataIndex = 0;

    while ( ReadNextBlock( db, buf ) ) { 
        if ( db.Type == PL_ExtEventType ) {
            if ( db.Channel == event_to_extract ) {
                *( pdEvData + dataIndex ) = db.Unit;
                *( pdTsData + dataIndex ) = GetTimestamp( db );
                dataIndex++;
            }
        }
    }
    return nCount;
}

std::size_t CPlexMethods::plx_ad( int channel_to_extract, mxArray **ppmxTs, mxArray **ppmxFn, mxArray **ppmxWf, Unit mV, int& thisChannelADFrequency )
{
    double* pdWfData = NULL;
    int header_num = -1;

    if ( !GetAnalogHeaderIndexFromAnalogChannelNumber( channel_to_extract, header_num ) ) {
        mexPrintf( "\nplx_ad: no header for the specified A/D channel.\n" );
        return 0;
    }

    // this routine relies on accurate tallies
    if ( !m_bTallied ) {
        if ( !tallyPlxFile() ){
            return 0;
        }
    }

    long long longCount  = GetAnalogChannelSampleCount( channel_to_extract );
    if ( longCount == 0 ) {
        mexPrintf( "\nplx_ad: data for channel %d not found.\n", channel_to_extract );
        return 0;
    }
    if ( longCount > MAX_NUMBER_OF_VALUES_IN_ONE_MATLAB_MATRIX || longCount > MAX_SIZE_T_VALUE ) {
        mexPrintf( "\nplx_ad: too many data values for channel %d.\n", channel_to_extract );
        return 0;
    }

    size_t nCount = ( size_t )longCount;

    int gain = m_AnalogChannelHeaders[header_num].Gain;
    thisChannelADFrequency = m_AnalogChannelHeaders[header_num].ADFreq;
    if ( thisChannelADFrequency <= 0 ) {
        mexPrintf( "\nplx_ad: no A/D frequency for channel %d.\n", channel_to_extract );
        return 0;
    }
    if ( gain <= 0 ) {
        mexPrintf( "\nplx_ad: no gain for channel %d.\n", channel_to_extract );
        return 0;
    }

    if ( !SeekToDataStart() ) {
        return 0;
    }

    if ( mV != NONE ) {
        *ppmxWf = mxCreateDoubleMatrix( nCount, 1, mxREAL );
        pdWfData = mxGetPr( *ppmxWf );
    }

    std::list<double>   ts_list;
    std::list<unsigned> fn_list;
    TimeStampType       current_ts = -1; // fragment timestamp
    unsigned            current_fn = 0; // number of data points in this fragment

    TransHelper<> to_mV( calculateRawToMilliVoltsCoeffForAnalogChannel( header_num ) );

    PL_DataBlockHeader db;
    std::vector<short> buf;
    buf.reserve( 512 );

    while ( ReadNextBlock( db, buf ) ) { 
        if ( db.Type != PL_ADDataType || db.Channel != channel_to_extract ) {
            continue;
        }

        TimeStampType ts = toTimeStampType( db.UpperByteOf5ByteTimestamp, db.TimeStamp );
        // here we check if the current timestamp is at the expected distance from the fragment start
        // the condition is: (ts-start_ts)/fh.Adfr == num_points/channel_fr
        // we then multiply by fh.Adfr*channel_fr
        // if the timestamp is not where we expect it for uninterrupted recording,
        // this means that we have a gap in the data and we add a new fragment
        if ( ( ts - current_ts ) * thisChannelADFrequency != ( TimeStampType )current_fn * m_FileHeader.ADFrequency ) {
            if ( current_fn > 0 ) {
                ts_list.push_back( ( ( double )current_ts ) / m_FileHeader.ADFrequency );
                fn_list.push_back( current_fn );
            }
            current_ts = ts;
            current_fn = 0;
        }

        current_fn += ( unsigned int )buf.size();

        switch ( mV ) {
case RAW:
    pdWfData = std::copy( buf.begin(), buf.end(), pdWfData );
    break;
case MV:
    pdWfData = std::transform( buf.begin(), buf.end(), pdWfData, to_mV );
    break;
        }
    }

    // save the last run
    if ( current_fn > 0 ) {
        ts_list.push_back( ( ( double )current_ts ) / m_FileHeader.ADFrequency );
        fn_list.push_back( current_fn );
    }

    // create proper ts and fn arrays
    *ppmxTs = mxCreateDoubleMatrix( ( int )ts_list.size(), 1, mxREAL );
    std::copy( ts_list.begin(), ts_list.end(), mxGetPr( *ppmxTs ) );
    *ppmxFn = mxCreateDoubleMatrix( ( int )fn_list.size(), 1, mxREAL );
    std::copy( fn_list.begin(), fn_list.end(), mxGetPr( *ppmxFn ) );

    return nCount;
}

size_t CPlexMethods::plx_ad_span( int channel_to_extract, long long startCount, long long endCount, mxArray **ppmxWf, bool mV, int& ADFrequency )
{
    if ( startCount <= 0 ) {
        mexPrintf( "\nplx_ad_span: startCount should be positive.\n" );
        return 0;
    }
    if ( endCount <= 0 ) {
        mexPrintf( "\nplx_ad_span: endCount should be positive.\n" );
        return 0;
    }
    size_t dataSize = 0;
    if ( endCount > startCount ) {
        if ( ( endCount - startCount ) > MAX_NUMBER_OF_VALUES_IN_ONE_MATLAB_MATRIX || ( endCount - startCount ) > MAX_SIZE_T_VALUE ) {
            mexPrintf( "\nplx_ad_span: too many data values for channel %d.\n", channel_to_extract );
            return 0;
        }
        dataSize = ( size_t )( endCount - startCount + 1 );
    } else {
        mexPrintf( "\nplx_ad_span: invalid startCount and endCount values.\n" );
        mexPrintf( "plx_ad_span: endCount value should be greater than startCount value.\n" );
        return 0;
    }

    int header_num = -1;
    if ( !GetAnalogHeaderIndexFromAnalogChannelNumber( channel_to_extract, header_num ) ) {
        mexPrintf( "\nplx_ad: no header for the specified A/D channel.\n" );
        return 0;
    }
    // this routine relies on accurate tallies
    if ( !m_bTallied ) {
        if ( !tallyPlxFile() ){
            return 0;
        }
    }
    long long nCount  = GetAnalogChannelSampleCount( channel_to_extract );
    if ( nCount == 0 ) {
        mexPrintf( "\nplx_ad_span: data for channel %d not found.\n", channel_to_extract );
        return 0;
    }
    if ( startCount > nCount || endCount > nCount ) {
        mexPrintf( "\nplx_ad_span: startCount and endCount should be less than or equal to the number of samples.\n" );
        return 0;
    }

    int gain = m_AnalogChannelHeaders[header_num].Gain;
    ADFrequency = m_AnalogChannelHeaders[header_num].ADFreq;
    if ( ADFrequency <= 0 || gain <= 0 ) {
        mexPrintf( "\nplx_ad_span: no A/D frequency or gain for the specified A/D channel.\n" );
        return 0;
    }

    if ( !SeekToDataStart() ) {
        return 0;
    }

    *ppmxWf = mxCreateDoubleMatrix( dataSize, 1, mxREAL );
    double* pdWfData = mxGetPr( *ppmxWf );

    double  to_mV = mV ? calculateRawToMilliVoltsCoeffForAnalogChannel( header_num ) : 1.0;

    long long k = 1;
    size_t counter = 0;

    PL_DataBlockHeader db;
    std::vector<short> buf;
    buf.reserve( 512 );

    while ( ReadNextBlock( db, buf ) ) { 
        if ( db.Type == PL_ADDataType ) {
            if ( db.Channel == channel_to_extract ) {
                for ( size_t i = 0; i < buf.size(); i++ ) {
                    if ( ( k >= startCount ) && ( k <= endCount ) ) {
                        *( pdWfData + counter ) = to_mV * buf[i];
                        counter++;
                    }
                    k++;
                }
            }
        }
        if ( k > endCount ) break;
    }
    return dataSize;
}

int CPlexMethods::plx_chan_gains( mxArray **ppmxGain )
{
    if ( m_SpikeChannelHeaders.size() == 0 ) {
        mexPrintf( "\nplx_chan_gains: file contains no spike channels.\n" );
        return 0;
    }

    *ppmxGain = mxCreateDoubleMatrix( m_SpikeChannelHeaders.size(), 1, mxREAL );
    double* pdGain = mxGetPr( *ppmxGain );

    int preAmpGain = 1000;
    if ( m_FileHeader.Version >= 105 ) {
        preAmpGain = m_FileHeader.SpikePreAmpGain;
    }

    for ( size_t j = 0; j < m_SpikeChannelHeaders.size(); j++ ) {
        int gain = m_SpikeChannelHeaders[j].Gain;
        if ( gain == 0 ) gain = 1;
        *( pdGain + j ) = gain * preAmpGain;
    }
    return (int)m_SpikeChannelHeaders.size();
}

int CPlexMethods::plx_chan_thresholds( mxArray **ppmxTresh )
{
    if ( m_SpikeChannelHeaders.size() == 0 ) {
        mexPrintf( "\nplx_chan_thresholds: file contains no spike channels.\n" );
        return 0;
    }
    *ppmxTresh = mxCreateDoubleMatrix( m_SpikeChannelHeaders.size(), 1, mxREAL );
    double* pdTresh = mxGetPr( *ppmxTresh );

    for ( size_t j = 0; j < m_SpikeChannelHeaders.size(); j++ ) {
        *( pdTresh + j ) = m_SpikeChannelHeaders[j].Threshold;
    }
    return (int)m_SpikeChannelHeaders.size();
}

int CPlexMethods::plx_chan_filters( mxArray **ppmxFilter )
{
    if ( m_SpikeChannelHeaders.size() == 0 ) {
        mexPrintf( "\nplx_chan_filters: file contains no spike channels.\n" );
        return 0;
    }

    *ppmxFilter = mxCreateDoubleMatrix( m_SpikeChannelHeaders.size(), 1, mxREAL );
    double* pdFilter = mxGetPr( *ppmxFilter );

    for ( size_t j = 0; j < m_SpikeChannelHeaders.size(); j++ ) {
        *( pdFilter + j ) = m_SpikeChannelHeaders[j].Filter;
    }
    return (int)m_SpikeChannelHeaders.size();
}

int CPlexMethods::plx_adchan_gains( mxArray **ppmxGain )
{
    if ( m_AnalogChannelHeaders.size() == 0 ) {
        mexPrintf( "\nplx_adchan_gain: file contains no A/D channels.\n" );
        return 0;
    }

    *ppmxGain = mxCreateDoubleMatrix( m_AnalogChannelHeaders.size(), 1, mxREAL );
    double* pdGain = mxGetPr( *ppmxGain );

    for ( size_t j = 0; j < m_AnalogChannelHeaders.size(); j++ ) {
        int gain = m_AnalogChannelHeaders[j].Gain * m_AnalogChannelHeaders[j].PreAmpGain;
        if ( gain == 0 ) gain = 1;
        *( pdGain + j ) = ( double )gain;
    }
    return (int)m_AnalogChannelHeaders.size();
}

int CPlexMethods::plx_adchan_freqs( mxArray **ppmxFreq )
{
    if ( m_AnalogChannelHeaders.size() == 0 ) {
        mexPrintf( "\nplx_adchan_freq: file contains no A/D channels.\n" );
        return 0;
    }

    *ppmxFreq = mxCreateDoubleMatrix( m_AnalogChannelHeaders.size(), 1, mxREAL );
    double* pdFreq = mxGetPr( *ppmxFreq );

    for ( size_t j = 0; j < m_AnalogChannelHeaders.size(); j++ ) {
        *( pdFreq + j ) = ( double )m_AnalogChannelHeaders[j].ADFreq;
    }
    return (int)m_AnalogChannelHeaders.size();

}

int CPlexMethods::plx_adchan_samplecounts( mxArray **ppmxSamples )
{
    if ( m_AnalogChannelHeaders.size() == 0 ) {
        mexPrintf( "\nplx_adchan_samplecounts: file contains no A/D channels.\n" );
        return 0;
    }

    if ( !m_bTallied ) {
        if ( !tallyPlxFile() ){
            return 0;
        }
    }

    *ppmxSamples = mxCreateDoubleMatrix( m_AnalogChannelHeaders.size(), 1, mxREAL );
    double* pdSamples = mxGetPr( *ppmxSamples );

    for ( size_t j = 0; j < m_AnalogChannelHeaders.size(); j++ ) {
        pdSamples[j] = ( double )GetAnalogChannelSampleCount( m_AnalogChannelHeaders[j].Channel );
    }
    return ( int )m_AnalogChannelHeaders.size();
}

void CPlexMethods::plx_information( int *Version, int* Freq, char** pszComment, int *Trodalness, int *NPW, int *PreTresh,
                                   int *SpikePeakV, int *SpikeADResBits, int *SlowPeakV,
                                   int *SlowADResBits, double *Duration, char* DateTime, size_t DateTimeLen )
{
    ( *Version ) = ( int )m_FileHeader.Version;
    ( *Freq ) = ( int )m_FileHeader.ADFrequency;
    ( *pszComment ) = m_FileHeader.Comment;
    ( *Trodalness ) = ( int )m_FileHeader.Trodalness;
    ( *NPW ) = m_FileHeader.NumPointsWave;
    ( *PreTresh ) = m_FileHeader.NumPointsPreThr;
    ( *SpikePeakV ) = m_FileHeader.SpikeMaxMagnitudeMV;
    ( *SpikeADResBits ) = m_FileHeader.BitsPerSpikeSample;
    ( *SlowPeakV ) = m_FileHeader.SlowMaxMagnitudeMV;
    ( *SlowADResBits ) = m_FileHeader.BitsPerSlowSample;
    ( *Duration )  = m_FileHeader.LastTimestamp / ( *Freq );
    sprintf( DateTime, "%2d/%2d/%4d %2d:%2d:%2d",
        m_FileHeader.Month, m_FileHeader.Day, m_FileHeader.Year,
        m_FileHeader.Hour, m_FileHeader.Minute, m_FileHeader.Second );
}

int CPlexMethods::plx_chan_names( mxArray **ppmxNames )
{
    if ( m_SpikeChannelHeaders.size() == 0 ) {
        mexPrintf( "\nplx_chan_names: file contains no spike channels.\n" );
        return 0;
    }

    std::vector<const char*> strings;
    strings.resize( m_SpikeChannelHeaders.size() );
    for ( size_t i = 0; i < m_SpikeChannelHeaders.size(); i++ ) {
        strings[i] = m_SpikeChannelHeaders[i].Name;
    }
    *ppmxNames = mxCreateCharMatrixFromStrings( m_SpikeChannelHeaders.size(), &strings[0] );
    return ( int )m_SpikeChannelHeaders.size();
}

int CPlexMethods::plx_adchan_names( mxArray **ppmxNames )
{
    if ( m_AnalogChannelHeaders.size() == 0 ) {
        mexPrintf( "\nplx_adchan_names: file contains no A/D channels.\n" );
        return 0;
    }

    std::vector<const char*> strings;
    strings.resize( m_AnalogChannelHeaders.size() );
    for ( size_t i = 0; i < m_AnalogChannelHeaders.size(); i++ ) {
        strings[i] = m_AnalogChannelHeaders[i].Name;
    }
    *ppmxNames = mxCreateCharMatrixFromStrings( m_AnalogChannelHeaders.size(), &strings[0] );
    return ( int )m_AnalogChannelHeaders.size();
}

int CPlexMethods::plx_event_names( mxArray **ppmxNames )
{
    if ( m_EventChannelHeaders.size() == 0 ) {
        mexPrintf( "\nplx_event_names: file contains no event channels.\n" );
        return 0;
    }

    std::vector<const char*> strings;
    strings.resize( m_EventChannelHeaders.size() );
    for ( size_t i = 0; i < m_EventChannelHeaders.size(); i++ ) {
        strings[i] = m_EventChannelHeaders[i].Name;
    }
    *ppmxNames = mxCreateCharMatrixFromStrings( m_EventChannelHeaders.size(), &strings[0] );
    return ( int )m_EventChannelHeaders.size();
}

// helpers
double  CPlexMethods::calculateRawToMilliVoltsCoeffForAnalogChannel( int analogChannelHeaderIndex ) const
{
    double  preamp_gain = m_FileHeader.Version > 101 ? m_AnalogChannelHeaders[analogChannelHeaderIndex].PreAmpGain : 1000.0;
    double  max_magnitude = m_FileHeader.Version > 102 ? m_FileHeader.SlowMaxMagnitudeMV : 5000.0;
    double  sample_size = m_FileHeader.Version > 102 ? 1 << ( m_FileHeader.BitsPerSlowSample - 1 ) : 2048.0;

    return  max_magnitude / ( sample_size * preamp_gain * m_AnalogChannelHeaders[analogChannelHeaderIndex].Gain );
}

double  CPlexMethods::calculateRawToMilliVoltsCoeffForSpikeChannel( int spikeChannelHeaderIndex ) const
{
    double  preamp_gain = m_FileHeader.Version > 104 ? m_FileHeader.SpikePreAmpGain : 1000.0;
    double  max_magnitude = m_FileHeader.Version > 102 ? m_FileHeader.SpikeMaxMagnitudeMV : 3000.0;
    double  sample_size = m_FileHeader.Version > 102 ? 1 << ( m_FileHeader.BitsPerSpikeSample - 1 ) : 2048.0;

    return  max_magnitude / ( sample_size * preamp_gain * m_SpikeChannelHeaders[spikeChannelHeaderIndex].Gain );
}

// VT interpreter
namespace
{
    enum VT_Type { BAD = 0, X1 = 1, Y1 = 2, X2 = 4, Y2 = 8, X3 = 16, Y3 = 32, CX1 = 64, CY1 = 128, CM = 256 };

    enum VT_Mode {
        UNKNOWN,
        CENTROID,               // 1 set of coordinates, no motion
        CENTROID_WITH_MOTION,   // 1 set of coordinates, with motion
        LED_1,                  // 1 set of coordinates
        LED_2,
        LED_3,
        LED_12,                 // 2 sets of coordinates
        LED_13,
        LED_23,
        LED_123,                // 3 sets of coordinates
    };

    static const unsigned   MAP_MODE_TO_SIZE[ LED_123 + 1 ] = { 0, 3, 4, 3, 3, 3, 5, 5, 5, 7 };

    struct VT_Data {
        static const unsigned short DATA_MASK = 0x03FF;
        static const unsigned short DIBT_MASK = 0x8000;
        static const unsigned short TYPE_MASK = 0x3C00;

        static const unsigned short CENTROID_X1 = 0x0000;
        static const unsigned short CENTROID_Y1 = 0x0400;
        static const unsigned short CENTROID_MOTION = 0x1000;
        static const unsigned short LED_X1 = 0x1400;
        static const unsigned short LED_Y1 = 0x1800;
        static const unsigned short LED_X2 = 0x1C00;
        static const unsigned short LED_Y2 = 0x2000;
        static const unsigned short LED_X3 = 0x2400;
        static const unsigned short LED_Y3 = 0x2800;

        VT_Data( double val ) {
            unsigned short  v = static_cast<unsigned short>( val );
            value = v & DATA_MASK;
            switch ( v & TYPE_MASK ) {
case CENTROID_X1:
    type = CX1;
    break;
case CENTROID_Y1:
    type = CY1;
    break;
case CENTROID_MOTION:
    type = CM;
    break;
case LED_X1:
    type = X1;
    break;
case LED_Y1:
    type = Y1;
    break;
case LED_X2:
    type = X2;
    break;
case LED_Y2:
    type = Y2;
    break;
case LED_X3:
    type = X3;
    break;
case LED_Y3:
    type = Y3;
    break;
default:
    type = BAD;
    break;
            }
        }

        VT_Type     type;
        unsigned    value;
    };

    struct VT_Acc {
        VT_Acc() : present( 0 ) {}

        double          ts;
        unsigned short  x1, y1, x2, y2, x3, y3, cx1, cy1, cm;
        unsigned        present;

        static const double MAX_DELAY;

        static const unsigned short LED_123_MASK = X1 + Y1 + X2 + Y2 + X3 + Y3;
        static const unsigned short LED_12_MASK = X1 + Y1 + X2 + Y2;
        static const unsigned short LED_13_MASK = X1 + Y1 + X3 + Y3;
        static const unsigned short LED_23_MASK = X2 + Y2 + X3 + Y3;
        static const unsigned short LED_1_MASK = X1 + Y1;
        static const unsigned short LED_2_MASK = X2 + Y2;
        static const unsigned short LED_3_MASK = X3 + Y3;
        static const unsigned short CENTROID_MASK = CX1 + CY1;
        static const unsigned short CENTROID_WITH_MOTION_MASK = CX1 + CY1 + CM;

        void    clear() { present = 0; }

        bool    accept( double timestamp, const VT_Data& data ) {
            if ( present && ( timestamp < ts || timestamp > ts + MAX_DELAY ) )   return  false;
            if ( data.type == BAD || ( present & data.type ) )   return  false;
            ts = timestamp;
            switch ( data.type ) {
case X1:
    x1 = data.value;
    break;
case Y1:
    y1 = data.value;
    break;
case X2:
    x2 = data.value;
    break;
case Y2:
    y2 = data.value;
    break;
case X3:
    x3 = data.value;
    break;
case Y3:
    y3 = data.value;
    break;
case CX1:
    cx1 = data.value;
    break;
case CY1:
    cy1 = data.value;
    break;
case CM:
    cm = data.value;
    break;
            }
            present |= data.type;
            return  true;
        }

        // helper
        static bool test( unsigned short value, unsigned short mask ) { return ( value & mask ) == mask; }

        VT_Mode mode() const {
            if ( test( present, LED_123_MASK ) ) return  LED_123;
            if ( test( present, LED_12_MASK ) )  return  LED_12;
            if ( test( present, LED_13_MASK ) )  return  LED_13;
            if ( test( present, LED_23_MASK ) )  return  LED_23;
            if ( test( present, CENTROID_WITH_MOTION_MASK ) )    return  CENTROID_WITH_MOTION;
            if ( test( present, CENTROID_MASK ) )                return  CENTROID;
            if ( test( present, LED_1_MASK ) )   return  LED_1;
            if ( test( present, LED_2_MASK ) )   return  LED_2;
            if ( test( present, LED_3_MASK ) )   return  LED_3;
            return  UNKNOWN;
        }
    };

    const double VT_Acc::MAX_DELAY = 1.0 / 105.0;
}

void CPlexMethods::vt_interpret(
                                const mxArray* ts,
                                const mxArray* sv,
                                mxArray** nCoords,
                                mxArray** nDim,
                                mxArray** vtMode,
                                mxArray** coords
                                )
{
    int nts = ( int )mxGetM( ts );
    const double* pts = mxGetPr( ts );
    const double* psv = mxGetPr( sv );

    // VT groups
    std::vector<VT_Acc> v;

    // mode counters
    unsigned    m[ LED_123 + 1 ];
    std::fill( m, m + LED_123 + 1, 0 );

    // build array of VT groups
    for ( int i = 0; i < nts; ++i ) {
        double  ts = pts[i];
        VT_Data val( psv[i] );
        if ( v.empty() || !v.back().accept( ts, val ) ) {
            if ( !v.empty() )    ++m[ v.back().mode() ];
            VT_Acc  acc;
            if ( acc.accept( ts, val ) ) v.push_back( acc );
        }
    }
    if ( !v.empty() )    ++m[ v.back().mode() ];

    // find the most populous mode
    VT_Mode     mode  = static_cast<VT_Mode>( std::max_element( m + 1, m + LED_123 + 1 ) - m );
    unsigned    count = m[mode];
    unsigned    size  = MAP_MODE_TO_SIZE[mode];

    // allocate output variables
    *nCoords = mxCreateDoubleMatrix( 1, 1, mxREAL );
    mxGetPr( *nCoords )[0] = count;
    *nDim = mxCreateDoubleMatrix( 1, 1, mxREAL );
    mxGetPr( *nDim )[0] = size;
    *vtMode = mxCreateDoubleMatrix( 1, 1, mxREAL );
    mxGetPr( *vtMode )[0] = mode;
    *coords = mxCreateDoubleMatrix( count, size, mxREAL );
    double* coordinates = mxGetPr( *coords );

    unsigned    j = 0;
    for ( std::vector<VT_Acc>::const_iterator i = v.begin(); i != v.end(); ++i ) {
        if ( i->mode() != mode ) continue;
        coordinates[ j + 0 * count ] = i->ts;
        switch ( mode ) {
case CENTROID_WITH_MOTION:
    coordinates[ j + 3 * count ] = i->cm;
case CENTROID:
    coordinates[ j + 1 * count ] = i->cx1;
    coordinates[ j + 2 * count ] = i->cy1;
    break;
case LED_123:
    coordinates[ j + 5 * count ] = i->x3;
    coordinates[ j + 6 * count ] = i->y3;
case LED_12:
    coordinates[ j + 3 * count ] = i->x2;
    coordinates[ j + 4 * count ] = i->y2;
case LED_1:
    coordinates[ j + 1 * count ] = i->x1;
    coordinates[ j + 2 * count ] = i->y1;
    break;
case LED_23:
    coordinates[ j + 3 * count ] = i->x3;
    coordinates[ j + 4 * count ] = i->y3;
case LED_2:
    coordinates[ j + 1 * count ] = i->x2;
    coordinates[ j + 2 * count ] = i->y2;
    break;
case LED_13:
    coordinates[ j + 1 * count ] = i->x1;
    coordinates[ j + 2 * count ] = i->y1;
    coordinates[ j + 3 * count ] = i->x3;
    coordinates[ j + 4 * count ] = i->y3;
    break;
case LED_3:
    coordinates[ j + 1 * count ] = i->x3;
    coordinates[ j + 2 * count ] = i->y3;
    break;
        }
        ++j;
    }
}

void CPlexMethods::ddt_read( const char* fileName, mxArray **ppmxDdt, bool mV, int& numberOfChannels, unsigned int& dataCount, int& ADfrequency )
{
    short buf[512];
    DigFileHeader dg;
    long long fileSize = 0;

    // note that we are using local file here since the file is not persisted between calls
    FILE* localFile = fopen( fileName, "rb" );
    if ( localFile == NULL ) {
        ReportError( "Cannot open ddt file", fileName );
        return;
    }
    fseek( localFile, 0, SEEK_END );
    fileSize = ftell( localFile );
    fseek( localFile, 0, SEEK_SET );
    if ( fread( &dg, sizeof( dg ), 1, localFile ) != 1 ) {
        mexErrMsgTxt( "Cannot read ddt file header." );
        fclose( localFile );
        return;
    }

    numberOfChannels = dg.NChannels;
    dataCount = 0;
    ADfrequency = ( int )dg.Freq;

    int dataOffset = dg.DataOffset;
    long long nbuf64 = ( fileSize - dataOffset ) / ( dg.NChannels * 2 );
    if ( nbuf64 >= UINT_MAX ) {
        mexErrMsgTxt( "ddt file contains more than 2^32 samples per channel\n" );
        fclose( localFile );
        return;
    }

    unsigned int nbuf = ( unsigned int )nbuf64;
    dataCount = nbuf;

    *ppmxDdt = mxCreateDoubleMatrix( dg.NChannels, nbuf, mxREAL );
    double* pdDdData = mxGetPr( *ppmxDdt );

    if ( mV ) {
        double to_mV[64];
        double mx = 1 << ( dg.Version >= 101 ? dg.BitsPerSample : 12 );
        if ( dg.Version < 102 ) {
            to_mV[0] = 10.0 / dg.Gain / mx;
            std::fill( to_mV + 1, to_mV + dg.NChannels, to_mV[0] );
        } else {
            int j = 0;
            double my = ( dg.Version >= 103 ? dg.MaxMagnitudeMV : 5000 ) * 2.0 / dg.Gain / mx;
            for ( int i = 0; i < 64; ++i ) {
                if ( dg.ChannelGain[i] != 255 ) {
                    to_mV[j++] = my / dg.ChannelGain[i];
                    if ( j >= dg.NChannels ) break;
                }
            }
        }
        while ( feof( localFile ) == 0 && fread( buf, sizeof( short ), dg.NChannels, localFile ) == dg.NChannels ) {
            pdDdData = std::transform( buf, buf + dg.NChannels, to_mV, pdDdData, std::multiplies<double>() );
        }
    } else {

        while ( feof( localFile ) == 0 && fread( buf, sizeof( short ), dg.NChannels, localFile ) == dg.NChannels ) {
            pdDdData = std::copy( buf, buf + dg.NChannels, pdDdData );
        }
    }

    fclose( localFile );
}

int CPlexMethods::ddt_write( const char* filename, int nch, int npoints, double freq, const mxArray** d )
{
    // note that we are using local file here
    FILE* localFile = fopen( filename, "wb" );
    // open the file
    if ( localFile == NULL ) {
        ReportError( "Cannot open ddt file", filename );
        return	0;
    }

    // create a header
    DigFileHeader	header;
    memset( &header, 0, sizeof( header ) );
    header.Version = 103;
    header.DataOffset = sizeof( header );
    header.Freq = freq;
    header.NChannels = nch;
    header.BitsPerSample = 16;
    header.MaxMagnitudeMV = 5000;

    strncpy( header.Comment, "Made with Plexon's Matlab Offline SDK", sizeof( header.Comment ) );
    memset( header.ChannelGain, 1, nch );
    memset( header.ChannelGain + nch, 255, sizeof( header.ChannelGain ) - nch );

    time_t	now = time( NULL );
    tm*	t;
    t = localtime( &now );
    header.Year   = t->tm_year + 1900;
    header.Month  = t->tm_mon  + 1;
    header.Day    = t->tm_mday;
    header.Hour   = t->tm_hour;
    header.Minute = t->tm_min;
    header.Second = t->tm_sec;

    double from_mV = 0x10000 / 10000.0;

    const double* p = mxGetPr( *d );

    // calculate scale
    double scale = 1;
    {
        // find minmax
        const double* x = p + nch * npoints - 1;
        double minimum = *x;
        double maximum = *x;
        if ( x != p ) {
            for ( --x; x != p; --x ) {
                if ( maximum < *x ) {
                    maximum = *x;
                } else if ( *x < minimum ) {
                    minimum = *x;
                }
            }
        }
        minimum = ( minimum < 0 ? -32768 / from_mV / minimum : 0 );
        maximum = ( maximum > 0 ? 32767 / from_mV / maximum : 0 );
        if ( minimum ) {
            scale = maximum ? min( minimum, maximum ) : minimum;
        } else {
            if ( maximum ) scale = maximum;
        }
    }

    // round off the scale and calculate the gain
    int	gain = int( scale );
    for ( int order = 1000000; order > 1; order /= 10 ) {
        if ( gain > order ) {
            gain -= gain % order;
            from_mV *= gain;
            break;
        }
    }

    header.Gain = gain;

    // write the header
    if ( !fwrite( &header, sizeof( header ), 1, localFile ) ) {
        fclose( localFile );
        ReportError( "Cannot write ddt file header", filename );
        return	0;
    }

    // write data
    short	buf[64];
    TransHelper<short>	mV( from_mV );
    for ( int i = 0; i < npoints; ++i ) {
        std::transform( p, p + nch, buf, mV );
        p += nch;
        if ( !fwrite( buf, sizeof( buf[0] ), nch, localFile ) ) {
            fclose( localFile );
            ReportError( "Cannot write ddt file data", filename );
            return	0;
        }
    }

    fclose( localFile );
    return	1;
}

///////////// Reader class /////////////////

Reader::Reader()
:m_pFile(NULL)
{
}

Reader::~Reader()
{
    Close();
}

bool Reader::Open( const char* fileName )
{
    m_pFile = fopenEx( fileName, "rb" );
    return ( m_pFile != NULL );
}

unsigned long long Reader::GetFileSize()
{
    if ( m_pFile == NULL ) {
        return 0;
    }
    if( fseekEx( m_pFile, 0, SEEK_END ) != 0 ) {
        return 0;
    }
    unsigned long long fileSize = ftellEx( m_pFile );
    fseekEx( m_pFile, 0, SEEK_SET );
    return fileSize;
}

bool Reader::Seek( unsigned long long seekPositionFromStartOfFile )
{
    if ( m_pFile == NULL ) {
        return false;
    }
    if ( fseekEx( m_pFile, seekPositionFromStartOfFile, SEEK_SET ) == 0 ){
        return true;
    }
    return false;
}

unsigned int Reader::Read( void* buffer, unsigned int numberOfBytesToRead )
{
    if ( m_pFile == NULL ) {
        return 0;
    }
    if ( fread( buffer, numberOfBytesToRead, 1, m_pFile ) == 1 ) {
        return numberOfBytesToRead;
    } else {
        return 0;
    }
}

bool Reader::IsAtEndOfFile()
{
    if ( m_pFile == NULL ) {
        return true;
    }
    return ( feof( m_pFile ) != 0 );
}

void Reader::Close()
{
    if ( m_pFile != NULL ) {
        fclose(m_pFile);
        m_pFile = NULL;
    }
}
