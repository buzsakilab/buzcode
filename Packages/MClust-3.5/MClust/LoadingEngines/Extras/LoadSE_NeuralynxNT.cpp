/*---------------------------------
* LoadST0_nt
* MEX file
*
* input:
*    fn = file name string
*    records_to_get = an array that is either a range of values
*    record_units = 
*         1: timestamp list.(a vector of timestamps to load (uses binsearch to find them in file))
*         2: record number list (a vector of records to load)
*	  3: range of timestamps (a vector with 2 elements: a start and an end timestamp)
*	  4: range of records (a vector with 2 elements: a start and an end record number)
*	  5: return the count of spikes (records_to_get should be [] in this case)
*	 if only fn is passed in then the entire file is opened.
*    if only fn is passed AND only t is provided as the output, then all 
*       of the timestamps for the entire file are returned.
*
* output:
*    [t, wv]
*    t = n x 1: timestamps of each spike in file
*    wv = n x 4 x 32 waveforms
*
* version 5.1
*
* Reads Neuralynx tetrode (TT) files, based on the Windows-family processors.
* by checking if a header exists (for sun TTfiles) or not (for NT-TTfiles)
* Checks for standard Neuralynx header if present (in Cheeath versions >= 1.3)
* and automatically skips header.
*
* ADR May 2002
* modified from LoadTT0_sun_and_nt
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
*--------------------------------*/

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

#ifdef __GNUC__
#define __int64 long long 
#endif
  
#ifdef __GNUC__
  const long long  TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFFLL;  //( = 2^52-1; the largest integer fitting
#else
  const __int64 TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFF;  //( = 2^52-1; the largest integer fitting
#endif

// needed for structure to align correctly.
#pragma pack(1) 

#define BY_TIMESTAMP 1
#define BY_RECORD 2
#define BY_TIMESTAMP_RANGE 3
#define BY_RECORD_RANGE 4
#define COUNTSPIKES 5

//___________________________________________________________________________________
int SkipCheetahNTHeader(FILE *fp)
{
  // Cheetah NT versions 1.3 or greater have a standard Neuralynx header of ascii strings of total
  // length of 16384 bytes in the beginning of the binary file. 
  // Check if Cheetah header present (in versions > 1.3) and skip header
  //
  // returns 1 if new Neuralynx header is present and was skipped
  //         0 if NO  header is present
  
  char headerflag[8];   // if first 8 bytes of AnyCheetahfile.dat contain "########" there is a header 
  // of length 16,384 bytes to skip (from the very beginning including the eight '#') 
  const int NHEADERBYTES = 16384;	
  
  fread(headerflag, sizeof(char), 8, fp);  
  //mexPrintf("Headerflag =  %8s \n", headerflag);
  fseek(fp,0,0);       // reset filepointer to the beginning of file
  if (strncmp(headerflag,"########",8) == 0){
    fseek(fp,NHEADERBYTES,0);  // set filepointer after byte NHEADERBYTES
    mexPrintf("NT-Header skipped (%d bytes)\n",NHEADERBYTES);
    return 1;
  } 
  return 0;
}

//___________________________________________________________________________________
int SkipHeader(FILE *fp)
{
  /* returns 0 if header present and skipped  (success) 
  **         1 if NO header present in file (indicates a new NT_cheetah TT file) */
  long curpos = ftell(fp);
  char headerline[81];
  
  fgets(headerline, 80, fp);
  if (strncmp(headerline, "%%BEGINHEADER", 13) == 0){
    while (strncmp(headerline, "%%ENDHEADER",11) != 0)
      fgets(headerline, 80, fp);
    return 0;
  } else {
    fseek(fp, curpos, SEEK_SET);
    return 1;
  }
}

//___________________________________________________________________________________
int bigendianMachine(void)
{
  /* returns 0 if is a littleendian machine, else returns nonzero */
  /* key is that it looks to see if short's second byte is the low order bits */
  short fullLoByte = 0xFF;
  char *byteOrder = (char *) &fullLoByte;
  return (byteOrder[1]);	
}


//___________________________________________________________________________________
inline short swapbytes(short ii)
  // swap byte order of a short: (0,1) -> (1,0)
{
  union {
    short s;
    char c[2];
  } tmp0,tmp;
  
  tmp.s = ii;
  tmp0.c[0] = tmp.c[1];					
  tmp0.c[1] = tmp.c[0];
  
  return tmp0.s;
}

//___________________________________________________________________________________
inline unsigned long swapbytes(unsigned long ii)
  // swap byte order of an unsigned long: (0,1,2,3) -> (3,2,1,0)
{
  union {
    unsigned long ul;
    char c[4];
  } tmp0,tmp;
  
  tmp.ul = ii;
  tmp0.c[0] = tmp.c[3];					
  tmp0.c[1] = tmp.c[2];
  tmp0.c[2] = tmp.c[1];					
  tmp0.c[3] = tmp.c[0];
  
  return tmp0.ul;
}

//___________________________________________________________________________________
inline long swapbytes(long ii)
  // swap byte order of a long: (0,1,2,3) -> (3,2,1,0)
{
  union {
    long l;
    char c[4];
  } tmp0,tmp;
  
  tmp.l = ii;
  tmp0.c[0] = tmp.c[3];					
  tmp0.c[1] = tmp.c[2];
  tmp0.c[2] = tmp.c[1];					
  tmp0.c[3] = tmp.c[0];
  
  return tmp0.l;
}

//___________________________________________________________________________________
inline __int64 swapbytes(__int64 ii)
  // swap byte order of a long long: (0,1,2,3,4,5,6,7) -> (7,6,5,4,3,2,1,0)
{
  union {
    __int64 ll;
    char c[8];
  } tmp0,tmp;
  
  tmp.ll = ii;
  tmp0.c[0] = tmp.c[7];					
  tmp0.c[1] = tmp.c[6];
  tmp0.c[2] = tmp.c[5];					
  tmp0.c[3] = tmp.c[4];
  tmp0.c[4] = tmp.c[3];					
  tmp0.c[5] = tmp.c[2];
  tmp0.c[6] = tmp.c[1];					
  tmp0.c[7] = tmp.c[0];
  
  return tmp0.ll;
}

//___________________________________________________________________________________
int binsearch(int n, double* data, double* key)
{
  // n = number of records in data.
  // data is an array of records to search
  // key is the value to find.
  // 
  // binsearch returns the index of the closest item in data
  //
  // from ADR adapted by cowen. 
  // fixed to return always the lower integreal timestamp (floor(*data)) if key is not integral (PL)
  // changed end = n-1 to end = n to allow retrieval of last record in file (ncst)
  
  int start = 0;
  int end = 0; 
  int mid = 0;
  double tmp;
  
  // binary search 
  start = 0;
  end = n;                    
  while (start < (end-1))
    {
      mid = (int) floor((double)(start + end)/2);
      tmp = floor(data[mid]);
      if ((*key) == tmp)
	start = end = mid;
      if ((*key) < tmp)
	end = mid;
      if ((*key) > tmp)
	start = mid;
    }
  
  return start;
}

//___________________________________________________________________________________
//===================================================================================
//-----------------------------------------------------------------------------------
// Loading Functions

const int recSize = 112;
const int recSizeAfterTimestamp = recSize - sizeof(__int64);

//___________________________________________________________________________________
void GetOneRecord(FILE *fp, double &T, double wv[128])
  // reads one record and goes to the next record
{
  int bigendianFlag = bigendianMachine();

  // NT TT record
  __int64 qwTimeStamp, qwTimeStamp0;
  long dwParams[10];
  short snData[64];
  
  fread(&qwTimeStamp0,  sizeof(char),  8, fp);
  fread(dwParams,      sizeof(char),  40, fp);
  fread(snData  ,      sizeof(char),  64, fp);
			
  if(bigendianFlag){				// convert from NT(little endian) to Sun (big endian)
    qwTimeStamp = swapbytes(qwTimeStamp0);
    if(qwTimeStamp > TIMESTAMP_MAX){
      mexPrintf(" ERROR: timestamp is too large to fit in a double!\n");
      mexPrintf(" Converted timestamps MAY or MAY NOT be valid - proceed with care! \n");
    }
    T = (double) qwTimeStamp/100.0;
    for (int j = 0; j<1; j++)
      for(int k = 0; k<32; k++)
	wv[j + 4*k] = (double) swapbytes(snData[k]);  
    for (int j = 0; j<4; j++)
      for(int k = 0; k<32; k++)
	wv[j + 4*k] = 0;  
    
  } else {
				// don't convert, just copy
    qwTimeStamp = qwTimeStamp0;
    if(qwTimeStamp > TIMESTAMP_MAX){
      mexPrintf(" ERROR: timestamp is too large to fit in a double!\n");
      mexPrintf(" Converted timestamps MAY or MAY NOT be valid - proceed with care! \n");
    }
    T = (double) qwTimeStamp/100.0;
    for (int j = 0; j<1; j++)
      for(int k = 0; k<32; k++)
	wv[j + 4*k] = (double) snData[k];  
    for (int j =0; j<4; j++)
      for(int k = 0; k<32; k++)
	wv[j + 4*k] = 0;  
  }
  
}

//___________________________________________________________________________________
int GetNumberOfSpikes(char* fn){
  
  // open file 
  FILE* fp = fopen(fn, "rb");
  if (!fp)
    mexErrMsgTxt("Could not open file.");
  
  //skip header and determine file record size
  int new_NT_format = SkipHeader(fp);  
  if (!new_NT_format)
    mexErrMsgTxt("Old SunOS formats are not supported by this loading engine.");
  SkipCheetahNTHeader(fp);     // Skip standard Neuralynx header if present (Cheetah versions >= 1.3)

  // get filesize
  int postHeaderPos = ftell(fp);     // beginnig of file after header (if any)
  fseek(fp,0,2);                     // goto end of file
  int nbytes = ftell(fp) - postHeaderPos;
  
  int nSpikes = nbytes/recSize; // no need to skip last record for NT_cheetah files
  mexPrintf("Reading file %s:\nRecordSize = %d,  %d spikes, %d bytes.\n",
	    fn, recSize, nSpikes, nbytes);
  
  // cleanup
  fclose(fp);	
  
  return nSpikes;
}


//___________________________________________________________________________________
void ReadTT(char* fn, int nSpikes, double *t, double *wv){
  
  double t0; 
  double wv0[128];
  
  // open file 
  FILE *fp = fopen(fn, "rb");
  if (!fp)
    mexErrMsgTxt("ERROR: Could not open file.");
  
  // skip header 
  int new_NT_format = SkipHeader(fp);  // flag for new NT_format TT files (0=old SUN, 1=new NT) 
  if (!new_NT_format) 
    mexErrMsgTxt("Old sun format not supported by this loading engine.");
  SkipCheetahNTHeader(fp);     // Skip standard Neuralynx header if present (Cheetah versions >= 1.3)
  long postHeaderPos = ftell(fp);
  
  // read records and convert all to double
  fseek(fp, postHeaderPos, SEEK_SET);
  for (int i = 0; i < nSpikes; i++){
    GetOneRecord(fp,t0,wv0);
    t[i] = t0;
    for (int j = 0; j<4; j++)
      for(int k = 0; k<32; k++)
	wv[i + nSpikes*j + nSpikes*4*k] = (double) wv0[j + 4*k];  
  }
    
  fclose(fp);	
}

//_________________________________________________________________________________________________
void ReadTTByRecord(char* fn, double *records_to_get, int n_records_to_get, double *t, double *wv)
  // Open the file and fseek to just those record number passed in 
  // in the array: records_to_get. The last record of records to get 
  // indicates the end of records. It's code is END_OF_RECORDS.
{
  int i = 0;
  double t0;
  double wv0[128];

  // open file 
  FILE *fp = fopen(fn, "rb");
  if (!fp)
    mexErrMsgTxt("ERROR: Could not open file.");
  
  // skip header 
  int new_NT_format = SkipHeader(fp);  // flag for new NT_format TT files (0=old SUN, 1=new NT) 
  if (!new_NT_format) 
    mexErrMsgTxt("Old sun format not supported by this loading engine.");
  SkipCheetahNTHeader(fp);     // Skip standard Neuralynx header if present (Cheetah versions >= 1.3)
  long postHeaderPos = ftell(fp);
	
  // read records and convert all to double
  while(i < n_records_to_get)
    {
      // Go directly to the record in question. Do not pass go. NO $200.
      fseek(fp, postHeaderPos+sizeof(char)*(recSize)*((long)records_to_get[i] - 1), SEEK_SET);
      GetOneRecord(fp,t0,wv0);
      t[i] = t0;
      for (int j = 0; j<4; j++)
	for(int k = 0; k<32; k++)
	  wv[i + n_records_to_get*j + n_records_to_get*4*k] = (double) wv0[j + 4*k];	
      i++;
    }
  fclose(fp);	
}	

//___________________________________________________________________________________
void ReadTT_timestamps(char* fn, int nSpikes, double *t)
{
  double t0;
  double wv0[128];
  
  // open file 
  FILE *fp = fopen(fn, "rb");
  if (!fp)
    mexErrMsgTxt("ERROR: Could not open file.");
  
  // skip header 
  int new_NT_format = SkipHeader(fp);  // flag for new NT_format TT files (0=old SUN, 1=new NT) 
  if (!new_NT_format) 
    mexErrMsgTxt("Old sun format not supported by this loading engine.");
  SkipCheetahNTHeader(fp);     // Skip standard Neuralynx header if present (Cheetah versions >= 1.3)
  long postHeaderPos = ftell(fp);
	
  // read records and convert all to double
  fseek(fp, postHeaderPos, SEEK_SET);
  for (int i = 0; i < nSpikes; i++)
    {
      GetOneRecord(fp,t0,wv0);
      t[i] = t0;
    }
  fclose(fp);	
}	

//___________________________________________________________________________________
void mexFunction(int nOUT, mxArray *pOUT[],
				 int nINP, const mxArray *pINP[])
{
	
	int i;
	double *records_to_get, *range_to_get,*t,*wv, *all_timestamps;
	int n_records_to_get = 0;
	int record_units = 0;
	int nSpikesInFile = 0;
	int length_records_to_get = 0;	
	int start_idx=0;
	int end_idx = 0;
	int idx = 0;	
	/* check number of arguments: expects 1 input */
	if (nINP != 3 && nINP != 1)
				mexErrMsgTxt("Call with fn or fn and array of recs to load(vector), and 1(timestamp) or 2(record number	) as inputs.");
	if (nOUT > 2)
				mexErrMsgTxt("Requires two outputs (t, wv).");

	/* read inputs */
	int fnlen = (mxGetM(pINP[0]) * mxGetN(pINP[0])) + 1;
	char *fn = (char *) mxCalloc(fnlen, sizeof(char)); 
	if (!fn)
		mexErrMsgTxt("Not enough heap space to hold converted string.");
	int errorstatus = mxGetString(pINP[0], fn,fnlen);    
	if (errorstatus)
		mexErrMsgTxt("Could not convert string data.");


	nSpikesInFile = GetNumberOfSpikes(fn);
	////////////////////////////////////////////////////////
	// If only one input is passed, assume the whole file is
	// to be loaded. If only one output is present, assume it is
	// only timestamps.
	////////////////////////////////////////////////////////

	if(nINP == 1 && nOUT == 1)
	{
		/* Return the timestamps of all of the records */
			// create outputs 
		pOUT[0] = mxCreateDoubleMatrix(nSpikesInFile, 1, mxREAL);
		t = mxGetPr(pOUT[0]);

		// load tt file fn into t and wv arrays 
		ReadTT_timestamps(fn,nSpikesInFile,t);

		// cleanup
		mxFree(fn);
		return;

	}else if(nINP == 1 && nOUT == 2)
	{
		////////////////////////////////////////////////////////
		// create outputs 
		////////////////////////////////////////////////////////
		mexPrintf("Getting %i records.\n",nSpikesInFile);
		pOUT[0] = mxCreateDoubleMatrix(nSpikesInFile , 1, mxREAL);
		t = mxGetPr(pOUT[0]);
		int wvDims[] = {nSpikesInFile , 4, 32};    // wv = nSpikes x 4 x 32   FORTRAN (column major) array
		pOUT[1] = mxCreateNumericArray(3, wvDims, mxDOUBLE_CLASS, mxREAL);
		wv = mxGetPr(pOUT[1]);

		////////////////////////////////////////////////////////
		// load tt file fn into t and wv arrays 
		////////////////////////////////////////////////////////

		ReadTT(fn,nSpikesInFile ,t,wv);

		////////////////////////////////////////////////////////
		// cleanup
		////////////////////////////////////////////////////////

		mxFree(fn);
		return;
	}
	////////////////////////////////////////////////////////
	// unpack inputs 
	////////////////////////////////////////////////////////
	length_records_to_get = mxGetM(pINP[1]) * mxGetN(pINP[1]);
	records_to_get = mxGetPr(pINP[1]);
	record_units = (int) mxGetScalar(pINP[2]);

	switch(record_units)
	{
		case BY_TIMESTAMP:
			////////////////////////////////////////////////////////
			// Convert the timestamps into record numbers. This will
			// make loading these records easy.
			////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////
			// Create a very large array of all of the timestamps
			////////////////////////////////////////////////////////
			n_records_to_get = length_records_to_get;
			all_timestamps = (double *)calloc( nSpikesInFile, sizeof( double ) );
			if (all_timestamps == NULL)
				mexErrMsgTxt("NOT ENOUGH MEMORY");
			range_to_get = (double *) calloc( n_records_to_get, sizeof( double ) );
			if (range_to_get == NULL)
				mexErrMsgTxt("NOT ENOUGH MEMORY");

			ReadTT_timestamps(fn,nSpikesInFile,all_timestamps);

			for (i = 0;i<n_records_to_get;i++)
			{
				idx = binsearch(nSpikesInFile, all_timestamps, &records_to_get[i]);
				range_to_get[i] = idx + 1; // Add one since records are assumed
											 // to be from 1 to end.
			}

			free(all_timestamps);
			break;

		case BY_RECORD:
			////////////////////////////////////////////////////////
			// Get records asked for. First, subract one since matlab
			// users typically think records as being ordered from 
			// 1 to ....
			////////////////////////////////////////////////////////

			n_records_to_get = length_records_to_get;
			range_to_get = records_to_get;
			break;

		case BY_TIMESTAMP_RANGE:

			////////////////////////////////////////////////////////
			// Error checking
			////////////////////////////////////////////////////////

			if(length_records_to_get != 2)
			{
				mexErrMsgTxt("Must pass in two arguements for parameter 2.");
				return;
			}

			////////////////////////////////////////////////////////
			// Create a very large array of all of the timestamps
			////////////////////////////////////////////////////////
			all_timestamps = (double *)calloc( nSpikesInFile, sizeof( double ) );
			if (all_timestamps == NULL)
				mexErrMsgTxt("NOT ENOUGH MEMORY");

			ReadTT_timestamps(fn,nSpikesInFile,all_timestamps);
			////////////////////////////////////////////////////////
			// Find the index in all_timestamps of the start record.
			////////////////////////////////////////////////////////

			start_idx = binsearch(nSpikesInFile, all_timestamps, &records_to_get[0]) + 1; // Add one since records are assumed
																					  	  // to be from 1 to end.

			////////////////////////////////////////////////////////
			// Find the index in all_timestamps of the end record.
			////////////////////////////////////////////////////////
			end_idx = binsearch(nSpikesInFile, all_timestamps, &records_to_get[1]) + 1;  // Add one since records are assumed
																				  	 // to be from 1 to end.


			free(all_timestamps);
			
			n_records_to_get = end_idx - start_idx + 1;
			////////////////////////////////////////////////////////
			// Allocate the space
			// Using ints would be much more efficient, but it would make it necessary
			// to rewrite the ReadTT program. I don't want to do that.
			////////////////////////////////////////////////////////

			range_to_get = (double *) calloc( n_records_to_get, sizeof( double ) );
			if (range_to_get == NULL)
			{
				mexErrMsgTxt("NOT ENOUGH MEMORY");
				return;
			}

			for (i = 0; i<n_records_to_get;i++)
				range_to_get[i] = start_idx + i;

			break;

		case BY_RECORD_RANGE:
			////////////////////////////////////////////////////////
			// This is easy. First check to make sure the user passed in
			// the proper number of arguements.
			////////////////////////////////////////////////////////
			if(length_records_to_get != 2)
			{
				mexErrMsgTxt("Must pass in a start and end record number for argument 2.");
				return;
			}
			start_idx = (int) records_to_get[0] ;
			end_idx   = (int) records_to_get[1] ;

			////////////////////////////////////////////////////////
			n_records_to_get = end_idx - start_idx + 1;
			////////////////////////////////////////////////////////
			// Allocate the space
			////////////////////////////////////////////////////////
			range_to_get = (double  *)calloc( n_records_to_get, sizeof( double  ) );
			if (range_to_get == NULL)
			{
				mexErrMsgTxt("NOT ENOUGH MEMORY");
				return;
			}
			for (i = 0; i<n_records_to_get;i++)
			{
				range_to_get[i] = start_idx + i;
			}
			break;
			
		case COUNTSPIKES: // ADR/NCST
			printf("Counting spikes.\n");
			pOUT[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			t = mxGetPr(pOUT[0]);
			t[0] = nSpikesInFile;		
			break;

		default:
			mexErrMsgTxt("Incorrect parameter 3.");
	}
	if (!(record_units==COUNTSPIKES))
	{
		////////////////////////////////////////////////////////
		// Allocate the space
		////////////////////////////////////////////////////////
		printf("Getting %i records.\n",n_records_to_get);
		pOUT[0] = mxCreateDoubleMatrix(n_records_to_get, 1, mxREAL);
		t = mxGetPr(pOUT[0]);
		int wvDims[] = {n_records_to_get, 4, 32};    // wv = nSpikes x 4 x 32   FORTRAN (column major) array
		pOUT[1] = mxCreateNumericArray(3, wvDims, mxDOUBLE_CLASS, mxREAL);
		wv = mxGetPr(pOUT[1]);
		////////////////////////////////////////////////////////
		// Get the data. Record_units will whether to get by 
		// timestamp or by record number.
		////////////////////////////////////////////////////////
		ReadTTByRecord(fn,range_to_get,n_records_to_get,t,wv);

		// Free this variable only if it was not a matlab variable.
		if (record_units == BY_TIMESTAMP_RANGE || record_units == BY_RECORD_RANGE )
			free(range_to_get);
	}
	// cleanup
	mxFree(fn);

}

