/*---------------------------------
* Load TT
* MEX file
*
* PL 1999 (based on ADR 1998)
*
* input:
*    fn = file name string
* 
* output:
*    [t, wv]
*    t = n x 1: timestamps of each spike in file
*    wv = n x 4 x 32 waveforms
*
* version 5.0
*
* Reads both sun TTfiles and NT-TTfiles and distinguishes them
* by checking if a header exists (for sun TTfiles) or not (for NT-TTfiles)
*
--------------------------------*/

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <matrix.h>

#ifdef __GNUC__
#define __int64 int64_t 
#endif

//___________________________________________________________________________________
int32_t SkipHeader(FILE *fp)
{
	/* returns 0 if header present and skipped  (success) 
	**         1 if NO header present in file (indicates a new NT_cheetah TT file) */
	int32_t curpos = ftell(fp);
	char headerline[81];

	fgets(headerline, 80, fp);
	if (strncmp(headerline, "%%BEGINHEADER", 13) == 0){
		while (strncmp(headerline, "%%ENDHEADER",11) != 0)
			fgets(headerline, 80, fp);
		return 0;
	} else {
		fseek(fp, curpos+16384, SEEK_SET);
		return 1;
	}
}

//___________________________________________________________________________________
int32_t bigendianMachine(void)
{
	/* returns 0 if is a littleendian machine, else returns nonzero */
	/* key is that it looks to see if int16_t's second byte is the low order bits */
	int16_t fullLoByte = 0xFF;
	char *byteOrder = (char *) &fullLoByte;
	return (byteOrder[1]);	
}


//___________________________________________________________________________________
inline int16_t swapbytes(int16_t ii)
// swap byte order of a int16_t: (0,1) -> (1,0)
{
	union {
		int16_t s;
		char c[2];
	} tmp0,tmp;
	
	tmp.s = ii;
	tmp0.c[0] = tmp.c[1];					
	tmp0.c[1] = tmp.c[0];

	return tmp0.s;
}

//___________________________________________________________________________________
inline uint32_t swapbytes(uint32_t ii)
// swap byte order of an uint32_t: (0,1,2,3) -> (3,2,1,0)
{
	union {
		uint32_t ul;
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
inline int32_t swapbytes(int32_t ii)
// swap byte order of a int32_t: (0,1,2,3) -> (3,2,1,0)
{
	union {
		int32_t l;
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
inline int64_t swapbytes(int64_t ii)
// swap byte order of a int64_t: (0,1,2,3,4,5,6,7) -> (7,6,5,4,3,2,1,0)
{
	union {
		int64_t ll;
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
int GetNumberOfSpikes(char* fn, int32_t nProbes){
	
	// open file 
	FILE* fp = fopen(fn, "rb");
	if (!fp)
		mexErrMsgTxt("Could not open file.");

	//skip header and determine file record size
	int32_t new_NT_format = SkipHeader(fp);
	int32_t recSize = 260;
	if (new_NT_format) recSize = 64 * nProbes + 48;; 
        
	// get filesize
	int32_t postHeaderPos = ftell(fp);     // beginnig of file after header (if any)
	fseek(fp,0,2);                     // goto end of file
	int32_t nbytes = ftell(fp) - postHeaderPos;

	int32_t nSpikes = nbytes/recSize - 1;    /* skip last record since it may be incomplete */
	mexPrintf("Reading file %s:\nRecordSize = %d,  %d spikes, %d bytes.\n",
		fn, recSize, nSpikes, nbytes);

	// cleanup
	fclose(fp);	

	return nSpikes;
}


//___________________________________________________________________________________
void ReadTT(char* fn, int32_t nSpikes, double *t, double *wv, int32_t nProbes){
	
#ifdef __GNUC__
	const int64_t  TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFFLL;  //( = 2^52-1; the largest integer fitting
#else
	const int64_t TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFF;  //( = 2^52-1; the largest integer fitting
#endif
	//   in a double IEEE mantissa (=7 bytes)
	int32_t bigendianFlag = bigendianMachine();
	
	// NT TT record
	int64_t qwTimeStamp, qwTimeStamp0;
	int32_t dwParams[10];
	int16_t snData[32*nProbes];
	
	// sun TT record
	int32_t tmpT;
	int16_t tmpWV[32*nProbes];
	
	// open file 
	FILE *fp = fopen(fn, "rb");
	if (!fp)
		mexErrMsgTxt("ERROR: Could not open file.");
	
	// skip header 
	int32_t new_NT_format = SkipHeader(fp);  // flag for new NT_format TT files (0=old SUN, 1=new NT) 
	int32_t postHeaderPos = ftell(fp);
	
	if (new_NT_format){ 
		
		// read records and convert all to double
		fseek(fp, postHeaderPos, SEEK_SET);
		for (int32_t i = 0; i < nSpikes; i++){
			
			fread(&qwTimeStamp0,  sizeof(char),   8, fp);
			fread(dwParams,      sizeof(char),  40, fp);
			fread(snData  ,      sizeof(char), 64*nProbes, fp);
			
			if(bigendianFlag){
				// convert from NT(little endian) to Sun (big endian)
				qwTimeStamp = swapbytes(qwTimeStamp0);
				if(qwTimeStamp > TIMESTAMP_MAX){
					mexPrintf(" ERROR: timestamp %d in file %s is too large to fit in a double!\n",i,fn);
				   mexPrintf(" Converted timestamps MAY or MAY NOT be valid - proceed with care! \n");
				}
				t[i] = (double) qwTimeStamp/100.0;
				for (int32_t j = 0; j<nProbes; j++)
					for(int32_t k = 0; k<32; k++)
						wv[i + nSpikes*j + nSpikes*nProbes*k] = (double) swapbytes(snData[j + nProbes*k]);  
				
			} else {
				// don't convert, just copy
				qwTimeStamp = qwTimeStamp0;
				if(qwTimeStamp > TIMESTAMP_MAX){
					mexPrintf(" ERROR: timestamp %d in file %s is too large to fit in a double!\n",i,fn);
				   mexPrintf(" Converted timestamps MAY or MAY NOT be valid - proceed with care! \n");
				}
				t[i] = (double) qwTimeStamp/100.0;
				for (int32_t j = 0; j<nProbes; j++)
					for(int32_t k = 0; k<32; k++)
						wv[i + nSpikes*j + nSpikes*nProbes*k] = (double) snData[j + nProbes*k];  
			}
			
		}
		
	} else {    /*  OLD sun cheetah format */
		
		// read t, wv 
		fseek(fp, postHeaderPos, SEEK_SET);
		for (int32_t i = 0; i < nSpikes; i++){
			
			fread(&tmpT,  sizeof(char),   4, fp);
			fread(tmpWV, sizeof(char), 256, fp);
			
			if (bigendianFlag){
				// don't convert, just copy into Fortran (column oriented) arrays t,wv
				t[i] = (double) tmpT;
				for (int32_t j = 0; j<4; j++)
					for(int32_t k = 0; k<32; k++)
				  	   wv[i + nSpikes*j + nSpikes*4*k] = (double) tmpWV[j + 4*k];
			} else { 
				// convert from Sun (big endian) to NT(little endian) 
				t[i] = (double) swapbytes(tmpT);
				for (int32_t j = 0; j<4; j++)
					for(int32_t k = 0; k<32; k++)
				  	   wv[i + nSpikes*j + nSpikes*4*k] = (double) swapbytes(tmpWV[j + 4*k]);
			}	
		}
	}
   fclose(fp);	
}	



//___________________________________________________________________________________
void mexFunction(int32_t nOUT, mxArray *pOUT[],
				 int32_t nINP, const mxArray *pINP[])
{
	
	int32_t nProbes = 4;
	/* check number of arguments: expects 1 input */
	if (nINP != 1 && nINP != 2)
				mexErrMsgTxt("Call with fn or fn, nProbes (default =4) as inputs.");
	if (nOUT != 2)
				mexErrMsgTxt("Requires two outputs (t, wv)");

	/* read inputs */
	int32_t fnlen = (mxGetM(pINP[0]) * mxGetN(pINP[0])) + 1;
	char *fn = (char *) mxCalloc(fnlen, sizeof(char)); 
	if (!fn)
		mexErrMsgTxt("Not enough heap space to hold converted string.");
	int32_t errorstatus = mxGetString(pINP[0], fn,fnlen);    
	if (errorstatus)
		mexErrMsgTxt("Could not convert string data.");
	
	if (nINP == 2)
		nProbes = (int32_t)mxGetScalar(pINP[1]);
	int32_t nSpikes = GetNumberOfSpikes(fn, nProbes);

	

	// create outputs 
	pOUT[0] = mxCreateDoubleMatrix(nSpikes, 1, mxREAL);
	double *t = mxGetPr(pOUT[0]);
	
	int32_t wvDims[] = {nSpikes, nProbes, 32};    // wv = nSpikes x 4 x 32   FORTRAN (column major) array
	pOUT[1] = mxCreateNumericArray(3, wvDims, mxDOUBLE_CLASS, mxREAL);
	double *wv = mxGetPr(pOUT[1]);

	// load tt file fn into t and wv arrays 
	ReadTT(fn,nSpikes,t,wv, nProbes);

	// cleanup
	mxFree(fn);


}

