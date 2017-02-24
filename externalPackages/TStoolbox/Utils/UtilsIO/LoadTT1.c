/*---------------------------------
* Load TT
* MEX file
*
* ADR 1998
*
* modified by batta 1999 to read the TT file in chunks
* modified by batta 1999 to get rid of bad timestamps, i.e. those
* generating negative ISIs.
*
* input:
*    fn = file name string
*    start = start from this spike
*    chunk_length = length of the lbock to be read
* 
* output:
*    [t, wv, tbad, wvbad]
*    t = n x 1: timestamps of each (good) spike in file
*    wv = n x 4 x 32 waveforms
*    tbad = nbad x 1,
*    wvbad = nbad x 4 x 32: timestamps and waveforms of "bad" spikes
*
* 
*
* v 5.0 modified to skip last incomplete records
* status: BETA



--------------------------------*/

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

#define CHECK_FOR_BAD_ISI

void SkipHeader(FILE *fp)
{
	long curpos = ftell(fp);
	char headerline[81];

	fgets(headerline, 80, fp);
	if (strncmp(headerline, "%%BEGINHEADER", 13) == 0)
		while (strncmp(headerline, "%%ENDHEADER",11) != 0)
			fgets(headerline, 80, fp);
	else
		fseek(fp, curpos, SEEK_SET);
}

int bigendianMachine(void)
{
	/* returns 0 if is a littleendian machine, else returns nonzero */
	/* key is that it looks to see if short's second byte is the low order bits */
	short fullLoByte = 0xFF;
	char *byteOrder = (char *) &fullLoByte;
	return (byteOrder[1]);	
}

#define RECSIZE ((sizeof(unsigned long)) + 4 * 32 * sizeof(short))

void mexFunction(int nOUT, mxArray *pOUT[],
				 int nINP, const mxArray *pINP[])
{
	int errorstatus;
	int bigendianFlag = bigendianMachine();
	long postHeaderPos, endPos;
	int fnlen;
	int start, chunk_length;
	int nGoodSpikes, nBadSpikes = 0;
	int outputBadSpikes = 0;
	int *BadSpikesIdx;
	

	char *fn;
	FILE *fp;
	int nSpikes;
	double *t;
	double *wv;
	double *tbad;
	double *wvbad;
	
	int i,j, k, b,q,trode, ch;  /* counters */

	/* check number of arguments: expects 1 or three inputs */
	if (nINP != 1)
	  {
	    if(nINP != 3)
	      mexErrMsgTxt("Requires one or three inputs");
	    else
	      {
		start = mxGetScalar(pINP[1]);
		chunk_length = mxGetScalar(pINP[2]);
	      }
	  }
	else
	  {
	    start = 0;
	    chunk_length = 0;
	  }
	
	if (nOUT != 2 && nOUT != 4)
		mexErrMsgTxt("Requires two outputs (t, wv) or four outputs (t, wv, tbad, wvbad).");
	if(nOUT == 4)
	  outputBadSpikes = 1;
	
	/* read inputs */
	fnlen = (mxGetM(pINP[0]) * mxGetN(pINP[0])) + 1;
	fn = (char *) mxCalloc(fnlen, sizeof(char)); 
	if (!fn)
	  mexErrMsgTxt("Not enough heap space to hold converted string.");
	errorstatus = mxGetString(pINP[0], fn,fnlen);    
	if (errorstatus)
		mexErrMsgTxt("Could not convert string data.");
		
	/* open file */
	fp = fopen(fn, "rb");
	mxFree(fn);
	if (!fp)
		mexErrMsgTxt("Could not open file.");

	/* skip header */
	SkipHeader(fp);
	postHeaderPos = ftell(fp);

	/* count number of spikes 
	for (nSpikes = 0; !feof(fp); nSpikes ++)
	 {
		char buf[4 + 256];			
		fread(buf, sizeof(char), 4 + 256, fp);
	 }
	*/

        fseek(fp, 0, SEEK_END);
	endPos = ftell(fp);
	nSpikes = floor((endPos - postHeaderPos) / RECSIZE);
	if(chunk_length == 0)
	  chunk_length = nSpikes;
	if(start > nSpikes-1)
	  {
	    int wvDims[] = {0, 4, 32};
	    mexPrintf("Reading %d spikes.\n", 0);

	    pOUT[0] = mxCreateDoubleMatrix(0, 1, mxREAL);
	    pOUT[1] = mxCreateNumericArray(3, wvDims, mxDOUBLE_CLASS, mxREAL);
	    return;
	  }
	



	if((nSpikes-start)<chunk_length)
	  chunk_length = nSpikes-start;
	

	mexPrintf("Reading %d spikes.\n", chunk_length);
	t = (double *) mxCalloc(chunk_length, sizeof(double));
	wv = (double *) mxCalloc(chunk_length * 128, sizeof(double));
	if (!t || !wv)
		mexErrMsgTxt("Not enough heap space for TT components.");

	/* read t, wv */
	fseek(fp, postHeaderPos, SEEK_SET);

	
	fseek(fp, (256+4)*start, SEEK_CUR);

	for (i = 0; i < chunk_length; i++)
	 {
		union {
			char c[4];
			unsigned long ul;
		} tmpT,tmpT0;

		union {
			char c[2];
			short s;
		} tmpWV[128], tmpWV0[128];

		fread(&tmpT,  sizeof(char),   4, fp);
		fread(tmpWV, sizeof(char), 256, fp);

		if (bigendianFlag)
		 {
			memcpy(&tmpT0, &tmpT, 4);
			memcpy(tmpWV0, tmpWV, 256);
		 }
		else
		 { 
			tmpT0.c[0] = tmpT.c[3];
			tmpT0.c[1] = tmpT.c[2];
			tmpT0.c[2] = tmpT.c[1];
			tmpT0.c[3] = tmpT.c[0];
			for (j=0; j < 128; j++)
			 { 
				tmpWV0[j].c[0] = tmpWV[j].c[1];					
				tmpWV0[j].c[1] = tmpWV[j].c[0];
			 }
		}

		t[i] = (double) tmpT0.ul;
		/*  mexPrintf("i = %d t[i] = %g\n", i, t[i]); */
		
		for (j = 0; j<128; j++)
			wv[j + 128 * i] = (double) tmpWV0[j].s;
	 }

#ifdef CHECK_FOR_BAD_ISI
	/* check for negative timestamps and throw the bad points away */

	nGoodSpikes = 0;
	nBadSpikes = 0;
	i = 0;
	j = 0;
	mexPrintf("t[0] = %g\n");
	
	while(i < chunk_length)
	  {
	    
	    nGoodSpikes++;
	    j = 1;	    
	    if (i+j >=chunk_length) break;
	    while(t[i+j] < t[i] && i+j <chunk_length)
	      {
		
		nBadSpikes++;
		j++;
	      }
	    if (i+j >=chunk_length) break;
	    i = i+j;
	  }
	if((nGoodSpikes + nBadSpikes) != chunk_length)
	  {
	    char ErrMsg[120];
	    sprintf(ErrMsg, "Problem with the Bad timestamps code, nGoodSpikes = %d, nBadSpikes = %d", nGoodSpikes, nBadSpikes);
	    mexErrMsgTxt(ErrMsg);
	  }
	mexPrintf("nGoodSpikes = %d, nBadSpikes = %d\n", nGoodSpikes,
		  nBadSpikes);
	
	BadSpikesIdx = (int *) mxCalloc(nBadSpikes, sizeof(int));
	if(outputBadSpikes)
	  {
	    tbad = (double *) mxCalloc(nBadSpikes, sizeof(double));
	    wvbad = (double *) mxCalloc(nBadSpikes * 128, sizeof(double));
	  }
	
	
	nGoodSpikes = 0;
	nBadSpikes = 0;
	i = 0;
	j = 0;
	
	while(i < chunk_length)
	  {

	    nGoodSpikes++;
	    j = 1;	    
	    if (i+j >=chunk_length) break;
	    while(t[i+j] < t[i] && i+j <chunk_length)
	      {
		BadSpikesIdx[nBadSpikes] = i+j;
		nBadSpikes++;
		j++;
	      }
	    if (i+j >=chunk_length) break;
	    i = i+j;
	  }

	i = 0;
	j = 1;
	q = chunk_length;
	
	while(i < nBadSpikes)
	  {
	    b = BadSpikesIdx[i]-i;
		 
	    if(outputBadSpikes)
	      {
		tbad[i] = t[b];
		memcpy(&wvbad[128*i], &wv[128*b], 
		       128*sizeof(double));
	      }
	    if((BadSpikesIdx[i+1] - BadSpikesIdx[i]) == 1)
	      {
		i++;
		j++;
		continue;
	      }
	    for(k = BadSpikesIdx[i-j+1]-(i-j+1);k+j<q;k++)
	      {
		t[k] = t[k+j];
		memcpy(&wv[128*k], &wv[128*(k+j)], 128*sizeof(double));
		
	      }
	    q -= j;
	    j=1;
	    i++;
	    mexPrintf("i = %d\n", i);
	    
	  }
	
	if(q != nGoodSpikes)
	  {	
	    char ErrMsg[120];
	    
	    sprintf(ErrMsg, "Problem with the bad spikes discarding code, q = %d", q);
	    
	    mexErrMsgTxt(ErrMsg);
	  }
	
	chunk_length = nGoodSpikes;
	
	mexPrintf("Now creating outputs\n");
	
#endif
	/* create outputs */
	pOUT[0] = mxCreateDoubleMatrix(chunk_length, 1, mxREAL);
	memcpy(mxGetPr(pOUT[0]), t, chunk_length * sizeof(double));
	mxFree(t);






	{
	 int wvDims[] = {chunk_length, 4, 32};
	 double *wvptr;
	 pOUT[1] = mxCreateNumericArray(3, wvDims, mxDOUBLE_CLASS, mxREAL);
	 wvptr = mxGetPr(pOUT[1]);
 	 for (j = 0, i = 0; i<chunk_length; i++)
		for (ch = 0; ch < 32; ch++)
		for (trode = 0; trode < 4; trode++, j++)
			 {
				int subs[] = {i, trode, ch};
				int index = mxCalcSingleSubscript(pOUT[1], 3, subs);  
				wvptr[index] = wv[j];
			}

	}
	mxFree(wv);

#ifdef CHECK_FOR_BAD_ISI
	if(outputBadSpikes)
	  {
	    pOUT[2] = mxCreateDoubleMatrix(nBadSpikes, 1, mxREAL);
	    memcpy(mxGetPr(pOUT[2]), tbad, nBadSpikes * sizeof(double));
	    mxFree(tbad);

	    {
	      int wvDims[] = {nBadSpikes, 4, 32};
	      double *wvptr;
	      pOUT[3] = mxCreateNumericArray(3, wvDims, mxDOUBLE_CLASS, mxREAL);
	      wvptr = mxGetPr(pOUT[3]);
	      for (j = 0, i = 0; i<nBadSpikes; i++)
		for (ch = 0; ch < 32; ch++)
		  for (trode = 0; trode < 4; trode++, j++)
		    {
		      int subs[] = {i, trode, ch};
		      int index = mxCalcSingleSubscript(pOUT[3], 3, subs);  
		      wvptr[index] = wvbad[j];
		    }

	    }


	    
	  }

#endif	
	fclose(fp);

}
