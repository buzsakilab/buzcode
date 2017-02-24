 /*-----------------------------------
 * n-dimensional histogram
 * MEX file
 * 
 * ADR 1998
 * 
 * input: samples (nD x nSamp) -- value of each dimension at each sample
 *        bins (nD x 1) -- how many bins for each dimension
 *        mins (nD x 1) -- min values for each dimension
 *        maxes (nD x 1) -- max values for each dimension
 * output: histogram (ndim multidimensional array) -- histogram
 *
 * version L1.0
 * status promoted
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <string.h>

char MsgBuffer[100];
/*--------------------
 * to bin : converts subscripts -> bin entry
 *--------------------*/
int toBin(int nD, 
	  double *bins, double *mins, double *maxes, 
	  double *subs)
{
  int d;
  int retval;
  int binToAddTo;

  retval = 0;
  for (d=nD-1; d>=0; d--)
  {
	if (subs[d] <= mins[d])
		binToAddTo = 0;               /* put in lowest bin i.e. subs = mins -> add 0 */
	else if (subs[d] >= maxes[d])
		binToAddTo = (int)bins[d]-1;
	else
		/*binToAddTo = ((int) floor(((subs[d] - mins[d])/(maxes[d] - mins[d])) * (bins[d]-1)));
		* changed, ncst 20 Aug 02 */
		binToAddTo = ((int) floor(((subs[d] - mins[d])/(maxes[d] - mins[d])) * (bins[d])));
		
	 /*mexPrintf("%d: %d -> %d\n", d, (int) subs[d], binToAddTo);
	 *mexPrintf("Max: %d Min: %d #Bins %d\n", (int) (maxes[d] - mins[d]), (int) (subs[d] - mins[d]), (int) ((subs[d] - mins[d])/(maxes[d] - mins[d])) * (bins[d]-1)*100);}
	* mexWarnMsgTxt(MsgBuffer); */
	retval = retval * (int) bins[d] + binToAddTo;
  };
  return retval;
}

/*--------------------
 * histogram function
 *--------------------*/

void nDHist(int nD, int nSamp, double *S, 
	    double *bins, double *mins, double *maxes, 
	    double *H)
{
  int iS;

  for (iS = 0; iS < nSamp; iS++)
    {      
      int b = toBin(nD, bins, mins, maxes, S + iS * nD);
      H[b] = H[b] + 1;
    }
}

/*--------------------
 * gateway function
 --------------------*/

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  int nD, nSamp;
  double *S, *bins, *mins, *maxes, *H;

  /* check number of arguments: expects 2 inputs, 1 output */
  if (nINP != 4)
    mexErrMsgTxt("Call with samples, bins, mins, maxes as inputs");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs */
  nD = mxGetM(pINP[0]);
  nSamp = mxGetN(pINP[0]);
  if ((mxGetM(pINP[1]) != nD) || (mxGetN(pINP[1]) != 1))
      mexErrMsgTxt("Wrong dimensionality of bins.");
  if ((mxGetM(pINP[2]) != nD) || (mxGetN(pINP[2]) != 1))
    mexErrMsgTxt("Wrong dimensionality of mins.");
  if ((mxGetM(pINP[3]) != nD) || (mxGetN(pINP[3]) != 1))
    mexErrMsgTxt("Wrong dimensionality of maxes.");
  
  /* unpack inputs */
  S = mxGetPr(pINP[0]);
  bins  = mxGetPr(pINP[1]);
  mins  = mxGetPr(pINP[2]);
  maxes = mxGetPr(pINP[3]);

  /* create outputs */
  {
    int i,d, totalsize=1;
	int *dim = (int*) calloc(sizeof(int), nD);
	if (!dim)
		mexErrMsgTxt("Cannot create dimension array. Probably out of memory.");
    for (d=0; d<nD; d++)
      {
		dim[d] = (int)floor(bins[d]);
		totalsize = totalsize * dim[d];
      }

    if (nD == 1)
      pOUT[0] = mxCreateDoubleMatrix(dim[0], 1, mxREAL);
    else if (nD == 2)
      pOUT[0] = mxCreateDoubleMatrix(dim[0], dim[1], mxREAL);
    else 
      pOUT[0] = mxCreateNumericArray(nD, dim, 
				     mxDOUBLE_CLASS, mxREAL); 
    if (!pOUT[0]) 
      mexErrMsgTxt("Cannot create output array. Probably out of memory.");  
    H = mxGetPr(pOUT[0]); 
    if (!H)
      mexErrMsgTxt("INTERNAL: H not allocated.  Please send email to adr@nsma.");
    
    for (i=0; i<totalsize; i++)
      H[i] = 0;

	free(dim);
  }

  /* calculate */
  nDHist(nD, nSamp, S, bins, mins, maxes, H);
}
