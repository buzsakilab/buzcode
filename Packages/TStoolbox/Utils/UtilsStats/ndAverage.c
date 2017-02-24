/* n-dimensional average map, points are binned accrodign to their
   n-dimensional location and the average of data is computed 
   for each bin 
   input: 
   samples (nD x nSamp) -- value of each dimension at each sample
   data: (1 x nSamp) -- the data to average 
 *        bins (nD x 1) -- how many bins for each dimension
 *        mins (nD x 1) -- min values for each dimension
 *        maxes (nD x 1) -- max values for each dimension


fpbatta 2005, modified from ADR code 

*/ 



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
		binToAddTo = ((int) floor((subs[d] - mins[d])/(maxes[d] - mins[d]) * (bins[d]-1)));
	/* sprintf(MsgBuffer, "%d: %d -> %d\n", d, (int) subs[d], binToAddTo);
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



void nDSum(int nD, int nSamp, double *S, double *data,  
	    double *bins, double *mins, double *maxes, 
	    double *A)
{
  int iS;

  for (iS = 0; iS < nSamp; iS++)
    {      
      int b = toBin(nD, bins, mins, maxes, S + iS * nD);
      A[b] = A[b] + data[iS];
    }
}

/*--------------------
 * gateway function
 --------------------*/

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  int nD, nSamp, i, totalsize=1;
  double *S, *data, *bins, *mins, *maxes, *H, *A;

  /* check number of arguments: expects 2 inputs, 1 output */
  if (nINP != 5)
    mexErrMsgTxt("Call with samples, data, bins, mins, maxes as inputs");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs */
  nD = mxGetM(pINP[0]);
  nSamp = mxGetN(pINP[0]);

  if ((mxGetM(pINP[1]) != 1) || (mxGetN(pINP[1]) != nSamp))
    mexErrMsgTxt("Wrong dimensionality of data");
  
  if ((mxGetM(pINP[2]) != nD) || (mxGetN(pINP[2]) != 1))
      mexErrMsgTxt("Wrong dimensionality of bins.");
  if ((mxGetM(pINP[3]) != nD) || (mxGetN(pINP[3]) != 1))
    mexErrMsgTxt("Wrong dimensionality of mins.");
  if ((mxGetM(pINP[4]) != nD) || (mxGetN(pINP[4]) != 1))
    mexErrMsgTxt("Wrong dimensionality of maxes.");
  
  /* unpack inputs */
  S = mxGetPr(pINP[0]);
  data = mxGetPr(pINP[1]);
  bins  = mxGetPr(pINP[2]);
  mins  = mxGetPr(pINP[3]);
  maxes = mxGetPr(pINP[4]);

  /* create outputs */
  {
    int d;
	int *dim = (int*) mxCalloc(sizeof(int), nD);
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

    /* the histogram count */

/*     H = mxGetPr(pOUT[0]);  */
  
    H = mxCalloc(totalsize, sizeof(double));
    
    /* the output is the binned average */

    A = mxGetPr(pOUT[0]);
    

    mxFree(dim);
  }

  /* calculate */
  nDHist(nD, nSamp, S, bins, mins, maxes, H);
  nDSum(nD, nSamp, S, data, bins, mins, maxes, A);

  for(i = 0; i < totalsize; i++)
    A[i] /= H[i];
  


  mxFree(H);
  

}
