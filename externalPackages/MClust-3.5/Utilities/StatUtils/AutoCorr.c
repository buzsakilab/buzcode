 /* AutoCorr
 * auto correlations
 * MEX file
 * 
 * batta 1999, lipa 2000
 * 
 * input: t1: a time series to auto correlate in 1/10000 sec
 *                (assumed to be sorted) 
 *        binsize: the binsize for the auto corr histogram in msec
 *        nbins: the number of bins
 * output: C the auto correlation histogram
 *         B (optional) a vector with the times corresponding to the bin centers
 *
 * version 1.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  
  double *t1;
  double binsize;
  double *C, *cross, *B;
  double w, lbound, rbound;
  
  int nbins, nt1, nt2;
  int i1 = 0, i2 = 0, l, j, k;
  
  /* check number of arguments: expects 3 inputs, 1 or 2 outputs */
  if (nINP != 3)
    mexErrMsgTxt("Call with t1,  binsize and nbins  as inputs.");
  if (nOUT != 1 && nOUT != 2)
    mexErrMsgTxt("Requires one or two outputs.");

  /* check validity of inputs */
  if (mxGetM(pINP[0]) != 1 && mxGetN(pINP[0]) != 1)
    mexErrMsgTxt("t1 must be a row or column vector");
  if (mxGetM(pINP[1]) * mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("binsize must be scalar");
  if (mxGetM(pINP[2]) * mxGetN(pINP[2]) != 1)
    mexErrMsgTxt("nbins must be scalar");

  /* unpack inputs */
  nt1 = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  t1 = mxGetPr(pINP[0]);
  nt2 = nt1;
  binsize = mxGetScalar(pINP[1]);
  nbins = (int)mxGetScalar(pINP[2]);
  




  

  pOUT[0] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
  C = mxGetPr(pOUT[0]);
  if(nOUT == 2)
  {
      double m;
      
      pOUT[1] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
      B =  mxGetPr(pOUT[1]);
      m = binsize/2.0;
      for(j = 0; j < nbins; j++)	B[j] = m + j * binsize;

  }
  
  binsize *= 10;                 /* convert to timestamp units (msec/10) */
  /* cross correlations */
  
  w = nbins * binsize;
   
  for(i1 = 0; i1 < nt1; i1++)
  {
      l = i1+1;
      rbound = t1[i1];
      
      for(j = 0; j < nbins; j++)
	  {
		  k = 0;
		  rbound += binsize;
		  while(t1[l] < rbound && l < nt2-1)
		  {
			  l++;
			  k++;
		  }
		  C[j] += k;
	  }
  }
  
  for(j = 0; j < nbins; j++)  C[j] /= nt1 * binsize / 10000;
  
      
}
  
  
		 
