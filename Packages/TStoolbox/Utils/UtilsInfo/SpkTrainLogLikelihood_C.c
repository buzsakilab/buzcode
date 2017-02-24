 /* CrossCorr
 * cross correlations
 * MEX file
 * 
 * batta 1999
 * 
 * input: t1, t2: two time series to cross correlate in 1/10000 sec
 *                (assumed to be sorted) 
 *        binsize: the binsize for the cross corr histogram in msec
 *        nbins: the number of bins
 * output: C the cross correlation histogram
 *         B (optional) a vector with the times corresponding to the bins
 *
 * version 1.0
 -----------------------------------*/

/* copyright (c) 1999 Francesco P. Battaglia
 * This software is released under the GNU GPL
 * www.gnu.org/copyleft/gpl.html
 */



#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  
  double *q;
  double *f;
  double *L;
  int i1, nt1, nt2;
  
  /* check number of arguments: expects 4 inputs, 1 or 2 outputs */
  if (nINP != 2)
    mexErrMsgTxt("CAll with q and f");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs 
  if (mxGetM(pINP[0]) != 1 && mxGetN(pINP[0]) != 1)
    mexErrMsgTxt("t1 must be a row or column vector");
  if (mxGetM(pINP[1]) != 1 && mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("t2 must be a row or column vector");*/

  /* unpack inputs */
  nt1 = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  q = mxGetPr(pINP[0]);
  /*nt2 = mxGetM(pINP[1]) * mxGetN(pINP[1]);*/
  f = mxGetPr(pINP[1]);
  /*if (nt1 != nt2)
      mexErrMsgTxt("tvectors must be of same length");*/

  pOUT[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  L = mxGetPr(pOUT[0]);
  
  for(i1 = 0; i1 < nt1; i1++)
    {
      if (f[i1]>0)
              *L += log(f[i1]) * q[i1] - f[i1];
  }
            
}
  
  
		 
