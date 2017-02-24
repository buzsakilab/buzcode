/* Restrict_intervals.c MEX file, optimization of the 3 argument form
 * for the REstrict function
 * INPUTS:
 * Data, tstmp1, tstmp2: vectors of timestamps, first is the vector to
 * restrict, tstmp1 and tstmp2 are, respetively the beginnings and the
 * ends of the intervals in which to Restrict
 * OUTPUT:
 * ix: vector of indices in Data for the elements that are in the intervals
 *
 * reimplemented in a mex file from adr 1999 code
 * batta 2001 */




#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <string.h>

int binsearch_past(double *data, double key, int min_idx, int max_idx)
{
  
  int mid, start = min_idx, end = max_idx-1, result;
  while (start < (end-1))
    {
      mid = floor((start + end)/2);
      if ((key) == data[mid])
	start = end = mid;
      if ((key) < data[mid])
	end = mid;
      if ((key) > data[mid])
	start = mid;
    }

  result = end;
  return result;
  
}


void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{

  double *data, *tstmp1, *tstmp2, *ix;
  int n_data, n_tstmp;
  int i, j, k;
  


  /* check number of arguments: expects 2 inputs, 1 output */
  if (nINP != 3)
    mexErrMsgTxt("Call with Data, tstart, tend as inputs.");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs, they must be row or column vectors */

  if(mxGetM(pINP[0]) != 1 && mxGetN(pINP[0]) != 1)
    mexErrMsgTxt("Data must be a row or column vector.");
  if(mxGetM(pINP[1]) != 1 && mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("Tstart must be a row or column vector.");
  if(mxGetM(pINP[2]) != 1 && mxGetN(pINP[2]) != 1)
    mexErrMsgTxt("Tend must be a row or column vector.");




  /* unpack inputs */
  n_data = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  data = (double *)mxGetPr(pINP[0]);
  n_tstmp = mxGetM(pINP[1]) * mxGetN(pINP[1]);
  tstmp1 = (double *)mxGetPr(pINP[1]);
  
  /* tstmp1 and tstmp2 must have the same length */
  if((mxGetM(pINP[2]) * mxGetN(pINP[2])) != n_tstmp)
    mexErrMsgTxt("Tstart and Tend must have the same length");
  tstmp2 = (double *)mxGetPr(pINP[2]);
  

  ix = mxCalloc(n_data, sizeof(double));
  
  k = 0;
  j = 0;
  
  for(i = 0; i < n_tstmp; i++)
    {
      j = binsearch_past(data, tstmp1[i], j, n_data);
      while(data[j] <= tstmp2[i] && j < n_data) {
	ix[k++] = (double)j + 1.;
	j++;
      } 
      
    }
  
  /* pack outputs */
  pOUT[0] = mxCreateDoubleMatrix(k, 1, mxREAL);
  memcpy(mxGetPr(pOUT[0]), ix, k * sizeof(double));
  


}
