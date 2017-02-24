/* Histogram_intervals MEX file C implementation for the Histogram function 
 * INPUTS: 
 * data: the data to count
 * t0, t1, arrays of values for start and stop of the intervals
 * OUTPUTS: 
 * c: vector of counts
 *
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

  if(data[start] > key)
    {
      result = start;
      goto alldone;
    }
  while (start < (end-1))
    {
      mid = floor((start + end)/2);
      if ((key) == data[mid])
	{
	  start = end = mid;
	  break;
	}
      if ((key) < data[mid])
	end = mid;
      if ((key) > data[mid])
	start = mid;
    }

  result = end;
  if(key == start)
    result = start;

 alldone:
/*   mexPrintf("debug: key = %g, result = %g\n", key, data[result]); */

  return result;
  
}


void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{

  double *data, *tstmp1, *tstmp2, *c;
  int n_data, n_tstmp;
  int i, j, j0, k;
  


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
  

  /* pack outputs */
  pOUT[0] = mxCreateDoubleMatrix(n_tstmp, 1, mxREAL);
  c = mxGetPr(pOUT[0]);
  

  if(n_data == 0)
    return;

  /* tstmp1 and tstmp2 must have the same length */
  if((mxGetM(pINP[2]) * mxGetN(pINP[2])) != n_tstmp)
    mexErrMsgTxt("Tstart and Tend must have the same length");
  tstmp2 = (double *)mxGetPr(pINP[2]);
  


  

  j = 0;
  j0 = 0;
  
  for(i = 0; i < n_tstmp; i++)
    {
      k = 0;
      j = binsearch_past(data, tstmp1[i], j0, n_data);
      j0 = j;
      while(j < n_data && data[j] < tstmp2[i]) {
	k++;
	j++;
      } 
      c[i] = k;
      
    }
  

  


}
