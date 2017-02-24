/* Restrict_align.c MEX file optimization of the findAlignment routine
 * INPUTS:
 * Data, tstmp: vectors of timestamps 
 * OUTPUTS:
 * ix: vector of indices in Data for the elements that are the closest
 * to those in  tstmp
 * 
 * reimplemented in a mex file from adr 1999 code
 * batta 2001 */




#include "mex.h"
#include <math.h>
#include <matrix.h>


int binsearch_closest(double *data, double key, int min_idx, int max_idx)
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
  if (((key) - data[start]) <= (data[end] - (key)))
    result = start + 1;
  else 
    result = end + 1;
  return result;
  
}


void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{

  double *data, *tstmp, *ix;
  int n_data, n_tstmp;
  int i, j, k;
  


  /* check number of arguments: expects 2 inputs, 1 output */
  if (nINP != 2)
    mexErrMsgTxt("Call with Data, key as inputs.");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs, they must be row or column vectors */

  if(mxGetM(pINP[0]) != 1 && mxGetN(pINP[0]) != 1)
    mexErrMsgTxt("Data must be a row or column vector.");
  if(mxGetM(pINP[1]) != 1 && mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("Tstmp must be a row or column vector.");
  
  /* unpack inputs */
  n_data = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  data = (double *)mxGetPr(pINP[0]);
  n_tstmp = mxGetM(pINP[1]) * mxGetN(pINP[1]);
  tstmp = (double *)mxGetPr(pINP[1]);
  
  /* pack outputs */
  pOUT[0] = mxCreateDoubleMatrix(1, n_tstmp, mxREAL);
  ix = (double *) mxGetPr(pOUT[0]);


  j  = binsearch_closest(data, tstmp[0], 1, n_data);
  ix[0] = (double)j;
  
  for(i=1; i < n_tstmp; i++)
    {
      if(mxIsFinite(tstmp[i]))
	{
	  while((data[j] < tstmp[i]) && (j < n_data-1))
	    j++;
	  k = ((tstmp[i] - data[j-1]) < (data[j] - tstmp[i]) ?
	       (j-1) : j);
	  ix[i] = (double)k;
	}
      else
	ix[i] = mxGetNaN();
    }
  


}
