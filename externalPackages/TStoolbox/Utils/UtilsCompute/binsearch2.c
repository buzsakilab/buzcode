 /*-----------------------------------
 * binary search
 * MEX file
 * 
 * ADR 1998
 * 
 * input: Data -- n x 1, assumed to be sorted
 *        key
 * output: index into Data
 *
 * version 1.0
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  int n;
  int start, end, mid;
  double *key;
  double *data;
  double *result;
  
  /* check number of arguments: expects 2 inputs, 1 output */
  if (nINP != 2)
    mexErrMsgTxt("Call with Data, key as inputs.");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs */
  if (mxGetM(pINP[1]) * mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("Can only use one key at a time.");

  /* unpack inputs */
  n = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  data = (double *)mxGetPr(pINP[0]);
  key = (double *)mxGetPr(pINP[1]);

  /* binary search */
  start = 0;
  end = n-1;
  while (start < (end-1))
    {
      mid = floor((start + end)/2);
      /*fprintf(stderr, "binsearch (%.0f): %d(%.0f) %d(%.0f) %d(%.0f)\n",
       *key, start, data[start], mid, data[mid], end, data[end]);*/
      if ((*key) == data[mid])
	start = end = mid;
      if ((*key) < data[mid])
	end = mid;
      if ((*key) > data[mid])
	start = mid;
    }

  if(((*key) - data[start]) > (data[start+1] - (*key)))
    start++;
  
  
  /* pack outputs */
  pOUT[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  result = (double *) mxGetPr(pOUT[0]);
  *(result) = (double) start + 1;
}
  
  
		 
