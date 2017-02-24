 /*SlidingNorm-----------------------------------
 * Norm of a sliding window of elements
 * MEX file
 * 
 * batta 1999
 * 
 * input: Data -- n x 1
 *        winsize
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
  int i,j, ws;
  
  double *winsize;
  double *data;
  double r, *norm;
  
  /* check number of arguments: expects 2 inputs, 1 output */
  if (nINP != 2)
    mexErrMsgTxt("Call with Data, winsize as inputs.");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs */
  if (mxGetM(pINP[1]) * mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("Winsize must be a scalar");

  /* unpack inputs */
  n = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  data = (double *)mxGetPr(pINP[0]);
  winsize = (double *)mxGetPr(pINP[1]);
  ws = (int)*winsize;
  mexPrintf("winsize = %g ws = %d\n", *winsize, ws);
  
  /* pack outputs */
  pOUT[0] = mxCreateDoubleMatrix(1, n-ws+1, mxREAL);
  norm = (double *) mxGetPr(pOUT[0]);

  
  for(i=0; i<n-ws;i++) 
    {
      r = 0.;
      for(j=0; j<ws; j++)
	r += data[i+j] * data[i+j];
      norm[i] = sqrt(r);
    }
  
      
  


}
  
  
		 
