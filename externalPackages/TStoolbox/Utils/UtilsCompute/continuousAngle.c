 /* thetac = continuosAngle(theta)

 takes a variable intepreted as an angle, and at each point adds a
 multiple of 2*PI in order to avoid jumps in variable

 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  int n, i;

  double m, *data_in, *data_out;
  
  
  /* check number of arguments: expects 1 input, 1 output */
  if (nINP != 1)
    mexErrMsgTxt("Call with theta");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs */
  if (mxGetM(pINP[0]) > 1 &&  mxGetN(pINP[0]) > 1)
    mexErrMsgTxt("Input must be row or column vector");

  n = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  data_in = mxGetPr(pINP[0]);
  

  /* prepare output */
  pOUT[0] = mxCreateDoubleMatrix(mxGetM(pINP[0]), mxGetN(pINP[0]), mxREAL);
  data_out = mxGetPr(pOUT[0]);
  
  
  data_out[0]= data_in[0];

  
  
  for(i= 1;i < n; i++)
    {
      m = nearbyint((data_in[i]- data_out[i-1])/ (2*M_PI));
      data_out[i] = data_in[i] - m * 2 * M_PI;
    }


}
  
  
		 
