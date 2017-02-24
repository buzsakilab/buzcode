/* Restrict_idx_iSet.c: MEX file that optimizes the job of the
   Restrict for the intervalSet case. private function in the tsd
   class

   ix = Restrict_idx_iSet(range_tsd, int_start, int_end);

   range_tsd, is the ts vector of the input tsd
   int_start, int_end contain the extremes of the restriciting
   intervalSet
   returns 1-based indices in ix


*/  



#include <mex.h>
#include <math.h>
#include <matrix.h>
#include <string.h>



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{

  int n_points, n_intervals;
  double *t, *t0, *t1, *ix;
  int i, i1, k, n;
  

  if(nrhs != 3)
    mexErrMsgTxt("Requires three inputs");
  
  if(nlhs != 1)
    mexErrMsgTxt("Requires one output");
  
  
  /* unpack inputs */

  /* input 0 is t */

  if(mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("t must be a row or column vector");
  
  n_points = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  
  t = mxGetPr(prhs[0]);
  
  /* input 1 is t0 */

  if(mxGetM(prhs[1]) != 1 && mxGetN(prhs[1]) != 1)
    mexErrMsgTxt("t0 must be a row or column vector");
  
  n_intervals = mxGetM(prhs[1]) * mxGetN(prhs[1]);

  t0 = mxGetPr(prhs[1]);
  
  /* input 2 is t1 */

  if(mxGetM(prhs[2]) != 1 && mxGetN(prhs[2]) != 1)
    mexErrMsgTxt("t1 must be a row or column vector");

  if(mxGetM(prhs[2]) != mxGetM(prhs[2]) || mxGetN(prhs[2]) != mxGetN(prhs[2]))
    mexErrMsgTxt("t1 and t0 must have same size");
  
  t1 = mxGetPr(prhs[2]);
  
  for(i = 0; i < n_intervals; i++)
    if(t1[i] < t0[i])
      mexErrMsgTxt("t1 must be greater than t0");
  

  /* prepare output */
  ix = mxCalloc(n_points, sizeof(double));
  



  i = 0;
  i1 = 0;
  k = 0;  
  n = 0;
  
  while(i < n_points)
    {
      i = i1;
      
      while(i < n_points && t[i] < t0[n])
	i++;
      i1 = i;
      
 
      while(i < n_points && t[i] < t1[n])
	{
	  ix[k] = i+1.;
	  
	  i++;
	  k++;
	}
      
      n++;
      if(n >= n_intervals)
	break;
      
    }

  plhs[0] = mxCreateNumericMatrix(k, 1, mxDOUBLE_CLASS, mxREAL);
  
  memcpy(mxGetPr(plhs[0]), ix, k * sizeof(double));
  
  mxFree(ix);
  
  
  
}
