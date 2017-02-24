/*  function c = Histogram_intervals(t, t0, t1) */
/*  % c = Histogram_intervals(t, t0, t1) (MEX file) */
/*  % */
/*  % returns a vector of counts of elements of t contained in each of the */
/*  % intervals defined by the extremes t0 (low extreme) and t1 (high */
/*  % extreme) */
/*  % INPUTS: */
/*  % t: a row or column array of time points, assumed to be sorted */
/*  % t0, t1: sorted arrays of interval extremes such that t0(i) < t1 (i) for */
/*  % all i */
/*  % OUTPUTS: */
/*  % x = an array of counts */
  
/*  % status: under construction */
/*  % batta 2002   */

/* copyright (c) 2004 Francesco P. Battaglia */
/* This software is released under the GNU GPL */
/* www.gnu.org/copyleft/gpl.html */

#include <matrix.h>
#include <mex.h>





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{

  int n_points, n_intervals;
  double *t, *t0, *t1, *c;
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

  plhs[0] = mxCreateNumericMatrix(1, n_intervals, mxDOUBLE_CLASS, mxREAL);
  
  c = mxGetPr(plhs[0]);
  
  
  if(n_points == 0)
    return;
  

  /* do the counting */

  i = 0;
  i1 = 0;
  
  n = 0;
  
  while(i < n_points)
    {
      i = i1;
      
      while(i < n_points && t[i] < t0[n])
	i++;
      i1 = i;
      
      k = 0;
      while(i < n_points && t[i] < t1[n])
	{
	  i++;
	  k++;
	}
      
      c[n] = k;
      n++;
      if(n >= n_intervals)
	break;
      
    }
  
}
