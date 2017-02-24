/* intervalMean_c.c C implementation of the intervalMean function.
   inputs are:
   1: array of timestamps for tsd to average 
   2: data array in a transposed format, so that time index is the
   last one (major one) 
   3: array of timestamps for starting times of intervals
   4: array of timestamps for end times of intervals
*/  

/* copyright (c) 2004 Francesco P. Battaglia */
/* This software is released under the GNU GPL */
/* www.gnu.org/copyleft/gpl.html */

#include <matrix.h>
#include <mex.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{

  int n_points, n_intervals;
  double *t, *t0, *t1, *c, *d;
  int i, i1, k, n, j;
  int ndim;
  const int *dims;
  int *tmp_dims;
  int *sub_dims; /* dimensionality for submatrix to be used for more
		    complex calculations */  
  int n_data; /* size of the data component for each time */
  
  int d_ix, o_ix;
  
  double *d_ptr, *o_ptr;
  

  if(nrhs != 4)
    mexErrMsgTxt("Requires four inputs");
  
  if(nlhs != 1)
    mexErrMsgTxt("Requires one output");
  
  
  /* unpack inputs */

  /* input 0 is t */

  if(mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("t must be a row or column vector");
  
  n_points = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  
  t = mxGetPr(prhs[0]);
  
  /* input 1 is d, the data */

  ndim = mxGetNumberOfDimensions(prhs[1]);
  dims = mxGetDimensions(prhs[1]);
  
  d = mxGetPr(prhs[1]);
  
  if(dims[ndim-1] != n_points)
    mexErrMsgTxt("Dimensionality mismatch between timestamps and data");
  
  



  /* input 2 is t0 */

  if(mxGetM(prhs[2]) != 1 && mxGetN(prhs[2]) != 1)
    mexErrMsgTxt("t0 must be a row or column vector");
  
   n_intervals = mxGetM(prhs[2]) * mxGetN(prhs[2]);

  t0 = mxGetPr(prhs[2]);
  
  /* input 3 is t1 */

  if(mxGetM(prhs[3]) != 1 && mxGetN(prhs[3]) != 1)
    mexErrMsgTxt("t1 must be a row or column vector");

  if(mxGetM(prhs[3]) != mxGetM(prhs[2]) || mxGetN(prhs[3]) != mxGetN(prhs[2]))
    mexErrMsgTxt("t1 and t0 must have same size");
  
  t1 = mxGetPr(prhs[3]);
  

  /* sanity check on the intervals */

  for(i = 0; i < n_intervals; i++)
    if(t1[i] < t0[i])
      mexErrMsgTxt("t1 must be greater than t0");
  

  /* prepare output */
  
  tmp_dims = mxCalloc(ndim, sizeof(int));
  
  n_data = 1;
  
  for(i = 0; i < ndim-1; i++)
    {
      tmp_dims[i] = dims[i];
      n_data *= dims[i];
      
    }
  
  tmp_dims[ndim-1] = n_intervals;
  




  plhs[0] = mxCreateNumericArray(ndim, tmp_dims, mxDOUBLE_CLASS, mxREAL);
  
  c = mxGetPr(plhs[0]);
  
  
  
  for(i = 0; i < ndim-1; i++)
    tmp_dims[i] = 0;
  /* do the computation */

  i = 0;
  i1 = 0;
  
  n = 0;
  
  while(i < n_points)
    {
      i = i1; /* initialize search for new interval */
      
      while(i < n_points && t[i] < t0[n]) /* search beginning of interval */
	i++;
      i1 = i;
      
      k = 0; /* number of points in the interval */
      
      tmp_dims[ndim-1] = n;
      o_ix = mxCalcSingleSubscript(plhs[0], ndim, tmp_dims);
      o_ptr = c + o_ix; /* a bit of pointer arithmetic to find the
			   point where to write */


      while(i < n_points && t[i] < t1[n])
	{
	  tmp_dims[ndim-1] = i;
	  d_ix = mxCalcSingleSubscript(prhs[1], ndim, tmp_dims);
	  d_ptr = d + d_ix; /* position yourself on the data for i-th
			       timestamps */
	  
	  for(j = 0; j < n_data; j++)
	    o_ptr[j] += d_ptr[j];
	  i++;
	  k++;

	}

      for(j = 0; j < n_data; j++) /* normalize mean */
	o_ptr[j] /= (double)k;

      n++;
      if(n >= n_intervals)
	break;
      
    }
  



  mxFree(tmp_dims);
  
}
