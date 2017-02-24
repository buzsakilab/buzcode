/* intervalFun_c.c C implementation of the intervalFun function.
   inputs are:
   1: array of timestamps for tsd to average 
   2: data array in a transposed format, so that time index is the
   last one (major one) 
   3: array of timestamps for starting times of intervals
   4: array of timestamps for end times of intervals
   5: string indicating which function to call 
*/  

/* copyright (c) 2004 Francesco P. Battaglia */
/* This software is released under the GNU GPL */
/* www.gnu.org/copyleft/gpl.html */

#include <string.h>
#include <math.h>

#include <matrix.h>
#include <mex.h>

typedef mxArray *( *intFun_t)(double *, double *, int, int *, double *, double *, int);




mxArray *doMean(double *t, double *d, int ndim, int *dims, double *t0, double *t1, int n_intervals)
{
  int i, i1, j, n, k, o_ix, d_ix;
  double *o_ptr, *d_ptr;
  
  int n_points;
  
  int *tmp_dims;
/*   int *sub_dims; */ /* dimensionality for submatrix to be used for more
		    complex calculations */  
  int n_data; /* size of the data component for each time */
  /* prepare output */
  
  mxArray *o;
  double *c;
  

  tmp_dims = mxCalloc(ndim, sizeof(int));
  
  n_data = 1;
  
  for(i = 0; i < ndim-1; i++)
    {
      tmp_dims[i] = dims[i];
      n_data *= dims[i];
      
    }
  
  tmp_dims[ndim-1] = n_intervals;
  n_points = dims[ndim-1];
  
  o = mxCreateNumericArray(ndim, tmp_dims, mxDOUBLE_CLASS, mxREAL);
  
  c = mxGetPr(o);
  
  
  
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
      o_ix = mxCalcSingleSubscript(o, ndim, tmp_dims);
      o_ptr = c + o_ix; /* a bit of pointer arithmetic to find the
			   point where to write */


      while(i < n_points && t[i] < t1[n])
	{
	  tmp_dims[ndim-1] = i;
	  d_ix = i * n_data;
	  
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

  return o;
  
}

mxArray *doMax(double *t, double *d, int ndim, int *dims, double *t0, double *t1, int n_intervals)
{
  int i, i1, j, n, k, o_ix, d_ix;
  double *o_ptr, *d_ptr;
  
  int n_points;
  
  int *tmp_dims;
/*   int *sub_dims; */ /* dimensionality for submatrix to be used for more
		    complex calculations */  
  int n_data; /* size of the data component for each time */
  /* prepare output */
  
  mxArray *o;
  double *c;
  

  tmp_dims = mxCalloc(ndim, sizeof(int));
  
  n_data = 1;
  
  for(i = 0; i < ndim-1; i++)
    {
      tmp_dims[i] = dims[i];
      n_data *= dims[i];
      
    }
  
  tmp_dims[ndim-1] = n_intervals;
  n_points = dims[ndim-1];
  
  o = mxCreateNumericArray(ndim, tmp_dims, mxDOUBLE_CLASS, mxREAL);
  
  c = mxGetPr(o);
  
  
  
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
      o_ix = mxCalcSingleSubscript(o, ndim, tmp_dims);
      o_ptr = c + o_ix; /* a bit of pointer arithmetic to find the
			   point where to write */
      d_ix = i * n_data;
      d_ptr = d + d_ix; 
      for(j = 0; j < n_data; j++)
	o_ptr[j] = d_ptr[j];
      

      while(i < n_points && t[i] < t1[n])
	{
	  d_ix = i * n_data;
	  
	  d_ptr = d + d_ix; /* position yourself on the data for i-th
			       timestamps */
	  
	  for(j = 0; j < n_data; j++)
	    if(d_ptr[j] > o_ptr[j])
	      o_ptr[j] = d_ptr[j];
	  i++;
	  k++;

	}


      n++;
      if(n >= n_intervals)
	break;
      
    }
  



  mxFree(tmp_dims);

  return o;
  
}

mxArray *doMin(double *t, double *d, int ndim, int *dims, double *t0, double *t1, int n_intervals)
{
  int i, i1, j, n, k, o_ix, d_ix;
  double *o_ptr, *d_ptr;
  
  int n_points;
  
  int *tmp_dims;
/*   int *sub_dims; */ /* dimensionality for submatrix to be used for more
		    complex calculations */  
  int n_data; /* size of the data component for each time */
  /* prepare output */
  
  mxArray *o;
  double *c;
  

  tmp_dims = mxCalloc(ndim, sizeof(int));
  
  n_data = 1;
  
  for(i = 0; i < ndim-1; i++)
    {
      tmp_dims[i] = dims[i];
      n_data *= dims[i];
      
    }
  
  tmp_dims[ndim-1] = n_intervals;
  n_points = dims[ndim-1];
  
  o = mxCreateNumericArray(ndim, tmp_dims, mxDOUBLE_CLASS, mxREAL);
  
  c = mxGetPr(o);
  
  
  
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
      o_ix = mxCalcSingleSubscript(o, ndim, tmp_dims);
      o_ptr = c + o_ix; /* a bit of pointer arithmetic to find the
			   point where to write */
      d_ix = i * n_data;
      d_ptr = d + d_ix; 
      for(j = 0; j < n_data; j++)
	o_ptr[j] = d_ptr[j];
      

      while(i < n_points && t[i] < t1[n])
	{
	  d_ix = i * n_data;
	  
	  d_ptr = d + d_ix; /* position yourself on the data for i-th
			       timestamps */
	  
	  for(j = 0; j < n_data; j++)
	    if(d_ptr[j] < o_ptr[j])
	      o_ptr[j] = d_ptr[j];
	  i++;
	  k++;

	}


      n++;
      if(n >= n_intervals)
	break;
      
    }
  



  mxFree(tmp_dims);

  return o;
  
}




mxArray *doVar(double *t, double *d, int ndim, int *dims, double *t0, double *t1, int n_intervals)
{
  int i, i1, j, n, k, o_ix, d_ix;
  double *o_ptr, *d_ptr;
  
  int n_points;
  
  int *tmp_dims;
/*   int *sub_dims; */ /* dimensionality for submatrix to be used for more
		    complex calculations */  
  int n_data; /* size of the data component for each time */
  /* prepare output */
  
  mxArray *o;
  double *c;
  
  double *m, *m2;
  

  tmp_dims = mxCalloc(ndim, sizeof(int));
  
  n_data = 1;
  
  for(i = 0; i < ndim-1; i++)
    {
      tmp_dims[i] = dims[i];
      n_data *= dims[i];
      
    }
  
  tmp_dims[ndim-1] = n_intervals;
  n_points = dims[ndim-1];
  

  m = (double *)mxCalloc(n_data, sizeof(double));
  m2 = (double *)mxCalloc(n_data, sizeof(double));
  
  o = mxCreateNumericArray(ndim, tmp_dims, mxDOUBLE_CLASS, mxREAL);
  
  c = mxGetPr(o);
  
  
  
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
      o_ix = mxCalcSingleSubscript(o, ndim, tmp_dims);
      o_ptr = c + o_ix; /* a bit of pointer arithmetic to find the
			   point where to write */
      for(j = 0; j < n_data; j++)
	{
	  m[j] = 0.;
	  m2[j] = 0.;
	}
      


      while(i < n_points && t[i] < t1[n])
	{
	  tmp_dims[ndim-1] = i;
	  d_ix = i * n_data;
	  
	  d_ptr = d + d_ix; /* position yourself on the data for i-th
			       timestamps */
	  
	  for(j = 0; j < n_data; j++)
	    {
	      m[j] += d_ptr[j];
	      m2[j] += d_ptr[j] * d_ptr[j];
	    }
	  
	  i++;
	  k++;

	}

      for(j = 0; j < n_data; j++) /* normalize mean */
	o_ptr[j] /= m2[j] / k - (m[j] / k) * (m[j] / k);
	  

      n++;
      if(n >= n_intervals)
	break;
      
    }
  



  mxFree(tmp_dims);
  mxFree(m);
  mxFree(m2);
  

  return o;
  
}



mxArray *doStd(double *t, double *d, int ndim, int *dims, double *t0, double *t1, int n_intervals)
{
  int i, i1, j, n, k, o_ix, d_ix;
  double *o_ptr, *d_ptr;
  
  int n_points;
  
  int *tmp_dims;
/*   int *sub_dims; */ /* dimensionality for submatrix to be used for more
		    complex calculations */  
  int n_data; /* size of the data component for each time */
  /* prepare output */
  
  mxArray *o;
  double *c;
  
  double *m, *m2;
  

  tmp_dims = mxCalloc(ndim, sizeof(int));
  
  n_data = 1;
  
  for(i = 0; i < ndim-1; i++)
    {
      tmp_dims[i] = dims[i];
      n_data *= dims[i];
      
    }
  
  tmp_dims[ndim-1] = n_intervals;
  n_points = dims[ndim-1];
  

  m = (double *)mxCalloc(n_data, sizeof(double));
  m2 = (double *)mxCalloc(n_data, sizeof(double));
  
  o = mxCreateNumericArray(ndim, tmp_dims, mxDOUBLE_CLASS, mxREAL);
  
  c = mxGetPr(o);
  
  
  
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
      o_ix = mxCalcSingleSubscript(o, ndim, tmp_dims);
      o_ptr = c + o_ix; /* a bit of pointer arithmetic to find the
			   point where to write */
      for(j = 0; j < n_data; j++)
	{
	  m[j] = 0.;
	  m2[j] = 0.;
	}
      


      while(i < n_points && t[i] < t1[n])
	{
	  tmp_dims[ndim-1] = i;
	  d_ix = i * n_data;
	  
	  d_ptr = d + d_ix; /* position yourself on the data for i-th
			       timestamps */
	  
	  for(j = 0; j < n_data; j++)
	    {
	      m[j] += d_ptr[j];
	      m2[j] += d_ptr[j] * d_ptr[j];
	    }
	  
	  i++;
	  k++;

	}

      for(j = 0; j < n_data; j++) /* normalize mean */
	o_ptr[j] /= sqrt(m2[j] / k - (m[j] / k) * (m[j] / k));
	  

      n++;
      if(n >= n_intervals)
	break;
      
    }
  



  mxFree(tmp_dims);
  mxFree(m);
  mxFree(m2);

  return o;
  
}






void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{

  int n_points, n_intervals;
  double *t, *t0, *t1,  *d;
  int i;
  int ndim;
  const int *dims;

  char func_str[80];
  intFun_t func;
  

  if(nrhs != 5)
    mexErrMsgTxt("Requires five inputs");
  
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
  
  /* input 4 is func */

  mxGetString(prhs[4], func_str, 79);
  

  /* sanity check on the intervals */

  for(i = 0; i < n_intervals; i++)
    if(t1[i] < t0[i])
      mexErrMsgTxt("t1 must be greater than t0");
  

  if(!strcmp(func_str, "mean")) {
    func = doMean;
  }
  else if(!strcmp(func_str, "var")) {
    func = doVar;
  }
  else if(!strcmp(func_str, "std")) {
    func = doStd;
  }
  else if(!strcmp(func_str, "min")) {
    func = doMin;
  }
  else if(!strcmp(func_str, "max")) {
    func = doMax;
  }
  else 
    {
      mexPrintf("#%s#\n", func_str);
      
      mexErrMsgTxt("Unknown option");
    }
  

      
  


  plhs[0] = func(t, d, ndim, dims, t0, t1, n_intervals);
    

  







  /***********************************************************************/


}
