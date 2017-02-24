/* intervalSplit_c.c C implementation of the intervalSplit function.
   inputs are:
   1: array of timestamps for tsd to average 
   2: array of timestamps for starting times of intervals
   3: array of timestamps for end times of intervals

   outputs are:
   1:first index of tsd points included in i-th interval
   2:last index of tsd points included in i-th interval
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
  double *t, *t0, *t1;
  
  int i, i1, n;
  double *ix_first, *ix_last;
  
  

  if(nrhs != 3)
    mexErrMsgTxt("Requires three inputs");
  
  if(nlhs != 2)
    mexErrMsgTxt("Requires two outputs");
  
  
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
  

  /* sanity check on the intervals */

  for(i = 0; i < n_intervals; i++)
    if(t1[i] < t0[i])
      mexErrMsgTxt("t1 must be greater than t0");
  

  /* prepare output */
  
  plhs[0] = mxCreateDoubleMatrix(n_intervals, 1, mxREAL);
  ix_first = mxGetPr(plhs[0]);
  
  plhs[1] = mxCreateDoubleMatrix(n_intervals, 1, mxREAL);
  ix_last = mxGetPr(plhs[1]);
  

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
      ix_first[n] = i + 1;
      
 
      while(i < n_points && t[i] < t1[n])
	{
	  i++;

	}
      ix_last[n] = i;
      

      n++;
      if(n >= n_intervals)
	break;
      
    }
  


  
}
