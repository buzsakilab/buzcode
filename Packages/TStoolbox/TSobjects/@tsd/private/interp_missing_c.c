/* interp_missing_c.c c implementation of part of the interp_missing function
   inputs are:
   1: array of timestamps for tsd to interpolate
   2: the guess on the assumed inter-event interval
   3: the guess to use to pre-allocate memory (should be smaller than previous
   outputs are:
   1: array of interpolated timestamps
   2: array of indices indicating the position of the original
   timestamps into the interpoalted array
*/  

/* copyright (c) 2004 Francesco P. Battaglia */
/* This software is released under the GNU GPL */
/* www.gnu.org/copyleft/gpl.html */

#include <string.h>
#include <math.h>

#include <matrix.h>
#include <mex.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{

  int n_points;
  double *t, *nt, *orig_ix;
  double mdt ; /* the desidered inter-event interval */
  double mdt_min; /* the inter-event interval to consider to pre-allocate memory */ 
  
  double curr_nt;
  
  int i, k, ok, n_max_nt;
  

  if(nrhs != 2 && nrhs != 3)
    mexErrMsgTxt("Requires two or three inputs");
  
  if(nlhs != 2)
    mexErrMsgTxt("Requires two output");
  
  
  /* unpack inputs */

  /* input 0 is t */

  if(mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("t must be a row or column vector");
  
  n_points = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  
  t = mxGetPr(prhs[0]);
  
  /* input 1 is d, the data */


  
  mdt = mxGetScalar(prhs[1]);
  mdt_min = mdt;
  
  if(nrhs == 3)
    mdt_min = mxGetScalar(prhs[2]);
  

  /* prepare output */
  
  n_max_nt = floor((t[n_points-1] - t[0]) / mdt_min) + 2;
  
  nt = mxCalloc(n_max_nt, sizeof(double));
  

  plhs[1] = mxCreateDoubleMatrix(n_points, 1, mxREAL);
 
  orig_ix = mxGetPr(plhs[1]);
  
 
  
  
  
  /* do the computation */

  i = 0;
  k = 0;
  ok = 0;
  
  nt[k] = t[i];
  orig_ix[ok] = i;
  curr_nt = t[i];

  i++;
  k++;
  ok++;
  
  
  while(i < n_points)
    {
      while(t[i] >= (curr_nt + 2. * mdt))
	{
	  curr_nt += mdt;
      if( k > n_max_nt-1)
          mexErrMsgTxt("k overflow 1");
      nt[k] = curr_nt;
      k++;
      }
      if( k > n_max_nt-1)
          mexErrMsgTxt("k overflow 2");
      
      nt[k] = t[i];
      curr_nt = t[i];

      if (ok > n_points-1)
        mexErrMsgTxt("ok overflow 1");
      orig_ix[ok] = i;
      i++;
      k++;
      ok++;
    }

 

  /* craft first output */
  
  plhs[0] = mxCreateDoubleMatrix(k, 1, mxREAL);
  memcpy(mxGetPr(plhs[0]), nt, k * sizeof(double));
  
  mxFree(nt);

  
}
