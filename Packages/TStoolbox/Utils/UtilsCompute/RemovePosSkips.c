/* gp = RemovePosSkips(x, y, t, skip_thresh_x, skip_thresh_y, skip_thresh_t) 

   Remove tracker jumping points

   INPUTS:
   x, y, t: tracker x, y, and time
   skip_thresh_x, skip_thresh_y: maximum jumping values allowed
   skip_thresh_t: only act when the time step is less than skip_thresh_t

   OUTPUT:
   gp: the indices of the "good points"
*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

#define INP_X_IDX 0
#define INP_Y_IDX 1
#define INP_T_IDX 2
#define INP_SKIP_THRESH_X 3
#define INP_SKIP_THRESH_Y 4
#define INP_SKIP_THRESH_T 5

#define OUT_GP_IDX 0

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  
  double *x, *y, *t;
  double *gp, *gpp;
  
  int n_points, n_good;
  int i, j;
  int n_mean = 3;
  
  double st_x, st_y, st_t;
  





  /* check number of arguments: expects 6 inputs, 1  outputs */
  if (nINP != 6)
    mexErrMsgTxt("Call with x, y, t, skip_thresh_x, skip_thresh_y, skip_thresh_t as inputs.");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one outputs.");

  x = mxGetPr(pINP[INP_X_IDX]);
  y = mxGetPr(pINP[INP_Y_IDX]);
  t = mxGetPr(pINP[INP_T_IDX]);
  
  st_x = mxGetScalar(pINP[INP_SKIP_THRESH_X]);
  st_y = mxGetScalar(pINP[INP_SKIP_THRESH_Y]);
  st_t = mxGetScalar(pINP[INP_SKIP_THRESH_T]);
  
  
  n_points = mxGetM(pINP[INP_X_IDX]) * mxGetN(pINP[INP_X_IDX]);
  gp = mxCalloc(n_points, sizeof(double));
  
  gp[0] = 0;
  n_good = 1;
  
  i = 0;


  
/*   while(i < n_points-1) */
/*     { */
/*       j = i+1; */
      
/*       while(j < (n_points-1) && (((x[j]-x[i]) > st_x) || ((y[j]-y[i]) > st_y)) /\* && */
/* 									    (t[j]-t[j-1] < st_t) *\/) */
/* 	j++; */
/*       i = j; */
/*       gp[n_good] = i; */
/*       n_good++; */
/*     } */
  
  for(i = 1; i < n_points; i++)
    {
      double mx, my;
      int i1, i2;
      i1 = i - n_mean;
      i1 = (i1 >= 0 ? i1 : 0);
      i2 = i + n_mean;
      
      i2 = (i2 < n_points ? i2 : 0);
      mx = 0.;
      my = 0.;
      
      for(j = i1; j < i2; j++)
	{
	  mx += x[j];
	  my += y[j];
	}
      
      mx /= (i2-i1);
      my /= (i2-i1);
      
      if( (t[i]-t[i-1] > st_t) || (fabs(x[i]-mx) < st_x && fabs(y[i]-my) < st_y))
	{
	  gp[n_good] = i;
	  n_good++;
	  
	}
    }
  




  pOUT[OUT_GP_IDX] = mxCreateDoubleMatrix(n_good, 1, mxREAL);
  gpp = mxGetPr(pOUT[OUT_GP_IDX]);
  
  for(i = 0; i < n_good; i++)
    gpp[i] = gp[i] + 1.;
  
  mxFree(gp);
  
}





	      
