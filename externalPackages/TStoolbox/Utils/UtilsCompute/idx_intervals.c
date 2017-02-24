/* ix_c = idx_intervals(t, istart, iend) */
/* intervals described by istart iend are considered to be in growing order of istart, but not necessariliy disjoint */

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  int n;
  int n_int;
  

  double *t;
  double *istart, *iend;



  
  /* check number of arguments: expects 3 inputs, 1 output */
  if (nINP != 3)
    mexErrMsgTxt("Call with t, istart, iend as inputs.");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");



  /* unpack inputs */
  n = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  t = (double *)mxGetPr(pINP[0]);
  n_int = mxGetM(pINP[1]) * mxGetN(pINP[1]);

  if(mxGetM(pINP[2]) * mxGetN(pINP[2]) != n_int)
    mexErrMsgTxt("istart and iend must have same size");
  
  istart = (double *)mxGetPr(pINP[1]);
  iend = (double *)mxGetPr(pINP[2]);
  
  pOUT[0] = mxCreateCellMatrix(n_int, 1);
  
  /* sequential search */
  
  {
    int i, j, k;
    int ii, ll;
    
    i = 0;
    j = 0;
    
    mxArray *mxa;
    double *pr;
    
    while (i < n && j < n_int)
      {
	while(t[i] < istart[j] && i < n) i++;
	
	k = i-1;
	
	while(t[k] < iend[j] && k < n) k++;
	
	mxa = mxCreateDoubleMatrix(1, k-i+1, mxREAL);
	pr = mxGetPr(mxa);
	
	ll = 0;
	
	for (ii = i; ii <= k; ii++)
	  {
	    pr[ll] = ii;
	    ll++;
	  }
	
	mxSetCell(pOUT[0], j, mxa);
	
	j++;
      }
    
    for(ii=j+1; ii< n_int; ii++)
      {
	mxa = mxCreateDoubleMatrix(1, 0, mxREAL);
	mxSetCell(pOUT[0], ii, mxa);
      }
  }
  
	
	 
	
	
	


 
  
  

}
  
  
		 
