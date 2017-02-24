/* [sPF, c_idx, pf_idx] = SplitPF(PF, n_bin_bound, bidir)
   MEX file
   Splits the place field histograms based on the criterion used in
   Mehta, PNAS,  1997.
   Place field boundaries are placed where there are two consecutive
   bins with f.r. less than 10% than peak.

   INPUTS:
   PF: a n_cells X n_bins matrix containing the place field histograms 
   n_bin_bound: number of bins to use for the boundary criterion
   bidi: flag indicating if place fields are "bidirectional". If set,
   the first n_bins/2 bins are the p.f. in one directions, second half
   is the other direction
   OUTPUTS:
   sPF: a n_pf x n_bins matrix containing separated histograms for
   each p.f.
   c_idx: a n_pf X 1 vector with indices indicating to which cell each
   p.f. belongs to
   pf_idx: n_pf X 1 vector with progressive indices for p.f. for each
   cell.

   batta 2002, initially written for CRAM p.f. analysis
*/


#include "mex.h"
#include <math.h>
#include <matrix.h>

inline double computeMax(double *x, int l)
{
  double norm;
  int i = 0, j = 0;
  while(isnan(x[j])) j++;
  norm = x[j];
  

  for(i = j; i < l; i++)
    if(!isnan(x[i]) && x[i] > norm)
      norm = x[i];
  
  
  return norm;
}



void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{

  
  double **pf, *pf_ptr, *c_ptr, *i_ptr;
  int n_bin_bound;
  int bidir = 0;
  double ***spf;
  
  int *bd_init, *bd_end;
  

  int *n_pf;
  int n_pf_tot = 0;
  
  int n_cells, n_bins;
  int i, j, k, l;
  int subs[2];
  
  
  
  /* check number of arguments: expects 4 inputs, 1 or 2 outputs */
  if (nINP != 2 && nINP != 3)
    mexErrMsgTxt("Call with PF, n_bin_bound,  bidir (optional) as inputs.");
  if (nOUT != 3)
    mexErrMsgTxt("Requires three outputs.");

  /* check validity of inputs */

  /* input 0 is a matrix*/ 
  n_cells = mxGetM(pINP[0]);
  n_bins = mxGetN(pINP[0]);
  
  /* input 1 is a scalar */
  if(mxGetM(pINP[1]) * mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("n_bin_bound must be scalar");
  
  /* input 2 is a scalar */
  if(nINP == 3)
    if(mxGetM(pINP[2]) * mxGetN(pINP[2]) != 1)
      mexErrMsgTxt("bidir must be scalar (flag)");
  

  /* unpack inputs */
  pf = (double **)mxCalloc(n_cells, sizeof(double *));
  for(i = 0; i < n_cells; i++)
    pf[i] = (double *)mxCalloc(n_bins, sizeof(double));
  
  pf_ptr = mxGetPr(pINP[0]);
  
  for(i = 0; i < n_cells; i++)
    for(j = 0; j < n_bins; j++)
      {
	subs[0] = i;
	subs[1] = j;
	pf[i][j] = 
	  pf_ptr[mxCalcSingleSubscript(pINP[0], 2, subs)];
      }
  

  n_bin_bound = mxGetScalar(pINP[1]);
  
  if(nINP == 3)
    bidir = (int)mxGetScalar(pINP[2]);
  


  /* prepare intermediatae variables */

  spf = (double ***)mxCalloc(n_cells, sizeof(double **));
  n_pf = (int *)mxCalloc(n_cells, sizeof(int));
  
  bd_init = (int *)mxCalloc(n_bins+1, sizeof(int));
  bd_end = (int *)mxCalloc(n_bins+1, sizeof(int));
  

  /* computation loop */

  for(i = 0; i< n_cells ; i++)
    {
      enum {IN, BOUNDARY, OUT
      } status;
      
      int n_bd = 0, n_bd_end = 0;
      int bd_cand = 0;
      
      double *p = pf[i];
      
      double thr;
      for(j = 0; j < n_bins+1; j++)
	{
	  bd_init[j] = 0;
	  bd_end[j] = 0;
	}
	      
      thr = computeMax(p, n_bins) / 10.;
      
      if(thr == 0)
	goto NO_PLACE_FIELD;
      
	       
      if(p[0] > thr)
	{
	  bd_init[0] = 1;
	  status = IN;
	}
      else
	{
	  status = OUT;
	}
      

      for(j = 0; j < n_bins; j++)
	{
	  switch(status)
	    {
	    case IN:
	      if(p[j] < thr)
		{
		  bd_cand = j;
		  status = BOUNDARY;
		}
	      break;

	    case OUT:
	      if(p[j] > thr)
		{
		  bd_init[j] = 1;
		  status = IN;
		}
	      break;
	      
	    case BOUNDARY:
	      if(p[j] > thr)
		status = IN;
	      else
		{
		  if(j-bd_cand+1==n_bin_bound)
		    {
		      bd_end[bd_cand] = 1;
		      status = OUT;
		    }
		}
	    }
	  if(bidir && j == (n_bins/2-1))
	    {
	      if(status==IN)
		bd_end[j+1] = 1;
	      if(status==BOUNDARY)
		bd_end[bd_cand] = 1;
	      status= OUT;
	    }
	}

      if(status==IN)
	bd_end[n_bins] = 1;
      if(status == BOUNDARY)
	bd_end[bd_cand] = 1;
      

      n_bd = 0;
      n_bd_end = 0;
      
      for(j = 0; j < (n_bins+1); j++)
	{
	  if(bd_init[j])
	    n_bd++;
	  if(bd_end[j])
	    n_bd_end++;
	}
      
      if(n_bd != n_bd_end)
	mexErrMsgTxt("Boundary mismatch, something weird going on");
      
      
      /* get rid of p.f. of same length or shorter than n_bin_bound */

      for(j = 0; j < (n_bins+1); j++)
	{
	  int j1;
	  if(bd_init[j])
	    {
	      j1 = j;
	      while(j < n_bins+1)
		{
		  j++;
		  if(bd_end[j])
		    {
		      if(j-j1 <= n_bin_bound)
			{
			  bd_init[j1]= 0;
			  bd_end[j] = 0;
			  n_bd--;
			}
		      j = j1;
		      
		      break;
		    }
		}
	    }
	}
      

    NO_PLACE_FIELD:
      
      n_pf[i] = n_bd;
      spf[i] = (double **)mxCalloc(n_pf[i], sizeof(double *));
      for(j = 0; j < n_pf[i]; j++)
	spf[i][j] = (double *)mxCalloc(n_bins, sizeof(double));
      
      {

	int j1, j2;
	j = 0;
	
	n_bd = 0;
	while(n_bd < n_pf[i])
	  {

	    while((!bd_init[j]) && (j < n_bins+1)) j++;
	    j1 = j;
	    j++;
	    
	    while((!bd_end[j]) && (j < n_bins+1)) j++;
	    j2 = j;
	    
	    for(k= j1; k < j2; k++)
	      spf[i][n_bd][k] = p[k];

	    n_bd++;
	  }
	

	
      }
      n_pf_tot += n_pf[i];
      
      
    }
  



  /* creating outputs */
	  
  pOUT[0] = mxCreateDoubleMatrix(n_pf_tot, n_bins, mxREAL);
  pf_ptr = mxGetPr(pOUT[0]);
  
  pOUT[1] = mxCreateDoubleMatrix(n_pf_tot, 1, mxREAL);
  c_ptr = mxGetPr(pOUT[1]);
	    
  pOUT[2] = mxCreateDoubleMatrix(n_pf_tot, 1, mxREAL);
  i_ptr = mxGetPr(pOUT[2]);
	    
  l = 0;
  
 
  for(i = 0; i < n_cells; i++)
    {
      for(j = 0; j < n_pf[i]; j++)
	{
	  for(k = 0; k < n_bins; k++)
	    {
	      subs[0] = l;
	      subs[1] = k;
	      pf_ptr[mxCalcSingleSubscript(pOUT[0], 2, subs)] = 
		spf[i][j][k];
	    }
    
	  c_ptr[l] = i+1;
	  i_ptr[l] = j+1;
	  l++;
	}
    }
  
	  


  /* cleanup */

  for(i = 0; i < n_cells; i++)
    mxFree((void *)pf[i]);
  mxFree((void *)pf);
  

  for(i = 0; i < n_cells; i++)
    {
      for(j = 0; j < n_pf[i]; j++)
	mxFree((void *)spf[i][j]);
      mxFree((void *)spf[i]);
    }
  mxFree((void *)spf);
  mxFree((void *)n_pf);
  
  mxFree((void *)bd_init);
  mxFree((void *)bd_end);
  

  
	

 
  

}



		 
