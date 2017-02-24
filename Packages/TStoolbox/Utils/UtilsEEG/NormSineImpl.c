/* routine to compute the specgram of a  */

#include <matrix.h>
#include <mex.h>
#include <math.h>
#include <string.h>



#ifndef PI
#define PI M_PI
#endif







#define MATLAB

#ifdef MATLAB
#define REAL(i, j) ( (sp_r)[(i) + ixmax * (j)] )
#define COMPLEX(i, j)  ( (sp_i)[(i) + ixmax * (j)] )

#endif


/* [ns, phase, s1, s2] = NormSineImpl(ds) */
/* take an oscillating function, such as the output of a theta filter,
   and normalizes it piecewise so that the peaks/troughs are 1/-1

   INPUTS:
   ds = a row/column vector with the oscillating function
   OUTPUTS:
   ns = the normalized function
   s1, s2, the first and last point (1-based) in ds corresponding to
   points in ns

   only for use within wrapper function NormSine
*/




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{
  

  double *ds;
  int n;
  int *peaks,  *troughs;
  
  int ipeak, itrough, npeaks, ntroughs, peak_block, peak_max, trough_max;
  
  double *ns, *phase;
  int i, j;
  int s1, s2;
  
  


  if(nrhs != 1)
    mexErrMsgTxt("Requires one input");
  
  if(nlhs != 4)
    mexErrMsgTxt("Requires  four outputs");


  /* unpack input */

  
  /* input 0 is ds */

  if(mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("ds must be a row or column vector");
  
  n= mxGetM(prhs[0]) * mxGetN(prhs[0]);


  ds = mxGetPr(prhs[0]);
  
  /* allocate workspace */

  ns = (double *)mxCalloc(n, sizeof(double));
  phase = (double *)mxCalloc(n, sizeof(double));

  peak_block = n/10;
  
  peaks = (int *)mxCalloc(peak_block, sizeof(int));
  troughs = (int *)mxCalloc(peak_block, sizeof(int));
  peak_max = peak_block;
  trough_max = peak_block;
  

  /* find peak/troughs */

  ipeak = -1;
  itrough = -1;
  

  for(i = 1; i < (n-1); i++)
    {
      if((ds[i]-ds[i-1]) > 0 && (ds[i+1]-ds[i]) < 0)
	{
	  ipeak++;
	  if(ipeak >= peak_max)
	    {
	      peak_max += peak_block;
	      mxRealloc(peaks, sizeof(int) * peak_max);
	    }
	  peaks[ipeak] = i;

	}
      
      if((ds[i]-ds[i-1]) < 0 && (ds[i+1]-ds[i]) > 0)
	{
	  itrough++;
	  if(itrough >= trough_max)
	    {
	      trough_max += peak_block;
	      mxRealloc(troughs, sizeof(int) * trough_max);
	    }
	  troughs[itrough] = i;
	  
	}
      
    }
      

  npeaks = ipeak+1;
  ntroughs = itrough+1;
  ipeak = 0;
  itrough = 0;
  
  /* align peaks/troughs */
  
  while((ipeak+1) < npeaks && peaks[ipeak+1] < troughs[0])
    ipeak++;
  
  while(itrough < ntroughs && troughs[itrough] < peaks[ipeak])
    itrough++;
  

  s1 =peaks[ipeak] + 1;
  j = 0;
  
  while(ipeak < npeaks)
    {
      double dp, dt;
      dp = ds[peaks[ipeak]];
      dt = ds[troughs[itrough]];

      for(i = peaks[ipeak]; i <troughs[itrough]; i++)
	{
	  ns[j] = 2 * (ds[i] - dt) / (dp-dt) - 1;
	  phase[j] = acos(ns[j]);
	  
	  j++;
	  
	}
      
      ipeak++;
      if(ipeak >=npeaks)
	{
	  s2 = troughs[itrough];
	  break;
	}
      
      dp = ds[peaks[ipeak]];
      
      for(i = troughs[itrough]; i < peaks[ipeak]; i++)
	{
	  ns[j] = 2 * (ds[i] - dt) / (dp-dt) - 1;
	  phase[j] = 2 * PI - acos(ns[j]);
	  
	  j++;
	}
	  
      itrough++;
      
      if(itrough >= ntroughs)
	{
	  s2 = peaks[ipeak];
	  break;
	}
      
    }
  

  


  /* prepare output */

  plhs[0] = mxCreateDoubleMatrix(1, j, mxREAL);
  memcpy(mxGetPr(plhs[0]), ns, j*sizeof(double));
  plhs[1] = mxCreateDoubleMatrix(1, j, mxREAL);
  memcpy(mxGetPr(plhs[1]), phase, j*sizeof(double));
  

  plhs[2] = mxCreateScalarDouble((double)s1);
  plhs[3] = mxCreateScalarDouble((double)s2);
  

  mxFree((void *)ns);
  mxFree((void *)peaks);
  mxFree((void *)troughs);
  



  
}
