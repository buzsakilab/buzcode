


#include <stdlib.h>
#include <math.h>
#include <mex.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf_gamma.h>







void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{



  int N, k, Ncomb;
  gsl_combination *c;
  size_t *d;
  
  int i, j;
  

  int subs[2];
  double *ptr;
  

  if(nrhs != 2)
    mexErrMsgTxt("2 inputs required");

  if(nlhs != 1)
    mexErrMsgTxt("Requires one output.");
  

  if(mxGetM(prhs[0]) * mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("N must be scalar");
  
  N = mxGetScalar(prhs[0]);
  
  if(mxGetM(prhs[1]) * mxGetN(prhs[1]) != 1)
    mexErrMsgTxt("k must be scalar");
  
  k = mxGetScalar(prhs[1]);
  

  Ncomb = (int) (gsl_sf_fact(N) / (gsl_sf_fact(k) * gsl_sf_fact(N-k))); 
/*    mexPrintf("%g, %g, %g, Ncomb = %d\n", gsl_sf_fact(N),  */
/*  	    gsl_sf_fact(k), gsl_sf_fact(N-k), Ncomb); */
  

  
  
  
  c = gsl_combination_calloc (N, k);
  gsl_combination_init_first(c);
  



  plhs[0] = mxCreateDoubleMatrix(Ncomb, k, mxREAL);
  ptr = mxGetPr(plhs[0]);
  i = 0;
  
  do
    {
      d = gsl_combination_data(c);
      for(j = 0; j < k; j++)
	{

	  subs[0] = i;
	  subs[1] = j;
	  ptr[mxCalcSingleSubscript(plhs[0], 2, subs)] = 
	    (double)d[j] + 1.;
	  
	}
     /*   gsl_combination_fprintf (stdout, c, " %u"); */
      /*  mexPrintf("\na\n"); */

      i++;
      
    }
  while (gsl_combination_next(c) == GSL_SUCCESS);

  gsl_combination_free(c);

}
