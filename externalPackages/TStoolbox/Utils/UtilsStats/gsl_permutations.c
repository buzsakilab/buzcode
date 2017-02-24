#include <stdlib.h>
#include <math.h>
#include <mex.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_gamma.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{


  int i, j;
  int N, Nperm;
  int *d;
  
  int subs[2];
  double *ptr;
  gsl_permutation *c;
  

 if(nrhs != 1)
    mexErrMsgTxt("1 input required");

  if(nlhs != 1)
    mexErrMsgTxt("Requires one output.");
  
  if(mxGetM(prhs[0]) * mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("N must be scalar");
  
  N = mxGetScalar(prhs[0]);

  Nperm = (int) (gsl_sf_fact(N));
  
  c = gsl_permutation_calloc (N);
  gsl_permutation_init(c);



  plhs[0] = mxCreateDoubleMatrix(Nperm, N, mxREAL);
  ptr = mxGetPr(plhs[0]);

  for(i = 0; i < Nperm; i++)
    {
      d = gsl_permutation_data(c);
      for(j = 0; j < N; j++)
	{

	  subs[0] = i;
	  subs[1] = j;
	  ptr[mxCalcSingleSubscript(plhs[0], 2, subs)] = 
	    (double)d[j] + 1.;
	}
      
      gsl_permutation_next(c);
    }
  


  gsl_permutation_free(c);
  
}

	   
