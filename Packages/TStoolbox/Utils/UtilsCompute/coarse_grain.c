

#include <string.h>

#include <math.h> 
#include <mex.h>
#include <matrix.h>





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
		 *prhs[])
{

  double *ptr;

  int mrows, ncols, ndata, ntimes, nbins, i, b, c;
  double *data, *times, *bincenters, *boundaries, *cg, d;
  
  
   

  if(nrhs != 3)
    mexErrMsgTxt("Three inputs required");

   if(nlhs != 1)
    mexErrMsgTxt("Requires one output");
  



  /* input 0 is the datavector */
  
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  

  if(!(mrows == 1 || ncols == 1))
    mexErrMsgTxt("Input 1 must be a row or column vector");
  if(mrows == 1)
    ndata = ncols;
  else
    ndata = mrows;
  
  data = mxGetPr(prhs[0]);

   /* input 1 is the times*/

  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  
  if(!(mrows == 1 || ncols == 1))
    mexErrMsgTxt("Input 2 must be a row or column vector");
  if(mrows == 1)
    ntimes = ncols;
  else 
    ntimes = mrows;
  
  if(ntimes != ndata)
    mexErrMsgTxt("Inputs 1 and 2 must have the same length");
  
  times = mxGetPr(prhs[1]);
  
  
  
 


/* input 2 is the time bins centers vector */

  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  
  if(!(mrows == 1 ||  ncols == 1))
    mexErrMsgTxt("Input 2 must be a row or column vector");
  
  if(mrows == 1)
    nbins = ncols;
  else 
    nbins = mrows;
  

  bincenters = mxGetPr(prhs[2]);
  
  boundaries = (double *)mxCalloc(nbins-1, sizeof(double));
  for(i = 0; i < nbins-1;i++)
    boundaries[i] = (bincenters[i] + bincenters[i+1]) / 2.;

  
  d = 0.;
  
  b = 0;
  i = 0;
  c = 0;
  
  cg = (double *)mxCalloc(nbins, sizeof(double));
  

  while(i < ndata)
    {
      if(times[i] > boundaries[b])
	{
	  /*  if(c > 0) */
	    cg[b] = d / c;
	  /*	  else
		  cg[b] = 0.; */
	  
	  d = 0.;
	  c = 0;
	  while(b < nbins-1 && times[i] > boundaries[b])
	    {
	      b++;
	      cg[b] = mxGetNaN();
	    }
	  
	  if(b >=nbins)
	    break;
	  
	}
      d += data[i];
      c++;
      i++;
    }
/*    if(c > 0) */
    cg[nbins-1] = d/c;
  /*  else
      cg[nbins-1] = 0.; */
  
  plhs[0] = mxCreateDoubleMatrix(1, nbins, mxREAL);
  
  ptr = mxGetPr(plhs[0]);
  memcpy(ptr, cg, nbins * sizeof(double));
  
  


      





}



