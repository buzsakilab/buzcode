#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <matrix.h>


inline double l2norm(double *x, int l)
{
  double norm = 0.;
  int i = 0;
  
  for(i = 0; i < l; i++)
    norm += x[i] * x[i];
  
  return sqrt(norm);
}

  


void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  
  double **x, *xo,  *z = NULL, *c, *s, **zp, **zpold, *maxes, *mins;
  double nm;
  
  int n_samples, n, D, K, i, j, iterating, k, iter, *ns, max_iter = 100;
  int subs[2];
  
  /* check number of arguments: expects 2 or 3 inputs,  2 outputs */
  if (nINP != 2 && nINP != 3)
    mexErrMsgTxt("Call with 2 or 3 inputs");
  if (nOUT != 2)
    mexErrMsgTxt("Requires two outputs.");

  /* check validity of inputs */


   

  n_samples = mxGetM(pINP[0]);
  D = mxGetN(pINP[0]);
  K = (int)mxGetScalar(pINP[1]);
  
  if (mxGetM(pINP[1]) != 1 || mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("K must be scalar");
  if (nINP > 2)
    if(mxGetM(pINP[2]) != K || mxGetN(pINP[2]) != D)
      mexErrMsgTxt("x0 must be K x D");


  /* unpack inputs */
  xo = mxGetPr(pINP[0]);


  x = mxCalloc(n_samples, sizeof(double *));

   

   
  for(i = 0; i < n_samples; i++)
    {
      x[i] = mxCalloc(D, sizeof(double));
      
      if(x[i]  == NULL)
	
	{

	  
	  mexErrMsgTxt("problem allocating memory");
	}
      for(j = 0; j < D; j++)
	x[i][j] = 0.;

      
	  
    }



  for(i = 0; i < n_samples; i++)
    {

      for(j = 0; j < D; j++)
	{
	  int  idx;
	  subs[0] = i;
	  subs[1] = j;


  	  idx = mxCalcSingleSubscript(pINP[0], 2, subs); 

	  x[i][j] = xo[idx];
	  
	  
	}
    }
  
    

  
  



  if(nINP > 2)
    z = mxGetPr(pINP[2]);
  
  /* prepare inputs */

  pOUT[0] = mxCreateDoubleMatrix(K, D, mxREAL);
  pOUT[1] = mxCreateDoubleMatrix(1, n_samples, mxREAL);
  
  c = mxGetPr(pOUT[0]);
  s = mxGetPr(pOUT[1]);
  
  zp = mxCalloc(K, sizeof(double *));
  zpold = mxCalloc(K, sizeof(double *));
  
  for(i = 0; i < K; i++)
    {
      zp[i] = mxCalloc(D, sizeof(double));
      zpold[i] = mxCalloc(D, sizeof(double));
    }
  
  ns = mxCalloc(K, sizeof(int));
  
  /* prepare initial condition (if not provided) */




  maxes = mxCalloc(D, sizeof(double));
  mins = mxCalloc(D, sizeof(double));
  
  {

      int ii, l;
      


      for(ii=0; ii < D; ii++)
	{
	  double el;
	  el = x[0][ii];
	  
	  
	  maxes[ii] = el;
	  mins[ii] = el;
	}
      

      for(l = 0; l < n_samples; l++)
	for(ii=0; ii < D; ii++)
	  {

	    double el;
	    
	    el = x[l][ii];
	    
	    if(el > maxes[ii])
	      maxes[ii] = el;
	    if(el < mins[ii])
	      mins[ii] = el;
	    
	  }
  }


  if(nINP == 2)
    {
      double maxes[D], mins[D];
      int ii, l;
      


      
      for(l = 0; l < K; l++)
	for(ii=0; ii < D; ii++)
	  {
	    int subs[] = {l, ii};

	    zp[l][ii] = 
	      (maxes[ii] - mins[ii]) * drand48() + mins[ii];
	  }
      
    }
  else
    {
      

      for(i = 0; i < K; i++)
	for(j = 0; j < D; j++)
	  {
	    int subs[] = {i,j};
	    zp[i][j] = z[mxCalcSingleSubscript(pINP[2], 2, subs)];
	  }
  
    }
  
      
  for(i=0; i < n_samples; i++)
	s[i] = 0.;
      
      iterating = 1;
      k = 1;
      iter = 0;
      n = n_samples;
      

      while(iterating)
	{
	  int l, ii, jj;
	  int ns[K];
	 
	  double dist[K], dt[D];
	  
	  iter++;
	  
	  for(l = 0; l < n; l++)
	    {
	      double mmin;
	      int ind;

	      
	      for(ii = 0; ii < K; ii++)
		{
		  for(jj = 0; jj < D; jj++)
		    {
		      dt[jj] = x[l][jj] -
		      zp[ii][jj];
		    }
		  
		  dist[ii] = l2norm(dt, D);
		}

	      
	      mmin = dist[0];
	      ind = 0;
	      
	      for(i = 0; i < K; i++)
		if(dist[i] < mmin) {
		  mmin = dist[i];
		  ind = i;
		}
	      s[l] = ind;
	    }
	  
	  for(i = 0; i <K; i++)
	    {
	      memcpy(zpold[i], zp[i], sizeof(double) * D);
	      bzero(zp[i], sizeof(double) * D);
	    }
	  
	      

	  
	  for(i = 0; i < K; i++)
	    ns[i] = 0;
	  
	  for(i = 0 ; i < n; i++)
	    {
	      int is, j;
	      is = s[i];
	      ns[is]++;
	      for(j = 0; j < D; j++)
		zp[is][j] += x[i][j];
	    }
	  
		  
	  for(i = 0; i < K; i++)
	    {
	      
	    if(ns[i])
	      for(j = 0; j < D; j++)
		zp[i][j] /= ns[i];
	    else
	      {
		;
		
		
		
	      }
	    }
	  

	  nm = 0.;
	  for(i = 0; i < K; i++)
	    {
	      for(j = 0; j < D; j++)
		zpold[i][j] -= zp[i][j];
	      
	      nm += l2norm(zpold[i], D);
	    }
	  
	  if (nm == 0 || iter >= max_iter) 
	    iterating  = 0;
/*  	  mexPrintf("iter = %d\tnm = %g\n", iter, nm); */
	}
      
      
      /* pack output */

      for(i = 0; i < n_samples; i++)
	s[i]++;
      
      
      for(i = 0; i < K; i++)
	for(j = 0; j < D; j++)
	  {
	    int subs[] = {i, j};
	    c[mxCalcSingleSubscript(pOUT[0], 2, subs)] = zp[i][j];
	  }
      
      
}

  
	   

	    



	      
      

      

