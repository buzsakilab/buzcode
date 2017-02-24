/* "spherical kmeans, as suggested for example in Dhillon, 2000. 
INPUTS: 
X:  a n_samples X D matrix where each row is a data point vector. Each
    data point vector is expected to be normalized to unit L^2 norm 
K:  the number of clusters
Zo: a K X D matrix with the initial condition for the cluster centroids

OUTPUTS:
C: a K x D matrix of the K centroids 
S: a 1 x n_samples array of cluster belonging


*/

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

inline double cosine_distance(double *x, double *y, int l)
{ 
  double  norm = 0.;
  int i = 0;

  
  
  for(i = 0; i < l; i++)
    norm += ((x[i]*y[i]));
  
  
  
 return (1. -norm);
  
  
}

inline double neg_overlap(double *x, double *y, int l)
{ 
  double mx = 0., my = 0., norm = 0.;
  int i = 0;
  
  for(i = 0; i < l; i++)
    {
      mx += x[i];
      my += y[i];
    }

  mx /= l;
  my /= l;
  
  
  for(i = 0; i < l; i++)
    norm += (x[i]-mx)*(y[i]-my);
  
  norm /= sqrt((mx - mx * mx) * (my - my*my));
  
 return ( -norm);
  
  
}





void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  
  double **x, *xo,  *z = NULL, *c, *s, **zp, **zpold;
  double nm;
  
  int n_samples, n, D, K, i, j, iterating, k, iter,  max_iter = 10000;
  int subs[2];
  
  /* check number of arguments: expects 2 or 3 inputs,  2 outputs */
  if ( nINP != 3  && nINP != 4)
    mexErrMsgTxt("Call with 3 or 4 inputs");
  if (nOUT != 2 && nOUT != 3)
    mexErrMsgTxt("Requires two or three outputs.");

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
  
  if(nINP > 3) 
    max_iter = (int)mxGetScalar(pINP[3]);
  
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
  

  
  /* prepare initial condition (if not provided) */







      for(i = 0; i < K; i++)
	for(j = 0; j < D; j++)
	  {
	    int subs[] = {i,j};
	    zp[i][j] = z[mxCalcSingleSubscript(pINP[2], 2, subs)];
	  }
  

  
      
  for(i=0; i < n_samples; i++)
	s[i] = 0.;
      
      iterating = 1;
      k = 1;
      iter = 0;
      n = n_samples;
      

      while(iterating)
	{
	  int l, ii;

	 
	  double dist[K];
	  
	  iter++;
	  
	  for(l = 0; l < n; l++)
	    {
	      double mmin;
	      int ind;

	      
	      for(ii = 0; ii < K; ii++)
		{
		  /*  for(jj = 0; jj < D; jj++) */
/*  		    { */
/*  		      dt[jj] = x[l][jj] - */
/*  		      zp[ii][jj]; */
/*  		    } */
		  
		  dist[ii] = cosine_distance(x[l], zp[ii], D);
/*  		  dist[ii] = l2norm(dt, D); */
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
	  
	      

	  
	  
	  for(i = 0 ; i < n; i++)
	    {
	      int is, j;
	      is = s[i];

	      for(j = 0; j < D; j++)
		zp[is][j] += x[i][j];
	    }
	  
		  
	  for(i = 0; i < K; i++)
	    {
	      
	      double nz;
	      nz = l2norm(zp[i], D);
	      
	      if (nz > 0.)
		for(j = 0; j < D; j++)
		  zp[i][j] /= nz;
	      
	    
	    }
	  

	  nm = 0.;
	  for(i = 0; i < K; i++)
	    {
	      for(j = 0; j < D; j++)
		zpold[i][j] -= zp[i][j];
	      
	      nm += l2norm(zpold[i], D);
	    }
	  
	  if (nm == 0 || iter >= max_iter)
	    {
	      if( iter >= max_iter)
		mexPrintf("iteractions maxed out\n");
	      
	      iterating  = 0;
	    }
	  
/*  	  mexPrintf("iter = %d\tnm = %g\n", iter, nm); */
	}
      
      

      if(nOUT > 2)
	pOUT[2] = mxCreateScalarDouble((double) iter);
      
      /* pack output */

      for(i = 0; i < n_samples; i++)
	s[i]++;
      
      
      for(i = 0; i < K; i++)
	for(j = 0; j < D; j++)
	  {
	    int subs[] = {i, j};
	    c[mxCalcSingleSubscript(pOUT[0], 2, subs)] = zp[i][j];
	  }
      
      for(i = 0; i < n_samples; i++)
	{
	  mxFree((void *)x[i]);
	}
      
      mxFree((void *)x);
      

      for(i = 0; i < K; i++)
	{
	  mxFree((void *)zp[i]);
	  mxFree((void *)zpold[i]);
	}
      

      mxFree((void *)zp);
      mxFree((void *)zpold);
      
      

      


}

  
	   

	    



	      
      

      

