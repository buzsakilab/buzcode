/* "spherical kmeans, as suggested for example in Dhillon, 2000. 
 with the "refinement" variation
INPUTS: 
X:  a n_samples X D matrix where each row is a data point vector. Each
    data point vector is expected to be normalized to unit L^2 norm 
K:  the number of clusters
Zo: a K X D matrix with the initial condition for the cluster centroids
f:  the iteration factor for the FV part of the algorithm
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

inline double cosine_similarity(double *x, double *y, int l)
{ 
  double  norm = 0.;
  int i = 0;

  
  
  for(i = 0; i < l; i++)
    norm += ((x[i]*y[i]));
  
  
  
 return (norm);
  
  
}


#define MATLAB

#ifdef MATLAB

/* inline double* getRealItem(const mxArray *p, int i, int j) */
/* { */
/*   double *pr; */
/*   int nRow; */
/*   pr = mxGetPr(p); */
/*   nRow = mxGetM(p); */
/*   return pr + (i + nRow *  j); */
/* } */


inline double* getRealItem(const mxArray *p, int i, int j)
{
  double *pr;
  int subs[2];
  subs[0] = i;
  subs[1] = j;

  pr = mxGetPr(p);
  return &(pr[mxCalcSingleSubscript(p, 2, subs)]);
  
}

#endif

void updateXC(const mxArray *x, double **c, double **xc, 
	      int K, int D, int n_samples, const int *update_flags)
{
  
  int i, j, k;
  

  for(i = 0; i < K; i++)
    if(update_flags[i])
      {
	for(j = 0; j < n_samples; j++)
	  {
	    double cs = 0.;
	    for(k = 0; k < D; k++)
	      cs += c[i][k] * *(getRealItem(x, j, k));
	    xc[i][j] = cs;
	  }
      }
  
}


void updatePfromC_KM(int *p, double **xc, int K, int n_samples)
{
  int i, j;
  int mi;
  double m;
  
  for(j = 0; j < n_samples; j++)
    {
      mi = 0;
      m = xc[0][j];
      for(i = 1; i < K; i++)
	{
	  if(xc[i][j] > m)
	    {
	      m = xc[i][j];
	      mi = i;
	    }
	}
      p[j] = mi;
    }
  
}


  

void updateCfromP(double **c, double *s, const mxArray *x, int *p, 
		  int K, int D, int n_samples, const int *update_flags)
{
  
  int i, j, k;
  

  for(i = 0; i < K; i++)
    if(update_flags[i])
      for(k = 0; k < D; k++)
	c[i][k] = 0.;
  
  for(j = 0; j < n_samples; j++)
    {
      i = p[j];
      if(update_flags[i])
	{
	  for(k = 0; k < D; k++)
	    c[i][k] += *(getRealItem(x, j, k));
	}
    }
  
  for(i = 0; i < K; i++)
    if(update_flags[i])
      {
	s[i] = l2norm(c[i], D);
	for(k = 0; k < D; k++)
	  c[i][k] /= s[i];
      }
  
}

      



double KMrun(const mxArray *x, double **c, double *s, int *p, double **xc, 
	     int K, int D, int n_samples, double tol)
{
  
  int iter;
  double q = 0., qOld = 0.;
  
  int done = 0;
  int i;
  
  int update_flags[K];
  
  for(i = 0; i < K; i++)
    update_flags[i] = 1;
  
  int max_iter = 1000;
  
  iter = 0;
  
  while(!done)
    {
      updateXC(x, c, xc, K, D, n_samples, update_flags);
      
      updatePfromC_KM(p, xc, K, n_samples);
      
      updateCfromP(c, s, x, p, K, D, n_samples, update_flags);
      
 

      qOld = q;
      
      q = 0;
      
      for(i = 0; i < K; i++)
	q += s[i];

      iter++;
      
      mexPrintf("KM deltaQ = %g\n", q-qOld);

      if(iter >= max_iter || ((iter>=1) && (q - qOld) < tol))
	done = 1;
      

    }
  
  updateXC(x, c, xc, K, D, n_samples, update_flags);

  return q - qOld;
  
}

   
void FVrun(const mxArray *x, double **c, double *s, int *p, double **xc, 
	   int K, int D, int n_samples, int *U)
{
  int i, j, i_from;
  
  int mj, m_from, m_to;
  double m = -1e20, q1, q2;
  int update_flags[K];
  
  for(i = 0; i < K; i++)
    update_flags[i] = 0;
  

  for(j = 0; j < n_samples; j++)
    if(U[j])
      {
	i_from = p[j];
	q1 = sqrt(s[i_from] * s[i_from] - 2 * s[i_from] * xc[i_from][j] + 1) -
	  s[i_from];
	
	for(i = 0; i < K; i++)
	  if(i != i_from)
	    {
	      q2 = sqrt(s[i]*s[i] + 2 * s[i] * xc[i][j] + 1) - s[i];
	      if(q1+q2 > m)
		{
		  m = q1+q2;
		  mj = j;
		  m_from = i_from;
		  m_to = i;
		}
	    }
	
      }
  
  p[mj] = m_to;
  U[mj] = 0;
  
  update_flags[m_from] = 1;
  update_flags[m_to] = 1;
  
  updateCfromP(c, s, x, p, K, D, n_samples, update_flags);

  updateXC(x, c, xc, K, D, n_samples, update_flags);
}


  

  


   
double KLFVrun(const mxArray *x, double **c, double *s, int *p, double **xc, 
	    int K, int D, int n_samples, int f, int **pv)
{
  
  int i, j, l;
  
  double objChange[f];
  
  int U[n_samples];
  int update_flags[K];
  
  double q, qOld;
  double m, ms;
  int mi;
  
  for(i = 0; i < K; i++)
    update_flags[i] = 1;
  

  q = 0.;
  for(i = 0; i < K; i++)
    q += s[i];
  


  for(i = 0; i < n_samples; i++)
    U[i] = 1;
  
  
  for(i = 0; i < f; i++)
    {
      for(j = 0; j < n_samples; j++)
	pv[i][j] = p[j];

      FVrun(x, c, s, p, xc, K, D, n_samples, U);

      qOld = q;
      q = 0.;
      for(l = 0; l < K; l++)
	q += s[l];
     
      objChange[i] = q - qOld;
    }
  
  ms = 0.;
  
  m = 0.;
  mi = 0;
  for(i = 0; i < f; i++)
    {
      ms += objChange[i];
      if(ms > m)
	{
	  m = ms;
	  mi = i+1;
	}
    }
  
  if(mi < f)
    for(j = 0; j < n_samples; j++)
      p[j] = pv[mi][j];

  updateCfromP(c, s, x, p, K, D, n_samples, update_flags);

  
  {
    double *pr;
    mxArray *rhs[1];
    rhs[0] = mxCreateDoubleMatrix(1, f, mxREAL);
    pr = mxGetPr(rhs[0]);
    for(i = 0; i < f; i++)
      pr[i] = objChange[i];
    mexCallMATLAB(0, NULL, 1, rhs, "plotFVQuality");
/*     mexCallMATLAB(0, NULL, 0, NULL, "keyboard"); */
  }
  
    
    

      


  return m;
  
}

      



    
  









void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  const mxArray *x, *z;
  
  double   **c; /* K x D */
  double   *s;  /* K  (quality of the cluster) */
  int   *p;  /* n_samples (the partition) */
  double   **xc; /* K x n_samples (the element similarity) */
  int   **pv;   /* f x n_samples (to keep the partitions in memory
		    during KLFV runs */
  int  **bb;
  double q;
  
  double *quality;
  



  double deltaQ_KM, deltaQ_FV;
  
  
  double  tol = 1.; /* the tolerance */

  
  
  int n_samples, D, K, i, j, l,  f, iter,  max_iter = 10, done = 0;

  
  /* check number of arguments: expects 2 or 3 inputs,  2 outputs */
  if ( nINP != 4  && nINP != 5)
    mexErrMsgTxt("Call with 4 or 5 inputs");
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
  x = pINP[0];
 
  z = pINP[2];

  
  f = (int)mxGetScalar(pINP[3]);
  

  if(nINP > 4) 
    max_iter = (int)mxGetScalar(pINP[4]);
  

  

  

  
  /* allocate memory */

  xc = (double **)mxCalloc(K, sizeof(double *));
  for(i = 0 ; i < K ; i++)
    xc[i] = mxCalloc(n_samples, sizeof(double));
  
  c = (double **)mxCalloc(K, sizeof(double *));
  for(i = 0 ; i < K ; i++)
    c[i] = mxCalloc(D, sizeof(double));
  
  pv = (int **)mxCalloc(f, sizeof(int *));
  for(i = 0; i < f; i++)
    pv[i] = mxCalloc(n_samples, sizeof(int));
  
  quality = (double *)mxCalloc(max_iter, sizeof(double));
  

  p = (int *)mxCalloc(n_samples, sizeof(int));
  
  s = (double *)mxCalloc(K, sizeof(double));
  
  
  
  /* initialize centroids with provided initial condition */
	    
  for(i = 0 ; i < K ; i++)
    for(j = 0; j < D; j++)
      c[i][j] = *(getRealItem(z, i, j));
  

  done = 0;
  

  iter = 0;
  
  while(!done)
    {

      
      deltaQ_KM = KMrun(x, c, s, p, xc, K, D, n_samples, tol);

     
      deltaQ_FV = KLFVrun(x, c, s, p, xc, K, D, n_samples, f, pv);

      q = 0.;
      for(l = 0; l < K; l++)
	q += s[l];

      quality[iter] = q;
      
      iter++;
      mexPrintf("deltaQ_FV = %g\n", deltaQ_FV);
      mexPrintf("deltaQ_KM = %g\n", deltaQ_KM);
      
      if(deltaQ_FV + deltaQ_KM < tol || iter == max_iter)
	 done = 1;



    }
  
      
  /* prepare outputs */

  pOUT[0] = mxCreateDoubleMatrix(K, D, mxREAL);
  

  for(i = 0; i < K; i++)
    for(j = 0; j < D; j++)
      *(getRealItem(pOUT[0], i, j)) = c[i][j];
  



  pOUT[1] = mxCreateDoubleMatrix(1, n_samples, mxREAL);
  
  for(i = 0; i < n_samples; i++)
    *(getRealItem(pOUT[1], 0, i)) = (double)p[i] + 1.;


  if(nOUT > 2)
    {
      
      pOUT[2] = mxCreateDoubleMatrix(1, iter, mxREAL);
  
      for(i = 0; i < iter; i++)
	*(getRealItem(pOUT[2], 0, i)) = quality[i];
    }
  

      
  /* clean up */

  for(i = 0; i < K; i++)
    {
      mxFree(c[i]);
      mxFree(xc[i]);
    }
  
  for(i = 0; i < f; i++)
    mxFree(pv[i]);
  

  mxFree(c);
  mxFree(xc);
  mxFree(pv);
  mxFree(s);
  mxFree(p);
  
  pv = (int **)mxCalloc(f, sizeof(int *));
  for(i = 0; i < f; i++)
    pv[i] = mxCalloc(n_samples, sizeof(int));

   for(i = 0; i < f; i++)
    mxFree(pv[i]);
/*         bb = mxCalloc(1000, sizeof(int *)); */
/*       for(i = 0; i < 1000; i++) */
/* 	bb[i] = mxCalloc(1000, sizeof(int)); */
/*       for(i = 0; i < 1000; i++) */
/* 	mxFree(bb[i]); */
/*       mxFree(bb); */

  
}

