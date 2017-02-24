/* [ov] = MakeReactOverlap(samples, templates) */
/* MEX file  
 * for each samples, computes the maximum overlap with on eof
 * the templates
 * INPUTS:
 * samples: a n_samples x n_cells matrix 
 * templates: a n_templates x n_cells matrix
 * ov: a n_samplesX1 vector of overlap (Pearson formula)
 * 
 * batta 2002 */




#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <matrix.h>


inline double computeVar(double *x, int l)
{
  double norm = 0.;
  int i = 0;
  
  for(i = 0; i < l; i++)
    norm += x[i] * x[i];
  
  return norm / (l-1);
}

inline double computeMean(double *x, int l)
{
  double norm = 0.;
  int i = 0;
  
  for(i = 0; i < l; i++)
    norm += x[i];
  
  return norm / (l);
}


inline double computeCov(double *x, double *y, int l)
{
  double norm = 0.;
  int i = 0;
  
  for(i = 0; i < l; i++)
    norm += x[i] * y[i];
  
  return norm / (l);
}

inline double computeDot(double *x, double *y, int l)
{
  double norm = 0.;
  int i = 0;
  
  for(i = 0; i < l; i++)
    norm += x[i] * y[i];
  
  return norm;
}

inline double computeMax(double *x, int l)
{
  double norm = x[0];
  int i = 0;
  
  for(i = 1; i < l; i++)
    if(x[i] > norm)
      norm = x[i];
  
  
  return norm;
}


void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  
  double *xo, **x, **templates, *ov;

  
  
  
  int n_samples, n_cells, n_templates, i, j, l;
  int subs[2];
  
  /* check number of arguments: expects 2 inputs,  1 outputs */
  if (nINP != 2)
    mexErrMsgTxt("Call with 2 inputs");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs */


   

  n_samples = mxGetM(pINP[0]);
  n_cells = mxGetN(pINP[0]);

  
  if(mxGetN(pINP[1]) != n_cells)
    mexErrMsgTxt("Samples and templates musts have same number of columns");
  
  n_templates = mxGetM(pINP[1]);
  
  
  



  /* unpack inputs */

  /* the samples */

  xo = mxGetPr(pINP[0]);

  
  x = mxCalloc(n_samples, sizeof(double *));

  for(i = 0; i < n_samples; i++)
    {
      x[i] = mxCalloc(n_cells, sizeof(double));
      if(x[i]  == NULL)
	mexErrMsgTxt("problem allocating memory");
    }



  for(i = 0; i < n_samples; i++)
    {
      for(j = 0; j < n_cells; j++)
	{
	  int  idx;
	  subs[0] = i;
	  subs[1] = j;
  	  idx = mxCalcSingleSubscript(pINP[0], 2, subs); 
	  x[i][j] = xo[idx];
	}
    }
  
    
  /* the templates */

  
  xo = mxGetPr(pINP[1]);
  templates = mxCalloc(n_templates, sizeof(double *));
  for(i = 0; i < n_templates; i++)
    {
      templates[i] = mxCalloc(n_cells, sizeof(double));
      if(templates[i]  == NULL)
	mexErrMsgTxt("problem allocating memory");
    }
  
  for(i = 0; i < n_templates; i++)
    {
      for(j = 0; j < n_cells; j++)
	{
	  int  idx;
	  subs[0] = i;
	  subs[1] = j;
  	  idx = mxCalcSingleSubscript(pINP[1], 2, subs); 
	  templates[i][j] = xo[idx];
	}
    }





 



  
  /* prepare outputs */

  pOUT[0] = mxCreateDoubleMatrix(n_samples, n_templates, mxREAL);

  ov = mxGetPr(pOUT[0]);
  


  /* compute Pearsons */

  
  for(i = 0; i < n_samples; i++)
    for(l = 0; l < n_templates ; l++)
      {
	subs[0] = i;
	subs[1] = l;
	ov[mxCalcSingleSubscript(pOUT[0], 2, subs)] = 
	  computeDot(x[i], templates[l], n_cells);

      }
  



  /* clean-up stuff */
  for(i = 0; i < n_templates; i++)
    mxFree((void *)templates[i]);
  mxFree((void *)templates);
  
  for(i = 0; i < n_samples; i++)
    mxFree((void *)x[i]);
  mxFree((void *)x);
  


  



  
  
  
}
