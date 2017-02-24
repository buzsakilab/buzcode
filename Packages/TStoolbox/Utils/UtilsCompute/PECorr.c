 /* PECorr Peri-Event correlation coefficient
 * MEX file
 * 
 * batta 2002
 * 
 * input: ev: an event time series 
 *        t1, t2: two time series  in 1/10000 sec
 *                (assumed to be sorted) the correlation coefficient
 *        as a function of  time gap from eventr times will be computed
 *        binsize: the binsize for the cross corr histogram in msec
 *        nbins: the number of bins
 * output: C the peri-event correlation coefficient 
 *         B (optional) a vector with the times corresponding to the bins
 *         
 * version 0.1
 -----------------------------------*/

#define EV_IDX 0
#define T1_IDX 1
#define T2_IDX 2
#define BINSIZE_IDX 3
#define NBINS_IDX 4

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  double *ev;
  
  double *t1;
  double *t2;
  double binsize;
  double *C, *B;
  double *C1, *C2, *C12, *C11, *C22;
  
  double w, lbound, rbound;
  
  int nbins, nt1, nt2, nev;
  int i1 = 0, i2 = 0, ie = 0, l1, l2, j, k1, k2;
  int isempty = 0;
  
  /* check number of arguments: expects 4 inputs, 1 or 2 outputs */
  if (nINP != 5)
    mexErrMsgTxt("Call with ev, t1, t2, binsize and nbins  as inputs.");
  if (nOUT != 1 && nOUT != 2)
    mexErrMsgTxt("Requires one or two outputs.");

  /* check validity of inputs */
  if (mxGetM(pINP[EV_IDX]) != 1 && mxGetN(pINP[EV_IDX]) != 1)
    mexErrMsgTxt("ev must be a row or column vector");
  if (mxGetM(pINP[BINSIZE_IDX]) * mxGetN(pINP[BINSIZE_IDX]) != 1)
    mexErrMsgTxt("binsize must be scalar");
  if (mxGetM(pINP[NBINS_IDX]) * mxGetN(pINP[NBINS_IDX]) != 1)
    mexErrMsgTxt("nbins must be scalar");


  if (mxGetM(pINP[T1_IDX]) * mxGetN(pINP[T1_IDX]) == 0)
    isempty = 1;
  if (mxGetM(pINP[T2_IDX]) * mxGetN(pINP[T2_IDX]) == 0)
    isempty = 1;
  

  if(!isempty) 
    {
      if (mxGetM(pINP[T1_IDX]) != 1 && mxGetN(pINP[T1_IDX]) != 1)
	mexErrMsgTxt("t1 must be a row or column vector");
      if (mxGetM(pINP[T2_IDX]) != 1 && mxGetN(pINP[T2_IDX]) != 1)
	mexErrMsgTxt("t2 must be a row or column vector");
    }
  
  /* unpack inputs */
  if(!isempty)
    {
      nev = mxGetM(pINP[EV_IDX]) * mxGetN(pINP[EV_IDX]);
      ev = mxGetPr(pINP[EV_IDX]);  
      nt1 = mxGetM(pINP[T1_IDX]) * mxGetN(pINP[T1_IDX]);
      t1 = mxGetPr(pINP[T1_IDX]);
      nt2 = mxGetM(pINP[T2_IDX]) * mxGetN(pINP[T2_IDX]);
      t2 = mxGetPr(pINP[T2_IDX]);
    }
  
  binsize = mxGetScalar(pINP[BINSIZE_IDX]);
  nbins = (int)mxGetScalar(pINP[NBINS_IDX]);
  



  /* we want nbins to be odd */
  if ((nbins / 2) * 2 == nbins)
    nbins++;

  

  pOUT[0] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
  if(nOUT >= 2)
    {
      double m;
      
      pOUT[1] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
      B =  mxGetPr(pOUT[1]);
      m = - binsize * ((double)nbins / 2);
      
      for(j = 0; j < nbins; j++)
	B[j] = m + j * binsize;

    }


  if(!isempty)
    {
      C = mxGetPr(pOUT[0]);

  
      C1 = (double *)mxCalloc(nbins, sizeof(double));
      C2 = (double *)mxCalloc(nbins, sizeof(double));
      C12 = (double *)mxCalloc(nbins, sizeof(double));
      C11 = (double *)mxCalloc(nbins, sizeof(double));
      C22 = (double *)mxCalloc(nbins, sizeof(double));
  




      binsize *= 10;  
      /* cross correlations */
  
      w = ((double)nbins / 2) * binsize;

  
      for(ie = 0; ie < nev; ie++)
	{
	  lbound = ev[ie] - w;

	  while(t1[i1] < lbound && i1 < nt1 -1)
	    i1++;
	  while(t1[i1-1] > lbound && i1 > 1)
	    i1--;
	  rbound = lbound;
	  l1 = i1;

	  while(t2[i2] < lbound && i2 < nt2 -1)
	    i2++;
	  while(t2[i2-1] > lbound && i2 > 1)
	    i2--;
	  rbound = lbound;
	  l2 = i2;
      
	  for(j = 0; j < nbins; j++)
	    {

	      rbound += binsize;

	      k1 = 0;
	      while(t1[l1] < rbound && l1 < nt1-1)
		{
		  l1++;
		  k1++;
		}

	      k2 = 0;
	      while(t2[l2] < rbound && l2 < nt2-1)
		{
		  l2++;
		  k2++;
		}

	      C1[j] += k1;
	      C2[j] += k2;
	      C12[j] += k1 * k2;
	      C11[j] += k1 * k1;
	      C22[j] += k2 * k2;
	
	  
	    }
	}
  
  
      for(j = 0; j < nbins; j++)
	{
	  C1[j] /= nev;
	  C2[j] /= nev;
	  C12[j] /= nev;
	  C11[j] /= nev;
	  C22[j] /= nev;
      
	  C[j] = (C12[j] - (C1[j]*C2[j])) / 
	    sqrt((C11[j] - C1[j] * C1[j]) * (C22[j] - C2[j] * C2[j]));
	}
  
      
  

  


      mxFree(C1);
      mxFree(C2);
      mxFree(C12);
      mxFree(C11);
      mxFree(C22);
    }
  
      
}
  
  
		 
