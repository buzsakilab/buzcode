 /* ETAverage event-triggered average
  * MEX file 
  * fpbattaglia 2007 
  * input  ev: an event time series 
  *        t1 time series  in 1/10000 sec
  *               (assumed to be sorted)
  *        d the values in the dataset 
  *        binsize: the size of the bin in ms
  *        nbins: the number of bins to compute 
  * output C the event triggered average 
 *         B (optional) a vector with the times corresponding to the bins
 */







 /* adepted from PECorr original header
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
#define D_IDX 2 
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
  
  double *t1, *D;
  double binsize;
  double *C, *B, *S;

  
  double dNow,dNow2,tmpVal, w, lbound, rbound;
  
  int nbins, nt1, nev, nNow, *nSamples;
  int i1, ie, l1,j;
  int isempty = 0;

  i1=0;
  ie=0;
  
    if (nINP != 5)
    mexErrMsgTxt("Call with ev, t1, d, binsize and nbins  as inputs.");
  
  if (mxGetM(pINP[EV_IDX]) != 1 && mxGetN(pINP[EV_IDX]) != 1)
    mexErrMsgTxt("ev must be a row or column vector");
  if (mxGetM(pINP[BINSIZE_IDX]) * mxGetN(pINP[BINSIZE_IDX]) != 1)
    mexErrMsgTxt("binsize must be scalar");
  if (mxGetM(pINP[NBINS_IDX]) * mxGetN(pINP[NBINS_IDX]) != 1)
    mexErrMsgTxt("nbins must be scalar");



  if (mxGetM(pINP[T1_IDX]) * mxGetN(pINP[T1_IDX]) == 0)
    isempty = 1;

  

  if(!isempty) 
    {
      if (mxGetM(pINP[T1_IDX]) != 1 && mxGetN(pINP[T1_IDX]) != 1)
	mexErrMsgTxt("t1 must be a row or column vector");
 
      if (mxGetM(pINP[D_IDX]) != 1 && mxGetN(pINP[D_IDX]) != 1)
	mexErrMsgTxt("t1 must be a row or column vector");

      if(mxGetM(pINP[T1_IDX]) != mxGetM(pINP[D_IDX]) ||
	 mxGetN(pINP[T1_IDX]) != mxGetN(pINP[D_IDX]))
	mexErrMsgTxt("t1 and d must have the same dims");
      
    }
  
  /* unpack inputs */
  if(!isempty)
    {
      nev = mxGetM(pINP[EV_IDX]) * mxGetN(pINP[EV_IDX]);
      ev = mxGetPr(pINP[EV_IDX]);  
      nt1 = mxGetM(pINP[T1_IDX]) * mxGetN(pINP[T1_IDX]);
      t1 = mxGetPr(pINP[T1_IDX]);
      D = mxGetPr(pINP[D_IDX]);
      
    }
  
  binsize = mxGetScalar(pINP[BINSIZE_IDX]);
  nbins = (int)mxGetScalar(pINP[NBINS_IDX]);
  

  


  /* we want nbins to be odd */
  if ((nbins / 2) * 2 == nbins)
    nbins++;

  pOUT[0] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
  pOUT[1] = mxCreateDoubleMatrix(nbins, 1, mxREAL);

  if(nOUT >= 3)
    {
      double m;
      
      pOUT[2] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
      B =  mxGetPr(pOUT[2]);
      m = - binsize * ((double)(nbins-1) / 2);
      
      for(j = 0; j < nbins; j++)
        B[j] = m + j * binsize;

    }
  
  
  if(!isempty)
    {
      C = mxGetPr(pOUT[0]);
	  S = mxGetPr(pOUT[1]);

      nSamples = (int *)mxCalloc(nbins, sizeof(int));
      
      binsize *= 10;  
      /* cross correlations */
  
      w = ((double)(nbins-1) / 2) * binsize;

  
      for(ie = 0; ie < nev; ie++)
	{
	  lbound = ev[ie] - w;

	  while(t1[i1] < lbound && i1 < nt1 -1)
	    i1++;
	  while(i1 > 1 && t1[i1-1] > lbound)
	    i1--;

	  
	  rbound = lbound;
	  l1 = i1;

	  rbound = lbound;
      
	  for(j = 0; j < nbins; j++)
	    {

	      rbound += binsize;

	      nNow = 0;
	      dNow = 0.;
	      dNow2 = 0.;
	      while(t1[l1] < rbound && l1 < nt1-1)
		{
		  nNow++;
		  tmpVal = D[l1];	
		  dNow += tmpVal;
		  tmpVal *= tmpVal; 	
		  dNow2 += tmpVal;
			
		  l1++;
		  
		}

	
	      if(nNow > 0)
		{
		  nSamples[j]++;
		  C[j] += dNow/nNow;
		  S[j] += dNow2/nNow;	
		}
	      
	    }
	  
	}
  
  
      for(j = 0; j < nbins; j++)
	{
	  C[j] /= nSamples[j];
	  S[j] /= nSamples[j];
		tmpVal = C[j];
		tmpVal *= tmpVal;
	  S[j] -= tmpVal ;
	}
  
      
      mxFree((void *)nSamples);
      
    }
  
  
      
}
  
  
		 
