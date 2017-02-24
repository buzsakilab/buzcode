 /* MakeCCPairs
 * cross correlations for all the specified pairs in a set of spike trains
 * MEX file
 * 
 * batta 2002 (derived by CrossCorr)
 * 
 * input: S: a cell array of n_cells series to cross correlate in 1/10000 sec
 *                (assumed to be sorted) 
 * cell_I, cell_J, vectors indicating the indices of hte pairs to be processed
 *        binsize: the binsize for the cross corr histogram in msec
 *        nbins: the number of bins
 * output: C a n_pairs *nbins matrix containing the 
 *         cross correlation histograms
 *         B (optional) a vector with the times corresponding to the bins
 *
 * version 1.0
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void CrossCorr(double *, double *, int , double *, int , 
	       int , int );



void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{

  
  double *t1;
  double *t2;
  double binsize;
  double *C, *CC, *B;
  

  int *cell_I, *cell_J;
  double *cI, *cJ;
  

  int nbins, nt1, nt2, n_cells, n_pairs;
  int l, j;
  
  /* check number of arguments: expects 4 inputs, 1 or 2 outputs */
  if (nINP != 5)
    mexErrMsgTxt("Call with ts, cell_I, cell_J binsize and nbins  as inputs.");
  if (nOUT != 1 && nOUT != 2)
    mexErrMsgTxt("Requires one or two outputs.");

  /* check validity of inputs */

  /* input 0 is a cell array containing the spike trains */ 
  if ((mxGetM(pINP[0]) != 1 && mxGetN(pINP[0]) != 1) ||
      (!mxIsCell(pINP[0])))
    mexErrMsgTxt("t1 must be a row or column cell array");
  n_cells = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  



  if (mxGetM(pINP[1]) != 1 && mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("cell_I must be a row or column vector");
  n_pairs = mxGetM(pINP[1]) * mxGetN(pINP[1]);

  mexPrintf("n_pairs = %d\n\n", n_pairs);
  
  
  if (mxGetM(pINP[2]) != 1 && mxGetN(pINP[2]) != 1)
    mexErrMsgTxt("cell_J must be a row or column vector");
  if(mxGetM(pINP[2]) * mxGetN(pINP[2]) != n_pairs)
    mexErrMsgTxt("cell_I and cell_j must have the same size");
  



  if (mxGetM(pINP[3]) * mxGetN(pINP[3]) != 1)
    mexErrMsgTxt("binsize must be scalar");
  if (mxGetM(pINP[4]) * mxGetN(pINP[4]) != 1)
    mexErrMsgTxt("nbins must be scalar");

  /* unpack inputs */

  cell_I = (int *)mxCalloc(n_pairs, sizeof(int));
  cell_J = (int *)mxCalloc(n_pairs, sizeof(int));  
  cI = mxGetPr(pINP[1]);
  cJ = mxGetPr(pINP[2]);

  for(l = 0; l < n_pairs; l++)
    {
      cell_I[l] = (int)(cI[l]) - 1;
      cell_J[l] = (int)(cJ[l]) - 1;
    }
  

  binsize = mxGetScalar(pINP[3]);
  nbins = (int)mxGetScalar(pINP[4]);
  /* we want nbins to be odd */
  if ((nbins / 2) * 2 == nbins)
    nbins++;
  

  nt1 = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  t1 = mxGetPr(pINP[0]);
  nt2 = mxGetM(pINP[1]) * mxGetN(pINP[1]);
  t2 = mxGetPr(pINP[1]);
  
  pOUT[0] = mxCreateDoubleMatrix(n_pairs, nbins, mxREAL);
  CC = mxGetPr(pOUT[0]);
  C = (double *)mxCalloc(nbins, sizeof(double));
  





  for(l = 0; l <n_pairs; l++)
    {
      mxArray *mxa;
      int subs[2];
      
      if(l % 1000 == 0)
	mexPrintf("l = %d\n", l);

      mxa = mxGetCell(pINP[0], cell_I[l]);
      nt1 = mxGetM(mxa) * mxGetN(mxa);
      t1 = mxGetPr(mxa);
      mxa = mxGetCell(pINP[0], cell_J[l]);
      nt2 = mxGetM(mxa) * mxGetN(mxa);
      t2 = mxGetPr(mxa);

      CrossCorr(C, t1, nt1, t2, nt2, nbins, binsize);
      
      for(j = 0; j < nbins; j++)
	{
	  subs[0] = l;
	  subs[1] = j;
	  CC[mxCalcSingleSubscript(pOUT[0], 2, subs)] = C[j];
	}
      
    }
  





  if(nOUT == 2)
    {
      double m;
      
      pOUT[1] = mxCreateDoubleMatrix(1, nbins, mxREAL);
      B =  mxGetPr(pOUT[1]);
      m = - binsize * (nbins / 2);
      
      for(j = 0; j < nbins; j++)
	B[j] = m + j * binsize;

    }





  mxFree((void *)cell_I);
  mxFree((void *)cell_J);
  mxFree((void *)C);
  

}


void CrossCorr(double *C, double *t1, int nt1, double *t2, int nt2, 
		  int nbins, int binsize)
		    
{
    

  int i1 = 0, i2 = 0, l, j, k;
  double w, lbound, rbound;
  binsize *= 10;  
  /* cross correlations */
  
  w = (nbins / 2) * binsize;

  
  for(j = 0; j < nbins; j++)
    C[j] = 0.;
  
  if((nt1 == 0) ||  (nt2 == 0))
    return;
  
  
  for(i1 = 0; i1 < nt1; i1++)
    {
      lbound = t1[i1] - w;
      while(t2[i2] < lbound && i2 < nt2 -1)
	i2++;
      while(t2[i2-1] > lbound && i2 > 1)
	i2--;
      rbound = lbound;
      l = i2;
      
      for(j = 0; j < nbins; j++)
	{
	  k = 0;
	  rbound += binsize;
	  while(t2[l] < rbound && l < nt2-1)
	    {
	      l++;
	      k++;
	    }

	  C[j] += k;
	}
    }
  
  
/*    for(j = 0; j < nbins; j++) */
/*      C[j] /= nt1 * binsize / 10000; */
  
      
}
  
  
		 
