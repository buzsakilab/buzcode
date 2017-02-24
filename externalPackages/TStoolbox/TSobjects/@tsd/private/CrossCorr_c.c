 /* CrossCorrSD_b
 * cross correlations
 * MEX file
 * 
 * batta 1999
 * 
 * input: t1, t2: two time series to cross correlate in 1/10000 sec
 *                (assumed to be sorted) 
 *        binsize: the binsize for the cross corr histogram in msec
 *        nbins: the number of bins
 * output: C the cross correlation histogram
 *         B (optional) a vector with the times corresponding to the bins
 *         SC (optional) the S.D. of the crosscorrelation (of the
 *          profile of unit 2 as a PETH on the train of unit 1
 * this version takes care of the boundary conditions, to prevent
 spurious decline  effects in the CC, by discardign the points in the triggering time series that are closer than binsize * (nbins/2) to the first or last point in the triggered time series 

 * version 1.0

readapted for the new tsd class 

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  
  double *t1;
  double *t2;
  double binsize;
  double *C, *B, *SC;
  double w, lbound, rbound;


  int fix_boundaries = 0;
  int sem = 0;
  
  int nbins, nt1, nt2, t1start;
  int i1 = 0, i2 = 0, l, j, k;
  
  /* check number of arguments: expects 4 or 5 inputs, 2 or 3 outputs */
  if (nINP != 4 && nINP != 5 && nINP != 6)
    mexErrMsgTxt("Call with t1, t2, binsize, nbins, [fix_boundaries], [sem] as inputs.");
  if (nOUT != 2 && nOUT != 3)
    mexErrMsgTxt("Requires two or three outputs.");

  /* check validity of inputs */
  if (mxGetM(pINP[0]) != 1 && mxGetN(pINP[0]) != 1)
    mexErrMsgTxt("t1 must be a row or column vector");
  if (mxGetM(pINP[1]) != 1 && mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("t2 must be a row or column vector");
  if (mxGetM(pINP[2]) * mxGetN(pINP[2]) != 1)
    mexErrMsgTxt("binsize must be scalar");
  if (mxGetM(pINP[3]) * mxGetN(pINP[3]) != 1)
    mexErrMsgTxt("nbins must be scalar");

  if(nINP >= 5) 
    fix_boundaries = (int) mxGetScalar(pINP[4]);
  
  if (nINP >=6)
    sem = (int) mxGetScalar(pINP[5]);
  



  /* unpack inputs */
  nt1 = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  t1 = mxGetPr(pINP[0]);
  nt2 = mxGetM(pINP[1]) * mxGetN(pINP[1]);
  t2 = mxGetPr(pINP[1]);
  binsize = mxGetScalar(pINP[2]);
  nbins = (int)mxGetScalar(pINP[3]);
  



  /* we want nbins to be odd */
  if ((nbins / 2) * 2 == nbins)
    nbins++;

  

  pOUT[0] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
  C = mxGetPr(pOUT[0]);
  if(nOUT >= 2)
    {
      double m;
      
      pOUT[1] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
      B =  mxGetPr(pOUT[1]);
      m = - binsize * (nbins / 2);
      
      for(j = 0; j < nbins; j++)
	B[j] = m + j * binsize;

    }
  

  if(nOUT == 3)
    {
      pOUT[2] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
      SC =  mxGetPr(pOUT[2]);

    }

  if(nt1 == 0 || nt2 == 0) 
    {
/*       mexErrMsgTxt("no data points"); */
      
/*       mexPrintf("We return here\n"); */

      
      return;
    }
  
  /* cross correlations */
  
  w = (nbins / 2) * binsize;
  

  i1 =  0;
  t1start = 0;
  
  if (fix_boundaries)
    {
      
      while(i1 < nt1 && t1[i1] < t2[0]+w)
	i1++;
      t1start = i1;
  
      i1 = nt1-1;
      while(i1 >= 0 && t1[i1] > t2[nt2-1]-w)
	i1--;
  
      nt1 = i1+1;
    }
  
  
  
  for(i1 = t1start; i1 < nt1; i1++)
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
	  if(nOUT == 3)
	    SC[j] += k*k;
	  
	}
    }
  
  
  for(j = 0; j < nbins; j++)
      C[j] /= (nt1 - t1start) * binsize;
  
  {
    double norm;
    
    norm = (nt1 - t1start) * binsize * binsize;
    
    if(nOUT == 3)
      for(j = 0; j < nbins; j++)
	{
	  SC[j] /= norm;;
	  SC[j] = sqrt(SC[j]);
	  if(sem)
	    SC[j] /= sqrt((double)nt1 - t1start);

	}
  
  }
  
      
}
  
  
		 
