 /* JPETH
 * Joint Peri Event Time Histogram
 * MEX file
 * inspired from CrossCorr by francesco battaglia (header below)
 * adrien peyrache 2007
 * 
 * input: tref, t1,t2: three time series to cross correlate in 1/10000 sec
 *                (assumed to be sorted). Tref is the common reference.
 *        binsize: the binsize for the cross corr histogram in msec
 *        nbins: the number of bins
 * output: C the cross correlation histogram
 *         B (optional) a vector with the times corresponding to the bins
 *
 * version 1.0
 -----------------------------------*/


 /* CrossCorr
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
 *
 * version 1.0
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
  double *t3;
  double binsize;
  double *C, *cross, *B;
  double w, lbound, rbound, rbound2;
  
  int nbins, nt1, nt2, nt3;
  int i1 = 0, i2 = 0, i3 = 0, l2, l3, j, k, valj, valk;
  
  /* check number of arguments: expects 4 inputs, 1 or 2 outputs */
  if (nINP != 5)
    mexErrMsgTxt("Call with tref, t2, t3, binsize and nbins  as inputs.");
  if (nOUT != 1 && nOUT != 2)
    mexErrMsgTxt("Requires one or two outputs.");

  /* check validity of inputs */
  if (mxGetM(pINP[0]) != 1 && mxGetN(pINP[0]) != 1)
    mexErrMsgTxt("t1 must be a row or column vector");
  if (mxGetM(pINP[1]) != 1 && mxGetN(pINP[1]) != 1)
    mexErrMsgTxt("t2 must be a row or column vector");
  if (mxGetM(pINP[2]) != 1 && mxGetN(pINP[2]) != 1)
    mexErrMsgTxt("t2 must be a row or column vector");
  if (mxGetM(pINP[3]) * mxGetN(pINP[3]) != 1)
    mexErrMsgTxt("binsize must be scalar");
  if (mxGetM(pINP[4]) * mxGetN(pINP[4]) != 1)
    mexErrMsgTxt("nbins must be scalar");

  /* unpack inputs */
  nt1 = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  t1 = mxGetPr(pINP[0]);
  nt2 = mxGetM(pINP[1]) * mxGetN(pINP[1]);
  t2 = mxGetPr(pINP[1]);
  nt3 = mxGetM(pINP[2]) * mxGetN(pINP[2]);
  t3 = mxGetPr(pINP[2]);

  binsize = mxGetScalar(pINP[3]);
  nbins = (int)mxGetScalar(pINP[4]);
  



  /* we want nbins to be odd */
  if ((nbins / 2) * 2 == nbins)
    nbins++;

  

  pOUT[0] = mxCreateDoubleMatrix(nbins, nbins, mxREAL);
  C = mxGetPr(pOUT[0]);
  if(nOUT == 2)
    {
      double m;
      
      pOUT[1] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
      B =  mxGetPr(pOUT[1]);
      m = - binsize * (nbins / 2);
      
      for(j = 0; j < nbins; j++)
	B[j] = m + j * binsize;

    }
  
  binsize *= 10;
  /* cross correlations */
  
  w = (nbins / 2) * binsize;

  for(i1 = 0; i1 < nt1; i1++)
    {
      lbound = t1[i1] - w;
      while(t2[i2] < lbound && i2 < nt2 -1)
	i2++;
      while(i2 > 1 && t2[i2-1] > lbound )
	i2--;
      rbound = lbound;

      while(t3[i3] < lbound && i3 < nt3 -1)
	i3++;
      while(t3[i3-1] > lbound && i3 > 1)
	i3--;

      rbound = lbound;
      l2 = i2;

      for(j = 0; j < nbins; j++)
	{
		valj = 0;
		rbound += binsize;

		while(t2[l2] < rbound && l2 < nt2-1)
		  {
			l2++;
			valj++;
		  }
	
		rbound2 = lbound;
		l3 = i3;

		for(k = 0; k < nbins; k++)
		  {
			valk = 0;
			rbound2 += binsize;
			while(t3[l3] < rbound2 && l3 < nt3-1)
			  {
				l3++;
				valk++;
			  }
			
			valk = valk*valj;	
			C[k+j*nbins] += valk;
		  }
	}
    }
  
  
  for(j = 0; j < nbins; j++)
	for(k = 0; k < nbins; k++)
		C[k+j*nbins] /= nt1 * binsize / 10000;
  
      
}
  
  
		 
