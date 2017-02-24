/* LC = LagStateCorr(Q, maxbin, intervals)

   INPUTS: 
   Q: a nBin x nCells matrix of spike counts
   maxbin: the number of lagged correlations to be computed
   intervals: a nIntervals x 2 matrix of indices representing the
              beginning and the end of intervals to be considered for
	      the analysis, if not present take the whole matrix

*/


#include "mex.h"
#include <math.h>
#include <matrix.h>

#define MATLAB

#ifdef MATLAB
#define MAT_ITEM(m, i, j) ( (m)[(i) + nRow * (j)] )
#endif

#define Q_IDX 0
#define MAXBIN_IDX 1
#define INTERVALS_IDX 2

#define LC_IDX 0

static int nBins, nCells, nIntervals, maxBin;
static double *intervals, *Q, *LC;


static double Pearson_r(int i, int j, double *mean, double *std)
{
  int k, nRow = nBins;
  double cov = 0.;
  
  for(k = 0; k < nCells; k++)
    cov += (MAT_ITEM(Q, i, k) - mean[i]) * (MAT_ITEM(Q, j, k) - mean[j]);
  
  cov /= nCells;
  
  if(std[i] * std[j] > 0.)
    return cov / (std[i] * std[j]);
  else return 0.;
  
}

  



static void lagStateCorr()
{
  int i, j, nRow, intv;
  double *std, *mean, *LC_n;
  double start_int, stop_int;
  
  std = mxCalloc(nBins, sizeof(double));
  mean = mxCalloc(nBins, sizeof(double));
  LC_n = mxCalloc(maxBin, sizeof(double));
  


  nRow = nBins;
  for(i = 0; i < nBins; i++)
    {
      double v = 0., m = 0.;;
      for(j = 0; j < nCells; j++)
	{
	  v += MAT_ITEM(Q, i, j) * MAT_ITEM(Q, i, j) ;
	  m += MAT_ITEM(Q, i, j);
	}
      m /= nCells;
      v /= nCells;
      mean[i] = m;
      std[i] = sqrt((nCells / (nCells-1)) * (v - m * m));
    }

  for(i = 0; i < maxBin; i++)
    LC[0] = 0.;
  
  i = 0;
  intv = 0;
  
  while(i < nBins)
    {
      nRow = nIntervals;
      start_int = MAT_ITEM(intervals, intv, 0);
      stop_int = MAT_ITEM(intervals, intv, 1);
      if(start_int > nBins)
	start_int = nBins;
      if(stop_int > nBins)
	stop_int = nBins;
      
      nRow = nBins;
      while(i < start_int)
	i++;
      while(i < stop_int)
	{
	  j = 0;
	  while(j < maxBin && i+j < stop_int)
	    {
	      LC[j] += Pearson_r(i, i+j, mean, std);
	      LC_n[j] += 1.;
	      j++;
	    }
	  i++;
	}
      intv++;
    }
  
  for(i = 0; i < maxBin; i++)
    LC[i] /= LC_n[i];
  

  mxFree(std);
  mxFree(mean);
  mxFree(LC_n);
}

	  
  




void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  int i ;
  int nRow;
  


  /* check number of arguments: expects 2 or 3 inputs, 1  output */
  if (nINP < 2 || nINP > 3)
    mexErrMsgTxt("Call with Q, maxbin and intervals (optional)  as inputs.");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");


  nBins = mxGetM(pINP[Q_IDX]);
  nCells = mxGetN(pINP[Q_IDX]);
  Q = mxGetPr(pINP[Q_IDX]);
  



  if(mxGetM(pINP[MAXBIN_IDX]) * mxGetN(pINP[MAXBIN_IDX]) != 1)
    mexErrMsgTxt("maxbin must be scalar.");
  maxBin = mxGetScalar(pINP[MAXBIN_IDX]);
  
  if(nINP > 2)
    {
      if(mxGetN(pINP[INTERVALS_IDX]) != 2)
	mexErrMsgTxt("intervals must be nIntervals x 2");
      intervals = mxGetPr(pINP[INTERVALS_IDX]);
      nIntervals = mxGetM(pINP[INTERVALS_IDX]);
      for(i = 0; i< nIntervals; i++)
	{
	  nRow = nIntervals;
	  if(MAT_ITEM(intervals, i, 0) < 0)
	    MAT_ITEM(intervals, i, 0) = 0;
	  if(MAT_ITEM(intervals, i, 0) > nIntervals-1)
	    MAT_ITEM(intervals, i, 0) = nIntervals-1;
	  if(MAT_ITEM(intervals, i, 1) < 0)
	    MAT_ITEM(intervals, i, 1) = 0;
	  if(MAT_ITEM(intervals, i, 1) > nIntervals-1)
	    MAT_ITEM(intervals, i, 1) = nIntervals-1;
	}
    }
  else
    {
      intervals = mxCalloc(2, sizeof(double));
      nIntervals = 1;
      intervals[0] = 0;
      intervals[1] = nBins-1;

      
    }
  
  pOUT[LC_IDX] = mxCreateDoubleMatrix(1, maxBin, mxREAL);
  LC = mxGetPr(pOUT[LC_IDX]);
  
  lagStateCorr();
  if(nINP==2)
    mxFree(intervals);
  

}
