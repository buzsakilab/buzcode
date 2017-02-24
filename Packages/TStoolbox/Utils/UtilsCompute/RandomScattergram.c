 /*-----------------------------------
 * Random Scattergram: takes a 2-D density as an input and returns a set
 * of random points distributed with that density 
 * MEX file
 * 
 * batta 1999
 * 
 * input: density(ndim multidimensional array) -- histogram
 *        (normalized)
 *        bins (2 x 1) -- number of bins per each dimension
 *        mins (2 x 1) -- min values for each dimension
 *        maxes (2 x 1) -- max values for each dimension
 *        npoints -- number of points to output
 *        threshold -- threshold for displaying one dot in one bin
 * output: points (2 x npoints) 
 *
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <string.h>
#include <stdlib.h>
char MsgBuffer[100];


int RandomScattergram(double *, double *, double *, double *, int, double , double **);

/*--------------------
 * gateway function
 --------------------*/

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  int nD = 2, npoints;
  double *S, *bins, *mins, *maxes, *H, *points, threshold;

  /* check number of arguments: expects 6 inputs, 1 output */
  if (nINP != 6)
    mexErrMsgTxt("Call with samples, bins, mins, maxes as inputs");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs */
  nD = mxGetNumberOfDimensions(pINP[0]);
  if(nD != 2)  
    mexErrMsgTxt("Wrong dimensionality of density");  
  if ((mxGetM(pINP[1]) != nD) || (mxGetN(pINP[1]) != 1))
      mexErrMsgTxt("Wrong dimensionality of bins.");
  if ((mxGetM(pINP[2]) != nD) || (mxGetN(pINP[2]) != 1))
    mexErrMsgTxt("Wrong dimensionality of mins.");
  if ((mxGetM(pINP[3]) != nD) || (mxGetN(pINP[3]) != 1))
    mexErrMsgTxt("Wrong dimensionality of maxes.");
  
  /* unpack inputs */
  H = mxGetPr(pINP[0]);
  bins  = mxGetPr(pINP[1]);
  mins  = mxGetPr(pINP[2]);
  maxes = mxGetPr(pINP[3]);
  npoints = (int)mxGetScalar(pINP[4]);
  threshold = mxGetScalar(pINP[5]);
  /* create outputs */
  

  

  /* calculate */
  npoints = RandomScattergram(H,  bins, mins, maxes, npoints,
			      threshold, &points);
  pOUT[0] = mxCreateDoubleMatrix(npoints, 2, mxREAL);
  memcpy(mxGetPr(pOUT[0]), points, 2 * npoints * sizeof(double));
  

}

  


int RandomScattergram(double *H, double *bins, double *mins, double *maxes, int
		  npoints, double threshold, double **points)
{
  int i, j, l;
  int nmaxpoints, np, ntotpoints = 0;
  double *px;
  double *py;
  double *pt;
  
  double qx = (maxes[0] - mins[0]) / bins[0];
  double qy = (maxes[1] - mins[1]) / bins[1];
  




  nmaxpoints = npoints + bins[0] * bins[1];
  px = (double *)mxCalloc(2*nmaxpoints, sizeof(double));
  py = (double *)mxCalloc(2*nmaxpoints, sizeof(double));
  for(i = 0; i < bins[0]; i++)
    for(j = 0; j < bins[1]; j++)
      {
	np = (int)floor(npoints * H[i + (int)bins[0] * j]);
	if(np == 0 && H[i + (int)bins[0] * j] > threshold)
	  np = 1;
	for(l = 0;l < np; l++)
	  {
	    px[l + ntotpoints] = mins[0] + (i + drand48()) * qx;
	    py[l + ntotpoints] = mins[1] + (j + drand48()) * qy;
	    
	  }
	ntotpoints += np;
      }
  
  pt = (double *)mxMalloc(2 * ntotpoints * sizeof(double));
  memcpy(pt, px, ntotpoints*sizeof(double));
  memcpy(&(pt[ntotpoints]),py, ntotpoints*sizeof(double));
  *points = pt;
  
  return(ntotpoints);
  
}

	    
  
