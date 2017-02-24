 /*-----------------------------------
 * convex hull
 * MEX file
 * 
 * ADR 1999
 * 
 * input: 
 *         x - n x 1 array of x coordinates
 *         y - n x 1 array of y coordinates
 *
 * output: 
 *         k - 1d array of indices of which points 
 *             are on the convex hull
 *
 * algorithm: graham's scan.  
 *
 * version 1.0
 * status: promoted
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

#define MATH_PI 3.14159265358979323846

/*-----------------------------------
 * DEBUG consts 
 ----------------------------------*/
#define mexDEBUG(str) 
/*mexPrintf(str)*/
#define mexDEBUG1(str, v) 
/*mexPrintf(str,v)*/

/*-----------------------------------
 * point struct
 ----------------------------------*/
struct point {
	double rho;
	int key;
}; 

/*-----------------------------------
 * min point --- returns index of minimum point
 ----------------------------------*/
int MinPoint(double *a, double *b, int np)
{
	int ip; /* counter */
	int minIndex;
	double minA = 1e300, minB = 1e300;

	for (ip = 0; ip < np; ip++)
	{
		if (a[ip] > minA)
			continue;
		if (a[ip] == minA && b[ip] > minB)
			continue;
		minA = a[ip];
		minB = b[ip];
		minIndex = ip;
	};	/* for all points */					
	return minIndex;	
}	

/*----------------------------------
 * compare --- for qsort
 ----------------------------------*/
int compare(const void *a, const void *b)
{
	if (((struct point *)a)->rho > ((struct point *) b)->rho)
		return 1;
	if (((struct point *)a)->rho < ((struct point *) b)->rho)
		return -1;
	return 0;
}

/*-----------------------------------
 * convex hull function
 *
 * inputs: 
 *         double *x, *y - arrays of coordinates
 *         int np - number of points in x,y
 *         int **k - [output] array of indices of points on hull
 *         int *nk - [output] size of k
 *
 * outputs: 
 *         none
 -----------------------------------*/ 
void ConvexHull(
			double *x, double *y, int np,
			int **k, int *nk)
{
	int startPoint;
	struct point *p;
	int ip;				/* counter */
	double rho0, rho1, drho;
	int m1, m2;
	int *kL; int nkL;
    double d0,dplus,dminus;

	/* set up k */
	kL = (int *) malloc(sizeof(int) * np);
	nkL = 0;
	mexDEBUG("k allocated\n");

	/* find min y, min x point */
	startPoint = MinPoint(y,x,np);
	mexDEBUG1("startPoint = %d\n", startPoint);
	
	/* find angles */
	p = (struct point *) malloc(sizeof(struct point) * (np-1));
	for (ip = 0; ip < np; ip++)
		if (ip < startPoint)
			{
			 p[ip].rho = atan2(y[ip] - y[startPoint], x[ip] - x[startPoint]);
			 p[ip].key = ip;
			}
		else if (ip > startPoint)
			{
			 p[ip-1].rho = atan2(y[ip] - y[startPoint], x[ip] - x[startPoint]);			
			 p[ip-1].key = ip;
			};
	mexDEBUG("p constructed.\n");

	/*for (ip = 0; ip<np-1; ip++)
		mexDEBUG1("ip.rho = %f \n", p[ip].rho);*/

	/* sort points by angles */
	qsort((void *)p, np-1, sizeof(struct point), compare);
	mexDEBUG("p sorted.\n");

	/*
     for (ip = 0; ip<np-1; ip++)
        if (p[ip].rho == p[ip+1].rho) { 
            mexDEBUG1("%d", ip); mexDEBUG1("<%.0f,", x[ip]); mexDEBUG1("%.0f>",y[ip]);
            mexDEBUG1("rho = %f \n", p[ip].rho);
        };
     */
        
	/* do scan */
	kL[nkL++] = startPoint; mexDEBUG1("start point: key = %d \n", startPoint+1);
	kL[nkL++] = p[0].key;   mexDEBUG1("ip = %d ", 0);mexDEBUG1("key = %d ", p[0].key+1);mexDEBUG1("rho = %f\n", p[0].rho); 
	kL[nkL++] = p[1].key;   mexDEBUG1("ip = %d ", 1);mexDEBUG1("key = %d ", p[1].key+1);mexDEBUG1("rho = %f\n", p[1].rho);
	for (ip = 2; ip<np-1; ip++)
	{
		do {                                
			m1 = kL[nkL-1];
			m2 = kL[nkL-2];
			rho0 = atan2(y[m1]         - y[m2], x[m1]         - x[m2]);	/* angle m-2 -> m-1  */
			rho1 = atan2(y[p[ip].key]  - y[m2], x[p[ip].key]  - x[m2]); /* angle m-2 -> p */
			drho = rho0 - rho1;             

            if (drho <= -MATH_PI || (drho >= 0 && drho <= MATH_PI))
                nkL--;
		} while (drho <= -MATH_PI || (drho >= 0 && drho <= MATH_PI));		

		kL[nkL++] = p[ip].key;
	} /* for all points */        
*k = kL;
*nk = nkL; 
}


/*-----------------------------------
 * main MEX function
 -----------------------------------*/
 void mexFunction(
			int nOUT, mxArray *pOUT[],
			int nINP, const mxArray *pINP[])
{
  int nx, ny, nk;	/* number of points */
  double *x, *y;	/* input arrays */
  double *result;	/* actual output */
  int *k;			/* output */ 
  int ik;			/* counter */
  
  /* check number of arguments: expects 2 inputs, 1 output */
  if (nINP != 2)
    mexErrMsgTxt("Call with x,y as inputs.");
  if (nOUT != 1)
    mexErrMsgTxt("Requires one output.");

  /* check validity of inputs */
  nx = mxGetM(pINP[0]) * mxGetN(pINP[0]);
  ny = mxGetM(pINP[1]) * mxGetN(pINP[1]);   
  if (nx != ny)
    mexErrMsgTxt("Inputs must be the same size.");

  /* unpack inputs */
  x = (double *)mxGetPr(pINP[0]);
  y = (double *)mxGetPr(pINP[1]);

  mexDEBUG("All inputs unpacked.\n");

  /* convex hull */
  ConvexHull(x, y, nx, &k, &nk);

  /* pack outputs */
  pOUT[0] = mxCreateDoubleMatrix(nk+1, 1, mxREAL);
  result = (double *) mxGetPr(pOUT[0]);
  for (ik=0; ik<nk; ik++)
	result[ik] = k[ik]+1;
  result[ik] = k[0]+1;
}
