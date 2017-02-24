 /*-----------------------------------
  * in polygon
  * MEX file
  *
  * ADR 1999
  *
  * input: px - nP x 1
  *        py - nP x 1
  *        cx - nC x 1
  *        cy - nC x 1
  * output:
  *        b - boolean array of 1 if px(i),py(i) is
  *            in the polygon defined by cx,cy
  *
  * version 2.0
  * ADR 11 APR 2008 : added case when point was a corner
 %
% Status: PROMOTED (Release version)
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
 -----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
int nOUT, mxArray *pOUT[],
int nINP, const mxArray *pINP[])
{
    int nP;
    int nC;
    double *px, *py, *cx, *cy;
    double *result;
    
  /* check number of arguments: expects 4 inputs, 1 output */
    if (nINP != 4)
        mexErrMsgTxt("Call with px,py,cx,cy as inputs.");
    if (nOUT != 1)
        mexErrMsgTxt("Requires one output.");
    
  /* unpack inputs */
    nP = mxGetM(pINP[0]) * mxGetN(pINP[0]);
    px = (double *)mxGetPr(pINP[0]);
    py = (double *)mxGetPr(pINP[1]);
    
    nC = mxGetM(pINP[2]) * mxGetN(pINP[2]);
    cx = (double *)mxGetPr(pINP[2]);
    cy = (double *)mxGetPr(pINP[3]);
    
  /* pack outputs */
    pOUT[0] = mxCreateDoubleMatrix(nP, 1, mxREAL);
    result = (double *) mxGetPr(pOUT[0]);
    
  /* calculate */
    {
        int iP, iC;
        double xmin, xmax, ymin, ymax;
        double nIntersect;
        
    /* first find min and max x and y */
        xmin = xmax = cx[0]; ymin = ymax = cy[0];
        for (iC=0; iC<nC; iC++)
        {
            if (cx[iC] < xmin) xmin = cx[iC];
            if (cx[iC] > xmax) xmax = cx[iC];
            if (cy[iC] < ymin) ymin = cy[iC];
            if (cy[iC] > ymax) ymax = cy[iC];
        }
        
    /* for each point is it in the polygon or not */
        for (iP=0; iP<nP; iP++)
        {
    /* is it even in the rectangle defined by the polygon? */
            if (px[iP] < xmin) {result[iP] = 0; continue;};
            if (px[iP] > xmax) {result[iP] = 0; continue;};
            if (py[iP] < ymin) {result[iP] = 0; continue;};
            if (py[iP] > ymax) {result[iP] = 0; continue;};
            
    /* count the number of intersections */
            nIntersect = 0;
            for (iC=0; iC<(nC-1); iC++)
            {
        /* does the line PQ intersect the line AB? */
                double ax = cx[iC];
                double ay = cy[iC];
                double bx = cx[iC+1];
                double by = cy[iC+1];
                double intersecty;
                
                if (ax == bx) continue;
        /* ensure order correct */
                if (ax > bx)
                {
                    double tmp;
                    tmp = ax; ax = bx; bx = tmp;
                    tmp = ay; ay = by; by = tmp;
                }
                
                if (ax > px[iP]) continue;
                if (bx < px[iP]) continue;
                intersecty = ay + (px[iP] - ax)/(bx - ax) * (by - ay);
		if (py[iP] > 0)
		  {
		    if ((intersecty < py[iP]) && ((ax==px[iP]) || (bx==px[iP])))
		      nIntersect = nIntersect+0.5;
		    else if (intersecty <= py[iP])
		      nIntersect++;
		  } else 
		  {
		    if ((intersecty > py[iP]) && ((ax==px[iP]) || (bx==px[iP])))
		      nIntersect = nIntersect+0.5;
		    else if (intersecty >= py[iP])
		      nIntersect++;
		  } 
            } /* for iC */
            
            if (nIntersect != floor(nIntersect))
                mexErrMsgTxt("Non-integer nIntersect.");
            
            result[iP] = (nIntersect - 2*floor(nIntersect/2));
            
    /* for each point, is it actually the corner itself? */
    /* ADDED ADR 2008 */
            for (iC = 0; iC<(nC-1); iC++)
                if (px[iP] == cx[iC] && py[iP] == cy[iC])
                    result[iP] = 1;
            
        } /* for iP */
        
        
    } /* calculate */
    
}





