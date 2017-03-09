/***************************************************************************
    copyright            : (C) 2013 by Michaël Zugaro
    email                : michael.zugaro@college-de-france.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "mex.h"
#include <string.h>

int m,n;
double *status;

void Contiguous(int i,int j)
{
	/* If out of bounds, ignore */
	if ( i < 0 || i >= m || j < 0 || j > n ) return;

	/* If this item is equal to 1, mark it so it will not be checked twice */
	if ( status[i+m*j] == 1 )
	{
		/* Mark as belonging to the patch by setting to 2 */
		status[i+m*j] = 2;
		/*  Examine neighbours */
		Contiguous(i-1,j);
		Contiguous(i+1,j);
		Contiguous(i,j-1);
		Contiguous(i,j+1);
	}
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	int i,j;
	double *s;

	/* Check input parameters */
	if ( nlhs == 0 ) return;
	if ( nrhs != 3 ) mexErrMsgTxt("Incorrect number of input parameters.");
	if ( mxGetClassID(prhs[0]) != mxDOUBLE_CLASS || mxGetClassID(prhs[1]) != mxDOUBLE_CLASS || mxGetClassID(prhs[2]) != mxDOUBLE_CLASS ) mexErrMsgTxt("Incorrect parameters.");
	
	/* Get input parameter values */
	s = mxGetPr(prhs[0]);
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	i = mxGetScalar(prhs[1]) - 1;
	j = mxGetScalar(prhs[2]) - 1;

	/* Create output matrix */
	plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
	status = mxGetPr(plhs[0]);
 	if ( status == NULL ) mexErrMsgTxt("Not enough memory.");
	
	/* Copy input matrix and find patch */
	memcpy(status,s,m*n*mxGetElementSize(prhs[0]));
	Contiguous(i,j);
}

