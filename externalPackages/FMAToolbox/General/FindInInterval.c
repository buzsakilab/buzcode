/***************************************************************************
                            FindInInterval.c
                            --------------------
    begin                : June 2008
    copyright            : (C) 2008 by Michaël Zugaro
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

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double t0,t1,start;
  double *x,*indices;
  int i,m;

	if (nrhs < 2 || nrhs > 3) mexErrMsgTxt("FindInInterval: wrong number of input arguments.");
	if (nlhs > 1) mexErrMsgTxt("FindInInterval: wrong number of output arguments.");

	x = mxGetPr(prhs[0]);
	m = mxGetM(prhs[0]);
	if ( m == 1 ) m = mxGetN(prhs[0]);

	t0 = mxGetPr(prhs[1])[0];
	t1 = mxGetPr(prhs[1])[1];
	if (nrhs == 2) start = 0; else start = *mxGetPr(prhs[2])-1;

	plhs[0] = mxCreateDoubleMatrix(2,1,mxREAL);
	indices = mxGetPr(plhs[0]);
	indices[0] = indices[1] = 0;

	for ( i = start ; i < m && x[i] < t0 ; i++ );
	if ( i >= m )
	{
		mxDestroyArray(plhs[0]);
		plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
	}
	else
	{
		indices[0] = i+1;

		for ( ; i < m && x[i] <= t1 ; i++ );
		indices[1] = i;
	}
}
