/***************************************************************************
                            CountInIntervals.c
                            --------------------
    begin                : December 2009-2011
    copyright            : (C) 2009-2011 by Michaël Zugaro
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
  double *x,*intervals,*counts;
  int i,j,k,m,n;

	if (nrhs != 2) mexErrMsgTxt("Incorrect number of parameters (type 'help <a href=\"matlab:help CountInIntervals\">CountInIntervals</a>' for details).");

	/* There are m intervals... */
	intervals = mxGetPr(prhs[1]);
	if ( mxGetN(prhs[1]) != 2 ) mexErrMsgTxt("Incorrect intervals (type 'help <a href=\"matlab:help CountInIntervals\">CountInIntervals</a>' for details).");
	m = mxGetM(prhs[1]);

	/* ... and n values x to test */
	x = mxGetPr(prhs[0]);
	n = mxGetM(prhs[0]);
	if ( n == 1 ) n = mxGetN(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);
	if ( m == 0 || n == 0 ) return;
	counts = mxGetPr(plhs[0]);

	j = 0;
	/* Loop through intervals */
	for ( i = 0; i < m ; ++i )
	{
		/* Find first x[j] that belongs to interval i */
		while ( x[j] < intervals[i] && j < n ) ++j;
		if ( j > n ) break;
		/* How many of x[j], x[j+1]... belong to this interval? */
		k = j;
		while ( x[k] <= intervals[m+i] && k < n ) { ++counts[i]; ++k; }
	}
}
