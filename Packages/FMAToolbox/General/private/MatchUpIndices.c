/***************************************************************************
                          MatchUpIndices.c  -  description
                             -------------------
    copyright            : (C) 2012 by Michaël Zugaro
    email                : mzugaro@andromeda.rutgers.edu
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
	double   *list1,*list2,*matched;
	int      i,j,m,n;

	if ( nrhs != 2 )
		mexErrMsgTxt("MatchUpIndices: wrong number of input arguments.");
	if ( nlhs > 1 )
		mexErrMsgTxt("MatchUpIndices: wrong number of output arguments.");

	list1 = mxGetPr(prhs[0]);
	m = mxGetM(prhs[0]);
	if ( m == 1 ) m = mxGetN(prhs[0]);
	list2 = mxGetPr(prhs[1]);
	n = mxGetM(prhs[1]);
	if ( n == 1 ) n = mxGetN(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);
	matched = mxGetPr(plhs[0]);

	j = 0;
	for ( i = 0 ; i < m ; ++i )
	{
		while ( list1[i] > list2[j] && j < n ) ++j;
		matched[i] = j+1;
	}

/*	for ( i = 0 ; i < m ; ++i )
	{
		mexPrintf("list1[%d] = %f\n",i,list1[i]);
		while ( list1[i] > list2[j] && j < n ) { mexPrintf(" %f > list1[%d] = %f\n",i,list1[i],j,list2[j]);++j; }
		mexPrintf(" match: list2[%d] = %f\n",j,list2[j]);
		mexPrintf(" matched[%d] = %d\n",i,j+1);
		matched[i] = j+1;
	}*/
}
