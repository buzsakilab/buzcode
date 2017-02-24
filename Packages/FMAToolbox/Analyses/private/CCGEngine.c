/***************************************************************************
                           CCGEngine.c  -  description
                             -------------------
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
#include <math.h>

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	int nTimes,nIDs,halfBins,nBins,i,size,index,current,previous,next,id1,id2,bin;
	double *times,*id,*ccg;
	double binSize,maxDt;
	double t1,t2;

	/* Check parameters */
	if ( nlhs == 0 ) return;
	if ( nrhs != 4 ) mxErrMsgTxt("Incorrect number of input parameters.");
	if ( mxGetClassID(prhs[0]) != mxDOUBLE_CLASS || mxGetClassID(prhs[1]) != mxDOUBLE_CLASS || mxGetClassID(prhs[2]) != mxDOUBLE_CLASS || mxGetClassID(prhs[3]) != mxDOUBLE_CLASS ) mxErrMsgTxt("Incorrect parameters.");
	times = mxGetPr(prhs[0]);
	nTimes = mxGetM(prhs[0]);
	id = mxGetPr(prhs[1]);
	nIDs = mxGetM(prhs[1]);
	binSize = mxGetScalar(prhs[2]);
	halfBins = mxGetScalar(prhs[3]);
	if ( nIDs != nTimes ) mxErrMsgTxt("Times and IDs have different lengths.");

	nBins = 2*halfBins + 1;
	maxDt = binSize * (halfBins + 0.5);

	/* Number of different IDs */
	nIDs = 0;
	for ( i = 0 ; i < nTimes ; ++i )
		if ( id[i] > nIDs ) nIDs = id[i];

	/* Output parameter */
	size = nIDs * nIDs * nBins;
	plhs[0] = mxCreateNumericMatrix(size,1,mxDOUBLE_CLASS,mxREAL);
	ccg = mxGetPr(plhs[0]);
	if ( ccg == NULL ) mxErrMsgTxt("Not enought memory.");

	/* Loop through events */
	for ( current = 0; current < nTimes ; ++current )
	{
		id1 = id[current];
		t1 = times[current];

		/* Previous events */
		for ( previous = current-1 ; previous >= 0; --previous )
		{
			id2 = id[previous];
			t2 = times[previous];

			/* Make sure the event falls in the time window */
			if ( t1 - t2 > maxDt ) break;

			/* Add an event pair in the CCG */
			bin = halfBins + floor(0.5+(t2-t1)/binSize);
			index = nBins * nIDs * (id1-1) + nBins * (id2-1) + bin;
			if ( index < 0 || index >= size) mxErrMsgTxt("Index out of bounds");
			ccg[index]++;
		}

		/* Next events */
		for ( next = current+1 ; next < nTimes ; next++ )
		{
			id2 = id[next];
			t2 = times[next];

			/* Make sure the event falls in the time window */
			if( t2 - t1 >= maxDt ) break;

			/* Add an event pair in the CCG */
			bin = halfBins + floor(0.5+(t2-t1)/binSize);
			index = nBins * nIDs * (id1-1) + nBins * (id2-1) + bin;
			if ( index < 0 || index >= size) mxErrMsgTxt("Index out of bounds");
			ccg[index]++;
		}
	}
}
