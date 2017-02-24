/***************************************************************************
                          FindField.c  -  description
                             -------------------
    begin                : Wed May 22 2002
    copyright            : (C) 2002-2012 by Michaël Zugaro
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

int      xMax,yMax;
double   threshold,*firingMap;
bool     circX,circY;

void FindField(double x,double y,bool *firingField,bool *visited)
{
	/*mexPrintf("before x %f y %f\n",x,y);*/

	if ( x < 1 )
		if ( circX ) x = xMax;
		else return;
	if ( x > xMax )
		if ( circX ) x = 1;
		else return;
	if ( y < 1 )
		if ( circY ) y = yMax;
		else return;
	if ( y > yMax )
		if ( circY ) y = 1;
		else return;

	int i = (int) (x-1)*yMax+(y-1);

	if ( visited[i] ) return;

	visited[i] = true;
	if (firingMap[i] >= threshold)
	{
		firingField[i] = true;
		FindField(x+1,y,firingField,visited);
		FindField(x-1,y,firingField,visited);
		FindField(x,y+1,firingField,visited);
		FindField(x,y-1,firingField,visited);
	}
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	double   x,y;
	bool     *firingField,*visited;
	int      i,j;

	if ( nrhs != 6 )
		mexErrMsgTxt("FindField: wrong number of input arguments.");
	if ( nlhs > 1 )
		mexErrMsgTxt("FindField: wrong number of output arguments.");

	firingMap = mxGetPr(prhs[0]);
	yMax = mxGetM(prhs[0]);
	xMax = mxGetN(prhs[0]);
	x = *mxGetPr(prhs[1]);
	y = *mxGetPr(prhs[2]);
	threshold = *mxGetPr(prhs[3]);
	circX = *mxGetLogicals(prhs[4]);
	circY = *mxGetLogicals(prhs[5]);
	plhs[0] = mxCreateLogicalMatrix(yMax,xMax);
	firingField = mxGetLogicals(plhs[0]);

	visited = (bool*) mxCalloc((xMax+1)*(yMax+1),sizeof(bool));
	for ( i = 0;i < xMax*yMax ; ++i ) visited[i] = false;

	FindField(x,y,firingField,visited);

	/*mxFree(visited);*/
}
