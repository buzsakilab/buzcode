/* findIdx_c.c
 * find compound events in time series
 * MEX file 
 * 
 *
 * inputs:
 * C: a compoundEvents object (all sanity checks are assumed as
 * already performed  
 * t: a time series

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

*********************************************/
 





#include <stdlib.h>
#include <string.h>
#include <mex.h>
#include <math.h>
#include <matrix.h>

#define MEX_FILE 1

#ifdef MEX_FILE
#define calloc mxCalloc
#endif


typedef enum {  LESS_THAN, EQUAL, GREATER_THAN } CondSign_t;

int cEventFindFirstIndex(int eventSize, 
			 double *condTime,
			 CondSign_t *condSign,
			 int *condSkip, 
			 double *t,
			 int n_points,
			 int *events)
{
  int n_events;
  int i, i0, i1, j;
  int found = 0;
  int eventLength;
  double interval;
  
  eventLength = 1;
  for(i = 0; i < eventSize; i++)
    eventLength += eventSize;
  
  




  n_events = 0;
  
  i0 = 0;
  while(i0 <= n_points-eventLength)
    {
      found = 1;
      i1 = i0;
      
      for(j = 0; j < eventSize; j++)
	{
	  interval = t[i1 + condSkip[j]] - t[i1];
/* 	  mexPrintf("interval = %g\n", interval); */
	  
	  i1 = i1 + condSkip[j];
	  switch(condSign[j])
	    {
	    case LESS_THAN:
	      if(interval >= condTime[j])
		found = 0;
	      break;
	    case GREATER_THAN:
	      if(interval <= condTime[j])
		found = 0;
	      break;
	    case EQUAL:
	      if(interval != condTime[j])
		found = 0;
	      break;
	    default:
 	      mexErrMsgTxt("What the hell?"); 
	      


	    }
	  if(found == 0)
	    break;
	}
      
      if(found == 1)
	{
/* 	  mexPrintf("i0 %d\n", i0); */
	  events[n_events] = i0;
	  n_events++;
	}
      


      i0++;
    }
  

  return n_events;

  
}

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  double *t;
  double *condTime, *cs;
  int eventSize, n_points, n_events;
  int i;
  


  CondSign_t *condSign;
  
  int *condSkip, *events;
  


  
  if(nINP != 2)
    mexErrMsgTxt("Call with C, t as inputs");
  

  if(nOUT != 1)
    mexErrMsgTxt("Call with ix as output");
  
  /* extract inputs */
  
  condTime = mxGetPr(mxGetField(pINP[0], 0, "condTime"));
  eventSize = mxGetM(mxGetField(pINP[0], 0, "condTime")) *
    mxGetN(mxGetField(pINP[0], 0, "condTime"));
  


  cs = mxGetPr(mxGetField(pINP[0], 0, "condSign"));
  condSign = mxCalloc(eventSize, sizeof(CondSign_t));
  for(i = 0; i < eventSize; i++)
    {
      if (cs[i] < 0)
	condSign[i] = LESS_THAN;
      else if (cs[i] > 0)
	condSign[i] = GREATER_THAN;
      else
	condSign[i] = EQUAL;
    }
  
  cs = mxGetPr(mxGetField(pINP[0], 0, "condSkip"));

  condSkip = mxCalloc(eventSize, sizeof(int));
  

  for(i = 0; i < eventSize; i++)
    condSkip[i] = (int)cs[i];
  
  t = mxGetPr(pINP[1]);
  n_points = mxGetM(pINP[1]) * mxGetN(pINP[1]);
  
  events = mxCalloc(n_points, sizeof(int));
  


  n_events = cEventFindFirstIndex(eventSize, condTime, condSign, condSkip, t, n_points, events);
  
/*   mexPrintf("n_events = %d\n", n_events); */
  
  pOUT[0] = mxCreateDoubleMatrix(n_events, 1, mxREAL);
  cs = mxGetPr(pOUT[0]);
  
  for(i = 0; i < n_events; i++)
    cs[i] = (double)events[i];
			     
  


  mxFree(condSkip);
  mxFree(condSign);
  mxFree(events);
  




  
}





