/* do_intersect.c MEX file, engine for the interval_set intersection
   function 
   INPUTS:
   tstart1, tend1, tstart2, tend2, ... tstartN, tendN, start and stop
   arrays for the objects to intersect
   OUTPUTs Ostart, Oend
   batta 2001 
   starting version
*/

/* modified fpbatta 2004/05/23 to allow intervals in each set to be
   overlapping */ 


#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <stdlib.h>

struct int_bound 
{
  double t;
  int flag;
};



int compare_int_bound(const struct int_bound *pa, const struct int_bound *pb)
{
  double a, b;
  
  a = pa->t;
  b = pb->t;
  if(a == b) 
    return (pa->flag < pb->flag ? -1 : 1);
  
  return(((a < b) ? -1 : 1));
  
}


void mexFunction(
		 int nOUT, mxArray *pOUT[],
		 int nINP, const mxArray *pINP[])
{
  int i, j, k;
  int n_intervals, n_sets;
  
  double *pr, *pr_start, *pr_stop;
  struct int_bound *ib;
  int *in_out, status, st_temp, n_intervals_out;
  

  /* checking inputs and outputs */ 
  if (nINP < 4 || nINP % 2 != 0)
    mexErrMsgTxt("Call with an even number of inputs (at least four)");
  if (nOUT != 2)
    mexErrMsgTxt("Call with two outputs.");
  
  n_sets = nINP/2;
  

  /* count elements to allocate memory */
  
  n_intervals = 0;
  for(i = 0; i < nINP; i += 2)
    n_intervals += mxGetM(pINP[i]);
      

  /* allocate memory */ 
  ib = (struct int_bound *)mxCalloc(n_intervals * 2, sizeof(struct
							    int_bound));
  k = 0;

  /* enter data from inputs args the t field in the int_bound struct
     contain the interval boundary timestamp, flag contains i for the
     beginning of the i-th intervalSet, and -i for the end of the i-th
     intervalSet */
  for(i = 0; i < nINP; i += 2)
    {
      pr = mxGetPr(pINP[i]);
      
      for(j = 0; j < mxGetM(pINP[i]); j++)
	{
	  ib[k].t = pr[j];
	  ib[k].flag = i/2+1;
	  k++;
	}
      
      pr = mxGetPr(pINP[i+1]);
      
      for(j = 0; j < mxGetM(pINP[i+1]); j++)
	{
	  ib[k].t = pr[j];
	  ib[k].flag = -(i/2+1);
	  k++;
	}

    }
  

  
  qsort(ib, n_intervals*2, sizeof(struct int_bound),
	compare_int_bound);
  
  /* prepare variables for the algorithm */

  in_out = (int *)mxCalloc(n_sets+1, sizeof(int));
  status = 0;
  k = 0;


  /* first sweep to compute n_intervals_out */
  for(i = 0; i < n_intervals * 2; i++)
    {
      if(ib[i].flag > 0)
	in_out[ib[i].flag] += 1;
      else
	in_out[-(ib[i].flag)] -= 1;
      
      st_temp = 1;
      for(j = 1; j <= n_sets; j++)
	st_temp *= (in_out[j] > 0);
      
      if ((!status) && st_temp) /* entering interval */
	{
	  /* insert update code here */
	  k++;
	}
      else if(status && (!st_temp))
	{
	  /* insert update code here */
	}
      status = st_temp;
    }
  
  n_intervals_out = k;
  
  pOUT[0] = mxCreateDoubleMatrix(n_intervals_out, 1, mxREAL);
  pr_start = mxGetPr(pOUT[0]);
  pOUT[1] = mxCreateDoubleMatrix(n_intervals_out, 1, mxREAL);
  pr_stop = mxGetPr(pOUT[1]);
  
  /* second sweep, actually get the result */
    
  for(j = 1; j <= n_sets; j++)
    in_out[j] = 0;
  
  status = 0;
  k = 0;
  /* first sweep to compute n_intervals_out */
  for(i = 0; i < n_intervals * 2; i++)
    {
      if(ib[i].flag > 0)
	in_out[ib[i].flag] += 1;
      else
	{
	  in_out[-(ib[i].flag)] -= 1;
	}
      
      st_temp = 1;
      for(j = 1; j <= n_sets; j++)
	st_temp *= (in_out[j] > 0);
      
      if ((!status) && st_temp) /* entering interval */
	{
	  pr_start[k] = ib[i].t;
	  
	}
      else if(status && (!st_temp))
	{
	  pr_stop[k] = ib[i].t;
	  k++;
	}
      status = st_temp;
    }


  
}
