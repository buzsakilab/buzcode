/* Restrict_idx_align.c: MEX file that deals with the align case of
   restrict: three possible behaviors are possible for this function:
   'next': for each key, the first element equal of greater of
   each key is returned
   'closest': for each key, the element closest to the key is returned
   'equal': if an element exactly equal to each key is not present,
   an error is returned 

   call like this

   ix = Restrict_idx_align(Range_tsd, st, align_type);

   Range_tsd is the Range of the tsd input object, st is the
   restricting alignment array, 
   align_type is a string that can be one of 'closest', 'next', 'prev',
   'equal'

   returns the 1=based indices of the restricted tsd in ix

   private function in the @tsd class 

   copyright (c) 2004 Francesco P. Battaglia
   This software is released under the GNU GPL
   www.gnu.org/copyleft/gpl.html
*/

   


#include <mex.h>
#include <math.h>
#include <matrix.h>
#include <string.h>




int binsearch_next(double *data, double key, int min_idx, int max_idx)
{
  
  int mid, start = min_idx, end = max_idx-1;
  while (start < (end-1))
    {
      mid = floor((start + end)/2);
      if ((key) == data[mid])
	start = end = mid;
      if ((key) < data[mid])
	end = mid;
      if ((key) > data[mid])
	start = mid;
    }

  if(key == data[start])
    return start;
  
  return end;
  
}

int binsearch_prev(double *data, double key, int min_idx, int max_idx)
{
  
  int mid, start = min_idx, end = max_idx-1;
  if (key > data[end])
    return end+1;
  
  
  while (start < (end-1))
    {
      mid = floor((start + end)/2);
      if ((key) == data[mid])
	start = end = mid;
      if ((key) < data[mid])
	end = mid;
      if ((key) > data[mid])
	start = mid;
    }

/*   if(key == data[end]) */
/*     return end+1; */
  
  return start;
  
}

int binsearch_closest(double *data, double key, int min_idx, int max_idx)
{
  
  int mid, start = min_idx, end = max_idx-1, result;
  while (start < (end-1))
    {
      mid = floor((start + end)/2);
      if ((key) == data[mid])
	start = end = mid;
      if ((key) < data[mid])
	end = mid;
      if ((key) > data[mid])
	start = mid;
    }
  if (((key) - data[start]) <= (data[end] - (key)))
    result = start;
  else 
    result = end ;
  return result;
  
}


typedef int (*binsearch_t)(double*, double, int, int);

typedef enum {NEXT, PREV, CLOSEST, EQUAL} align_type_t;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{
 int n_points, n_keys;
  double *t, *t0, *ix;
  int i, j;
  char align_type_str[20];
  
  align_type_t align_type = CLOSEST;
  
  

  if(nrhs != 3)
    mexErrMsgTxt("Requires three inputs");
  
  if(nlhs != 1)
    mexErrMsgTxt("Requires one output");
  
  
  /* unpack inputs */

  /* input 0 is t, time to restrict*/

  if(mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("t must be a row or column vector");
  
  n_points = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  
  t = mxGetPr(prhs[0]);
  
  /* input 1 is t0 restricting keys  */

  if(mxGetM(prhs[1]) != 1 && mxGetN(prhs[1]) != 1)
    mexErrMsgTxt("t0 must be a row or column vector");
  
  n_keys = mxGetM(prhs[1]) * mxGetN(prhs[1]);

  t0 = mxGetPr(prhs[1]);
  

  /* input 2 is the align_type string */

  if(!mxIsChar(prhs[2])) 
    mexErrMsgTxt("align_type must be string");
  
  mxGetString(prhs[2], align_type_str, 19);
  


  /* prepare output */

  plhs[0] = mxCreateNumericMatrix(n_keys, 1, mxDOUBLE_CLASS, mxREAL);
      
  ix = mxGetPr(plhs[0]);
  
  
  if(!strcmp(align_type_str, "next") ) 
    align_type = NEXT;
  else if(!strcmp(align_type_str, "prev"))
    align_type = PREV;
  else if(!strcmp(align_type_str, "closest"))
    align_type = CLOSEST;
  else if(!strcmp(align_type_str, "equal"))
    align_type = EQUAL;
  else
    mexErrMsgTxt("Unrecognized option");
  


  
  j = binsearch_prev(t, t0[0], 0, n_points-1);
/*   mexPrintf("t[j] = %g t0[0] = %g\n", t[j], t0[j]); */
  
  j = (j > 0 ? j : 0);
  


  switch(align_type) 
    {
    case PREV:
      for(i = 0; i < n_keys; i++)
	{
	  if(mxIsFinite(t0[i]))
	    {
	      while((j < n_points) && (t[j] <= t0[i])  )
		j++;
	      if(j==0)
		j = 1;
	      
	      ix[i] = (double)j;
	    }
	  else
	    {
	      ix[i] = mxGetNaN();
	    }
	}
      break;
    case NEXT:
      
      
      for(i = 0; i < n_keys; i++)
	{
	  if(mxIsFinite(t0[i]))
	    {
	      while((t[j] < t0[i]) && (j < n_points-1))
		j++;
	      ix[i] = (double)j;
	      ix[i]++;
	      
/* 	      ix[i] = (ix[i] < n_points-1 ? ix[i]+1 : ix[i]); */
	      
	    }
	  else
	    {
	      ix[i] = mxGetNaN();
	    }
	}
      break;
      
    case CLOSEST:
    case EQUAL:
      for(i = 0; i < n_keys; i++)
	{
	  if(mxIsFinite(t0[i]))
	    {
	      while((t[j] < t0[i]) && (j < n_points-1))
		j++;
	      ix[i] = (double)j;
	      ix[i] = ((t0[i] - t[j-1]) <= (t[j] - t0[i]) ?
		       ix[i] :
		       ix[i] + 1);
	      if (ix[i] == 0) ix[i] = 1;
	      
	      
	    }
	  else
	    {
	      ix[i] = mxGetNaN();
	    }
	}
      break;
    }
  

  if(align_type == EQUAL)
    {
      for(i = 0; i < n_keys; i++)
	{
	  if(t[(int)(ix[i])-1] != t0[i])
	    mexErrMsgTxt("Alignment failure");
	}
    }
  
}

  



  
