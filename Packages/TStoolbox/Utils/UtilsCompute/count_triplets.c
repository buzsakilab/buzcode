/* count_triplets.c: MEX file, "engine" for the spike_triplets
   function, counts the triples of spikes  occurring within a time
   span
   INPUTS:
   t: a n_spikes x 1 array of timestamps
   idxs: a n_spikes x 1 array of cell indices corresponding to the
   timestamps
   n_cells: scalar indicating the number of cells
   time_span 
   OUTPUTS:
   triplets_freq: a n_comb(N_cells, 3) X 6 matrix of frequencies, rows
   a are cell triplets columns possible orderings, according to 
   1: 1 2 3
   2: 1 3 2 
   3: 2 1 3 
   4: 2 3 1 
   5: 3 1 2 
   6: 3 2 1 
   triples_table: a
   n_comb(N_cells, 3) X 3 matrix, a lookup table for the row dimension in
   triplets_freq, i-th row is the set of cell indices for the i-th
   triple, in increasing order
   
   this is really meant to be called by spike_triplets m
   function. Many assumptions are made that are not assured to hold in
   standalone use

   code by batta 2002, supporting a collaboration with A.Treves

*/


#include <math.h>
#include "mex.h"
#include "matrix.h"





/* this implements the lookup table for the permutations */

inline int order_code(int *ii)
{
  int cd;
  
  int i1, i2, i3;
  
  i1 = ii[0];
  i2 = ii[1];
  i3 = ii[2];
  

  if(i1 == i2 || i1 == i3 || i2 == i3)
    return -1;
  
  if(i1 < i2 && i1 < i3)
    cd = 0;
  else if(i1 < i3 || i1 < i2)
    cd = 2;
  else
    cd = 4;
  
  if(i2 > i3)
    cd++;
  
  return cd;
}


inline void order_tuple(int *tp, int n) /* order a n-tuple */
{
  int tmp;
  

  if (n != 3)
    mexErrMsgTxt("order_tuple specialized for triplets!");
  
  if(tp[0] > tp[1])
    {
      tmp = tp[1];
      tp[1] = tp[0];
      tp[0] = tmp;
    }
  
  if(tp[1] > tp[2])
    if(tp[0] > tp[2])
      {
	tmp = tp[2];
	tp[2] = tp[1];
	tp[1] = tp[0];
	tp[0] = tmp;
      }
    else
      {
	tmp = tp[2];
	tp[2] = tp[1];
	tp[1] = tmp;
      }
}

  
int n_cells; /* number of cells */

inline int lookup_idx(int i, int j, int k)
{
  return i + n_cells * j + n_cells * n_cells * k;
}


  

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{
  
  double *tstamp; /* spike timestamps (in 1/10000 sec) */
  int *idxs;
  int cidx[3];
  
  int midx[2];

  int n_spikes; /* total number of spike */

  int n_combs; /* number of combinations */
  
  int i,j,k, i_last;
  

  int *triple_lookup;
  double *triplets_freq;
  double *triples_table;
  int *triples_lookup; /* direct lookup table */

  
  double time_span;
  

  enum {TSTAMP=0, IDXS, NCELLS, TSPAN
  } in_arg_order;
  
  enum {TRIPLETS_FREQ, TRIPLES_TABLE 
  } out_arg_order;
  
  
  if(nrhs != 4 )
    mexErrMsgTxt("Four inputs required");
  

  if(nlhs != 2 )
    mexErrMsgTxt("Two outputs required");

  

  /* MEX related stuff, not really related with algorithms */

  /* extracting timestamps */

  
  n_spikes = mxGetM(prhs[TSTAMP]);
  if(mxGetN(prhs[TSTAMP]) != 1) 
    mexErrMsgTxt("Timestamp must be a column vector");
  
  tstamp = mxGetPr(prhs[TSTAMP]);
  
  /* extracting cell indices */

  if(mxGetM(prhs[IDXS]) != n_spikes)
    mexErrMsgTxt("indices and timestamps must have same length");
  
  {
    int i;
    double *dixs;
    
    dixs = mxGetPr(prhs[IDXS]);
    
    idxs = (int *) mxCalloc(n_spikes, sizeof(int));
    for(i = 0; i < n_spikes; i++)
      idxs[i] = (int)(dixs[i])-1;
  }
  


  /* extracting n_cells */

  if(mxGetM(prhs[NCELLS]) * mxGetN(prhs[NCELLS]) != 1)
    mexErrMsgTxt("n_cells must be a scalar");
  n_cells = (int)mxGetScalar(prhs[NCELLS]);
  

  /* extracting time_span*/

  if(mxGetM(prhs[TSPAN]) * mxGetN(prhs[TSPAN]) != 1)
    mexErrMsgTxt("time_span must be a scalar");
  time_span = (int)mxGetScalar(prhs[TSPAN]);
  


  /* create outputs */ 
  n_combs = n_cells * (n_cells-1) * (n_cells-2) / 6;
  
  plhs[TRIPLETS_FREQ] = mxCreateNumericMatrix(n_combs, 6,
					      mxDOUBLE_CLASS,
					      mxREAL);

  triplets_freq = mxGetPr(plhs[TRIPLETS_FREQ]);
  

  plhs[TRIPLES_TABLE] = mxCreateNumericMatrix(n_combs, 3,
					      mxDOUBLE_CLASS,
					      mxREAL);
  
  triples_table = mxGetPr(plhs[TRIPLES_TABLE]);
  

  /* create lookup tables direct and inverse */

  triples_lookup = (int *) mxCalloc(n_cells*n_cells*n_cells, sizeof(int));

  

  

  
  {
    int ix = 0;
    
    for(i = 0; i < n_cells; i++)
      for(j = i+1; j < n_cells; j++)
	for(k = j+1; k < n_cells; k++)
	  {
	    /*  	  int ix = lookup3[i] + (n_cells - i) *(j - i - 1)  */
	    /*  	    - (j - i - 1) * (j -i) /2  + (k-j-1); */
	    midx[0] = ix;
	    midx[1] = 0;
	    triples_table[mxCalcSingleSubscript(plhs[TRIPLES_TABLE], 2, midx)] =
	      i+1;
	    midx[1] = 1;
	    triples_table[mxCalcSingleSubscript(plhs[TRIPLES_TABLE], 2, midx)] =
	      j+1;
	    midx[1] = 2;
	    triples_table[mxCalcSingleSubscript(plhs[TRIPLES_TABLE], 2, midx)] =
	      k+1;
	    triples_lookup[lookup_idx(i,j,k)] = ix;
	    ix++;
	    
	  }
  }
  

  /* now complete the direct lookup table */

  for(i = 0; i < n_cells; i++)
    for(j = 0; j < n_cells; j++)
      for(k = 0; k < n_cells; k++)
	{
	  cidx[0] = i;
	  cidx[1] = j;
	  cidx[2] = k;
	  order_tuple(cidx, 3);
	  triples_lookup[lookup_idx(i,j,k)] = 
	    triples_lookup[lookup_idx(cidx[0], cidx[1], cidx[2])];
	}
  

  /* and now compute the frequency table */
  i_last = 0;
  
  for(i = 0; i < n_spikes-2; i++)
    {
      
      double t1 = tstamp[i];
      double tl = t1 + time_span;
      int ix = idxs[i];
      int jx, kx, cx, px;
 
     if(! (i % 10000))
	mexPrintf("i: %d\n", i);
      
      while(i_last < n_spikes && tstamp[i_last] < tl) i_last++;
      
      for(j = i+1; j < i_last; j++)
	if(idxs[j] != ix)
	  {
	    jx = idxs[j];
	    for(k = j+1; k < i_last; k++)
	      if(idxs[k] != ix && idxs[k] != jx)
		{
		  kx = idxs[k];
		  cidx[0] = ix;
		  cidx[1] = jx;
		  cidx[2] = kx;
		  px = order_code(cidx);
		  
		  order_tuple(cidx, 3);
		  cx = triples_lookup[lookup_idx(cidx[0], cidx[1], cidx[2])];
		  midx[0] = cx;
		  midx[1] = px;
		  triplets_freq[mxCalcSingleSubscript(plhs[TRIPLETS_FREQ], 2, midx)]++;
		}
	  }
    }
  
		
	  
      


  mxFree((void *)triples_lookup);

  

  
}
