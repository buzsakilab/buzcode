/* count_pairs.c: MEX file, "engine" for the spike_pairs
   function, counts the pairs of spikes  occurring within a time
   span, also outputs, in two cell arrays, the timestamps at which
   each ordered  pair occurred, and the locations where they occurred
   INPUTS:
   t: a n_spikes x 1 array of timestamps
   idxs: a n_spikes x 1 array of cell indices corresponding to the
   timestamps
   n_cells: scalar indicating the number of cells
   time_span 
   OUTPUTS:
   pairs_freq: a n_comb(N_cells, 2) X 2 matrix of frequencies, rows
   a are cell pairs columns possible orderings, according to 
   1: 1 2 
   2: 2 1 

   pairs_table: a
   n_comb(N_cells, 2) X 2 matrix, a lookup table for the row dimension in
   pairs_freq, i-th row is the set of cell indices for the i-th
   pair, in increasing order
   pairs_ts: a n_combs X 2 cell array of arrays, containing the
   timestamps at which each pair occurred

   this is really meant to be called by spike_pairs m
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
  
  int i1, i2;
  
  i1 = ii[0];
  i2 = ii[1];

  

  if(i1 == i2 )
    return -1;
  
  cd = 1;
  
  if(i1 < i2)
    cd = 0;
  
  
  return cd;
}


inline void order_tuple(int *tp, int n) /* order a n-tuple */
{
  int tmp;
  

  if (n != 2)
    mexErrMsgTxt("order_tuple specialized for pairs!");
  
  if(tp[0] > tp[1])
    {
      tmp = tp[1];
      tp[1] = tp[0];
      tp[0] = tmp;
    }
  
 
}

  
int n_cells; /* number of cells */

inline int lookup_idx(int i, int j)
{
  return i + n_cells * j;
}


  

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{
  
  double *tstamp; /* spike timestamps (in 1/10000 sec) */

  double *on_same_tet;
  
/*    double ***pairs_ts_ptr; */
  

/*    int **pairs_count; */
  int **pairs_last_ts;
  
  int *idxs;
  int cidx[2];
  
  int midx[2];

  int n_spikes; /* total number of spike */
  int n_pairs_true = 0 /*  n_pairs_chk = 0 */; /* diagnostics */
  
  int n_combs; /* number of combinations */
  
  int i,j, i_last;
  
  int omit_double = 0; /* if TRUE, then a spike can't be member of
			  more than one pair  of the same kind,
			  this helps avoiding many doubles that would
			  make the interpretation hard */
  

  double *pairs_freq;
  double *pairs_table;
  int *pairs_lookup; /* direct lookup table */

  
  double time_span;
  

  enum {TSTAMP=0, IDXS, NCELLS, TSPAN, ON_SAME_TET
  } in_arg_order;
  
  enum {PAIRS_FREQ, PAIRS_TABLE, PAIRS_TS  
  } out_arg_order;
  
  
  if(nrhs != 5 )
    mexErrMsgTxt("Five inputs required");
  

  if(nlhs != 3 )
    mexErrMsgTxt("Three outputs required");

  

  /* MEX related stuff, not really related with algorithms */

  /* extracting timestamps */

  mexPrintf("class of tstamp: %d\n", mxGetClassID(prhs[TSTAMP])); /* DEBUG */
  
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
  
  /* extracting on_same_tet */

  if(mxGetM(prhs[ON_SAME_TET]) != n_cells || 
     mxGetN(prhs[ON_SAME_TET]) != n_cells)
    mexErrMsgTxt("on_same_tet must be a n_cells X n_cells matrix");
  
  on_same_tet = mxGetPr(prhs[ON_SAME_TET]);

  
  /* create outputs */ 

  n_combs = 0;
  
  for(i = 0; i < n_cells; i++)
    for(j = i+1; j < n_cells; j++)
      {
	/*  	  int ix = lookup3[i] + (n_cells - i) *(j - i - 1)  */
	/*  	    - (j - i - 1) * (j -i) /2  + (k-j-1); */
	int i1;
	midx[0] = i;
	midx[1] = j;
	i1 = on_same_tet[mxCalcSingleSubscript(prhs[ON_SAME_TET], 2, midx)];

	if(!(i1))
	  n_combs++;
	    
	    
      }
  
  mexPrintf("n_combs = %d\n", n_combs);
  

  
  plhs[PAIRS_FREQ] = mxCreateNumericMatrix(n_combs, 2,
					      mxDOUBLE_CLASS,
					      mxREAL);

  pairs_freq = mxGetPr(plhs[PAIRS_FREQ]);
  

  plhs[PAIRS_TABLE] = mxCreateNumericMatrix(n_combs, 2,
					      mxDOUBLE_CLASS,
					      mxREAL);
  
  pairs_table = mxGetPr(plhs[PAIRS_TABLE]);
  

  plhs[PAIRS_TS] = mxCreateCellMatrix(n_combs, 2);





  /* create lookup tables direct and inverse */

  pairs_lookup = (int *) mxCalloc(n_cells*n_cells, sizeof(int));
  for(i = 0; i < n_cells*n_cells; i++)
    pairs_lookup[i] = -1;
  
  
  pairs_last_ts = (int **)mxCalloc(2, sizeof(int *));
  for(i = 0; i < 2; i++)
    pairs_last_ts[i] = (int *)mxCalloc(n_combs, sizeof(int));
  
  

  
  {
    int ix = 0;
    
    for(i = 0; i < n_cells; i++)
      for(j = i+1; j < n_cells; j++)
	{
	  /*  	  int ix = lookup3[i] + (n_cells - i) *(j - i - 1)  */
	  /*  	    - (j - i - 1) * (j -i) /2  + (k-j-1); */
	  int i1;
	  midx[0] = i;
	  midx[1] = j;
	  i1 = on_same_tet[mxCalcSingleSubscript(prhs[ON_SAME_TET], 2, midx)];
	
	    
	      
	  if(!(i1))
	    {
		
	      midx[0] = ix;
	      midx[1] = 0;
	      pairs_table[mxCalcSingleSubscript(plhs[PAIRS_TABLE], 2, midx)] =
		i+1;
	      midx[1] = 1;
	      pairs_table[mxCalcSingleSubscript(plhs[PAIRS_TABLE], 2, midx)] =
		j+1;

	      pairs_lookup[lookup_idx(i,j)] = ix;
	      ix++;
	    }
	    
	}
  }
  

  /* now complete the direct lookup table */

  for(i = 0; i < n_cells; i++)
    for(j = 0; j < n_cells; j++)
      {
	cidx[0] = i;
	cidx[1] = j;

	order_tuple(cidx, 2);
	pairs_lookup[lookup_idx(i,j)] = 
	  pairs_lookup[lookup_idx(cidx[0], cidx[1])];
      }
  

  /* and now compute the frequency table */
  i_last = 0;
  
  for(i = 0; i < n_spikes-2; i++)
    {
      
      double t1 = tstamp[i];
      double tl = t1 + time_span;
      int ix = idxs[i];
      int jx, cx, px;
 
     if(! (i % 10000))
	mexPrintf("i: %d\n", i);
      
      while(i_last < n_spikes && tstamp[i_last] < tl) i_last++;
      
      for(j = i+1; j < i_last; j++)
	if(idxs[j] != ix)
	  {
	    jx = idxs[j];
		{

		  

		  cidx[0] = ix;
		  cidx[1] = jx;

		  px = order_code(cidx);
		  
		  order_tuple(cidx, 2);
		  cx = pairs_lookup[lookup_idx(cidx[0], cidx[1])];
		  if(cx >= 0 && (!omit_double ||
		      (pairs_last_ts[0][cx] != tstamp[i] &&
		       pairs_last_ts[1][cx] != tstamp[j])))
		    {
		      pairs_last_ts[0][cx] = tstamp[i];
		      pairs_last_ts[1][cx] = tstamp[j];

		      
		      n_pairs_true++;
		      midx[0] = cx;
		      midx[1] = px;
		      pairs_freq[mxCalcSingleSubscript(plhs[PAIRS_FREQ], 2, midx)]++;
		    }
		  
		}
	  }
    }
  
  
	

   /* now allocate memory for timestamps */
  /*
  pairs_ts_ptr = (double ***)mxCalloc(2, sizeof(double **));
  pairs_count = (int **)mxCalloc(2, sizeof(int *));
  
  for(i = 0; i < 2; i++)
    {
      pairs_ts_ptr[i] = (double **)mxCalloc(n_combs, sizeof(double **));
      pairs_count[i] = (int *)mxCalloc(n_combs, sizeof(int));
    }
  
  

  for(i = 0; i < n_combs; i++)
    for(j = 0; j < 2; j++)
      {
	int n_pairs;
	mxArray *tmp;
	
	midx[0] = i;
	midx[1] = j;
	n_pairs = 
	  pairs_freq[mxCalcSingleSubscript(plhs[PAIRS_FREQ], 2,
	  midx)]; */
/*  	tmp = mxCreateNumericMatrix(n_pairs, 1, mxDOUBLE_CLASS, */
/*  				    mxREAL); */
/*  	pairs_ts_ptr[j][i] = mxGetPr(tmp); */
	
/*  	mxSetCell(plhs[PAIRS_TS],  */
/*  		  mxCalcSingleSubscript(plhs[PAIRS_TS], 2, midx), */
/*  		  tmp); */

/*  	n_pairs_chk += n_pairs; */
	


	


  /*      } */
  
	
/*    if (n_pairs_chk != n_pairs_true) */
/*      { */
/*        char errmsg[80]; */
/*        sprintf(errmsg, "n_pairs = %d, n_pairs_chk = %d!!!", n_pairs_true,  */
/*  	      n_pairs_chk); */
      
/*  	mexErrMsgTxt(errmsg); */
/*      } */
  

  /* now fill timestamp and position arrays */
	  

/*    for (i = 0; i < 2; i++) */
/*      for(j = 0; j < n_combs; j++) */
/*        pairs_last_ts[i][j] = 0; */
  
  
/*      i_last = 0; */
  
/*    for(i = 0; i < n_spikes-2; i++) */
/*      { */
      
/*        double t1 = tstamp[i]; */
/*        double tl = t1 + time_span; */
/*        int ix = idxs[i]; */
/*        int jx, cx, px; */
 
/*       if(! (i % 10000)) */
/*  	mexPrintf("i: %d\n", i); */
      
/*        while(i_last < n_spikes && tstamp[i_last] < tl) i_last++; */
      
/*        for(j = i+1; j < i_last; j++) */
/*  	if(idxs[j] != ix) */
/*  	  { */
/*  	    jx = idxs[j]; */
/*  		{ */
/*  		  cidx[0] = ix; */
/*  		  cidx[1] = jx; */
/*  		  px = order_code(cidx); */
		  
/*  		  order_tuple(cidx, 2); */
/*  		  cx = pairs_lookup[lookup_idx(cidx[0], cidx[1])]; */
/*  		  if(cx >= 0 && (!omit_double || */
/*  		      (pairs_last_ts[0][cx] != tstamp[i] && */
/*  		       pairs_last_ts[1][cx] != tstamp[j]))) */
/*  		    { */
/*  		      pairs_last_ts[0][cx] = tstamp[i]; */
/*  		      pairs_last_ts[1][cx] = tstamp[j]; */

/*  		      pairs_ts_ptr[px][cx][pairs_count[px][cx]] =  */
/*  			tstamp[i]; */
		  

/*  		      pairs_count[px][cx]++; */
/*  		    } */
		  

		  
/*  		} */
/*  	  } */
/*      } */
/*    */ 

/*    mxFree((void *)pairs_lookup); */

/*    for(i = 0; i < 2; i++) */
/*      { */

/*        mxFree((void *)pairs_ts_ptr[i]); */
/*        mxFree((void *)pairs_count[i]); */
      
/*      } */

/*    for(i = 0; i < 2; i++) */
/*      mxFree((void *)pairs_last_ts[i]); */
  
      
/*    mxFree((void *)pairs_ts_ptr); */
/*    mxFree((void *)pairs_count); */
/*    mxFree((void *)pairs_last_ts); */
  
/*    mexCallMATLAB(0, NULL, 1, &plhs[PAIRS_PHI], "disp"); */
  

  
}
