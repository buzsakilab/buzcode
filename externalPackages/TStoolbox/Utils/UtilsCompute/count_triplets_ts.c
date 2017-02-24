/* count_triplets.c: MEX file, "engine" for the spike_triplets
   function, counts the triples of spikes  occurring within a time
   span, also outputs, in two cell arrays, the timestamps at which
   each ordered  triplet occurred, and the locations where they occurred
   INPUTS:
   t: a n_spikes x 1 array of timestamps
   idxs: a n_spikes x 1 array of cell indices corresponding to the
   x_pos:
   y_pos: 
   phi_pos: 3 n_spikes x 1 arrays of position coordinates, for each
   spike (x, y, and angular)
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
   triplets_ts: a n_combs X 6 cell array of arrays, containing the
   timestamps at which each triplet occurred
   triplets_x
   triplets_y 
   triplets_phi: 3 n_combs X 6 cell array of arrays, containing the
   locations 
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
  double *x_pos, *y_pos, *phi_pos;
  double *on_same_tet;
  
  double ***triplets_ts_ptr, ***triplets_x_ptr, ***triplets_y_ptr,
    ***triplets_phi_ptr;
  int **triplets_count;
  int **triplets_last_ts;
  
  int *idxs;
  int cidx[3];
  
  int midx[2];

  int n_spikes; /* total number of spike */
  int n_triplets_true = 0, n_triplets_chk = 0; /* diagnostics */
  
  int n_combs; /* number of combinations */
  
  int i,j,k, i_last;
  
  int omit_double = 1; /* if TRUE, then a spike can't be member of
			  more than one triplet  of the same kind,
			  this helps avoiding many doubles that would
			  make the interpretation hard */
  

  double *triplets_freq;
  double *triples_table;
  int *triples_lookup; /* direct lookup table */

  
  double time_span;
  

  enum {TSTAMP=0, IDXS, X_POS, Y_POS, PHI_POS, NCELLS, TSPAN, ON_SAME_TET
  } in_arg_order;
  
  enum {TRIPLETS_FREQ, TRIPLES_TABLE, TRIPLETS_TS, TRIPLETS_X,
	TRIPLETS_Y, TRIPLETS_PHI
  } out_arg_order;
  
  
  if(nrhs != 8 )
    mexErrMsgTxt("Eight inputs required");
  

  if(nlhs != 6 )
    mexErrMsgTxt("Six outputs required");

  

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
  

  /* extracting location */

  mexPrintf("class of x: %d\n", mxGetClassID(prhs[X_POS])); /* DEBUG */

  if(mxGetN(prhs[X_POS]) != 1) 
    mexErrMsgTxt("x_pos must be a column vector");
  
  x_pos = mxGetPr(prhs[X_POS]);

  if(mxGetN(prhs[Y_POS]) != 1) 
    mexErrMsgTxt("y_pos must be a column vector");
  
  y_pos = mxGetPr(prhs[Y_POS]);

  if(mxGetN(prhs[PHI_POS]) != 1) 
    mexErrMsgTxt("phi_pos must be a column vector");
  
  phi_pos = mxGetPr(prhs[PHI_POS]);

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
  
  on_same_tet = (int)mxGetPr(prhs[ON_SAME_TET]);

  
  /* create outputs */ 

  n_combs = 0;
  
  for(i = 0; i < n_cells; i++)
    for(j = i+1; j < n_cells; j++)
      for(k = j+1; k < n_cells; k++)
	  {
	    /*  	  int ix = lookup3[i] + (n_cells - i) *(j - i - 1)  */
	    /*  	    - (j - i - 1) * (j -i) /2  + (k-j-1); */
	    int i1, i2, i3;
	    midx[0] = i;
	    midx[1] = j;
	    i1 = on_same_tet[mxCalcSingleSubscript(prhs[ON_SAME_TET], 2, midx)];
	    midx[0] = i;
	    midx[1] = k;
	    i2 = on_same_tet[mxCalcSingleSubscript(prhs[ON_SAME_TET], 2, midx)];
	    midx[0] = j;
	    midx[1] = k;
	    i3 = on_same_tet[mxCalcSingleSubscript(prhs[ON_SAME_TET], 2, midx)];
	    if(!(i1 || i2 || i3))
	      n_combs++;
	    
	    
	  }
  
  mexPrintf("n_combs = %d\n", n_combs);
  

  
  plhs[TRIPLETS_FREQ] = mxCreateNumericMatrix(n_combs, 6,
					      mxDOUBLE_CLASS,
					      mxREAL);

  triplets_freq = mxGetPr(plhs[TRIPLETS_FREQ]);
  

  plhs[TRIPLES_TABLE] = mxCreateNumericMatrix(n_combs, 3,
					      mxDOUBLE_CLASS,
					      mxREAL);
  
  triples_table = mxGetPr(plhs[TRIPLES_TABLE]);
  

  plhs[TRIPLETS_TS] = mxCreateCellMatrix(n_combs, 6);

  plhs[TRIPLETS_X] = mxCreateCellMatrix(n_combs, 6);

  plhs[TRIPLETS_Y] = mxCreateCellMatrix(n_combs, 6);

  plhs[TRIPLETS_PHI] = mxCreateCellMatrix(n_combs, 6);
  





  /* create lookup tables direct and inverse */

  triples_lookup = (int *) mxCalloc(n_cells*n_cells*n_cells, sizeof(int));
  for(i = 0; i < n_cells*n_cells*n_cells; i++)
    triples_lookup[i] = -1;
  
  
  triplets_last_ts = (int **)mxCalloc(3, sizeof(int *));
  for(i = 0; i < 3; i++)
    triplets_last_ts[i] = (int *)mxCalloc(n_combs, sizeof(int));
  
  

  
  {
    int ix = 0;
    
    for(i = 0; i < n_cells; i++)
      for(j = i+1; j < n_cells; j++)
	for(k = j+1; k < n_cells; k++)
	  {
	    /*  	  int ix = lookup3[i] + (n_cells - i) *(j - i - 1)  */
	    /*  	    - (j - i - 1) * (j -i) /2  + (k-j-1); */
	    int i1, i2, i3;
	    midx[0] = i;
	    midx[1] = j;
	    i1 = on_same_tet[mxCalcSingleSubscript(prhs[ON_SAME_TET], 2, midx)];
	    midx[0] = i;
	    midx[1] = k;
	    i2 = on_same_tet[mxCalcSingleSubscript(prhs[ON_SAME_TET], 2, midx)];
	    midx[0] = j;
	    midx[1] = k;
	    i3 = on_same_tet[mxCalcSingleSubscript(prhs[ON_SAME_TET], 2, midx)];
	    
	      
	    if(!(i1 || i2 || i3))
	      {
		
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
		  if(cx >= 0 && (!omit_double ||
		      (triplets_last_ts[0][cx] != tstamp[i] &&
		       triplets_last_ts[1][cx] != tstamp[j] &&
		       triplets_last_ts[2][cx] != tstamp[k])))
		    {
		      triplets_last_ts[0][cx] = tstamp[i];
		      triplets_last_ts[1][cx] = tstamp[j];
		      triplets_last_ts[2][cx] = tstamp[k];
		      
		      n_triplets_true++;
		      midx[0] = cx;
		      midx[1] = px;
		      triplets_freq[mxCalcSingleSubscript(plhs[TRIPLETS_FREQ], 2, midx)]++;
		    }
		  
		}
	  }
    }
  
	

  /* now allocate memory for timestamps, x, y */

  triplets_ts_ptr = (double ***)mxCalloc(6, sizeof(double **));
  triplets_x_ptr = (double ***)mxCalloc(6, sizeof(double **));
  triplets_y_ptr = (double ***)mxCalloc(6, sizeof(double **));
  triplets_phi_ptr = (double ***)mxCalloc(6, sizeof(double **));
  triplets_count = (int **)mxCalloc(6, sizeof(int *));
  
  for(i = 0; i < 6; i++)
    {
      triplets_ts_ptr[i] = (double **)mxCalloc(n_combs, sizeof(double **));
      triplets_x_ptr[i] = (double **)mxCalloc(n_combs, sizeof(double *));
      triplets_y_ptr[i] = (double **)mxCalloc(n_combs, sizeof(double *));
      triplets_phi_ptr[i] = (double **)mxCalloc(n_combs, sizeof(double *));
      triplets_count[i] = (int *)mxCalloc(n_combs, sizeof(int));
    }
  
  

  for(i = 0; i < n_combs; i++)
    for(j = 0; j < 6; j++)
      {
	int n_triplets;
	mxArray *tmp;
	
	midx[0] = i;
	midx[1] = j;
	n_triplets = 
	  triplets_freq[mxCalcSingleSubscript(plhs[TRIPLETS_FREQ], 2,
					      midx)];
	tmp = mxCreateNumericMatrix(n_triplets, 1, mxDOUBLE_CLASS,
				    mxREAL);
	triplets_ts_ptr[j][i] = mxGetPr(tmp);
	
	mxSetCell(plhs[TRIPLETS_TS], 
		  mxCalcSingleSubscript(plhs[TRIPLETS_TS], 2, midx),
		  tmp);

	n_triplets_chk += n_triplets;
	
	tmp = mxCreateNumericMatrix(n_triplets, 1, mxDOUBLE_CLASS,
				    mxREAL);
	triplets_x_ptr[j][i] = mxGetPr(tmp);
	
	mxSetCell(plhs[TRIPLETS_X], 
		  mxCalcSingleSubscript(plhs[TRIPLETS_X], 2, midx),
		  tmp);


	tmp = mxCreateNumericMatrix(n_triplets, 1, mxDOUBLE_CLASS,
				    mxREAL);
	triplets_y_ptr[j][i] = mxGetPr(tmp);
	
	mxSetCell(plhs[TRIPLETS_Y], 
		  mxCalcSingleSubscript(plhs[TRIPLETS_Y], 2, midx),
		  tmp);
	


	tmp = mxCreateNumericMatrix(n_triplets, 1, mxDOUBLE_CLASS,
				    mxREAL);
	triplets_phi_ptr[j][i] = mxGetPr(tmp);
	
	mxSetCell(plhs[TRIPLETS_PHI], 
		  mxCalcSingleSubscript(plhs[TRIPLETS_PHI], 2, midx),
		  tmp);
      }
  
	
  if (n_triplets_chk != n_triplets_true)
    {
      char errmsg[80];
      sprintf(errmsg, "n_triplets = %d, n_triplets_chk = %d!!!", n_triplets_true, 
	      n_triplets_chk);
      
	mexErrMsgTxt(errmsg);
    }
  

  /* now fill timestamp and position arrays */
	  

  for (i = 0; i < 3; i++)
    for(j = 0; j < n_combs; j++)
      triplets_last_ts[i][j] = 0;
  
  
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
		  if(cx >= 0 && (!omit_double ||
		      (triplets_last_ts[0][cx] != tstamp[i] &&
		       triplets_last_ts[1][cx] != tstamp[j] &&
		       triplets_last_ts[2][cx] != tstamp[k])))
		    {
		      triplets_last_ts[0][cx] = tstamp[i];
		      triplets_last_ts[1][cx] = tstamp[j];
		      triplets_last_ts[2][cx] = tstamp[k];
		      triplets_ts_ptr[px][cx][triplets_count[px][cx]] = 
			tstamp[i];
		  
		      triplets_x_ptr[px][cx][triplets_count[px][cx]] = 
			x_pos[i];

		      triplets_y_ptr[px][cx][triplets_count[px][cx]] = 
			y_pos[i];

		      triplets_phi_ptr[px][cx][triplets_count[px][cx]] = 
			phi_pos[i];

		      triplets_count[px][cx]++;
		    }
		  

		  
		}
	  }
    }


  mxFree((void *)triples_lookup);

  for(i = 0; i < 6; i++)
    {

      mxFree((void *)triplets_ts_ptr[i]);
      mxFree((void *)triplets_x_ptr[i]);
      mxFree((void *)triplets_y_ptr[i]);
      mxFree((void *)triplets_phi_ptr[i]);
      mxFree((void *)triplets_count[i]);
      
    }

  for(i = 0; i < 3; i++)
    mxFree((void *)triplets_last_ts[i]);
  
      
  mxFree((void *)triplets_ts_ptr);
  mxFree((void *)triplets_x_ptr);
  mxFree((void *)triplets_y_ptr);
  mxFree((void *)triplets_phi_ptr);
  mxFree((void *)triplets_count);
  mxFree((void *)triplets_last_ts);
  
/*    mexCallMATLAB(0, NULL, 1, &plhs[TRIPLETS_PHI], "disp"); */
  

  
}
