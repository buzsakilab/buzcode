

#include <math.h>
#include "mex.h"
#include "matrix.h"


#define ValidInterval(x,y) ((idx[y]-idx[x]) > minT && (idx[y]-idx[x])  < maxT)



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
		 *prhs[])
{
  int i, j, k, l, n;
  int mrows, ncols;
  int ndata, dimdata;
  int nidx;
  int ndims;
  
  int subs[2];
  int osubs[3];
  int nseq;
  
  int sz;
  int cv, pv;
  
  double *ptr;
  double **data;
  double *idx;
  double minT, maxT;
  
  mxArray *M;
  
  

  if(nrhs != 5 )
    mexErrMsgTxt("Five inputs required");
  
  if(nlhs < 1)
    mexErrMsgTxt("At least one output required");
  
  if(nlhs > 2) 
    mexErrMsgTxt("Too many output arguments");
  

  /* input 1 is the data matrix */

  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  ndata = mrows;
  dimdata = ncols;

  
  

  ptr = mxGetPr(prhs[0]);


  data = (double **)mxCalloc(mrows, sizeof(double *));
  for(i = 0; i < mrows; i++)
    data[i] = mxCalloc(ncols, sizeof(double));
  
  
  for(i = 0; i < mrows; i++)
    for(j = 0; j < ncols; j++)
      {
	subs[0] = i;
	subs[1] = j;
	data[i][j] = ptr[mxCalcSingleSubscript(prhs[0], 2, subs)];
      }
  
  
  /* input 2 is the index array */

  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  

  if(mrows != 1 && ncols != 1)
    mexErrMsgTxt("input 2 must be a vector.");
  
  idx = mxGetPr(prhs[1]);
  nidx = mrows * ncols;
  


  /* input 3 is the sequence size */

  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  
  if(mrows * ncols != 1)
    mexErrMsgTxt("input 3 must be a scalar.");
  
  sz = mxGetScalar(prhs[2]);
  
    /* input 4 is the minimum ISI for c.s. */

  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  
  if(mrows * ncols != 1)
    mexErrMsgTxt("input 3 must be a scalar.");
  
  minT = mxGetScalar(prhs[3]);

  /* input 5 is the maximum ISI for c.s. */

  mrows = mxGetM(prhs[4]);
  ncols = mxGetN(prhs[4]);
  
  if(mrows * ncols != 1)
    mexErrMsgTxt("input 3 must be a scalar.");
  
  maxT = mxGetScalar(prhs[4]);
  
  /* compute the number of sequences */

     nseq = 0;
     pv = 0;
     
     for(i = 1; i < nidx; i++)
       {
      	 if((idx[i] - idx[i-1]) < 0)
	   mexErrMsgTxt("idx must be non decreasing.");
	 cv = ValidInterval(i-1, i);
	 
	 if(cv && (!pv))
	   nseq++;
	 pv = cv;
	 
       }
     


     osubs[0] = nseq;
     osubs[1] = sz;
     osubs[2] = dimdata;
     ndims = 3;
     if(dimdata == 1)
       ndims = 2;
     
     
     M = mxCreateNumericArray(ndims, osubs, mxDOUBLE_CLASS, mxREAL);
     
     ptr = mxGetPr(M);
     
     i = 0;
     j = 0;
     n = 0;
     
     pv = 0;
     
     while(i < nidx) 
       {
	 if(ValidInterval(i-1,i))
	   {
	     
	     osubs[0] = n;
	     osubs[1] = j;
	     for(k = 0; k < dimdata; k++) 
	       {
		 osubs[2] = k;
		 ptr[mxCalcSingleSubscript(M, ndims, osubs)] = 
		   data[i-1][k];
	       }
	     j++;
	     if(j >= sz)
	       mexErrMsgTxt("A sequence exceeds sz.");
	   
	 

	     if(!ValidInterval(i,i+1))
	       {
		 osubs[1] = j;
		 for(k = 0; k < dimdata; k++) 
		   {
		     osubs[2] = k;
		     ptr[mxCalcSingleSubscript(M, ndims, osubs)] = 
		       data[i][k];
		   }
		 j++;
		 if(j >= sz)
		   mexErrMsgTxt("A sequence exceeds sz.");

		 for(l=j; l < sz;l++)
		   {
		     osubs[1] = l;
		     for(k = 0; k < dimdata; k++) 
		       {
			 osubs[2] = k;
			 ptr[mxCalcSingleSubscript(M, ndims, osubs)] = 0./0.;
		       }
		 
		   
		   }
	     
		 j = 0;
		 n++;
	       }

	   }
	 i++;
	 
       }
     
     
/*       for(l=j; l < sz;l++) */
/*         { */
/*  	 osubs[1] = l; */
/*  	 for(k = 0; k < dimdata; k++)  */
/*  	   { */
/*  	     osubs[2] = k; */
/*  	     ptr[mxCalcSingleSubscript(M, ndims, osubs)] = 0./0.; */
/*  	   } */
		 
     /* } */
     
     plhs[0] = M;
     for(i = 0; i < ndata; i++)
       mxFree((void *)data[i]);
     
     mxFree((void *)data);
     


	
     
}


  



    
  
