/*=================================================================
 *
 * YPRIME.C	Sample .MEX file corresponding to YPRIME.M
 *	        Solves simple 3 body orbit problem 
 *
 * The calling syntax is:
 *
 *		[yp] = yprime(t, y)
 *
 *  You may also want to look at the corresponding M-code, yprime.m.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2011 The MathWorks, Inc.
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	T_IN	prhs[0]
#define	Y_IN	prhs[1]


/* Output Arguments */

#define	YP_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

static	double	mu = 1/82.45;
static	double	mus = 1 - 1/82.45;


static void yprime(
		   double	yp[],
		   double	*t,
 		   double	y[]
		   )
{
    double	r1,r2;
    
    (void) t;     /* unused parameter */

    r1 = sqrt((y[0]+mu)*(y[0]+mu) + y[2]*y[2]); 
    r2 = sqrt((y[0]-mus)*(y[0]-mus) + y[2]*y[2]);

    /* Print warning if dividing by zero. */    
    if (r1 == 0.0 || r2 == 0.0 ){
        mexWarnMsgIdAndTxt( "MATLAB:yprime:divideByZero", 
                "Division by zero!\n");
    }
    
    yp[0] = y[1];
    yp[1] = 2*y[3]+y[0]-mus*(y[0]+mu)/(r1*r1*r1)-mu*(y[0]-mus)/(r2*r2*r2);
    yp[2] = y[3];
    yp[3] = -2*y[1] + y[2] - mus*y[2]/(r1*r1*r1) - mu*y[2]/(r2*r2*r2);
    return;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *yp; 
    double *t,*y; 
    size_t m,n; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 2) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidNumInputs",
                "Two input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:yprime:maxlhs",
                "Too many output arguments."); 
    } 
    
    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
    
    m = mxGetM(Y_IN); 
    n = mxGetN(Y_IN);
    if (!mxIsDouble(Y_IN) || mxIsComplex(Y_IN) || 
	(MAX(m,n) != 4) || (MIN(m,n) != 1)) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidY",
                "YPRIME requires that Y be a 4 x 1 vector."); 
    } 
    
    /* Create a matrix for the return argument */ 
    YP_OUT = mxCreateDoubleMatrix( (mwSize)m, (mwSize)n, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    yp = mxGetPr(YP_OUT);
    
    t = mxGetPr(T_IN); 
    y = mxGetPr(Y_IN);
        
    /* Do the actual computations in a subroutine */
    yprime(yp,t,y); 
    return;
    
}
