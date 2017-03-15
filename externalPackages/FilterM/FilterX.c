// FilterX.c
// FilterX - Fast C-Mex filter
// [Y, Z] = FilterX(b, a, X, Z, Reverse)
// INPUT:
//   b, a: Filter parameters as DOUBLE vectors. If the vectors have different
//      lengths, the shorter one is padded with zeros.
//   X: Signal as DOUBLE or SINGLE vector or array. The signal is filtered along
//      the first dimension (!even if X is a row vector!).
//   Z: Initial conditions as DOUBLE or SINGLE array. The size must be:
//        [(Order) - 1, SIZE(X,2), ..., SIZE(X, NDIMS(X))]
//      Optional, default: Zeros.
//   Reverse: The signal is processed in reverse order, if this is TRUE or
//      'reverse'. b, a and Z are not affected - see examples.
//      This is supports a faster FILTFILT operation.
//      Optional, default: FALSE.
//
// OUTPUT:
//   Y: Filtered signal with the same size and type as X.
//   Z: Final conditions as DOUBLE (!) array.
//
// NOTES:
//   - The output equals the output of FILTER exactly for DOUBLEs.
//   - For signals in SINGLE format, all intermediate values are DOUBLEs to
//     increase the accuracy. Therefore the output Z is a DOUBLE also.
//   - The parameters [a] and [b] are normalized to the 1st element of a, if
//     a(1) differs from 1.0.
//   - This function filters along the 1st dimension only. Use FilterM as
//     wrapper to process other dimensions also.
//
// COMPILATION:
//   mex -O FilterX.c
// Consider C99 comments on Linux with GCC:
//   mex -O CFLAGS="\$CFLAGS -std=c99" FilterX.c
// Pre-compiled Mex file can be downloaded at: http://www.n-simon.de/mex
// Run the unit-test uTest_FilterX after compiling.
//
// Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32bit
//         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008
// Assumed Compatibility: higher Matlab versions, Mac, Linux, 64bit
// Author: Jan Simon, Heidelberg, (C) 2011 matlab.THISYEAR(a)nMINUSsimon.de

/*
% $JRev: R-m V:013 Sum:I9sHyNUP/PnE Date:20-Jul-2011 01:12:31 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_DGradient $
% $File: Tools\Mex\Source\FilterX.c $
% History:
% 011: 17-Jul-2011 22:06, First stable version.
%      This is a naive implementation of the time domain difference equations.
%      I cannot imagine why this is faster than Matlab's FILTER.
%      The multi-threaded FILTER of Matlab 2011a might be faster for > 3 cores,
%      but this works only if the signal has > 3 columns also.
*/

// Headers: --------------------------------------------------------------------
#include "mex.h"
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>

// Definitions: ----------------------------------------------------------------
// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T               // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Limit number of dimensions of the input - this saves 2% computing time if
// the signal is tiny (e.g. [16 x 1]):
#define MAX_NDIMS 32

// Disable the /fp:precise flag to increase the speed on MSVC compiler:
#ifdef _MSC_VER
#pragma float_control(except, off)    // disable exception semantics
#pragma float_control(precise, off)   // disable precise semantics
#pragma fp_contract(on)               // enable contractions
// #pragma fenv_access(off)           // disable fpu environment sensitivity
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_HEAD "*** FilterX[mex]: "
#define ERR_ID   "JSimon:FilterX:"

// Some macros to reduce the confusion:
#define b_in   prhs[0]
#define a_in   prhs[1]
#define X_in   prhs[2]
#define Z_in   prhs[3]
#define Rev_in prhs[4]
#define Y_out  plhs[0]
#define Z_out  plhs[1]

// Prototypes: -----------------------------------------------------------------
void CoreDoubleN(double *X, mwSize MX, mwSize NX, double *a, double *b,
        mwSize order, double *Z, double *Y);
void CoreDouble2(double *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, double *Y);
void CoreDouble3(double *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, double *Y);
void CoreDouble4(double *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, double *Y);
void CoreDouble5(double *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, double *Y);
void CoreDouble6(double *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, double *Y);
void CoreDouble7(double *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, double *Y);

void CoreSingleN(float *X, mwSize MX, mwSize NX, double *a, double *b,
        mwSize order, double *Z, float *Y);
void CoreSingle2(float *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, float *Y);
void CoreSingle3(float *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, float *Y);
void CoreSingle4(float *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, float *Y);
void CoreSingle5(float *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, float *Y);
void CoreSingle6(float *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, float *Y);
void CoreSingle7(float *X, mwSize MX, mwSize NX, double *a, double *b,
        double *Z, float *Y);

void CoreDoubleNR(double *X, mwSize MX, mwSize NX, double *a, double *b,
        mwSize order, double *Z, double *Y);
void CoreSingleNR(float *X, mwSize MX, mwSize NX, double *a, double *b,
        mwSize order, double *Z, float *Y);

void CopySingleToDouble(double *Z, float *Zf, mwSize N);
void NormalizeBA(double *ab, mwSize nParam);
        
// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *a, *b, *Z, *X, *Y, a0;
  float  *Xf, *Yf;
  mwSize na, nb, order, nParam, MX, NX, Xndims, Zdims[MAX_NDIMS];
  bool   allocate_ba = false, allocate_Z = false, forward = true, hasZInput;
  char   Rev[2];
  const mwSize *Xdims;
  
  // Check number of inputs and outputs:
  if (nrhs < 3 || nrhs > 5) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNInput",
                       ERR_HEAD "3 to 5 inputs required.");
  }
  if (nlhs > 2) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNOutput",
                       ERR_HEAD "2 outputs allowed.");
  }
  
  // Check type of inputs:
  if (!mxIsDouble(b_in) || !mxIsDouble(a_in)) {
     mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput1_2",
                       ERR_HEAD "Filter parameters must be DOUBLES.");
  }
  
  // Get dimensions of inputs:
  nb = mxGetNumberOfElements(b_in);
  na = mxGetNumberOfElements(a_in);
  MX = mxGetM(X_in);    // First dimension
  NX = mxGetN(X_in);    // Product of trailing dimensions
  
  // Reply empty array if the parameters or the signal is empty:
  if (na * nb * MX * NX == 0) {
     plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
     if (nrhs == 2) {
        plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
     }
     return;
  }
  
  // Get a and b as vectors of the same length:
  if (na == nb) {       // Use input vectors directly, if possible:
     b      = mxGetPr(b_in);
     a      = mxGetPr(a_in);
     nParam = nb;
     
     allocate_ba = (bool) (a[0] != 1.0);  // Normalization is needed
     
  } else {              // na differs from nb:
     nParam      = na > nb ? na : nb;
     allocate_ba = true;
  }
  order = nParam - 1;
  
  if (allocate_ba) {    // Create local copy of expanded [b] and [a]:
     // It is slightly cheaper to allocate one array only:
     if ((b = mxCalloc(2 * nParam, sizeof(double))) == NULL) {
        mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                          ERR_HEAD "Cannot get memory for parameters.");
     }
     a = b + nParam;    // Use 2nd half of vector [b] as [a]
     memcpy(b, mxGetPr(b_in), nb * sizeof(double));
     memcpy(a, mxGetPr(a_in), na * sizeof(double));
     
     // Normalize if 1st element of [a] does not equal 1.0:
     if (a[0] != 1.0) {
        NormalizeBA(b, nParam);
     }
  }
  
  // Create array for final conditions, insert value of initial conditions:
  Xndims = mxGetNumberOfDimensions(X_in);
  Xdims  = mxGetDimensions(X_in);
  if (Xndims > MAX_NDIMS) {
     mexErrMsgIdAndTxt(ERR_ID   "BadSizeSignal",
                       ERR_HEAD "Signal cannot have more than %d dimensions.",
                       MAX_NDIMS);
  }
  
  // Z has dimensions [order, DIMS(X, 2:end)]:
  memcpy(Zdims, Xdims, Xndims * sizeof(mwSize));
  Zdims[0] = order;
  
  // 4th input is existing and not empty:
  hasZInput = false;
  if (nrhs >= 4) {
     hasZInput = (bool) (!mxIsEmpty(Z_in));
  }
  
  if (hasZInput) { // Initial conditions provided as input:
     // Check dimensions of Z:
     if (mxGetM(Z_in) != order || mxGetN(Z_in) != NX) {
        mexErrMsgIdAndTxt(ERR_ID "BadSizeZ",
                   ERR_HEAD
                   "Dimensions of initial conditions do not match the signal.");
     }
     
     if (nlhs == 2) {  // Final conditions wanted:
        Z_out = mxCreateNumericArray(Xndims, Zdims, mxDOUBLE_CLASS, mxREAL);
        Z     = mxGetPr(Z_out);
     } else {          // Final conditions not caught in the caller:
        allocate_Z = true;
        if ((Z = mxCalloc(order * NX, sizeof(double))) == NULL) {
           mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                             ERR_HEAD "No memory for initial conditions.");
        }
     }
     
     // Copy value from input with a conversion to DOUBLE on demand:
     if (mxIsDouble(Z_in)) {
        memcpy(Z, mxGetPr(Z_in), order * NX * sizeof(double));
     } else if (mxIsSingle(Z_in)) {
        CopySingleToDouble(Z, (float *) mxGetData(Z_in), order * NX);
     } else {
        mexErrMsgIdAndTxt(ERR_ID "BadTypeInput4",
               ERR_HEAD "Initial conditions must be a DOUBLE or SINGLE array.");
     }
     
  } else {             // Initial conditions not provided:
     if (nlhs == 2) {  // Final conditions wanted:
        Z_out = mxCreateNumericArray(Xndims, Zdims, mxDOUBLE_CLASS, mxREAL);
        Z     = mxGetPr(Z_out);
     } else {          // Cheaper to to use malloc than creating a Matlab array:
        allocate_Z = true;
        if ((Z = mxCalloc(order * NX, sizeof(double))) == NULL) {
           mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                             ERR_HEAD "No memory for initial conditions.");
        }
     }
  }

  // Flag for forward processing: ----------------------------------------------
  // 'Reverse', 1, TRUE: Process signal in reverse order:
  if (nrhs == 5) {
     if (!mxIsEmpty(Rev_in)) {
        if (mxIsChar(Rev_in)) {
           mxGetString(Rev_in, Rev, 2);
           forward = (bool) ((*Rev != 'r') && (*Rev != 'R'));
        } else if (mxIsLogical(Rev_in)) {
           forward = (bool) (mxGetScalar(Rev_in) == 0);
        }
     } else {
        // Most likely the user tries to define the dimension as in FILTER:
        mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput5",
                          ERR_HEAD "Reverse flag must be a string or LOGICAL.");
     }
  }
       
  // Call the calulator: -------------------------------------------------------
  if (mxIsDouble(X_in)) {
     // Create the output array:
     Y_out = mxCreateNumericArray(Xndims, Xdims, mxDOUBLE_CLASS, mxREAL);
     Y     = mxGetPr(Y_out);
     X     = mxGetPr(X_in);
     
     if (forward) {
#    if defined(__BORLAND__)
        // BCC 5.5 runs the unrolled loops with the half speed!
        CoreDoubleN(X, MX, NX, a, b, order, Z, Y);
#    else
        // Unrolled loops for common filter lengths:
        switch (order) {
           case 1:   CoreDouble2(X, MX, NX, a, b, Z, Y);  break;
           case 2:   CoreDouble3(X, MX, NX, a, b, Z, Y);  break;
           case 3:   CoreDouble4(X, MX, NX, a, b, Z, Y);  break;
           case 4:   CoreDouble5(X, MX, NX, a, b, Z, Y);  break;
           case 5:   CoreDouble6(X, MX, NX, a, b, Z, Y);  break;
           case 6:   CoreDouble7(X, MX, NX, a, b, Z, Y);  break;
           default:  CoreDoubleN(X, MX, NX, a, b, order, Z, Y);
        }
#    endif
     } else {  // Reverse order:
        CoreDoubleNR(X, MX, NX, a, b, order, Z, Y);
     }
     
  } else if (mxIsSingle(X_in)) {
     // Create the output array:
     Y_out = mxCreateNumericArray(Xndims, Xdims, mxSINGLE_CLASS, mxREAL);
     Yf    = (float *) mxGetData(Y_out);
     Xf    = (float *) mxGetData(X_in);
     
     if (forward) {
#    if defined(__BORLAND__)
        // BCC 5.5 runs the unrolled loops with the half speed!
        CoreSingleN(X, MX, NX, a, b, order, Z, Y);
#    else
        // Unrolled loops for common filter lengths:
        switch (order) {
           case 1:   CoreSingle2(Xf, MX, NX, a, b, Z, Yf);  break;
           case 2:   CoreSingle3(Xf, MX, NX, a, b, Z, Yf);  break;
           case 3:   CoreSingle4(Xf, MX, NX, a, b, Z, Yf);  break;
           case 4:   CoreSingle5(Xf, MX, NX, a, b, Z, Yf);  break;
           case 5:   CoreSingle6(Xf, MX, NX, a, b, Z, Yf);  break;
           case 6:   CoreSingle7(Xf, MX, NX, a, b, Z, Yf);  break;
           default:  CoreSingleN(Xf, MX, NX, a, b, order, Z, Yf);
        }
#    endif
     } else {  // Reverse order:
        CoreSingleNR(Xf, MX, NX, a, b, order, Z, Yf);
     }
     
  } else {  // Signal is neither a DOUBLE nor a SINGLE:
     mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput3",
                       ERR_HEAD "Signal must be a DOUBLE or SINGLE array.");
  }
  
  // Cleanup:
  if (allocate_Z) {
     mxFree(Z);
  }
  if (allocate_ba) {
     mxFree(b);       // Frees [a] implicitely!
  }
  
  return;
}

// =============================================================================
void CopySingleToDouble(double *Z, float *Zf, mwSize N)
{
  // Copy value of SINGLE array to DOUBLE array.
  mwSize i;
  for (i = 0; i < N; i++) {
     Z[i] = (double) Zf[i];
  }
  
  return;
}

// =============================================================================
void NormalizeBA(double *ba, mwSize nParam)
{
  // Normalize filter parameters such that a[0] is 1.0.
  double a0 = ba[nParam];
  mwSize i = 0, f = 2 * nParam;
  
  // Catch division by zero as error:
  if (a0 == 0.0) {
     mexErrMsgIdAndTxt(ERR_ID   "BadValueA",
                       ERR_HEAD "1st element of A cannot be 0.");
  }
        
  while (i < f) {
     ba[i++] /= a0;
  }
  
  return;
}

//  ****************************************************************************
//  ***                               DOUBLE                                 ***
//  ****************************************************************************

void CoreDoubleN(double *X, mwSize MX, mwSize NX, double *a, double *b,
                 mwSize order, double *Z, double *Y)
{
  // Direct form II transposed method for general filter length.
  // Implemented as time domain difference equations.
  // INPUT:
  //   X:  Double array. Operation happens of 1st dimension.
  //   MX: Number of elements in the 1st dimension
  //   NX: Number of columns, considers mutliple dimensions.
  //   a, b: Double vector of filter parameters. Both have nParam elements.
  //       The first element of a is 1.
  //   Z:  DOUBLE array, initial conditions.
  //   nParam: Number of filter parameters, order of filter + 1.
  // OUTPUT:
  //   Z:  DOUBLE array, final conditions.
  //   Y:  Double array, allocated by the caller.
   
  double Xi, Yi;
  mwSize i, j, R;
  
  i = 0;
  while (NX--) {                         // Next slice
     R = i + MX;                         // End of the column
     while (i < R) {
        Xi = X[i];                       // Get signal
        Yi = b[0] * Xi + Z[0];           // Filtered value
        for (j = 1; j < order; j++) {    // Update conditions
           Z[j - 1] = b[j] * Xi + Z[j] - a[j] * Yi;
        }
        Z[order - 1] = b[order] * Xi - a[order] * Yi;
        
        Y[i++] = Yi;                      // Write to output
     }
     Z += order;                          // Next condition vector
  }
  
  return;
}

// =============================================================================
void CoreDouble2(double *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, double *Y)
{
  // Filter with loop unrolled for 2 parameters (filter order 1).
  // Same input as the CoreDoubleN, but ommited [order], because it is 1.
   
  double Xi, Yi, z0, a1 = a[1];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi - a1 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
  }

  return;
}

// =============================================================================
void CoreDouble3(double *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, double *Y)
{
  double Xi, Yi, z0, z1, a1 = a[1], a2 = a[2];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi      - a2 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
  }

  return;
}

// =============================================================================
void CoreDouble4(double *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, double *Y)
{
  double Xi, Yi, z0, z1, z2, a1 = a[1], a2 = a[2], a3 = a[3];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi      - a3 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
  }

  return;
}

// =============================================================================
void CoreDouble5(double *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, double *Y)
{
  double Xi, Yi, z0, z1, z2, z3, a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     z3 = Z[3];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi + z3 - a3 * Yi;
        z3 = b[4] * Xi      - a4 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
     *Z++ = z3;
  }

  return;
}

// =============================================================================
void CoreDouble6(double *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, double *Y)
{
  double Xi, Yi, z0, z1, z2, z3, z4,
         a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     z3 = Z[3];
     z4 = Z[4];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi + z3 - a3 * Yi;
        z3 = b[4] * Xi + z4 - a4 * Yi;
        z4 = b[5] * Xi      - a5 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
     *Z++ = z3;
     *Z++ = z4;
  }

  return;
}

// =============================================================================
void CoreDouble7(double *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, double *Y)
{
  // Still 33% faster than the loop method.
  double Xi, Yi, z0, z1, z2, z3, z4, z5,
         a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5], a6 = a[6];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     z3 = Z[3];
     z4 = Z[4];
     z5 = Z[5];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi + z3 - a3 * Yi;
        z3 = b[4] * Xi + z4 - a4 * Yi;
        z4 = b[5] * Xi + z5 - a5 * Yi;
        z5 = b[6] * Xi      - a6 * Yi;
        Y[i++] = Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
     *Z++ = z3;
     *Z++ = z4;
     *Z++ = z5;
  }
  
  return;
}

// *****************************************************************************
// ***                             DOUBLE REVERSE                            ***
// *****************************************************************************

void CoreDoubleNR(double *X, mwSize MX, mwSize NX, double *a, double *b,
                  mwSize order, double *Z, double *Y)
{
  // Method for general filter length.
  // Signal X is process backwards, but a, b, and  Z have standard direction.
   
  double Xi, Yi;
  mwSize i, j, R;
  
  R = 0;
  while (NX--) {
     i = R + MX - 1;
     while (i >= R) {
        Xi = X[i];
        Yi = b[0] * Xi + Z[0];
        for (j = 1; j < order; j++) {
           Z[j - 1] = b[j] * Xi + Z[j] - a[j] * Yi;
        }
        Z[order - 1] = b[order] * Xi - a[order] * Yi;
        
        Y[i--] = Yi;
     }
     Z += order;
     R += MX;
  }
  
  return;
}


//  ****************************************************************************
//  ***                               SINGLE                                 ***
//  ****************************************************************************

void CoreSingleN(float *X, mwSize MX, mwSize NX, double *a, double *b,
                 mwSize order, double *Z, float *Y)
{
  // Method for general filter length.
  // INPUT:
  //   X:  SINGLE array. Operation happens of 1st dimension.
  //   MX: Number of elements in the 1st dimension
  //   NX: Number of columns, considers mutliple dimensions.
  //   a, b: DOUBLE vector of filter parameters. Both have (order+1) elements.
  //       The first element of a is 1.
  //   Z:  DOUBLE (!) array, initial conditions.
  //   order: Number of filter parameters.
  // OUTPUT:
  //   Z:  DOUBLE (!) array, final conditions.
  //   Y:  SINGLE array, allocated by the caller.
  //
  // All intermediate values are DOUBLEs to increase the accuracy.
  // This is 30% faster than the DOUBLE method.
   
  double Xi, Yi;
  mwSize i, j, R;
  
  i = 0;
  while (NX--) {
     R = i + MX;
     while (i < R) {
        Xi = X[i];
        Yi = b[0] * Xi + Z[0];
        for (j = 1; j < order; j++) {
           Z[j - 1] = b[j] * Xi + Z[j] - a[j] * Yi;
        }
        Z[order - 1] = b[order] * Xi - a[order] * Yi;
        
        Y[i++] = (float) Yi;
     }
     Z += order;
  }
  
  return;
}

// =============================================================================
void CoreSingle2(float *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, float *Y)
{
  // Filter with loop unrolled for 2 parameters (filter order 1).
  // Same input as the CoreSingleN, but ommited order, because it is 1.
   
  double Xi, Yi;
  double z0, a1 = a[1];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi - a1 * Yi;
        Y[i++] = (float) Yi;
     }
     *Z++ = z0;
  }

  return;
}

// =============================================================================
void CoreSingle3(float *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, float *Y)
{
  double Xi, Yi;
  double z0, z1, a1 = a[1], a2 = a[2];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi      - a2 * Yi;
        Y[i++] = (float) Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
  }

  return;
}

// =============================================================================
void CoreSingle4(float *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, float *Y)
{
  double Xi, Yi;
  double z0, z1, z2, a1 = a[1], a2 = a[2], a3 = a[3];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi      - a3 * Yi;
        Y[i++] = (float) Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
  }

  return;
}

// =============================================================================
void CoreSingle5(float *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, float *Y)
{
  double Xi, Yi;
  double z0, z1, z2, z3, a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     z3 = Z[3];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi + z3 - a3 * Yi;
        z3 = b[4] * Xi      - a4 * Yi;
        Y[i++] = (float) Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
     *Z++ = z3;
  }

  return;
}

// =============================================================================
void CoreSingle6(float *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, float *Y)
{
  double Xi, Yi;
  double z0, z1, z2, z3, z4,
         a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     z3 = Z[3];
     z4 = Z[4];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi + z3 - a3 * Yi;
        z3 = b[4] * Xi + z4 - a4 * Yi;
        z4 = b[5] * Xi      - a5 * Yi;
        Y[i++] = (float) Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
     *Z++ = z3;
     *Z++ = z4;
  }

  return;
}


// =============================================================================
void CoreSingle7(float *X, mwSize MX, mwSize NX, double *a, double *b,
                 double *Z, float *Y)
{
  // Still 33% faster than the loop method.
  double Xi, Yi;
  double z0, z1, z2, z3, z4, z5,
         a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5], a6 = a[6];
  mwSize i = 0, C;
  
  while (NX--) {
     z0 = Z[0];
     z1 = Z[1];
     z2 = Z[2];
     z3 = Z[3];
     z4 = Z[4];
     z5 = Z[5];
     C  = i + MX;
     while (i < C) {
        Xi = X[i];
        Yi = b[0] * Xi + z0;
        z0 = b[1] * Xi + z1 - a1 * Yi;
        z1 = b[2] * Xi + z2 - a2 * Yi;
        z2 = b[3] * Xi + z3 - a3 * Yi;
        z3 = b[4] * Xi + z4 - a4 * Yi;
        z4 = b[5] * Xi + z5 - a5 * Yi;
        z5 = b[6] * Xi      - a6 * Yi;
        Y[i++] = (float) Yi;
     }
     *Z++ = z0;
     *Z++ = z1;
     *Z++ = z2;
     *Z++ = z3;
     *Z++ = z4;
     *Z++ = z5;
  }

  return;
}

// *****************************************************************************
// ***                             SINGLE REVERSE                            ***
// *****************************************************************************

void CoreSingleNR(float *X, mwSize MX, mwSize NX, double *a, double *b,
                  mwSize order, double *Z, float *Y)
{
  // Method for general filter length.
  // Signal X is process backwards, but a, b, and  Z have standard direction.
   
  double Xi, Yi;
  mwSize i, j, R;
  
  R = 0;
  while (NX--) {
     i = R + MX - 1;
     while (i >= R) {
        Xi = X[i];
        Yi = b[0] * Xi + Z[0];
        for (j = 1; j < order; j++) {
           Z[j - 1] = b[j] * Xi + Z[j] - a[j] * Yi;
        }
        Z[order - 1] = b[order] * Xi - a[order] * Yi;
        
        Y[i--] = (float) Yi;
     }
     Z += order;
     R += MX;
  }
  
  return;
}
