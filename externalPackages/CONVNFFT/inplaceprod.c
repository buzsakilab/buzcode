/*************************************************************************
 * MATLAB MEX ROUTINE inplaceprod.c
 *
 * > inplaceprod(A,B) 
 * performs inplace (i.e., without allocate array) complex array product
 * > A(:) = A(:).*B(:);
 *
 * User must make sure A is not shared by other array.
 * B must have the same size as A, no check will be carried out
 *
 * A, B must be of class double or single (mixing is allowed)
 *
 * >> mex -O inplaceprod.c
 * For more recent MATLAB (R2018a) use
 * >> mex -O -R2018a inplaceprod.c
 *
 * Author: Bruno Luong <brunoluong@yahoo.com>
 * History
 *  Original: 16/Sept/2009
 ************************************************************************/

#include "mex.h"
#include "matrix.h"

#define A (mxArray*)(prhs[0])
#define B prhs[1]

#if MX_HAS_INTERLEAVED_COMPLEX
#define REC_SHIFT 2
#define ENGINE_RACB(Ar, Ai, Br, Bi, temp, Atype, Btype) { \
mexErrMsgTxt("INPLACEPROD: A must be complex."); \
}
#else
#define REC_SHIFT 1
#define ENGINE_RACB(Ar, Ai, Br, Bi, temp, Atype, Btype) { \
Ai = (Atype*)mxMalloc(n*sizeof(Atype)); \
mxSetImagData(A, (void*)Ai); \
if (Ai==NULL) \
    mexErrMsgTxt("INPLACEPROD: cannot allocate memory."); \
for (i=0; i<n; i++) { \
    Ai[i] = (Atype)(Ar[i]*Bi[i]); \
    Ar[i] *= (Atype)Br[i]; \
} \
}
#endif

/* Product engine */
#define ENGINE(Ar, Ai, Br, Bi, temp, Atype, Btype) { \
    if (Ai != NULL && Bi != NULL) { \
        for (i=0; i<n; i++) { \
            temp = *Ar; \
            *Ar = (Atype)(temp**Br - *Ai**Bi); \
            *Ai = (Atype)(temp**Bi + *Ai**Br); \
            Ar += REC_SHIFT; \
            Ai += REC_SHIFT; \
            Br += REC_SHIFT; \
            Bi += REC_SHIFT; \
        } \
    } \
    else if (Ai != NULL) { \
        for (i=0; i<n; i++) { \
            *Ar *= (Atype)*Br; \
            *Ai *= (Atype)*Br; \
            Ar += REC_SHIFT; \
            Ai += REC_SHIFT; \
            Br++; \
        } \
    } \
    else if (Bi != NULL) { \
           ENGINE_RACB(Ar, Ai, Br, Bi, temp, Atype, Btype); \
    } \
    else { \
        for (i=0; i<n; i++) { \
            Ar[i] *= (Atype)Br[i]; \
        } \
    } \
}
        
/* Gateway of inplaceprod */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    size_t i, n;
#if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble *Azd, *Bzd;
    mxComplexSingle *Azs, *Bzs;    
#endif
    double  *Ard, *Aid, *Brd, *Bid;
    float   *Ars, *Ais, *Brs, *Bis;    
    double  tempd;
    float   temps;
    mxClassID ClassA, ClassB;
    
    /* Get the classes of A and B */
    ClassA = mxGetClassID(A);
    ClassB = mxGetClassID(B);
    
    /* n = numel(A) */
    n = mxGetNumberOfElements(A);
    
#if MX_HAS_INTERLEAVED_COMPLEX
    if (ClassA == mxDOUBLE_CLASS)
    {
        if (mxIsComplex(A))
        {
            Azd = mxGetComplexDoubles(A);
            Ard = (double*)Azd;
            Aid = Ard+1;
        }
        else
        {
            Ard = mxGetDoubles(A);
            Aid = NULL;
        }
    }
    else
    {
        if (mxIsComplex(A))
        {
            Azs = mxGetComplexSingles(A);
            Ars = (float*)Azs;
            Ais = Ars+1;
        }
        else
        {
            Ars = mxGetSingles(A);
            Ais = NULL;
        }
    }
    if (ClassB == mxDOUBLE_CLASS) 
    {
        if (mxIsComplex(B))
        {
            Bzd = mxGetComplexDoubles(B);
            Brd = (double*)Bzd;
            Bid = Brd+1;
        }
        else
        {
            Brd = mxGetDoubles(B);
            Bid = NULL;
        }
    }
    else
    {
        if (mxIsComplex(B))
        {
            Bzs = mxGetComplexSingles(B);
            Brs = (float*)Bzs;
            Bis = Brs+1;
        }
        else
        {
            Brs = mxGetSingles(B);
            Bis = NULL;
        }
    }    
#else
    if (ClassA == mxDOUBLE_CLASS) {
        Aid = (double*)mxGetImagData(A);
    } else {
        Ais = (float*)mxGetImagData(A);
    }
    if (ClassB == mxDOUBLE_CLASS) {
        Bid = (double*)mxGetImagData(B);
    } else {
        Bis = (float*)mxGetImagData(B);
    }
#endif
    
    /* Call the macros depending of the classes of A and B */
    if ((ClassA == mxDOUBLE_CLASS) && (ClassB == mxDOUBLE_CLASS))
    {
        ENGINE(Ard, Aid, Brd, Bid, tempd, double, double);
    } else if (ClassA == mxDOUBLE_CLASS)
    {
        ENGINE(Ard, Aid, Brs, Bis, tempd, double, float);
    } else if (ClassB == mxDOUBLE_CLASS)
    {
        ENGINE(Ars, Ais, Brd, Bid, temps, float, double);
    } else
    {
        ENGINE(Ars, Ais, Brs, Bis, temps, float, float);
    }
    
    return;
} /* of inplaceprod */