/*
 * Calculate cross-correlograms with a range of options.
 * Called from chXcorr2.m.
 *
 * Copyright 2016 University of Surrey.
 * 
 */

#include "math.h"
#include "mex.h"
#include "matrix.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

/* structure for cross-correlation info */
struct xcorr
{
    int maxlag;
    double ic;
    double *c;
    double *nc;
    double *aL;
    double *aR;
};

/* cross-correlation functions */
void xcorr(double L[], double R[], int frame_length, int numsamples, struct xcorr * data);
void initXcorr(unsigned int maxlag, unsigned int length, struct xcorr * data);
void freeXcorr(struct xcorr * data);

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* ====================== DECLARATIONS ====================== */
	
    double
        *hc_L = mxGetPr(prhs[0]), /* input left signal */
        *hc_R = mxGetPr(prhs[1]), /* input right signal */
        *cout; /* output */
    
    int
        frame_count = *mxGetPr(prhs[2]), /* Input number of frames in input */
        frame_length = *mxGetPr(prhs[3]), /* Frame length */
        maxlag = *mxGetPr(prhs[4]), /* Maximum lag for cross-correlation */
        hop = *mxGetPr(prhs[5]), /* hop size */
        norm_flag = *mxGetPr(prhs[6]), /* normalisation flag */
        lengthCCG = (2*maxlag)+1; /* Cross-correlation length */
    
    mwSize
        numsamples = mxGetM(prhs[0]), /* Number of samples in HC data */
        numchans = mxGetN(prhs[0]); /* Number of frequency channels in input HC data */
    
    /* Indices */
    mwIndex
        m, /* lag index */
        i, /* Frequency channel index */
        j, /* Frame index */
        sample, /* sample index */
        index_ccg, /* Index to CCG */
        index_ic; /* index into interaural coherence */
        
    /* cross-correlation data */
    struct xcorr ccdata;
    initXcorr(maxlag, lengthCCG, &ccdata);
    
    /* choose output according to normalisation */
    if (norm_flag == 0) {
        cout = ccdata.c;
    }
    else {
        cout = ccdata.nc;
    }
    
    /* ====================== OUTPUTS ====================== */

    mwSize ccg_dims[3] = {(mwSize)lengthCCG, numchans, (mwSize)frame_count}; /* CCG dimensions */
    plhs[0] = mxCreateNumericArray(3,ccg_dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(frame_count,numchans,mxREAL);
    double *ccg_out = mxGetPr(plhs[0]),
        *ic_out = mxGetPr(plhs[1]);

    /* ====================== CROSS-CORRELATE ====================== */
    for ( j = 0; j < frame_count; j++ ) {
        for ( i = 0; i < numchans; i++ ) {
            /* indices */
            sample = (i*(mwIndex)numsamples)+(j*(mwIndex)hop);
            index_ccg = (i*lengthCCG)+(j*lengthCCG*numchans);
            index_ic = (i*(mwIndex)frame_count)+j;
            /* cross-correlate */
            xcorr(hc_L+sample,  hc_R+sample, frame_length, numsamples, &ccdata);
            /* outputs */
            ic_out[index_ic] = ccdata.ic;
            for ( m = 0; m < lengthCCG; m++ ) {
                ccg_out[index_ccg+m] = cout[m];
            }
        } /* end frequency loop */
    } /* end frame loop */
    /* Destroy arrays */
    freeXcorr(&ccdata);
    return;
}

void xcorr(double L[], double R[], int frame_length, int numsamples, struct xcorr * data)
{
    
    int m, /* lag index */
        n, /* sample index */
        maxlag = data->maxlag; /* maximum cross-correlation lag */
    double
        Rshift, /* shifted right signal */
        mL = 0.0, /* mean of left */
        mR = 0.0, /* mean of right */
        denom, /* denominator */
        *c = data->c, /* cross-correlation */
        *nc = data->nc, /* normalised cross-correlation */
        *aL = data->aL, /* left variance */
        *aR = data->aR, /* right variance */
        *ic = &(data->ic); /* interaural coherence */

    /* reset arrays to 0 */
    for (m = 0; m < 2*maxlag+1; m++) {
        c[m] = 0.0;
        nc[m] = 0.0;
        aL[m] = 0.0;
        aR[m] = 0.0;
    }
    
    /* calculate mean of time series */
    for (n = 0; n < frame_length; n++) {
        mL += L[n]/frame_length;
        mR += R[n]/frame_length;
    }
    
    /* calculate cross-correlations of time series */
    for (m = -maxlag; m <= maxlag; m++) {
        for (n = 0; n < frame_length; n++) {
            /* do not wrap, but set out-of-bounds to zero */
            if (((n-m) < 0) || ((n-m) > frame_length) || ((n-m) > numsamples)) {
                Rshift = 0.0;
            }
            else {
                Rshift = R[n-m];
            }
            /* cross-correlate */
            c[m+maxlag] += (L[n]-mL) * (Rshift-mR); /* cross-correlation */
            aL[m+maxlag] += pow(L[n]-mL,2.0); /* variance of left */
            aR[m+maxlag] += pow(Rshift-mR,2.0); /* variance of right */
        }
    }
    
    /* normalise the cross-correlation */
    data->ic = 0.0; /* reset to 0 */
    for (m = 0; m < 2*maxlag+1; m++) {
        denom = sqrt(aL[m])*sqrt(aR[m]); /* denominator */
        if (denom > 0.0) {
            nc[m] = c[m]/denom;
        }
        else { /* prevent divide by zero */
            nc[m] = 0.0;
        }
        /* calculate coherence */
        *ic = max(*ic,nc[m]);
    }
    
    return;
    
}

void initXcorr(unsigned int maxlag, unsigned int length, struct xcorr * data)
{
    /* initialise and allocate data */
    data->maxlag = maxlag;
    data->ic = 0.0;
    data->c = calloc(length,sizeof(double));
    data->nc = calloc(length,sizeof(double));
    data->aL = calloc(length,sizeof(double));
    data->aR = calloc(length,sizeof(double));
    if (length != 0 && (!data->c || !data->nc || !data->aL || !data->aR)) {
        freeXcorr(data);
    }
}

void freeXcorr(struct xcorr * data)
{
    /* free memory */
    if (data->c)
        free(data->c);
    if (data->nc)
        free(data->nc);
    if (data->aL)
        free(data->aL);
    if (data->aR)
        free(data->aR);
}
