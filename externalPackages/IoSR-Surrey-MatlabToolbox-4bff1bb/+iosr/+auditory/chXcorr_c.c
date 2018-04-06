/*
 * Calculate cross-correlograms with a wide range of options.
 * Called from chXcorr.m.
 *
 * Copyright 2016 University of Surrey.
 * 
 */

#include "math.h"
#include "mex.h"
#include "matrix.h"

#define INHIB_MULTIPLY 1
#define INHIB_SUBTRACT 2

#define MAX(A,B) ((A)>(B) ? (A) : (B))

#define A_TEMP_ARRAY mxCreateDoubleMatrix((mwSize)lengthCCG,numchans,mxREAL)
#define CCG_TEMP_ARRAY mxCreateNumericArray((mwSize)3,ccg_r_dims,mxDOUBLE_CLASS,mxREAL)

double xcorr(double L, double R, double prev, double tc) {
    return ((1.0/tc)*L*R) + ((1.0-(1.0/tc))*prev);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* ====================== DECLARATIONS ====================== */
	
    double *hc_L = mxGetPr(prhs[0]), /* input left signal */
        *hc_R = mxGetPr(prhs[1]); /* input right signal */
    
    int frameCount = *mxGetPr(prhs[2]), /* Input number of frames in input */
        frame_length = *mxGetPr(prhs[3]), /* Frame length */
        noverlap = *mxGetPr(prhs[4]), /* Number of frames to overlap in CCG calculations */
        maxlag = *mxGetPr(prhs[5]); /* Maximum lag for cross-correlation */
    
    double tau = *mxGetPr(prhs[6]), /* Time constant for cross-correlation calculation */
        *inhib = mxGetPr(prhs[7]), /* Inhibition data */
        ic_t = *mxGetPr(prhs[8]); /* Interaural coherence threshold */
    
    int norm_flag = *mxGetPr(prhs[9]), /* Normalisation flag: 0 = no normalisation, any other value indicates normalisation */
        inhib_mode_ID = *mxGetPr(prhs[10]); /* Inhibition mode: 1 = multiplication, 2 = subtraction */
    
    mwSize  numsamples = mxGetM(prhs[0]), /* Number of samples in HC data */
        numchans = mxGetN(prhs[0]); /* Number of frequency channels in input HC data */
	
    /* The IC for the cross-correlation at the current sample */
    double ic;
    
    int ic_count = 0, /* Counts number of samples over IC_T */
        lengthCCG = (2*maxlag)+1, /* Cross-correlation length */
        lookahead = 0; /* Number of frames to look ahead for CCG calculation */
    
    /* Indices */
    mwIndex i, /* Frequency channel index */
        j, /* Frame index */
        n, /* Sample index */
        m, /* Lag index */
        sample, /* sample index */
        ic_sample, /* IC sample index */
        index_ccg_r, /* Running cross-correlation index */
        index_ccg_r_ahead, /* Running cross-correlation index for look ahead frame */
        index_a, /* Index to auto- and cross-correlations */
        index_ccg, /* Index to CCG */
        leftn, /* Left index for calculating cross-correlations */
        rightn; /* Right index for calculating cross-correlations */
	
    /* Dimensions */
    mwSize ccg_dims[3] = {(mwSize)lengthCCG, numchans, (mwSize)frameCount}, /* CCG dimensions */
        ccg_r_dims[3] = {(mwSize)lengthCCG, (mwSize)frame_length*noverlap, numchans}; /* Running cross-correlation dimensions */
	
    /* ====================== TEMP ARRAYS ====================== */
    
    /* Running auto- and cross-correlations for previous sample */
    mxArray *a_LL_prev_mx,*a_RR_prev_mx,*a_LR_prev_mx;
    a_LL_prev_mx = A_TEMP_ARRAY;
    a_RR_prev_mx = A_TEMP_ARRAY;
    a_LR_prev_mx = A_TEMP_ARRAY;
    double *a_LL_prev = mxGetPr(a_LL_prev_mx),
        *a_RR_prev = mxGetPr(a_RR_prev_mx),
        *a_LR_prev = mxGetPr(a_LR_prev_mx);
    
    /* Running auto- and cross-correlations */
    mxArray *a_LL_mx,*a_RR_mx,*a_LR_mx;
    a_LL_mx = A_TEMP_ARRAY;
    a_RR_mx = A_TEMP_ARRAY;
    a_LR_mx = A_TEMP_ARRAY;
    double  *a_LL = mxGetPr(a_LL_mx),
        *a_RR = mxGetPr(a_RR_mx),
        *a_LR = mxGetPr(a_LR_mx);
    
    /* Cross-correlations */
    mxArray *ccg_norm_mx, *ccg_r_mx, *ccg_ic_mx;
    ccg_norm_mx = CCG_TEMP_ARRAY;
    ccg_r_mx = CCG_TEMP_ARRAY;
    ccg_ic_mx = mxCreateDoubleMatrix((mwSize)lengthCCG,(mwSize)1,mxREAL);
    double ccg_ic_temp,
        *ccg_norm = mxGetPr(ccg_norm_mx),
        *ccg_r = mxGetPr(ccg_r_mx),
        *ccg_ic = mxGetPr(ccg_ic_mx);	
    
    /* ====================== OUTPUTS ====================== */

    plhs[0] = mxCreateNumericArray(3,ccg_dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(numsamples,numchans,mxREAL);
    double *ccg_out = mxGetPr(plhs[0]),
        *ic_out = mxGetPr(plhs[1]);

    /* ====================== CROSS-CORRELATE ====================== */
    for ( j = 0; j < frameCount; j++ ) {
        for ( i = 0; i < numchans; i++ ) {
            index_ccg = (i*lengthCCG)+(j*lengthCCG*numchans);
            sample = (i*(mwIndex)numsamples)+((j+lookahead)*(mwIndex)frame_length);
            ic_sample = (i*(mwIndex)numsamples)+(j*(mwIndex)frame_length);
            for ( m = 0; m < lengthCCG; m++ ) {/* reset ccg for T-F unit */
                ccg_ic[m] = 0.0;
            }
            /* Calculate raw cross-correlations */
            for ( n = 0; n < frame_length*(noverlap-lookahead); n++ ) {
                index_a = (i*lengthCCG);
                index_ccg_r = ((n+(lookahead*frame_length))+i*(noverlap*frame_length))*lengthCCG;
                ic = 0.0; /* reset the IC value for the sample */
                for ( m = 0; m < lengthCCG; m++ ) { /* calculate cross-correlations */
                    leftn = sample+n+( (m>maxlag) ? (m-maxlag) : (0) );
                    rightn = sample+n+( (m<maxlag) ? (maxlag-m) : (0) );
                    a_LR[index_a+m] = xcorr(hc_L[leftn],hc_R[rightn],a_LR_prev[index_a+m],tau);
                    a_LL[index_a+m] = xcorr(hc_L[leftn],hc_L[leftn],a_LL_prev[index_a+m],tau);
                    a_RR[index_a+m] = xcorr(hc_R[rightn],hc_R[rightn],a_RR_prev[index_a+m],tau);
                    ccg_norm[index_ccg_r+m] = a_LR[index_a+m]/sqrt(a_LL[index_a+m]*a_RR[index_a+m]);
                    if ( (a_LL[index_a+m]*a_RR[index_a+m])==0.0 ) { /* prevent divide by zero */
                        ccg_norm[index_ccg_r+m] = 0.0;
                    }
                    ic = MAX(ic,ccg_norm[index_ccg_r+m]); /* calculate IC for sample (as max) */
                    if ( norm_flag == 0 ) {
                        /* Un-normalised */
                        ccg_r[index_ccg_r+m] = a_LR[index_a+m];
                    }
                    else {
                        /* Normalised */
                        ccg_r[index_ccg_r+m] = ccg_norm[index_ccg_r+m];
                    }
                    /* Set cross-correlations for what will become the previous sample */
                    a_LR_prev[index_a+m] = a_LR[index_a+m];
                    a_LL_prev[index_a+m] = a_LL[index_a+m];
                    a_RR_prev[index_a+m] = a_RR[index_a+m];
                }
                ic_out[sample+n] = ic; /* Write IC to output */
            }
            ic_count = 0; /* reset the count of samples with IC values over the IC threshold */
            /* Integrate cross-correlations */
            for ( n = 0; n < frame_length*noverlap; n++ ) {
                index_ccg_r = (n+i*(frame_length*noverlap))*lengthCCG;
                if ( ic_out[ic_sample+n] >= ic_t ) { /* check if IC exceeds IC threshold */
                    ic_count++; /* count the number of samples that exceed the IC threshold */
                    for ( m = 0; m < lengthCCG; m++ ) { /* inhibit (multiplication) */
                        switch (inhib_mode_ID) {
                            case INHIB_MULTIPLY:
                                ccg_ic[m] += ccg_r[index_ccg_r+m]*inhib[sample+n];
                                break;
                            case INHIB_SUBTRACT:
                                ccg_ic_temp = ccg_r[index_ccg_r+m]-((1.0/tau)*inhib[sample+n]);
                                ccg_ic[m] += MAX(ccg_ic_temp,0.0);
                                break;
                            default:
                                    mexErrMsgTxt("Unknown inhib_mode_ID");
                        }
                    }
                }
                else {
                    ic_out[ic_sample+n] = 0;
                }
            }
            if ( ic_count == 0 ) { /* write zeros to output if no samples in frame have high enough IC */
                for ( m = 0; m < lengthCCG; m++ ) {
                    ccg_out[index_ccg+m] = 0.0;
                }
            }
            else { /* average CCGs with high enough IC */
                for ( m = 0; m < lengthCCG; m++ ) {
                    ccg_out[index_ccg+m] = ccg_ic[m]/(double)ic_count; /* Write CCG to output */
                }
            }
            /* move ccg_r's backwards to append for subsequent frames */
            for ( n = 0; n < frame_length*(noverlap-1); n++ ) {
                index_ccg_r = (n+i*(frame_length*noverlap))*lengthCCG;
                index_ccg_r_ahead = ((n+frame_length)+i*(frame_length*noverlap))*lengthCCG;
                for ( m = 0; m < lengthCCG; m++ ) {
                    ccg_r[index_ccg_r+m] = ccg_r[index_ccg_r_ahead+m];
                }
            }
        } /* end frequency loop */
        lookahead = noverlap-1; /* Set value */
    } /* end frame loop */
    /* Destroy mx arrays */
    mxDestroyArray(ccg_norm_mx);
    mxDestroyArray(ccg_r_mx);
    mxDestroyArray(ccg_ic_mx);
    mxDestroyArray(a_LL_mx);
    mxDestroyArray(a_RR_mx);
    mxDestroyArray(a_LR_mx);
    mxDestroyArray(a_LL_prev_mx);
    mxDestroyArray(a_RR_prev_mx);
    mxDestroyArray(a_LR_prev_mx);
    return;
}
