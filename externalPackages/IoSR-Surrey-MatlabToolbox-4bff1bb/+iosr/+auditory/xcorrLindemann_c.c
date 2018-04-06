/* Lindemann's precedence model.
 * 
 * Copyright 2016 University of Surrey.
 * 
 */

#include "math.h"
#include "mex.h"

#define TEMP_ARRAY mxCreateDoubleMatrix((mwSize)length_c,(mwSize)1,mxREAL)
#define IN_ARRAY mxCreateDoubleMatrix((mwSize)n_sample,(mwSize)1,mxREAL)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* ====================== INPUTS & SCALARS ====================== */
    
    double *L = mxGetPr(prhs[0]), /* pointer to left signal input */
        *R = mxGetPr(prhs[1]), /* pointer to right signal input */
        fs = *mxGetPr(prhs[2]); /* samle frequency */
    
    mwSize maxlag = *mxGetPr(prhs[3]); /* maximum lag */
    
    double tinh = (*mxGetPr(prhs[4]))/1000.0, /* fade-off time constant */
        tint = (*mxGetPr(prhs[5]))/1000.0, /* integration time constant */
        td = 0.5*(1.0/fs), /* half sample period */
        alpha = exp(-td/tinh), /* fade-off factor */
        alpha2 = exp(-td/tint), /* integration factor */
        Mf = 6.0, /* fading constant for monaural sensitivity */
        wf = 0.035; /* monaural sensitivity of a correlator to the signal at the end of a delay line */
    
    /* sizes */
    mwSize numsamples = mxGetM(prhs[0]), /* number of samples */
        length_c = (2*maxlag)+1, /* length of cross-correlation (2*maxlag)+1 */
        in_length = length_c+numsamples-1, /* length of input to cross-correlation */
        n_sample = 2*in_length-1; /* length of inhibited input to cross-correlation */; 
    
    /* indices */
    mwIndex n, /* sample */
        m, /* lag */
        leftn, /* left samples to be cross-correlated */
        rightn; /* right samples to be cross-correlated */
    
    /* ====================== TEMP ARRAYS ====================== */
        
    /* monaural sensitivity functions */
    mxArray *wl_mx,*wr_mx;
    wl_mx = TEMP_ARRAY;
    wr_mx = TEMP_ARRAY;
    double *wl = mxGetPr(wl_mx),
        *wr = mxGetPr(wr_mx);
    for ( m = 0; m < length_c; m++ ) {
        wl[m] = wf*exp(-((double)m)/Mf);
        wr[m] = wf*exp(((double)m-(double)length_c+1.0)/Mf);
    }
    
    /* running cross-correlation */
    mxArray *k_mx;
    k_mx = TEMP_ARRAY;
    double *k = mxGetPr(k_mx);
    
    /* inhibited input to the cross-correlation */
    mxArray *lc_mx,*rc_mx;
    lc_mx = TEMP_ARRAY;
    rc_mx = TEMP_ARRAY;
    double *lc = mxGetPr(lc_mx),
        *rc = mxGetPr(rc_mx);
    
    /* input to cross-correlation */
    mxArray *lin_mx,*rin_mx;
    lin_mx = IN_ARRAY;
    rin_mx = IN_ARRAY;
    double *lin = mxGetPr(lin_mx),
        *rin = mxGetPr(rin_mx);
    
    /* inhibitory components */
    mxArray *il_mx,*ir_mx;
    il_mx = TEMP_ARRAY;
    ir_mx = TEMP_ARRAY;
    double *il = mxGetPr(il_mx),
        *ir = mxGetPr(ir_mx);
    
    /* dynamic inhibitory component for current and previous samples */
    mxArray *phin_mx,*phic_mx;
    phin_mx = TEMP_ARRAY;
    phic_mx = TEMP_ARRAY;
    double *phin = mxGetPr(phin_mx),
        *phic = mxGetPr(phic_mx);
    
    /* inhibited cross-correlation for current and previous samples */
    mxArray *sumn_mx,*sumc_mx;
    sumn_mx = TEMP_ARRAY;
    sumc_mx = TEMP_ARRAY;
    double *sumn = mxGetPr(sumn_mx),
        *sumc = mxGetPr(sumc_mx);
        
    /* ====================== OUTPUT ====================== */

    plhs[0] = TEMP_ARRAY;
    double *c_out = mxGetPr(plhs[0]);
        
    /* ====================== CROSS-CORRELATE ====================== */
        
    for ( n = 0 ; n < numsamples ; n++ ) {
        /* input to cross-correlation */
        lin[(2*n)] = L[n];
        rin[(2*n)] = R[n];
    }
    for ( n = 0; n < n_sample; n++ ) {
        for ( m = 0; m < 2*maxlag; m++ ) {
            /* inhibit input to cross-correlation */
            lc[m] = lc[m+1]*il[m+1];
            rc[(2*maxlag)-m] = rc[(2*maxlag)-m-1]*ir[(2*maxlag)-m-1];
        }
        lc[length_c-1] = lin[n];
        rc[0] = rin[n];
        /* cross-correlate */
        for ( m = 0; m < length_c; m++ ) {
            k[m] = (wl[m]+(1.0-wl[m])*rc[m])*(wr[m]+(1.0-wr[m])*lc[m]);
            phin[m] = k[m]+alpha*phic[m]*(1-k[m]);
            phic[m]=phin[m];
            il[m]=(1.0-rc[m])*(1-phic[m]);
            ir[m]=(1.0-lc[m])*(1-phic[m]);
            sumn[m]=sumc[m]*alpha2+(1.0-alpha2)*k[m];
            sumc[m]=sumn[m];
            c_out[m] += sumn[m];
        }
    }
    /* Destroy mx arrays */
    mxDestroyArray(wl_mx);
    mxDestroyArray(wr_mx);
    mxDestroyArray(k_mx);
    mxDestroyArray(lc_mx);
    mxDestroyArray(rc_mx);
    mxDestroyArray(il_mx);
    mxDestroyArray(ir_mx);
    mxDestroyArray(phin_mx);
    mxDestroyArray(phic_mx);
    mxDestroyArray(sumn_mx);
    mxDestroyArray(sumc_mx);
    mxDestroyArray(lin_mx);
    mxDestroyArray(rin_mx);
    return;
}
