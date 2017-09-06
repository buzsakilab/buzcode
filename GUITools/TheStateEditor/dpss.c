/*
 * MATLAB Compiler: 2.0.1
 * Date: Fri Feb  4 11:20:20 2000
 * Arguments: "-mw" "mtcsd" 
 */
#include "dpss.h"
#include "dpssdir.h"
#include "dpssload.h"
#include "fftfilt.h"
#include "sinc.h"
#include "signal_private_tridieig.h"
#include "signal_private_tridisolve.h"

static double __Array0_r[1] = { 0.0 };

static double __Array1_r[2] = { 2.0, 3.0 };

static mxArray * Mdpss_dpsscalc(mxArray * * V,
                                int nargout_,
                                mxArray * N,
                                mxArray * NW,
                                mxArray * k);
static mxArray * mlfDpss_dpsscalc(mxArray * * V,
                                  mxArray * N,
                                  mxArray * NW,
                                  mxArray * k);
static void mlxDpss_dpsscalc(int nlhs,
                             mxArray * plhs[],
                             int nrhs,
                             mxArray * prhs[]);
static mxArray * Mdpss_dpssint(mxArray * * V,
                               int nargout_,
                               mxArray * N,
                               mxArray * NW,
                               mxArray * M,
                               mxArray * int0,
                               mxArray * E,
                               mxArray * V_);
static mxArray * mlfDpss_dpssint(mxArray * * V,
                                 mxArray * N,
                                 mxArray * NW,
                                 mxArray * M,
                                 mxArray * int0,
                                 mxArray * E,
                                 mxArray * V_);
static void mlxDpss_dpssint(int nlhs,
                            mxArray * plhs[],
                            int nrhs,
                            mxArray * prhs[]);
static mxArray * Mdpss_parseinputs(mxArray * * k,
                                   mxArray * * Ni,
                                   mxArray * * traceFlag,
                                   mxArray * * err_msg,
                                   int nargout_,
                                   mxArray * N,
                                   mxArray * NW,
                                   mxArray * varargin);
static mxArray * mlfDpss_parseinputs(mxArray * * k,
                                     mxArray * * Ni,
                                     mxArray * * traceFlag,
                                     mxArray * * err_msg,
                                     mxArray * N,
                                     mxArray * NW,
                                     ...);
static void mlxDpss_parseinputs(int nlhs,
                                mxArray * plhs[],
                                int nrhs,
                                mxArray * prhs[]);
static mxArray * Mdpss_parseStrOpts(mxArray * * traceFlag,
                                    mxArray * * err_msg,
                                    int nargout_,
                                    mxArray * inputStr,
                                    mxArray * method_,
                                    mxArray * traceFlag_,
                                    mxArray * nopts,
                                    mxArray * opt_indx);
static mxArray * mlfDpss_parseStrOpts(mxArray * * traceFlag,
                                      mxArray * * err_msg,
                                      mxArray * inputStr,
                                      mxArray * method_,
                                      mxArray * traceFlag_,
                                      mxArray * nopts,
                                      mxArray * opt_indx);
static void mlxDpss_parseStrOpts(int nlhs,
                                 mxArray * plhs[],
                                 int nrhs,
                                 mxArray * prhs[]);

/*
 * The function "Mdpss" is the implementation version of the "dpss" M-function
 * from file "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 1-89). It
 * contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */
/*
 * function [E,V]=dpss(N, NW, varargin)
 */
static mxArray * Mdpss(mxArray * * V,
                       int nargout_,
                       mxArray * N,
                       mxArray * NW,
                       mxArray * varargin) {
    mxArray * E = mclGetUninitializedArray();
    mxArray * Ni = mclGetUninitializedArray();
    mxArray * Nlist = mclGetUninitializedArray();
    mxArray * dum = mclGetUninitializedArray();
    mxArray * err_msg = mclGetUninitializedArray();
    mxArray * i = mclGetUninitializedArray();
    mxArray * ind = mclGetUninitializedArray();
    mxArray * k = mclGetUninitializedArray();
    mxArray * method = mclGetUninitializedArray();
    mxArray * nargin_ = mclGetUninitializedArray();
    mxArray * traceFlag = mclGetUninitializedArray();
    mlfAssign(&nargin_, mlfNargin(1, N, NW, varargin, NULL));
    mclValidateInputs("dpss", 2, &N, &NW);
    /*
     * %DPSS   Discrete prolate spheroidal sequences (Slepian sequences).
     * %   [E,V] = DPSS(N,NW) are the first 2*NW discrete prolate spheroidal sequences
     * %   (DPSSs, or Slepian sequences) of length N (in the columns of E) and 
     * %   their corresponding concentrations (in vector V) in the frequency band 
     * %   |w|<=(2*pi*W) (where  W = NW/N is the half-bandwidth and w is in 
     * %   radians/sample).  E(:,1) is the length N signal most concentrated in the 
     * %   frequency band |w|<=(2*pi*W) radians/sample, E(:,2) is the signal 
     * %   orthogonal to E(:,1) which is most concentrated in this band, E(:,3) is the
     * %   signal orthogonal to both E(:,1) and E(:,2) which is most concentrated in 
     * %   this band, etc.  
     * %
     * %   For multi-taper spectral analysis, typical choices for NW are 2, 5/2, 3, 
     * %   7/2, or 4.
     * %
     * %   [E,V] = DPSS(N,NW,K) are the K most band-limited discrete prolate spheroidal
     * %   sequences.  [E,V] = DPSS(N,NW,[K1 K2]) returns the K1-th through the 
     * %   K2-th sequences.
     * %
     * %   [E,V] = DPSS(N,NW,'spline') uses spline interpolation to compute the DPSSs 
     * %   from existing DPSSs in the DPSS database with length closest to N.
     * %   [E,V] = DPSS(N,NW,'spline',Ni) interpolates from existing length Ni DPSSs.
     * %   DPSS(N,NW,'linear') and DPSS(N,NW,'linear',Ni) use linear interpolation, 
     * %   which is much faster but less accurate than spline interpolation.  
     * %   'linear' requires Ni > N. [E,V] = DPSS(N,NW,'calc') uses the direct 
     * %   algorithm (default).
     * %
     * %   Use a trailing 'trace' argument to find out which method DPSS uses, e.g.,
     * %   DPSS(...,'trace').
     * %
     * %   See also PMTM, DPSSLOAD, DPSSDIR, DPSSSAVE, DPSSCLEAR.
     * 
     * %   References: 
     * %     [1] Percival, D.B. and Walden, A.T., "Spectral Analysis For Physical
     * %         Applications", Cambridge University Press, 1993. 
     * 
     * %   Author: Eric Breitenberger, 10/3/95
     * %   Copyright (c) 1988-1999 The MathWorks, Inc. All Rights Reserved.
     * %   $Revision: 1.6 $  $Date: 1999/05/28 18:45:21 $
     * 
     * % Input parsing and validation
     * error(nargchk(2,6,nargin));
     */
    mlfError(mlfNargchk(mlfScalar(2.0), mlfScalar(6.0), nargin_));
    /*
     * [method,k,Ni,traceFlag,err_msg] = parseinputs(N,NW,varargin{:});
     */
    mlfAssign(
      &method,
      mlfDpss_parseinputs(
        &k,
        &Ni,
        &traceFlag,
        &err_msg,
        N,
        NW,
        mlfIndexRef(varargin, "{?}", mlfCreateColonIndex()),
        NULL));
    /*
     * error(err_msg);
     */
    mlfError(err_msg);
    /*
     * 
     * switch method
     * 
     * case 'calc'
     */
    if (mclSwitchCompare(method, mxCreateString("calc"))) {
        /*
         * if traceFlag,
         */
        if (mlfTobool(traceFlag)) {
            /*
             * disp('Computing the DPSS using direct algorithm...')
             */
            mlfDisp(
              mxCreateString("Computing the DPSS using direct algorithm..."));
        /*
         * end
         */
        }
        /*
         * [E,V] = dpsscalc(N,NW,k);
         */
        mlfAssign(&E, mlfDpss_dpsscalc(V, N, NW, k));
    /*
     * 
     * case {'spline','linear'}
     */
    } else if (mclSwitchCompare(
                 method,
                 mlfCellhcat(
                   mxCreateString("spline"), mxCreateString("linear"), NULL))) {
        /*
         * err_msg = '';
         */
        mlfAssign(&err_msg, mxCreateString(""));
        /*
         * if isempty(Ni)
         */
        if (mlfTobool(mlfIsempty(Ni))) {
            /*
             * ind = dpssdir(NW,'NW');
             */
            mlfAssign(&ind, mlfNDpssdir(1, NW, mxCreateString("NW")));
            /*
             * if ~isempty(ind)
             */
            if (mlfTobool(mlfNot(mlfIsempty(ind)))) {
                /*
                 * Nlist = [ind.N];
                 */
                mlfAssign(&Nlist, mlfHorzcat(mlfIndexRef(ind, ".N"), NULL));
                /*
                 * % find closest length and use that one
                 * [dum,i] = min(abs(N-Nlist));
                 */
                mlfAssign(
                  &dum, mlfMin(&i, mlfAbs(mlfMinus(N, Nlist)), NULL, NULL));
                /*
                 * Ni = Nlist(i);
                 */
                mlfAssign(&Ni, mlfIndexRef(Nlist, "(?)", i));
                /*
                 * if strcmp(method,'linear') & Ni<N
                 */
                {
                    mxArray * a_ = mclInitialize(
                                     mlfStrcmp(
                                       method, mxCreateString("linear")));
                    if (mlfTobool(a_) && mlfTobool(mlfAnd(a_, mlfLt(Ni, N)))) {
                        mxDestroyArray(a_);
                        /*
                         * if i<length(Nlist)
                         */
                        if (mlfTobool(mlfLt(i, mlfLength(Nlist)))) {
                            /*
                             * Ni = Nlist(i+1);
                             */
                            mlfAssign(
                              &Ni,
                              mlfIndexRef(
                                Nlist, "(?)", mlfPlus(i, mlfScalar(1.0))));
                        /*
                         * else
                         */
                        } else {
                            /*
                             * Ni = [];
                             */
                            mlfAssign(&Ni, mclCreateEmptyArray());
                            /*
                             * err_msg = sprintf('There is no DPSS with time-bandwidth product NW = %g and N > %g in database.',NW,N);
                             */
                            mlfAssign(
                              &err_msg,
                              mlfSprintf(
                                NULL,
                                mxCreateString(
                                  "There is no DPSS with time-bandwidth pr"
                                  "oduct NW = %g and N > %g in database."),
                                NW,
                                N,
                                NULL));
                        /*
                         * end
                         */
                        }
                    } else {
                        mxDestroyArray(a_);
                    }
                /*
                 * end
                 */
                }
            /*
             * else
             */
            } else {
                /*
                 * err_msg = sprintf('There is no DPSS with time-bandwidth product NW = %g in database.',NW);
                 */
                mlfAssign(
                  &err_msg,
                  mlfSprintf(
                    NULL,
                    mxCreateString(
                      "There is no DPSS with time-bandwid"
                      "th product NW = %g in database."),
                    NW,
                    NULL));
            /*
             * end
             */
            }
        /*
         * end
         */
        }
        /*
         * error(err_msg)
         */
        mlfError(err_msg);
        /*
         * 
         * if traceFlag,
         */
        if (mlfTobool(traceFlag)) {
            /*
             * disp(['Computing DPSS using ' method ' interpolation from length ' int2str(Ni) '...'])
             */
            mlfDisp(
              mlfHorzcat(
                mxCreateString("Computing DPSS using "),
                method,
                mxCreateString(" interpolation from length "),
                mlfInt2str(Ni),
                mxCreateString("..."),
                NULL));
        /*
         * end
         */
        }
        /*
         * 
         * [E,V]=dpssint(N,NW,Ni,method);
         */
        mlfAssign(&E, mlfDpss_dpssint(V, N, NW, Ni, method, NULL, NULL));
    /*
     * 
     * otherwise
     */
    } else {
        /*
         * error('Method string should be ''calc'', ''spline'' or ''linear''.')
         */
        mlfError(
          mxCreateString(
            "Method string should be 'calc', 'spline' or 'linear'."));
    /*
     * 
     * end
     */
    }
    mclValidateOutputs("dpss", 2, nargout_, &E, V);
    mxDestroyArray(Ni);
    mxDestroyArray(Nlist);
    mxDestroyArray(dum);
    mxDestroyArray(err_msg);
    mxDestroyArray(i);
    mxDestroyArray(ind);
    mxDestroyArray(k);
    mxDestroyArray(method);
    mxDestroyArray(nargin_);
    mxDestroyArray(traceFlag);
    /*
     * 
     * %----------------------------------------------------------------------
     * function [E,V] = dpsscalc(N,NW,k)
     * %DPSSCALC Calculate slepian sequences.
     * %   [E,V] = dpsscalc(N,NW,k) uses tridieig() to get eigenvalues 1:k if k is 
     * %   a scalar, and k(1):k(2) if k is a matrix, of the sparse tridiagonal matrix.
     * %   It then uses inverse iteration using the exact eigenvalues on a starting 
     * %   vector with approximate shape, to get the eigenvectors required.  It then 
     * %   computes the eigenvalues V of the Toeplitz sinc matrix using a fast 
     * %   autocorrelation technique.
     * 
     * %   Authors: T. Krauss, C. Moler, E. Breitenberger
     * W=NW/N;
     * if nargin < 3
     * k = min(round(2*N*W),N);
     * k = max(k,1);
     * end
     * if length(k) == 1
     * k = [1 k];
     * end
     * 
     * % Generate the diagonals
     * d=((N-1-2*(0:N-1)').^2)*.25*cos(2*pi*W);  % diagonal of B
     * ee=(1:N-1)'.*(N-1:-1:1)'/2;               % super diagonal of B
     * 
     * % Get the eigenvalues of B.
     * v = tridieig(d,[0; ee],N-k(2)+1,N-k(1)+1);
     * v = v(end:-1:1);
     * Lv = length(v);
     * 
     * %B = spdiags([[ee;0] d [0;ee]],[-1 0 1],N,N);
     * %I = speye(N,N);
     * 
     * % Compute the eigenvectors by inverse iteration with
     * % starting vectors of roughly the right shape.
     * E = zeros(N,k(2)-k(1)+1);
     * t = (0:N-1)'/(N-1)*pi;
     * warn_save = warning;  
     * warning('off')      % Turn off warnings in case tridisolve encounters
     * % an exactly singular matrix.
     * for j = 1:Lv
     * e = sin((j+k(1)-1)*t);
     * e = tridisolve(ee,d-v(j),e,N);
     * e = tridisolve(ee,d-v(j),e/norm(e),N);
     * e = tridisolve(ee,d-v(j),e/norm(e),N);
     * E(:,j) = e/norm(e);
     * end
     * warning(warn_save)
     * 
     * d=mean(E);
     * for i=k(1):k(2)
     * if rem(i,2)  % i is odd
     * % Polarize symmetric dpss
     * if d(i-k(1)+1)<0, E(:,i-k(1)+1)=-E(:,i-k(1)+1); end
     * else         % i is even
     * % Polarize anti-symmetric dpss
     * if E(2,i-k(1)+1)<0, E(:,i-k(1)+1)=-E(:,i-k(1)+1); end
     * end
     * end
     * 
     * % get eigenvalues of sinc matrix
     * %  Reference: [1] Percival & Walden, Exercise 8.1, p.390
     * s = [2*W; 4*W*sinc(2*W*(1:N-1)')];
     * q = zeros(size(E));
     * blksz = Lv;  % <-- set this to some small number if OUT OF MEMORY!!!
     * for i=1:blksz:Lv
     * blkind = i:min(i+blksz-1,Lv);
     * q(:,blkind) = fftfilt(E(N:-1:1,blkind),E(:,blkind));
     * end
     * V = q'*flipud(s);
     * 
     * % return 1 for any eigenvalues greater than 1 due to finite precision errors
     * V = min(V,1);
     * % return 0 for any eigenvalues less than 0 due to finite precision errors
     * V = max(V,0);
     * 
     * %---------------------------------------------------------------------
     */
    return E;
}

/*
 * The function "mlfDpss" contains the normal interface for the "dpss"
 * M-function from file "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines
 * 1-89). This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
mxArray * mlfDpss(mxArray * * V, mxArray * N, mxArray * NW, ...) {
    mxArray * varargin = mclUnassigned();
    int nargout = 1;
    mxArray * E = mclGetUninitializedArray();
    mxArray * V__ = mclGetUninitializedArray();
    mlfVarargin(&varargin, NW, 0);
    mlfEnterNewContext(1, -3, V, N, NW, varargin);
    if (V != NULL) {
        ++nargout;
    }
    E = Mdpss(&V__, nargout, N, NW, varargin);
    mlfRestorePreviousContext(1, 2, V, N, NW);
    mxDestroyArray(varargin);
    if (V != NULL) {
        mclCopyOutputArg(V, V__);
    } else {
        mxDestroyArray(V__);
    }
    return mlfReturnValue(E);
}

/*
 * The function "mlxDpss" contains the feval interface for the "dpss"
 * M-function from file "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines
 * 1-89). The feval function calls the implementation version of dpss through
 * this function. This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
void mlxDpss(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: dpss Line: 1 Column: 0 The function \"dpss\""
            " was called with more than the declared number of outputs (2)"));
    }
    for (i = 0; i < 2; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mprhs[2] = NULL;
    mlfAssign(&mprhs[2], mclCreateVararginCell(nrhs - 2, prhs + 2));
    mplhs[0] = Mdpss(&mplhs[1], nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 2 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 2; ++i) {
        mxDestroyArray(mplhs[i]);
    }
    mxDestroyArray(mprhs[2]);
}

/*
 * The function "Mdpss_dpsscalc" is the implementation version of the
 * "dpss/dpsscalc" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 89-164). It contains
 * the actual compiled code for that M-function. It is a static function and
 * must only be called from one of the interface functions, appearing below.
 */
/*
 * function [E,V] = dpsscalc(N,NW,k)
 */
static mxArray * Mdpss_dpsscalc(mxArray * * V,
                                int nargout_,
                                mxArray * N,
                                mxArray * NW,
                                mxArray * k) {
    mxArray * E = mclGetUninitializedArray();
    mxArray * Lv = mclGetUninitializedArray();
    mxArray * W = mclGetUninitializedArray();
    mxArray * ans = mclInitializeAns();
    mxArray * blkind = mclGetUninitializedArray();
    mxArray * blksz = mclGetUninitializedArray();
    mxArray * d = mclGetUninitializedArray();
    mxArray * e = mclGetUninitializedArray();
    mxArray * ee = mclGetUninitializedArray();
    mxArray * i = mclGetUninitializedArray();
    mclForLoopIterator iterator_0;
    mxArray * j = mclGetUninitializedArray();
    mxArray * nargin_ = mclGetUninitializedArray();
    mxArray * q = mclGetUninitializedArray();
    mxArray * s = mclGetUninitializedArray();
    mxArray * t = mclGetUninitializedArray();
    mxArray * v = mclGetUninitializedArray();
    mxArray * warn_save = mclGetUninitializedArray();
    mlfAssign(&nargin_, mlfNargin(0, N, NW, k, NULL));
    mclValidateInputs("dpss/dpsscalc", 3, &N, &NW, &k);
    mclCopyArray(&k);
    /*
     * %DPSSCALC Calculate slepian sequences.
     * %   [E,V] = dpsscalc(N,NW,k) uses tridieig() to get eigenvalues 1:k if k is 
     * %   a scalar, and k(1):k(2) if k is a matrix, of the sparse tridiagonal matrix.
     * %   It then uses inverse iteration using the exact eigenvalues on a starting 
     * %   vector with approximate shape, to get the eigenvectors required.  It then 
     * %   computes the eigenvalues V of the Toeplitz sinc matrix using a fast 
     * %   autocorrelation technique.
     * 
     * %   Authors: T. Krauss, C. Moler, E. Breitenberger
     * W=NW/N;
     */
    mlfAssign(&W, mlfMrdivide(NW, N));
    /*
     * if nargin < 3
     */
    if (mlfTobool(mlfLt(nargin_, mlfScalar(3.0)))) {
        /*
         * k = min(round(2*N*W),N);
         */
        mlfAssign(
          &k,
          mlfMin(
            NULL,
            mlfRound(mlfMtimes(mlfMtimes(mlfScalar(2.0), N), W)),
            N,
            NULL));
        /*
         * k = max(k,1);
         */
        mlfAssign(&k, mlfMax(NULL, k, mlfScalar(1.0), NULL));
    /*
     * end
     */
    }
    /*
     * if length(k) == 1
     */
    if (mlfTobool(mlfEq(mlfLength(k), mlfScalar(1.0)))) {
        /*
         * k = [1 k];
         */
        mlfAssign(&k, mlfHorzcat(mlfScalar(1.0), k, NULL));
    /*
     * end
     */
    }
    /*
     * 
     * % Generate the diagonals
     * d=((N-1-2*(0:N-1)').^2)*.25*cos(2*pi*W);  % diagonal of B
     */
    mlfAssign(
      &d,
      mlfMtimes(
        mlfMtimes(
          mlfPower(
            mlfMinus(
              mlfMinus(N, mlfScalar(1.0)),
              mlfMtimes(
                mlfScalar(2.0),
                mlfCtranspose(
                  mlfColon(
                    mlfScalar(0.0), mlfMinus(N, mlfScalar(1.0)), NULL)))),
            mlfScalar(2.0)),
          mlfScalar(.25)),
        mlfCos(mlfMtimes(mlfMtimes(mlfScalar(2.0), mlfPi()), W))));
    /*
     * ee=(1:N-1)'.*(N-1:-1:1)'/2;               % super diagonal of B
     */
    mlfAssign(
      &ee,
      mlfMrdivide(
        mlfTimes(
          mlfCtranspose(
            mlfColon(mlfScalar(1.0), mlfMinus(N, mlfScalar(1.0)), NULL)),
          mlfCtranspose(
            mlfColon(
              mlfMinus(N, mlfScalar(1.0)), mlfScalar(-1.0), mlfScalar(1.0)))),
        mlfScalar(2.0)));
    /*
     * 
     * % Get the eigenvalues of B.
     * v = tridieig(d,[0; ee],N-k(2)+1,N-k(1)+1);
     */
    mlfAssign(
      &v,
      mlfSignal_private_tridieig(
        d,
        mlfVertcat(
          mlfDoubleMatrix(1, 1, __Array0_r, NULL), mlfHorzcat(ee, NULL), NULL),
        mlfPlus(
          mlfMinus(N, mlfIndexRef(k, "(?)", mlfScalar(2.0))), mlfScalar(1.0)),
        mlfPlus(
          mlfMinus(N, mlfIndexRef(k, "(?)", mlfScalar(1.0))), mlfScalar(1.0)),
        NULL));
    /*
     * v = v(end:-1:1);
     */
    mlfAssign(
      &v,
      mlfIndexRef(
        v,
        "(?)",
        mlfColon(
          mlfEnd(v, mlfScalar(1), mlfScalar(1)),
          mlfScalar(-1.0),
          mlfScalar(1.0))));
    /*
     * Lv = length(v);
     */
    mlfAssign(&Lv, mlfLength(v));
    /*
     * 
     * %B = spdiags([[ee;0] d [0;ee]],[-1 0 1],N,N);
     * %I = speye(N,N);
     * 
     * % Compute the eigenvectors by inverse iteration with
     * % starting vectors of roughly the right shape.
     * E = zeros(N,k(2)-k(1)+1);
     */
    mlfAssign(
      &E,
      mlfZeros(
        N,
        mlfPlus(
          mlfMinus(
            mlfIndexRef(k, "(?)", mlfScalar(2.0)),
            mlfIndexRef(k, "(?)", mlfScalar(1.0))),
          mlfScalar(1.0)),
        NULL));
    /*
     * t = (0:N-1)'/(N-1)*pi;
     */
    mlfAssign(
      &t,
      mlfMtimes(
        mlfMrdivide(
          mlfCtranspose(
            mlfColon(mlfScalar(0.0), mlfMinus(N, mlfScalar(1.0)), NULL)),
          mlfMinus(N, mlfScalar(1.0))),
        mlfPi()));
    /*
     * warn_save = warning;  
     */
    mlfAssign(&warn_save, mlfWarning(NULL, NULL));
    /*
     * warning('off')      % Turn off warnings in case tridisolve encounters
     */
    mclPrintAns(&ans, mlfWarning(NULL, mxCreateString("off")));
    /*
     * % an exactly singular matrix.
     * for j = 1:Lv
     */
    for (mclForStart(&iterator_0, mlfScalar(1.0), Lv, NULL);
         mclForNext(&iterator_0, &j);
         ) {
        /*
         * e = sin((j+k(1)-1)*t);
         */
        mlfAssign(
          &e,
          mlfSin(
            mlfMtimes(
              mlfMinus(
                mlfPlus(j, mlfIndexRef(k, "(?)", mlfScalar(1.0))),
                mlfScalar(1.0)),
              t)));
        /*
         * e = tridisolve(ee,d-v(j),e,N);
         */
        mlfAssign(
          &e,
          mlfNSignal_private_tridisolve(
            0,
            mclValueVarargout(),
            ee,
            mlfMinus(d, mlfIndexRef(v, "(?)", j)),
            e,
            N,
            NULL));
        /*
         * e = tridisolve(ee,d-v(j),e/norm(e),N);
         */
        mlfAssign(
          &e,
          mlfNSignal_private_tridisolve(
            0,
            mclValueVarargout(),
            ee,
            mlfMinus(d, mlfIndexRef(v, "(?)", j)),
            mlfMrdivide(e, mlfNorm(e, NULL)),
            N,
            NULL));
        /*
         * e = tridisolve(ee,d-v(j),e/norm(e),N);
         */
        mlfAssign(
          &e,
          mlfNSignal_private_tridisolve(
            0,
            mclValueVarargout(),
            ee,
            mlfMinus(d, mlfIndexRef(v, "(?)", j)),
            mlfMrdivide(e, mlfNorm(e, NULL)),
            N,
            NULL));
        /*
         * E(:,j) = e/norm(e);
         */
        mlfIndexAssign(
          &E,
          "(?,?)",
          mlfCreateColonIndex(),
          j,
          mlfMrdivide(e, mlfNorm(e, NULL)));
    /*
     * end
     */
    }
    /*
     * warning(warn_save)
     */
    mclPrintAns(&ans, mlfWarning(NULL, warn_save));
    /*
     * 
     * d=mean(E);
     */
    mlfAssign(&d, mlfMean(E, NULL));
    /*
     * for i=k(1):k(2)
     */
    for (mclForStart(
           &iterator_0,
           mlfIndexRef(k, "(?)", mlfScalar(1.0)),
           mlfIndexRef(k, "(?)", mlfScalar(2.0)),
           NULL);
         mclForNext(&iterator_0, &i);
         ) {
        /*
         * if rem(i,2)  % i is odd
         */
        if (mlfTobool(mlfRem(i, mlfScalar(2.0)))) {
            /*
             * % Polarize symmetric dpss
             * if d(i-k(1)+1)<0, E(:,i-k(1)+1)=-E(:,i-k(1)+1); end
             */
            if (mlfTobool(
                  mlfLt(
                    mlfIndexRef(
                      d,
                      "(?)",
                      mlfPlus(
                        mlfMinus(i, mlfIndexRef(k, "(?)", mlfScalar(1.0))),
                        mlfScalar(1.0))),
                    mlfScalar(0.0)))) {
                mlfIndexAssign(
                  &E,
                  "(?,?)",
                  mlfCreateColonIndex(),
                  mlfPlus(
                    mlfMinus(i, mlfIndexRef(k, "(?)", mlfScalar(1.0))),
                    mlfScalar(1.0)),
                  mlfUminus(
                    mlfIndexRef(
                      E,
                      "(?,?)",
                      mlfCreateColonIndex(),
                      mlfPlus(
                        mlfMinus(i, mlfIndexRef(k, "(?)", mlfScalar(1.0))),
                        mlfScalar(1.0)))));
            }
        /*
         * else         % i is even
         */
        } else {
            /*
             * % Polarize anti-symmetric dpss
             * if E(2,i-k(1)+1)<0, E(:,i-k(1)+1)=-E(:,i-k(1)+1); end
             */
            if (mlfTobool(
                  mlfLt(
                    mlfIndexRef(
                      E,
                      "(?,?)",
                      mlfScalar(2.0),
                      mlfPlus(
                        mlfMinus(i, mlfIndexRef(k, "(?)", mlfScalar(1.0))),
                        mlfScalar(1.0))),
                    mlfScalar(0.0)))) {
                mlfIndexAssign(
                  &E,
                  "(?,?)",
                  mlfCreateColonIndex(),
                  mlfPlus(
                    mlfMinus(i, mlfIndexRef(k, "(?)", mlfScalar(1.0))),
                    mlfScalar(1.0)),
                  mlfUminus(
                    mlfIndexRef(
                      E,
                      "(?,?)",
                      mlfCreateColonIndex(),
                      mlfPlus(
                        mlfMinus(i, mlfIndexRef(k, "(?)", mlfScalar(1.0))),
                        mlfScalar(1.0)))));
            }
        /*
         * end
         */
        }
    /*
     * end
     */
    }
    /*
     * 
     * % get eigenvalues of sinc matrix
     * %  Reference: [1] Percival & Walden, Exercise 8.1, p.390
     * s = [2*W; 4*W*sinc(2*W*(1:N-1)')];
     */
    mlfAssign(
      &s,
      mlfVertcat(
        mlfHorzcat(mlfMtimes(mlfScalar(2.0), W), NULL),
        mlfHorzcat(
          mlfMtimes(
            mlfMtimes(mlfScalar(4.0), W),
            mlfSinc(
              mlfMtimes(
                mlfMtimes(mlfScalar(2.0), W),
                mlfCtranspose(
                  mlfColon(
                    mlfScalar(1.0), mlfMinus(N, mlfScalar(1.0)), NULL))))),
          NULL),
        NULL));
    /*
     * q = zeros(size(E));
     */
    mlfAssign(&q, mlfZeros(mlfSize(mclValueVarargout(), E, NULL), NULL));
    /*
     * blksz = Lv;  % <-- set this to some small number if OUT OF MEMORY!!!
     */
    mlfAssign(&blksz, Lv);
    /*
     * for i=1:blksz:Lv
     */
    for (mclForStart(&iterator_0, mlfScalar(1.0), blksz, Lv);
         mclForNext(&iterator_0, &i);
         ) {
        /*
         * blkind = i:min(i+blksz-1,Lv);
         */
        mlfAssign(
          &blkind,
          mlfColon(
            i,
            mlfMin(NULL, mlfMinus(mlfPlus(i, blksz), mlfScalar(1.0)), Lv, NULL),
            NULL));
        /*
         * q(:,blkind) = fftfilt(E(N:-1:1,blkind),E(:,blkind));
         */
        mlfIndexAssign(
          &q,
          "(?,?)",
          mlfCreateColonIndex(),
          blkind,
          mlfFftfilt(
            mlfIndexRef(
              E, "(?,?)", mlfColon(N, mlfScalar(-1.0), mlfScalar(1.0)), blkind),
            mlfIndexRef(E, "(?,?)", mlfCreateColonIndex(), blkind),
            NULL));
    /*
     * end
     */
    }
    /*
     * V = q'*flipud(s);
     */
    mlfAssign(V, mlfMtimes(mlfCtranspose(q), mlfFlipud(s)));
    /*
     * 
     * % return 1 for any eigenvalues greater than 1 due to finite precision errors
     * V = min(V,1);
     */
    mlfAssign(V, mlfMin(NULL, *V, mlfScalar(1.0), NULL));
    /*
     * % return 0 for any eigenvalues less than 0 due to finite precision errors
     * V = max(V,0);
     */
    mlfAssign(V, mlfMax(NULL, *V, mlfScalar(0.0), NULL));
    mclValidateOutputs("dpss/dpsscalc", 2, nargout_, &E, V);
    mxDestroyArray(Lv);
    mxDestroyArray(W);
    mxDestroyArray(ans);
    mxDestroyArray(blkind);
    mxDestroyArray(blksz);
    mxDestroyArray(d);
    mxDestroyArray(e);
    mxDestroyArray(ee);
    mxDestroyArray(i);
    mxDestroyArray(j);
    mxDestroyArray(k);
    mxDestroyArray(nargin_);
    mxDestroyArray(q);
    mxDestroyArray(s);
    mxDestroyArray(t);
    mxDestroyArray(v);
    mxDestroyArray(warn_save);
    /*
     * 
     * %---------------------------------------------------------------------
     * function [En,V] = dpssint(N, NW, M, int, E,V)
     * % Syntax: [En,V]=dpssint(N,NW); [En,V]=dpssint(N,NW,M,'spline');
     * %  Dpssint calculates discrete prolate spheroidal
     * %  sequences for the parameters N and NW. Note that
     * %  NW is normally 2, 5/2, 3, 7/2, or 4 - not i/N. The 
     * %  dpss are interpolated from previously calculated 
     * %  dpss of order M (128, 256, 512, or 1024). 256 is the 
     * %  default for M. The interpolation can be 'linear' 
     * %  or 'spline'. 'Linear' is faster, 'spline' the default.
     * %  Linear interpolation can only be used for M>N. 
     * %  Returns:
     * %              E: matrix of dpss (N by 2NW)
     * %              V: eigenvalue vector (2NW)
     * % 
     * % Errors in the interpolated dpss are very small but should be 
     * % checked if possible. The differences between interpolated
     * % values and values from dpsscalc are generally of order
     * % 10ee-5 or better. Spline interpolation is generally more
     * % accurate. Fractional errors can be very large near
     * % the zero-crossings but this does not seriously affect
     * % windowing calculations. The error can be reduced by using
     * % values for M which are close to N.
     * %
     * % Written by Eric Breitenberger, version date 10/3/95.
     * % Please send comments and suggestions to eric@gi.alaska.edu
     * %
     * 
     * W = NW/N;
     * 
     * if     nargin==2,
     * M=256; int='spline';
     * elseif nargin==3,
     * if isstr(M), int=M; M=256; 
     * else, int='spline';, end
     * end
     * 
     * if int=='linear' & N>M
     * error('Linear interpolation cannot be used for N>Ni. Use splining instead.')
     * end
     * 
     * if nargin<=4
     * [E,V] = dpssload(M,NW);
     * else
     * if size(E,1)~=M
     * error('Ni and row size of E don''t match.')
     * end
     * end
     * 
     * k=min(round(2*NW),N); % Return only first k values
     * k = max(k,1);
     * E=E(:,1:k);
     * V=V(1:k);
     * x=1:M;
     * 
     * % The scaling for the interpolation:
     * % This is not necessarily optimal, and 
     * % changing s can improve accuracy.
     * 
     * s=M/N;
     * midm=(M+1)/2;
     * midn=(N+1)/2;
     * delta=midm-s*midn;
     * xi=linspace(1-delta, M+delta, N);
     * 
     * % Interpolate from M values to N
     * % Spline interpolation is a bit better,
     * % but takes about twice as long.
     * % Errors from linear interpolation are 
     * % usually smaller than errors from scaling.
     * 
     * En=interp1(x,E,xi,['*' int]);
     * 
     * % Re-normalize the eigenvectors
     * En=En./(ones(N,1)*sqrt(sum(En.*En)));
     * 
     * %----------------------------------------------------------------------
     */
    return E;
}

/*
 * The function "mlfDpss_dpsscalc" contains the normal interface for the
 * "dpss/dpsscalc" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 89-164). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
static mxArray * mlfDpss_dpsscalc(mxArray * * V,
                                  mxArray * N,
                                  mxArray * NW,
                                  mxArray * k) {
    int nargout = 1;
    mxArray * E = mclGetUninitializedArray();
    mxArray * V__ = mclGetUninitializedArray();
    mlfEnterNewContext(1, 3, V, N, NW, k);
    if (V != NULL) {
        ++nargout;
    }
    E = Mdpss_dpsscalc(&V__, nargout, N, NW, k);
    mlfRestorePreviousContext(1, 3, V, N, NW, k);
    if (V != NULL) {
        mclCopyOutputArg(V, V__);
    } else {
        mxDestroyArray(V__);
    }
    return mlfReturnValue(E);
}

/*
 * The function "mlxDpss_dpsscalc" contains the feval interface for the
 * "dpss/dpsscalc" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 89-164). The feval
 * function calls the implementation version of dpss/dpsscalc through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
static void mlxDpss_dpsscalc(int nlhs,
                             mxArray * plhs[],
                             int nrhs,
                             mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: dpss/dpsscalc Line: 89 Colum"
            "n: 0 The function \"dpss/dpsscalc\" was called wit"
            "h more than the declared number of outputs (2)"));
    }
    if (nrhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: dpss/dpsscalc Line: 89 Colum"
            "n: 0 The function \"dpss/dpsscalc\" was called wit"
            "h more than the declared number of inputs (3)"));
    }
    for (i = 0; i < 2; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 3 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 3; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    mplhs[0] = Mdpss_dpsscalc(&mplhs[1], nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 2 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 2; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "Mdpss_dpssint" is the implementation version of the
 * "dpss/dpssint" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 164-240). It contains
 * the actual compiled code for that M-function. It is a static function and
 * must only be called from one of the interface functions, appearing below.
 */
/*
 * function [En,V] = dpssint(N, NW, M, int, E,V)
 */
static mxArray * Mdpss_dpssint(mxArray * * V,
                               int nargout_,
                               mxArray * N,
                               mxArray * NW,
                               mxArray * M,
                               mxArray * int0,
                               mxArray * E,
                               mxArray * V_) {
    mxArray * En = mclGetUninitializedArray();
    mxArray * W = mclGetUninitializedArray();
    mxArray * delta = mclGetUninitializedArray();
    mxArray * k = mclGetUninitializedArray();
    mxArray * midm = mclGetUninitializedArray();
    mxArray * midn = mclGetUninitializedArray();
    mxArray * nargin_ = mclGetUninitializedArray();
    mxArray * s = mclGetUninitializedArray();
    mxArray * x = mclGetUninitializedArray();
    mxArray * xi = mclGetUninitializedArray();
    mlfAssign(&nargin_, mlfNargin(0, N, NW, M, int0, E, V_, NULL));
    mclValidateInputs("dpss/dpssint", 6, &N, &NW, &M, &int0, &E, &V_);
    mclCopyArray(&M);
    mclCopyArray(&int0);
    mclCopyArray(&E);
    mclCopyInputArg(V, V_);
    /*
     * % Syntax: [En,V]=dpssint(N,NW); [En,V]=dpssint(N,NW,M,'spline');
     * %  Dpssint calculates discrete prolate spheroidal
     * %  sequences for the parameters N and NW. Note that
     * %  NW is normally 2, 5/2, 3, 7/2, or 4 - not i/N. The 
     * %  dpss are interpolated from previously calculated 
     * %  dpss of order M (128, 256, 512, or 1024). 256 is the 
     * %  default for M. The interpolation can be 'linear' 
     * %  or 'spline'. 'Linear' is faster, 'spline' the default.
     * %  Linear interpolation can only be used for M>N. 
     * %  Returns:
     * %              E: matrix of dpss (N by 2NW)
     * %              V: eigenvalue vector (2NW)
     * % 
     * % Errors in the interpolated dpss are very small but should be 
     * % checked if possible. The differences between interpolated
     * % values and values from dpsscalc are generally of order
     * % 10ee-5 or better. Spline interpolation is generally more
     * % accurate. Fractional errors can be very large near
     * % the zero-crossings but this does not seriously affect
     * % windowing calculations. The error can be reduced by using
     * % values for M which are close to N.
     * %
     * % Written by Eric Breitenberger, version date 10/3/95.
     * % Please send comments and suggestions to eric@gi.alaska.edu
     * %
     * 
     * W = NW/N;
     */
    mlfAssign(&W, mlfMrdivide(NW, N));
    /*
     * 
     * if     nargin==2,
     */
    if (mlfTobool(mlfEq(nargin_, mlfScalar(2.0)))) {
        /*
         * M=256; int='spline';
         */
        mlfAssign(&M, mlfScalar(256.0));
        mlfAssign(&int0, mxCreateString("spline"));
    /*
     * elseif nargin==3,
     */
    } else if (mlfTobool(mlfEq(nargin_, mlfScalar(3.0)))) {
        /*
         * if isstr(M), int=M; M=256; 
         */
        if (mlfTobool(mlfIsstr(M))) {
            mlfAssign(&int0, M);
            mlfAssign(&M, mlfScalar(256.0));
        /*
         * else, int='spline';, end
         */
        } else {
            mlfAssign(&int0, mxCreateString("spline"));
        }
    /*
     * end
     */
    }
    /*
     * 
     * if int=='linear' & N>M
     */
    {
        mxArray * a_ = mclInitialize(mlfEq(int0, mxCreateString("linear")));
        if (mlfTobool(a_) && mlfTobool(mlfAnd(a_, mlfGt(N, M)))) {
            mxDestroyArray(a_);
            /*
             * error('Linear interpolation cannot be used for N>Ni. Use splining instead.')
             */
            mlfError(
              mxCreateString(
                "Linear interpolation cannot be used"
                " for N>Ni. Use splining instead."));
        } else {
            mxDestroyArray(a_);
        }
    /*
     * end
     */
    }
    /*
     * 
     * if nargin<=4
     */
    if (mlfTobool(mlfLe(nargin_, mlfScalar(4.0)))) {
        /*
         * [E,V] = dpssload(M,NW);
         */
        mlfAssign(&E, mlfDpssload(V, M, NW));
    /*
     * else
     */
    } else {
        /*
         * if size(E,1)~=M
         */
        if (mlfTobool(
              mlfNe(mlfSize(mclValueVarargout(), E, mlfScalar(1.0)), M))) {
            /*
             * error('Ni and row size of E don''t match.')
             */
            mlfError(mxCreateString("Ni and row size of E don't match."));
        /*
         * end
         */
        }
    /*
     * end
     */
    }
    /*
     * 
     * k=min(round(2*NW),N); % Return only first k values
     */
    mlfAssign(
      &k, mlfMin(NULL, mlfRound(mlfMtimes(mlfScalar(2.0), NW)), N, NULL));
    /*
     * k = max(k,1);
     */
    mlfAssign(&k, mlfMax(NULL, k, mlfScalar(1.0), NULL));
    /*
     * E=E(:,1:k);
     */
    mlfAssign(
      &E,
      mlfIndexRef(
        E, "(?,?)", mlfCreateColonIndex(), mlfColon(mlfScalar(1.0), k, NULL)));
    /*
     * V=V(1:k);
     */
    mlfAssign(V, mlfIndexRef(*V, "(?)", mlfColon(mlfScalar(1.0), k, NULL)));
    /*
     * x=1:M;
     */
    mlfAssign(&x, mlfColon(mlfScalar(1.0), M, NULL));
    /*
     * 
     * % The scaling for the interpolation:
     * % This is not necessarily optimal, and 
     * % changing s can improve accuracy.
     * 
     * s=M/N;
     */
    mlfAssign(&s, mlfMrdivide(M, N));
    /*
     * midm=(M+1)/2;
     */
    mlfAssign(&midm, mlfMrdivide(mlfPlus(M, mlfScalar(1.0)), mlfScalar(2.0)));
    /*
     * midn=(N+1)/2;
     */
    mlfAssign(&midn, mlfMrdivide(mlfPlus(N, mlfScalar(1.0)), mlfScalar(2.0)));
    /*
     * delta=midm-s*midn;
     */
    mlfAssign(&delta, mlfMinus(midm, mlfMtimes(s, midn)));
    /*
     * xi=linspace(1-delta, M+delta, N);
     */
    mlfAssign(
      &xi, mlfLinspace(mlfMinus(mlfScalar(1.0), delta), mlfPlus(M, delta), N));
    /*
     * 
     * % Interpolate from M values to N
     * % Spline interpolation is a bit better,
     * % but takes about twice as long.
     * % Errors from linear interpolation are 
     * % usually smaller than errors from scaling.
     * 
     * En=interp1(x,E,xi,['*' int]);
     */
    mlfAssign(
      &En,
      mlfInterp1(x, E, xi, mlfHorzcat(mxCreateString("*"), int0, NULL), NULL));
    /*
     * 
     * % Re-normalize the eigenvectors
     * En=En./(ones(N,1)*sqrt(sum(En.*En)));
     */
    mlfAssign(
      &En,
      mlfRdivide(
        En,
        mlfMtimes(
          mlfOnes(N, mlfScalar(1.0), NULL),
          mlfSqrt(mlfSum(mlfTimes(En, En), NULL)))));
    mclValidateOutputs("dpss/dpssint", 2, nargout_, &En, V);
    mxDestroyArray(E);
    mxDestroyArray(M);
    mxDestroyArray(W);
    mxDestroyArray(delta);
    mxDestroyArray(int0);
    mxDestroyArray(k);
    mxDestroyArray(midm);
    mxDestroyArray(midn);
    mxDestroyArray(nargin_);
    mxDestroyArray(s);
    mxDestroyArray(x);
    mxDestroyArray(xi);
    /*
     * 
     * %----------------------------------------------------------------------
     * function [method,k,Ni,traceFlag,err_msg] = parseinputs(N,NW,varargin);
     * % PARSEINPUTS Parses the inputs to the DPSS function and validates it.
     * %
     * % Inputs:   -  The inputs to this function are the same as the ones 
     * %              passed to DPSS.  See the help for DPSS.
     * %
     * % Outputs:
     * %   method    - method updated with value entered by the user
     * %   k         - number of most band-limited DPSSs
     * %   Ni        - length of DPSSs to be interpolated
     * %   traceFlag - a boolean flag indicating if 'trace' was resquested
     * %   err_msg   - error message if an error occurred
     * %
     * 
     * % Here are all possible input combinations in varargin (after N,NW)...
     * %
     * % 1 Option specified
     * % (N,NW,k), (N,NW,method) or (N,NW,'trace')
     * %
     * % 2 Options Specified.
     * % (N,NW,k,method) or (N,NW,k,'trace') or  (N,NW,method,'trace')
     * % or (N,NW,method,Ni)
     * %
     * % 3 Options Specified.
     * % (N,NW,k,method,Ni), (N,NW,k,method,'trace') or (N,NW,'method',Ni,'trace')
     * %
     * % 4 Options Specified.
     * % (N,NW,k,method,Ni,'trace')
     * 
     * % Defined defaults
     * % If user didn't specify method or traceFlag set them to defaults.
     * method    = 'calc';
     * k         = [];
     * Ni        = [];
     * traceFlag = [];  % It gets set to 0 if user does not specified it.
     * err_msg   = '';
     * 
     * % Validate input arguments N and NW.
     * if length(N)>1,
     * err_msg = 'The length of the sequence, N, must be a scalar.';
     * return;
     * end
     * if length(NW)>1, 
     * err_msg = 'The Time-bandwidth product NW must be a scalar.';
     * return;
     * end
     * 
     * if NW >= N/2,
     * err_msg = 'Time-bandwidth product NW must be less than N/2.';
     * return;
     * end
     * 
     * % Default number of sequences to return
     * k = min(round(2*NW),N);
     * k = max(k,1);
     * 
     * % Validate and parse the optional input arguments
     * nopts = length(varargin);
     * 
     * for opt_indx = 1:nopts,
     * arg = varargin{opt_indx};
     * 
     * % Parse strings
     * if isstr(arg),
     * [method,traceFlag,err_msg] = ...
     * parseStrOpts(arg,method,traceFlag,nopts,opt_indx);
     * if err_msg, return; end
     * 
     * else % Parse numerics
     * 
     * % 1 Option Specified.
     * if opt_indx == 1,
     * k = arg;
     * if isempty(k) | any(k~=round(abs(k))) | any(k>N),
     * err_msg = 'K must be a positive integer in the range 1:N';
     * elseif length(k)>2 | prod(size(k))>2,
     * err_msg = 'K must be either a scalar or a two element vector.';
     * end
     * if err_msg, return; end
     * 
     * % 2 or 3 Options Specified
     * elseif any(opt_indx==[2,3]),
     * Ni = arg; 
     * if length(Ni)>1 | isempty(Ni),
     * err_msg = 'Ni must be a positive integer.';
     * return;
     * end      
     * end
     * 
     * end 
     * 
     * end % for-loop
     * 
     * if isempty(traceFlag),  % If user didn't specify it set it to 0 (default).
     * traceFlag = 0;
     * end
     * 
     * 
     * %----------------------------------------------------------------------
     */
    return En;
}

/*
 * The function "mlfDpss_dpssint" contains the normal interface for the
 * "dpss/dpssint" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 164-240). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
static mxArray * mlfDpss_dpssint(mxArray * * V,
                                 mxArray * N,
                                 mxArray * NW,
                                 mxArray * M,
                                 mxArray * int0,
                                 mxArray * E,
                                 mxArray * V_) {
    int nargout = 1;
    mxArray * En = mclGetUninitializedArray();
    mxArray * V__ = mclGetUninitializedArray();
    mlfEnterNewContext(1, 6, V, N, NW, M, int0, E, V_);
    if (V != NULL) {
        ++nargout;
    }
    En = Mdpss_dpssint(&V__, nargout, N, NW, M, int0, E, V_);
    mlfRestorePreviousContext(1, 6, V, N, NW, M, int0, E, V_);
    if (V != NULL) {
        mclCopyOutputArg(V, V__);
    } else {
        mxDestroyArray(V__);
    }
    return mlfReturnValue(En);
}

/*
 * The function "mlxDpss_dpssint" contains the feval interface for the
 * "dpss/dpssint" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 164-240). The feval
 * function calls the implementation version of dpss/dpssint through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
static void mlxDpss_dpssint(int nlhs,
                            mxArray * plhs[],
                            int nrhs,
                            mxArray * prhs[]) {
    mxArray * mprhs[6];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: dpss/dpssint Line: 164 Colum"
            "n: 0 The function \"dpss/dpssint\" was called with"
            " more than the declared number of outputs (2)"));
    }
    if (nrhs > 6) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: dpss/dpssint Line: 164 Colu"
            "mn: 0 The function \"dpss/dpssint\" was called wi"
            "th more than the declared number of inputs (6)"));
    }
    for (i = 0; i < 2; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 6 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 6; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(
      0, 6, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4], mprhs[5]);
    mplhs[0]
      = Mdpss_dpssint(
          &mplhs[1],
          nlhs,
          mprhs[0],
          mprhs[1],
          mprhs[2],
          mprhs[3],
          mprhs[4],
          mprhs[5]);
    mlfRestorePreviousContext(
      0, 6, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4], mprhs[5]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 2 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 2; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "Mdpss_parseinputs" is the implementation version of the
 * "dpss/parseinputs" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 240-339). It contains
 * the actual compiled code for that M-function. It is a static function and
 * must only be called from one of the interface functions, appearing below.
 */
/*
 * function [method,k,Ni,traceFlag,err_msg] = parseinputs(N,NW,varargin);
 */
static mxArray * Mdpss_parseinputs(mxArray * * k,
                                   mxArray * * Ni,
                                   mxArray * * traceFlag,
                                   mxArray * * err_msg,
                                   int nargout_,
                                   mxArray * N,
                                   mxArray * NW,
                                   mxArray * varargin) {
    mxArray * method = mclGetUninitializedArray();
    mxArray * arg = mclGetUninitializedArray();
    mclForLoopIterator iterator_0;
    mxArray * nopts = mclGetUninitializedArray();
    mxArray * opt_indx = mclGetUninitializedArray();
    mclValidateInputs("dpss/parseinputs", 2, &N, &NW);
    /*
     * % PARSEINPUTS Parses the inputs to the DPSS function and validates it.
     * %
     * % Inputs:   -  The inputs to this function are the same as the ones 
     * %              passed to DPSS.  See the help for DPSS.
     * %
     * % Outputs:
     * %   method    - method updated with value entered by the user
     * %   k         - number of most band-limited DPSSs
     * %   Ni        - length of DPSSs to be interpolated
     * %   traceFlag - a boolean flag indicating if 'trace' was resquested
     * %   err_msg   - error message if an error occurred
     * %
     * 
     * % Here are all possible input combinations in varargin (after N,NW)...
     * %
     * % 1 Option specified
     * % (N,NW,k), (N,NW,method) or (N,NW,'trace')
     * %
     * % 2 Options Specified.
     * % (N,NW,k,method) or (N,NW,k,'trace') or  (N,NW,method,'trace')
     * % or (N,NW,method,Ni)
     * %
     * % 3 Options Specified.
     * % (N,NW,k,method,Ni), (N,NW,k,method,'trace') or (N,NW,'method',Ni,'trace')
     * %
     * % 4 Options Specified.
     * % (N,NW,k,method,Ni,'trace')
     * 
     * % Defined defaults
     * % If user didn't specify method or traceFlag set them to defaults.
     * method    = 'calc';
     */
    mlfAssign(&method, mxCreateString("calc"));
    /*
     * k         = [];
     */
    mlfAssign(k, mclCreateEmptyArray());
    /*
     * Ni        = [];
     */
    mlfAssign(Ni, mclCreateEmptyArray());
    /*
     * traceFlag = [];  % It gets set to 0 if user does not specified it.
     */
    mlfAssign(traceFlag, mclCreateEmptyArray());
    /*
     * err_msg   = '';
     */
    mlfAssign(err_msg, mxCreateString(""));
    /*
     * 
     * % Validate input arguments N and NW.
     * if length(N)>1,
     */
    if (mlfTobool(mlfGt(mlfLength(N), mlfScalar(1.0)))) {
        /*
         * err_msg = 'The length of the sequence, N, must be a scalar.';
         */
        mlfAssign(
          err_msg,
          mxCreateString("The length of the sequence, N, must be a scalar."));
        /*
         * return;
         */
        goto return_;
    /*
     * end
     */
    }
    /*
     * if length(NW)>1, 
     */
    if (mlfTobool(mlfGt(mlfLength(NW), mlfScalar(1.0)))) {
        /*
         * err_msg = 'The Time-bandwidth product NW must be a scalar.';
         */
        mlfAssign(
          err_msg,
          mxCreateString("The Time-bandwidth product NW must be a scalar."));
        /*
         * return;
         */
        goto return_;
    /*
     * end
     */
    }
    /*
     * 
     * if NW >= N/2,
     */
    if (mlfTobool(mlfGe(NW, mlfMrdivide(N, mlfScalar(2.0))))) {
        /*
         * err_msg = 'Time-bandwidth product NW must be less than N/2.';
         */
        mlfAssign(
          err_msg,
          mxCreateString("Time-bandwidth product NW must be less than N/2."));
        /*
         * return;
         */
        goto return_;
    /*
     * end
     */
    }
    /*
     * 
     * % Default number of sequences to return
     * k = min(round(2*NW),N);
     */
    mlfAssign(
      k, mlfMin(NULL, mlfRound(mlfMtimes(mlfScalar(2.0), NW)), N, NULL));
    /*
     * k = max(k,1);
     */
    mlfAssign(k, mlfMax(NULL, *k, mlfScalar(1.0), NULL));
    /*
     * 
     * % Validate and parse the optional input arguments
     * nopts = length(varargin);
     */
    mlfAssign(&nopts, mlfLength(varargin));
    /*
     * 
     * for opt_indx = 1:nopts,
     */
    for (mclForStart(&iterator_0, mlfScalar(1.0), nopts, NULL);
         mclForNext(&iterator_0, &opt_indx);
         ) {
        /*
         * arg = varargin{opt_indx};
         */
        mlfAssign(&arg, mlfIndexRef(varargin, "{?}", opt_indx));
        /*
         * 
         * % Parse strings
         * if isstr(arg),
         */
        if (mlfTobool(mlfIsstr(arg))) {
            /*
             * [method,traceFlag,err_msg] = ...
             * parseStrOpts(arg,method,traceFlag,nopts,opt_indx);
             */
            mlfAssign(
              &method,
              mlfDpss_parseStrOpts(
                traceFlag, err_msg, arg, method, *traceFlag, nopts, opt_indx));
            /*
             * if err_msg, return; end
             */
            if (mlfTobool(*err_msg)) {
                goto return_;
            }
        /*
         * 
         * else % Parse numerics
         */
        } else {
            /*
             * 
             * % 1 Option Specified.
             * if opt_indx == 1,
             */
            if (mlfTobool(mlfEq(opt_indx, mlfScalar(1.0)))) {
                /*
                 * k = arg;
                 */
                mlfAssign(k, arg);
                /*
                 * if isempty(k) | any(k~=round(abs(k))) | any(k>N),
                 */
                {
                    mxArray * a_ = mclInitialize(mlfIsempty(*k));
                    if (mlfTobool(a_)) {
                        mlfAssign(&a_, mlfScalar(1));
                    } else {
                        mlfAssign(
                          &a_,
                          mlfOr(
                            a_, mlfAny(mlfNe(*k, mlfRound(mlfAbs(*k))), NULL)));
                    }
                    if (mlfTobool(a_)
                        || mlfTobool(mlfOr(a_, mlfAny(mlfGt(*k, N), NULL)))) {
                        mxDestroyArray(a_);
                        /*
                         * err_msg = 'K must be a positive integer in the range 1:N';
                         */
                        mlfAssign(
                          err_msg,
                          mxCreateString(
                            "K must be a positive integer in the range 1:N"));
                    /*
                     * elseif length(k)>2 | prod(size(k))>2,
                     */
                    } else {
                        mxDestroyArray(a_);
                        {
                            mxArray * a_ = mclInitialize(
                                             mlfGt(
                                               mlfLength(*k), mlfScalar(2.0)));
                            if (mlfTobool(a_)
                                || mlfTobool(
                                     mlfOr(
                                       a_,
                                       mlfGt(
                                         mlfProd(
                                           mlfSize(
                                             mclValueVarargout(), *k, NULL),
                                           NULL),
                                         mlfScalar(2.0))))) {
                                mxDestroyArray(a_);
                                /*
                                 * err_msg = 'K must be either a scalar or a two element vector.';
                                 */
                                mlfAssign(
                                  err_msg,
                                  mxCreateString(
                                    "K must be either a scalar o"
                                    "r a two element vector."));
                            } else {
                                mxDestroyArray(a_);
                            }
                        }
                    }
                /*
                 * end
                 */
                }
                /*
                 * if err_msg, return; end
                 */
                if (mlfTobool(*err_msg)) {
                    goto return_;
                }
            /*
             * 
             * % 2 or 3 Options Specified
             * elseif any(opt_indx==[2,3]),
             */
            } else if (mlfTobool(
                         mlfAny(
                           mlfEq(
                             opt_indx, mlfDoubleMatrix(1, 2, __Array1_r, NULL)),
                           NULL))) {
                /*
                 * Ni = arg; 
                 */
                mlfAssign(Ni, arg);
                /*
                 * if length(Ni)>1 | isempty(Ni),
                 */
                {
                    mxArray * a_ = mclInitialize(
                                     mlfGt(mlfLength(*Ni), mlfScalar(1.0)));
                    if (mlfTobool(a_)
                        || mlfTobool(mlfOr(a_, mlfIsempty(*Ni)))) {
                        mxDestroyArray(a_);
                        /*
                         * err_msg = 'Ni must be a positive integer.';
                         */
                        mlfAssign(
                          err_msg,
                          mxCreateString("Ni must be a positive integer."));
                        /*
                         * return;
                         */
                        goto return_;
                    } else {
                        mxDestroyArray(a_);
                    }
                /*
                 * end      
                 */
                }
            /*
             * end
             */
            }
        /*
         * 
         * end 
         */
        }
    /*
     * 
     * end % for-loop
     */
    }
    /*
     * 
     * if isempty(traceFlag),  % If user didn't specify it set it to 0 (default).
     */
    if (mlfTobool(mlfIsempty(*traceFlag))) {
        /*
         * traceFlag = 0;
         */
        mlfAssign(traceFlag, mlfScalar(0.0));
    /*
     * end
     */
    }
    /*
     * 
     * 
     * %----------------------------------------------------------------------
     * function [method,traceFlag,err_msg] = ...
     */
    return_:
    mclValidateOutputs(
      "dpss/parseinputs", 5, nargout_, &method, k, Ni, traceFlag, err_msg);
    mxDestroyArray(arg);
    mxDestroyArray(nopts);
    mxDestroyArray(opt_indx);
    /*
     * parseStrOpts(inputStr,method,traceFlag,nopts,opt_indx),
     * % PARSESTROPTS Parses the input string options passed to dpss.
     * %
     * % Inputs:
     * %   inputStr  - input string entered by the user
     * %   method    - interpolation method
     * %   traceFlag - a scalar with values 1 = 'trace' or 0 = 'notrace
     * %   nopts     - number of optional input arguments specified
     * %   opt_indx  - index of the input option being parsed
     * %
     * % Outputs:
     * %   method    - interpolation method specified by the user
     * %   traceFlag - a scalar with values 1 = 'trace' or 0 = 'notrace
     * %   err_msg   - error message if an error occurred
     * 
     * % Define output variable in case of early return.
     * err_msg = '';
     * 
     * % Define all possible string options; lower case, no punctuation:
     * strOpts  = {'calc','spline','linear',...
     * 'trace','notrace'};
     * 
     * i = strmatch(lower(inputStr), strOpts);
     * if isempty(i),
     * str1 = 'Invalid string option.  Valid strings are ''calc'', ''spline'' and';
     * str2 = '''linear'' for method, or ''trace'' to show which method is being used.';
     * err_msg = sprintf('%s\n%s\n',str1,str2);
     * return
     * else
     * % Determine option set to which this string applies:
     * switch i
     * case {1,2,3},
     * method = strOpts{i};
     * case 4,
     * traceFlag = 1;
     * case 5,         
     * traceFlag = 0;
     * end
     * end
     * 
     * if ~isempty(traceFlag) & nopts > opt_indx,
     * err_msg = 'The trailing string ''trace'' must be last in the input parameter list.'; 
     * end
     * 
     * % [EOF] dpss.m
     */
    return method;
}

/*
 * The function "mlfDpss_parseinputs" contains the normal interface for the
 * "dpss/parseinputs" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 240-339). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
static mxArray * mlfDpss_parseinputs(mxArray * * k,
                                     mxArray * * Ni,
                                     mxArray * * traceFlag,
                                     mxArray * * err_msg,
                                     mxArray * N,
                                     mxArray * NW,
                                     ...) {
    mxArray * varargin = mclUnassigned();
    int nargout = 1;
    mxArray * method = mclGetUninitializedArray();
    mxArray * k__ = mclGetUninitializedArray();
    mxArray * Ni__ = mclGetUninitializedArray();
    mxArray * traceFlag__ = mclGetUninitializedArray();
    mxArray * err_msg__ = mclGetUninitializedArray();
    mlfVarargin(&varargin, NW, 0);
    mlfEnterNewContext(4, -3, k, Ni, traceFlag, err_msg, N, NW, varargin);
    if (k != NULL) {
        ++nargout;
    }
    if (Ni != NULL) {
        ++nargout;
    }
    if (traceFlag != NULL) {
        ++nargout;
    }
    if (err_msg != NULL) {
        ++nargout;
    }
    method
      = Mdpss_parseinputs(
          &k__, &Ni__, &traceFlag__, &err_msg__, nargout, N, NW, varargin);
    mlfRestorePreviousContext(4, 2, k, Ni, traceFlag, err_msg, N, NW);
    mxDestroyArray(varargin);
    if (k != NULL) {
        mclCopyOutputArg(k, k__);
    } else {
        mxDestroyArray(k__);
    }
    if (Ni != NULL) {
        mclCopyOutputArg(Ni, Ni__);
    } else {
        mxDestroyArray(Ni__);
    }
    if (traceFlag != NULL) {
        mclCopyOutputArg(traceFlag, traceFlag__);
    } else {
        mxDestroyArray(traceFlag__);
    }
    if (err_msg != NULL) {
        mclCopyOutputArg(err_msg, err_msg__);
    } else {
        mxDestroyArray(err_msg__);
    }
    return mlfReturnValue(method);
}

/*
 * The function "mlxDpss_parseinputs" contains the feval interface for the
 * "dpss/parseinputs" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 240-339). The feval
 * function calls the implementation version of dpss/parseinputs through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
static void mlxDpss_parseinputs(int nlhs,
                                mxArray * plhs[],
                                int nrhs,
                                mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[5];
    int i;
    if (nlhs > 5) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: dpss/parseinputs Line: 240 Col"
            "umn: 0 The function \"dpss/parseinputs\" was called "
            "with more than the declared number of outputs (5)"));
    }
    for (i = 0; i < 5; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mprhs[2] = NULL;
    mlfAssign(&mprhs[2], mclCreateVararginCell(nrhs - 2, prhs + 2));
    mplhs[0]
      = Mdpss_parseinputs(
          &mplhs[1],
          &mplhs[2],
          &mplhs[3],
          &mplhs[4],
          nlhs,
          mprhs[0],
          mprhs[1],
          mprhs[2]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 5 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 5; ++i) {
        mxDestroyArray(mplhs[i]);
    }
    mxDestroyArray(mprhs[2]);
}

/*
 * The function "Mdpss_parseStrOpts" is the implementation version of the
 * "dpss/parseStrOpts" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 339-385). It contains
 * the actual compiled code for that M-function. It is a static function and
 * must only be called from one of the interface functions, appearing below.
 */
/*
 * function [method,traceFlag,err_msg] = ...
 */
static mxArray * Mdpss_parseStrOpts(mxArray * * traceFlag,
                                    mxArray * * err_msg,
                                    int nargout_,
                                    mxArray * inputStr,
                                    mxArray * method_,
                                    mxArray * traceFlag_,
                                    mxArray * nopts,
                                    mxArray * opt_indx) {
    mxArray * method = mclGetUninitializedArray();
    mxArray * i = mclGetUninitializedArray();
    mxArray * str1 = mclGetUninitializedArray();
    mxArray * str2 = mclGetUninitializedArray();
    mxArray * strOpts = mclGetUninitializedArray();
    mclValidateInputs(
      "dpss/parseStrOpts",
      5,
      &inputStr,
      &method_,
      &traceFlag_,
      &nopts,
      &opt_indx);
    mclCopyInputArg(&method, method_);
    mclCopyInputArg(traceFlag, traceFlag_);
    /*
     * parseStrOpts(inputStr,method,traceFlag,nopts,opt_indx),
     * % PARSESTROPTS Parses the input string options passed to dpss.
     * %
     * % Inputs:
     * %   inputStr  - input string entered by the user
     * %   method    - interpolation method
     * %   traceFlag - a scalar with values 1 = 'trace' or 0 = 'notrace
     * %   nopts     - number of optional input arguments specified
     * %   opt_indx  - index of the input option being parsed
     * %
     * % Outputs:
     * %   method    - interpolation method specified by the user
     * %   traceFlag - a scalar with values 1 = 'trace' or 0 = 'notrace
     * %   err_msg   - error message if an error occurred
     * 
     * % Define output variable in case of early return.
     * err_msg = '';
     */
    mlfAssign(err_msg, mxCreateString(""));
    /*
     * 
     * % Define all possible string options; lower case, no punctuation:
     * strOpts  = {'calc','spline','linear',...
     */
    mlfAssign(
      &strOpts,
      mlfCellhcat(
        mxCreateString("calc"),
        mxCreateString("spline"),
        mxCreateString("linear"),
        mxCreateString("trace"),
        mxCreateString("notrace"),
        NULL));
    /*
     * 'trace','notrace'};
     * 
     * i = strmatch(lower(inputStr), strOpts);
     */
    mlfAssign(&i, mlfStrmatch(mlfLower(inputStr), strOpts, NULL));
    /*
     * if isempty(i),
     */
    if (mlfTobool(mlfIsempty(i))) {
        /*
         * str1 = 'Invalid string option.  Valid strings are ''calc'', ''spline'' and';
         */
        mlfAssign(
          &str1,
          mxCreateString(
            "Invalid string option.  Valid strings are 'calc', 'spline' and"));
        /*
         * str2 = '''linear'' for method, or ''trace'' to show which method is being used.';
         */
        mlfAssign(
          &str2,
          mxCreateString(
            "'linear' for method, or 'trace' to "
            "show which method is being used."));
        /*
         * err_msg = sprintf('%s\n%s\n',str1,str2);
         */
        mlfAssign(
          err_msg,
          mlfSprintf(NULL, mxCreateString("%s\\n%s\\n"), str1, str2, NULL));
        /*
         * return
         */
        goto return_;
    /*
     * else
     */
    } else {
        /*
         * % Determine option set to which this string applies:
         * switch i
         * case {1,2,3},
         */
        if (mclSwitchCompare(
              i,
              mlfCellhcat(
                mlfScalar(1.0), mlfScalar(2.0), mlfScalar(3.0), NULL))) {
            /*
             * method = strOpts{i};
             */
            mlfAssign(&method, mlfIndexRef(strOpts, "{?}", i));
        /*
         * case 4,
         */
        } else if (mclSwitchCompare(i, mlfScalar(4.0))) {
            /*
             * traceFlag = 1;
             */
            mlfAssign(traceFlag, mlfScalar(1.0));
        /*
         * case 5,         
         */
        } else if (mclSwitchCompare(i, mlfScalar(5.0))) {
            /*
             * traceFlag = 0;
             */
            mlfAssign(traceFlag, mlfScalar(0.0));
        /*
         * end
         */
        }
    /*
     * end
     */
    }
    /*
     * 
     * if ~isempty(traceFlag) & nopts > opt_indx,
     */
    {
        mxArray * a_ = mclInitialize(mlfNot(mlfIsempty(*traceFlag)));
        if (mlfTobool(a_) && mlfTobool(mlfAnd(a_, mlfGt(nopts, opt_indx)))) {
            mxDestroyArray(a_);
            /*
             * err_msg = 'The trailing string ''trace'' must be last in the input parameter list.'; 
             */
            mlfAssign(
              err_msg,
              mxCreateString(
                "The trailing string 'trace' must be "
                "last in the input parameter list."));
        } else {
            mxDestroyArray(a_);
        }
    /*
     * end
     */
    }
    /*
     * 
     * % [EOF] dpss.m
     */
    return_:
    mclValidateOutputs(
      "dpss/parseStrOpts", 3, nargout_, &method, traceFlag, err_msg);
    mxDestroyArray(i);
    mxDestroyArray(str1);
    mxDestroyArray(str2);
    mxDestroyArray(strOpts);
    return method;
}

/*
 * The function "mlfDpss_parseStrOpts" contains the normal interface for the
 * "dpss/parseStrOpts" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 339-385). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
static mxArray * mlfDpss_parseStrOpts(mxArray * * traceFlag,
                                      mxArray * * err_msg,
                                      mxArray * inputStr,
                                      mxArray * method_,
                                      mxArray * traceFlag_,
                                      mxArray * nopts,
                                      mxArray * opt_indx) {
    int nargout = 1;
    mxArray * method = mclGetUninitializedArray();
    mxArray * traceFlag__ = mclGetUninitializedArray();
    mxArray * err_msg__ = mclGetUninitializedArray();
    mlfEnterNewContext(
      2, 5, traceFlag, err_msg, inputStr, method_, traceFlag_, nopts, opt_indx);
    if (traceFlag != NULL) {
        ++nargout;
    }
    if (err_msg != NULL) {
        ++nargout;
    }
    method
      = Mdpss_parseStrOpts(
          &traceFlag__,
          &err_msg__,
          nargout,
          inputStr,
          method_,
          traceFlag_,
          nopts,
          opt_indx);
    mlfRestorePreviousContext(
      2, 5, traceFlag, err_msg, inputStr, method_, traceFlag_, nopts, opt_indx);
    if (traceFlag != NULL) {
        mclCopyOutputArg(traceFlag, traceFlag__);
    } else {
        mxDestroyArray(traceFlag__);
    }
    if (err_msg != NULL) {
        mclCopyOutputArg(err_msg, err_msg__);
    } else {
        mxDestroyArray(err_msg__);
    }
    return mlfReturnValue(method);
}

/*
 * The function "mlxDpss_parseStrOpts" contains the feval interface for the
 * "dpss/parseStrOpts" M-function from file
 * "/u4/local/matlab/toolbox/signal/signal/dpss.m" (lines 339-385). The feval
 * function calls the implementation version of dpss/parseStrOpts through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
static void mlxDpss_parseStrOpts(int nlhs,
                                 mxArray * plhs[],
                                 int nrhs,
                                 mxArray * prhs[]) {
    mxArray * mprhs[5];
    mxArray * mplhs[3];
    int i;
    if (nlhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: dpss/parseStrOpts Line: 339 Col"
            "umn: 0 The function \"dpss/parseStrOpts\" was called "
            "with more than the declared number of outputs (3)"));
    }
    if (nrhs > 5) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: dpss/parseStrOpts Line: 339 Col"
            "umn: 0 The function \"dpss/parseStrOpts\" was called "
            "with more than the declared number of inputs (5)"));
    }
    for (i = 0; i < 3; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 5 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 5; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 5, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    mplhs[0]
      = Mdpss_parseStrOpts(
          &mplhs[1],
          &mplhs[2],
          nlhs,
          mprhs[0],
          mprhs[1],
          mprhs[2],
          mprhs[3],
          mprhs[4]);
    mlfRestorePreviousContext(
      0, 5, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 3 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 3; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}
