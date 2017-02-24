#include <math.h>
#include <mex.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>


/* [xo, to] = filterFFT(t, x, lowCut, lowAlpha, highCut, highAlpha) 
   band-pass filter timestamped data
   INPUTS:
     t: the timestampes for input data (assumed in 1/10000 sec
     x: the data
     lowCut, lowAlpha: indicates the high-pass frequency (in Hz)
                       and high-pass rolloff
                       a lowAlpha = 0 means a high-pass filter
     highCut, highAlpha: indicates the low-pass	frequency (in Hz)
                         and low-pass rolloff
			 a highAlpha = 0 means a low-pass filter
   OUTPUTS:
     to: timestamps for output data
     xo: filtered data

batta 2002 initial version
*/

void makeFilterCoeffs(double *filterCoeffs, 
		      int nfilter,
		      double lowCut,
		      double lowAlpha,
		      double highCut,
		      double highAlpha);



void FilterFFT(double *z, int n, double *filterCoeffs, 
	       int nfilter, int padding,
	       double *zo);

void filterBuffer(double *buf, double *filterCoeffs, int nfilter);


#undef DEAL_WITH_GAPS


#define T_IIDX 0
#define X_IIDX 1
#define LOWCUT_IIDX 2
#define LOWALPHA_IIDX 3
#define HIGHCUT_IIDX 4
#define HIGHALPHA_IIDX 5

#define TO_OIDX 1
#define XO_OIDX 0




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{

  double *t, *z; /* pointers to input */

  int n; /* data length */
   
  double lowCut, lowAlpha, highCut, highAlpha; /* filter parameters */
  

  double *filterCoeffs;
  
  int nfilter = 65536; /* length of the filter segment, a power of two */
  double tolerance = 1.5;
  int padding;
  

  double timestep;
  
  double *zo, *to;
  
  int data_is_continuous;
  
  
  int i;
  
  /* check and unpack inputs */

  if(nrhs != 6)
    mexErrMsgTxt("6 inputs required");

  if(nlhs != 2 && nlhs != 1)
    mexErrMsgTxt("Requires 1 or 2 outputs.");

  if(mxGetM(prhs[T_IIDX]) != 1 && mxGetN(prhs[T_IIDX]) != 1)
    mexErrMsgTxt("t must be row or column");

  n = mxGetM(prhs[T_IIDX]) * mxGetN(prhs[T_IIDX]);
  
  t = mxGetPr(prhs[T_IIDX]);
  
  if(mxGetM(prhs[X_IIDX]) != 1 && mxGetN(prhs[X_IIDX]) != 1)
    mexErrMsgTxt("x must be row or column");
  
  if(n != mxGetM(prhs[X_IIDX])  * mxGetN(prhs[X_IIDX]))
    mexErrMsgTxt("t and x must have same length");

  z = mxGetPr(prhs[X_IIDX]);


  if( mxGetM(prhs[LOWCUT_IIDX])  * mxGetN(prhs[LOWCUT_IIDX]) != 1)
    mexErrMsgTxt("lowCut must be scalar");
  lowCut = mxGetScalar(prhs[LOWCUT_IIDX]);
  
  if( mxGetM(prhs[LOWALPHA_IIDX])  * mxGetN(prhs[LOWALPHA_IIDX]) != 1)
    mexErrMsgTxt("lowAlpha must be scalar");
  lowAlpha = mxGetScalar(prhs[LOWALPHA_IIDX]);
  
  if( mxGetM(prhs[HIGHCUT_IIDX])  * mxGetN(prhs[HIGHCUT_IIDX]) != 1)
    mexErrMsgTxt("highCut must be scalar");
  highCut = mxGetScalar(prhs[HIGHCUT_IIDX]);
  
  if( mxGetM(prhs[HIGHALPHA_IIDX])  * mxGetN(prhs[HIGHALPHA_IIDX]) != 1)
    mexErrMsgTxt("highAlpha must be scalar");
  highAlpha = mxGetScalar(prhs[HIGHALPHA_IIDX]);
  
  if(lowCut > highCut)
    mexErrMsgTxt("What are you trying to tell me, lowCut > highCut ???");
  


  /* prepare outputs */
  
  if(nlhs > 1)
    {
      plhs[TO_OIDX] = mxCreateDoubleMatrix(n, 1, mxREAL);
      to = mxGetPr(plhs[TO_OIDX]);
      memcpy(to, t, n * sizeof(double));
    }
  
  plhs[XO_OIDX] = mxCreateDoubleMatrix(n, 1, mxREAL);
  zo = mxGetPr(plhs[XO_OIDX]);
  

  while(n < nfilter)
    nfilter /= 2;
  

  

  timestep = (t[n-1] - t[0])/(n-1);


  data_is_continuous = 1;
  
  for(i=0; i < n-1 ; i++)
    {
      if(t[i+1] - t[i] > timestep * tolerance)
#ifndef DEAL_WITH_GAPS
	mexErrMsgTxt("Data has gaps, we're not equipped to deal with that yet");
#else 
      data_is_continuous = 0;
#endif      
    }

  timestep /= 10000.; /* express in sec */  
  lowCut *= nfilter * timestep;
  lowAlpha /= nfilter * timestep;
  highCut *= nfilter *timestep;
  highAlpha /= nfilter * timestep;

  padding = (int) floor(10. * nfilter / lowCut );


  filterCoeffs = (double *) mxCalloc(nfilter, sizeof(double));
  makeFilterCoeffs(filterCoeffs, nfilter, lowCut, lowAlpha, highCut, highAlpha);
  
  
  
  FilterFFT(z, n, filterCoeffs, nfilter, padding, zo);
  


  
  

  



  mxFree((void *)filterCoeffs);
  

  
}




void makeFilterCoeffs(double *fc, 
		      int nfilter, 
		      double lowCut,
		      double lowAlpha,
		      double highCut,
		      double highAlpha)
{
  
  int i;
  
  fc[0] = (1.+tanh(lowAlpha *(-lowCut))) * (1.+ tanh(highAlpha*highCut)) /4.;
  
  for(i=1; i <= nfilter/2; i++)
    {
      fc[i] = 
	(1.+tanh(lowAlpha *(i-lowCut))) * 
	(1.+ tanh(highAlpha*(highCut-i))) /4.;
      fc[nfilter-i] = 
	(1.+tanh(lowAlpha *(i-lowCut))) * 
	(1.+ tanh(highAlpha*(highCut-i))) /4.;
    }
  
}

  


void FilterFFT(double *z, int n, double *filterCoeffs, int nfilter, int padding, double *zo)
{
  int  offset = 0;
  int copy_size = nfilter - (padding * 2);
  

  double *buf;
  
  buf = (double *)mxCalloc(nfilter, sizeof(double));
  

  
  /* initial chunk */

  memcpy(buf, z, nfilter * sizeof(double));
  
  filterBuffer(buf, filterCoeffs, nfilter);
  
  memcpy(zo+offset, buf, (padding+copy_size) * sizeof(double));
  
  offset += copy_size;
  

  /* do intermediate chunks */
  do 
    {
      memcpy(buf, z+offset, nfilter * sizeof(double));
      
      filterBuffer(buf, filterCoeffs, nfilter);
      
      memcpy(zo+(offset+padding), buf+padding, copy_size * sizeof(double));
      
      offset += copy_size;

    }
  while ((offset+nfilter) < n);
  
  
  /* final chunk */

  memcpy(buf, z + (n - nfilter), nfilter * sizeof(double));

  filterBuffer(buf, filterCoeffs, nfilter);

  memcpy(zo+(offset+padding), buf + (nfilter-(n-(offset+padding))) , (n-(offset+padding))*sizeof(double));
  

  mxFree((void *)buf);
  
  
}


void filterBuffer(double *buf, double *fc, int nfilter)
{

  int i;
  
  gsl_fft_real_radix2_transform(buf, 1, nfilter);
  
  for(i = 0; i < nfilter ; i++)
    buf[i] *= fc[i];
  
  gsl_fft_halfcomplex_radix2_inverse(buf, 1, nfilter);
}
