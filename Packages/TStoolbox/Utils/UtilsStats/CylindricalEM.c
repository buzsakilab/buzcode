/*  Cylindrical EM: 
    Expectation-Maximization algorithm for bivariate linear-circular data. 
    The distribution is gaussian in the linear dimension, Wrapped
    Normal for the circular dimension. The approximations performed
    here assume that the spread of each cluster along the circular
    dimension is << 2*PI. 
    Useful e.g. for clustering of theta-phase precession data, or any
    other kind of time, space dependent phase data


    batta 2003
    UNDER CONSTRUCTION

*/


#include <math.h>

#define MEXFILE

#ifdef MEXFILE

#include <matrix.h>
#include <mex.h>

#define calloc mxCalloc
#define free mxFree
#endif

double min_det = 1.;


/* parameters for a gaussian component */ 
typedef struct 
{
  double Pj;  /* the priori probability for the component */
  double muTheta; /* the mean on the circular dimension */
  double muY; /* the mean on teh linear direction */
  double S[3]; /* the covariance matrix:
	       S[0] = \sigma_{11}
	       S[1] = \sigma_{22}
	       S[2] = \sigma_{12}
	       */
  
  /* back up copy of the parameters, for convergence criterion
     purposes */
  double Pj_old;
  double muTheta_old;
  double muY_old;
  double S_old[3];
  double convergence;
  

  double detS; /* the determinant of S */
  double K; /* prefactor in the probability distribution */
  double A[3]; /* inverse matrix of S:
	       A[0] = a_{11}
	       A[1] = a_{22}
	       A[2] = a_{12}
	  */
  double thSlope; /* the slope of the component:
		     E(theta|y) = muTheta + thSlope (y - muY) */
  double r; /* correlation coefficient for the component */
   
   
} Component;  

 

/* update all the auxiliary members of the comp struct */
void updateComps(Component *comp_ptr)
{
  double convergence, ca;

  

  

  comp_ptr->detS = comp_ptr->S[0] * comp_ptr->S[1] - comp_ptr->S[2] * comp_ptr->S[2];
  
  /* avoid the gaussians becoming too small */


  if(comp_ptr->detS == 0.)
    {
/*       mexPrintf("Det = 0"); */
      comp_ptr->S[0] = 1.;
      comp_ptr->S[1] = 1.;
      comp_ptr->S[2] = 0.;
    }
  



  if(comp_ptr->detS < min_det)
    {
      double y1= 1;
      double th1 = 2 * M_PI / 180;
      
/*       mexPrintf("detS = %g\n", comp_ptr->detS); */
      
/*       comp_ptr->S[0] /= sqrt(comp_ptr->detS/min_det); */
/*       comp_ptr->S[1] /= sqrt(comp_ptr->detS/min_det); */
/*       comp_ptr->S[2] /= sqrt(comp_ptr->detS/min_det); */
      comp_ptr->S[0] = th1 /sqrt(min_det * y1 * th1);
      comp_ptr->S[1] = y1/sqrt(min_det * y1 * th1);
      comp_ptr->S[2] = 0.;

      comp_ptr->detS = min_det;
    }
  

  comp_ptr->K = 1. / ( 2. * M_PI * sqrt(comp_ptr->detS));
  comp_ptr->A[0] = comp_ptr->S[1] / comp_ptr->detS;
  comp_ptr->A[1] = comp_ptr->S[0] / comp_ptr->detS;
  comp_ptr->A[2] = - comp_ptr->S[2] / comp_ptr->detS;
  comp_ptr->thSlope = comp_ptr->S[2]/comp_ptr->S[1];
  comp_ptr->r = comp_ptr->S[2] / sqrt(comp_ptr->S[0] * comp_ptr->S[1]);
  
  
  convergence = 0.;

  ca = fabs(comp_ptr->Pj - comp_ptr->Pj_old) / comp_ptr->Pj;
  if(ca > convergence)
    convergence = ca;

  ca = fabs(comp_ptr->muTheta - comp_ptr->muTheta) / 2 * M_PI;
  if(ca > convergence)
    convergence = ca;
  
  ca = fabs(comp_ptr->muY - comp_ptr->muY) / comp_ptr->muY;
  if(ca > convergence)
    convergence = ca;
  
  ca = fabs(comp_ptr->S[0] - comp_ptr->S[0]) / comp_ptr->S[0];
  if(ca > convergence)
    convergence = ca;
 
  ca = fabs(comp_ptr->S[1] - comp_ptr->S[1]) / comp_ptr->S[1];
  if(ca > convergence)
    convergence = ca;
 
  ca = fabs(comp_ptr->S[2] - comp_ptr->S[2]) / comp_ptr->S[2];
  if(ca > convergence)
    convergence = ca;
 
  comp_ptr->convergence = convergence;
  

  comp_ptr->Pj_old = comp_ptr->Pj;
  comp_ptr->muTheta_old = comp_ptr->muTheta;
  comp_ptr->muY_old = comp_ptr->muY;
  comp_ptr->S_old[0] = comp_ptr->S[0];
  comp_ptr->S_old[1] = comp_ptr->S[1];
  comp_ptr->S_old[2] = comp_ptr->S[2];

}


inline double modPi(double th)
{
  double th1;
  
  th1 = fmod(th, 2 * M_PI);
  if(th1 > M_PI)
    th1 = 2 * M_PI - th1;
  
  return th1;
  

}

/* pick up the "right value" for the theta for probability computation
   */
double thStar(double th, double y, Component comp)
{
  double th_star = comp.muTheta + comp.thSlope * (y - comp.muY);
  
  return th_star + modPi(th-th_star);
}



double prob_comp(double th, double y, Component comp)
{
  double th0;
  double p;
  
  th0 = thStar(th, y, comp);
  

  p = comp.K * exp(- 0.5 * (comp.A[0] *(th0-comp.muTheta) * (th0-comp.muTheta) +
		       comp.A[1] *(y - comp.muY) * (y - comp.muY) + 
		       2. * comp.A[2] * (y - comp.muY) *(th0-comp.muTheta)));
  return p;
  
}





/* performs the EM
   N: number of data points
   theta: circular data
   y: linear data
   M: number of components
   comps: component parameters, supposed to be initialized 

   returns the final value of the likelihood
*/
double doEM(int N, double *theta, double *y, int M, Component *comps)
{
  
  double *P_x; /* probability of data points, give the "old" model */
  double **P_jx; /* posterior probability of j-th component, in the "old"
		 model, conditional to data point */
  int converged;
  double convergence;
  double criterion;
  double L;
  
  criterion = 0.01;
  
  int maxIter = 1000; /* max number fo iterations */
  int i, j;
  
  for(i = 0; i < M; i++) updateComps(&(comps[i]));

  P_x = (double *) calloc(N, sizeof(double));
  
  P_jx = (double **) calloc(M, sizeof(double *));
  for(i = 0; i < M; i++) P_jx[i] = (double *)calloc(N, sizeof(double));
  
   
  converged = 0;
  

  


  int iter = 0;

  while((!converged) && iter < maxIter)
    {
      double muTheta, muY;
      double S[3];
      double Pj;
      
      convergence = 0.;
      iter++;
      

      for(i = 0; i < N; i++)
	{
	  P_x[i] = 0.;
	  for(j = 0; j < M; j++)
	    P_x[i] += comps[j].Pj * prob_comp(theta[i], y[i], comps[j]);
/* 	  if(P_x[i] <= 0.) */
/* 	    mexPrintf("P_x[%d] = %g\n", i, P_x[i]); */
	  
	}
      
      
      for(i = 0; i < N; i++)
	for(j = 0; j < M; j++)
	  {
	    
	    P_jx[j][i] = prob_comp(theta[i], y[i], comps[j]) * comps[j].Pj /
	      P_x[i];
/* 	    if(isnan(P_jx[j][i])) */
/* 	      { */
		
/* 		mexPrintf("P_jx[%d][%d] is nan\n", i, j); */
/* 		mexPrintf("muTheta = %g muY = %g S[0] = %g S[1] = %g S[2] = %g pj = %g detS = %g\n", */
/* 			  comps[j].muTheta,  */
/* 			  comps[j].muY, */
/* 			  comps[j].S[0], */
/* 			  comps[j].S[1], */
/* 			  comps[j].S[2], */
/* 			  comps[j].Pj, */
/* 			  comps[j].detS); */

/* 	      } */
	    


	  }
      
      

      for(j = 0; j < M; j++)
	{
	  S[0] = 0.;
	  S[1] = 0.;
	  S[2] = 0.;
	  muTheta = 0.;
	  muY = 0.;
	  Pj = 0.;

	  for(i = 0; i < N; i++)
	    Pj += P_jx[j][i];
	  
	  for(i = 0; i < N; i++)	  
	    {

	      muTheta += theta[i] * P_jx[j][i];
	      muY += y[i] * P_jx[j][i];
	    }
	  
	  muTheta /= Pj;
	  muY /= Pj;
	  
	  comps[j].muTheta = muTheta;
	  comps[j].muY = muY;
	  

	  for(i = 0; i < N; i++)
	    {
	      double th0;
	      
	      th0 = thStar(theta[i], y[i], comps[j]);
	      S[0] += (th0-muTheta) * (th0-muTheta) * P_jx[j][i];
	      S[1] += (y[i]-muY) * (y[i]-muY)  * P_jx[j][i];
	      S[2] += (th0-muTheta) * (y[i]-muY) * P_jx[j][i];
	    }
	  
	  comps[j].S[0] = S[0] / Pj;
	  comps[j].S[1] = S[1] / Pj;
	  comps[j].S[2] = S[2] / Pj;
	  
	  comps[j].Pj = Pj / N;
	  
	  updateComps( &(comps[j]) );
	  
	  if(isnan(comps[j].muTheta))
	    {
	      L = 0. / 0.;
	      mexPrintf("L NaN\n");
	      mexPrintf("iter = %d S = %g %g %g\n", iter, comps[j].S[0], 
			comps[j].S[1], 
			comps[j].S[2]);
	      
	      goto the_end;
	    }
	  
	  
	  if(comps[j].convergence > convergence)
	    convergence = comps[j].convergence;
	  
	  
	}
      
      if(convergence > criterion)
	converged = 0;
      else
	converged = 1;
      
      
    }
  
  L = 0.;
  
  for(i = 0; i < N; i++)
    {
      P_x[i] = 0.;
      for(j = 0; j < M; j++)
	P_x[i] += comps[j].Pj * prob_comp(theta[i], y[i], comps[j]);
      L += log(P_x[i]);
      
    }
  

 the_end:

  free(P_x);
  for(i = 0; i < M; i++) free(P_jx[i]);
  free(P_jx);
  



  return (-L);
  
}



#define THETA_INP 0
#define Y_INP 1
#define MU_THETA_INP 2
#define MU_Y_INP 3
#define S_INP 4
#define PJ_INP 5

#define MU_THETA_OUT 0
#define MU_Y_OUT 1
#define S_OUT 2
#define PJ_OUT 3
#define L_OUT 4

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{


  double *theta, *y, *pr, L;
  int M, N;
  int i;
  const int *dims;
  
  Component *comps;
  
  /* check number of inputs and outputs */

  if(nrhs != 2 && nrhs != 6)
    mexErrMsgTxt("Must have 2 or 6 inputs");
  
  if(nlhs != 5)
    mexErrMsgTxt("Must have 5 outputs");
  
  


  /* unpacking inputs */

  if(mxGetM(prhs[THETA_INP]) != 1 && mxGetN(prhs[THETA_INP]) != 1)
    mexErrMsgTxt("Theta must be row or column vector");

  
  if(mxGetM(prhs[Y_INP]) != 1 && mxGetN(prhs[Y_INP]) != 1)
    mexErrMsgTxt("Y must be row or column vector");
  
  theta = mxGetPr(prhs[THETA_INP]);
  N = mxGetM(prhs[THETA_INP]) * mxGetN(prhs[THETA_INP]);
  
  y = mxGetPr(prhs[Y_INP]);
  if(mxGetM(prhs[Y_INP]) * mxGetN(prhs[Y_INP]) != N)
    mexErrMsgTxt("theta and y must have same length");
  


  if(nrhs > 2)
    {
      if(mxGetM(prhs[MU_THETA_INP]) != 1 && mxGetN(prhs[MU_THETA_INP]) != 1)
	mexErrMsgTxt("Mu_theta must be row or column vector");

      M = mxGetM(prhs[MU_THETA_INP]) * mxGetN(prhs[MU_THETA_INP]);
      
      comps = (Component *)calloc(M, sizeof(Component));
      pr = mxGetPr(prhs[MU_THETA_INP]);
      
      for(i = 0; i < M; i++)
	comps[i].muTheta = pr[i];
      
      if(mxGetM(prhs[MU_Y_INP]) != 1 && mxGetN(prhs[MU_Y_INP]) != 1)
	mexErrMsgTxt("Mu_y must be row or column vector");

      if(mxGetM(prhs[MU_Y_INP]) * mxGetN(prhs[MU_Y_INP]) != M)
	mexErrMsgTxt("all parameters should reflect the same number of components, M");
      
      pr = mxGetPr(prhs[MU_Y_INP]);
      
      for(i = 0; i < M; i++)
	comps[i].muY = pr[i];
     
  
      if(mxGetNumberOfDimensions(prhs[S_INP]) != 3)
	mexErrMsgTxt("S must be a 3-D array");
      
      dims = mxGetDimensions(prhs[S_INP]);
      if(dims[0] != M || dims[1] != 2 || dims[2] != 2)
	mexErrMsgTxt("S must be M x 2 x 2");
      
      
      pr = mxGetPr(prhs[S_INP]);
      
      for(i = 0; i < M; i++)
	{
	  int ix[3];
	  
	  ix[0] = i;
	  

	  ix[1] = 0;
	  ix[2] = 0;
	  comps[i].S[0] = pr[mxCalcSingleSubscript(prhs[S_INP], 3, ix)];
	  
	  ix[1] = 1;
	  ix[2] = 1;
	  comps[i].S[1] = pr[mxCalcSingleSubscript(prhs[S_INP], 3, ix)];
	  
	  ix[1] = 0;
	  ix[2] = 1;
	  comps[i].S[2] = pr[mxCalcSingleSubscript(prhs[S_INP], 3, ix)];
	  
	}
      
      

     if(mxGetM(prhs[PJ_INP]) != 1 && mxGetN(prhs[PJ_INP]) != 1)
	mexErrMsgTxt("Pj must be row or column vector");

      if(mxGetM(prhs[PJ_INP]) * mxGetN(prhs[PJ_INP]) != M)
	mexErrMsgTxt("all parameters should reflect the same number of components, M");
      
      pr = mxGetPr(prhs[PJ_INP]);
      
      for(i = 0; i < M; i++)
	comps[i].Pj = pr[i];
      
    }
  else
    mexErrMsgTxt("Not equipped yet to do independent initalization of components");
  
  
  L = doEM(N, theta, y, M, comps);
  


  /* pack outputs */

  plhs[MU_THETA_OUT] = mxCreateDoubleMatrix(M, 1, mxREAL);
  pr = mxGetPr(plhs[MU_THETA_OUT]);
  
  for(i = 0; i < M; i++)
    pr[i] = comps[i].muTheta;
  
  plhs[MU_Y_OUT] = mxCreateDoubleMatrix(M, 1, mxREAL);
  pr = mxGetPr(plhs[MU_Y_OUT]);
  
  for(i = 0; i < M; i++)
    pr[i] = comps[i].muY;
  
  {
    int ix[3];
    
    ix[0] = M;
    ix[1] = 2;
    ix[2] = 2;
    plhs[S_OUT] = mxCreateNumericArray(3, ix, mxDOUBLE_CLASS, mxREAL);
    pr = mxGetPr(plhs[S_OUT]);

  }
  

  for(i = 0; i < M; i++)
    {
      int ix[3];
      
      ix[0] = i;
      

      ix[1] = 0;
      ix[2] = 0;
      pr[mxCalcSingleSubscript(plhs[S_OUT], 3, ix)] = comps[i].S[0];
      
      ix[1] = 1;
      ix[2] = 1;
      pr[mxCalcSingleSubscript(plhs[S_OUT], 3, ix)] = comps[i].S[1];
      
      ix[1] = 0;
      ix[2] = 1;
      pr[mxCalcSingleSubscript(plhs[S_OUT], 3, ix)] = comps[i].S[2];
      
      ix[1] = 1;
      ix[2] = 0;
      pr[mxCalcSingleSubscript(plhs[S_OUT], 3, ix)] = comps[i].S[2];
      
    }
  

  plhs[PJ_OUT] = mxCreateDoubleMatrix(M, 1, mxREAL);
  pr = mxGetPr(plhs[PJ_OUT]);
  
  for(i = 0; i < M; i++)
    pr[i] = comps[i].Pj;

  plhs[L_OUT] = mxCreateScalarDouble(L);
  
  free(comps);
  


  
}

  
