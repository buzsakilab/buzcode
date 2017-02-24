
/* [A, lambda, b] = exp_fit_gsl(y, sigma)
   tries an exponential fit of the type 
   y(i) = A * exp(- lambda * t) + b on data y with sdev sigma
*/


#include <mex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>



struct data {
  size_t n;
  double * y;
  double * sigma;
};



int
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  mexPrintf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2), 
          gsl_blas_dnrm2 (s->f));
  return 0;
  
}

int
expb_f (const gsl_vector * x, void *params, 
        gsl_vector * f)
{
  size_t n = ((struct data *)params)->n;
  double *y = ((struct data *)params)->y;
  double *sigma = ((struct data *) params)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  double b = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = A * exp(-lambda * i) + b */
      double t = i;
      double Yi = A * exp (-lambda * t) + b;
      gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
    }

  return GSL_SUCCESS;
}

int
expb_df (const gsl_vector * x, void *params, 
         gsl_matrix * J)
{
  size_t n = ((struct data *)params)->n;
  double *sigma = ((struct data *) params)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      double t = i;
      double s = sigma[i];
      double e = exp(-lambda * t);
      gsl_matrix_set (J, i, 0, e/s); 
      gsl_matrix_set (J, i, 1, -t * A * e/s);
      gsl_matrix_set (J, i, 2, 1/s);

    }
  return GSL_SUCCESS;
}

int
expb_fdf (const gsl_vector * x, void *params,
          gsl_vector * f, gsl_matrix * J)
{
  expb_f (x, params, f);
  expb_df (x, params, J);

  return GSL_SUCCESS;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
                 *prhs[])
{

  double *y, *sigma;
  size_t n;
  size_t p = 3;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  struct data d;
  gsl_multifit_function_fdf f;
  double x_init[3] = { 1.0, 0.05, 0.0 };

  gsl_vector_view x = gsl_vector_view_array (x_init, p);



  int status;
  size_t  iter = 0;

  if(nrhs != 2)
    mexErrMsgTxt("2 inputs required");

  if(nlhs != 3)
    mexErrMsgTxt("Requires three output.");
  

  if(mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("y must be row or column");
  
  y = mxGetPr(prhs[0]);

  n = mxGetM(prhs[0])  * mxGetN(prhs[0]);
  
  if(mxGetM(prhs[1]) != 1 && mxGetN(prhs[1]) != 1)
    mexErrMsgTxt("sigma must be row or column");
  
  if(n != mxGetM(prhs[0])  * mxGetN(prhs[0]))
    mexErrMsgTxt("y and sigma must have same length");

  sigma = mxGetPr(prhs[1]);
  

  d.n = n;
  d.y = y;
  d.sigma = sigma;
  
  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);


  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      mexPrintf ("status = %s\n", gsl_strerror (status));

      print_state (iter, s);

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);
  gsl_multifit_covar (s->J, 0.0, covar);


#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  plhs[0] = mxCreateScalarDouble( (FIT(0)));
  plhs[1] = mxCreateScalarDouble( (FIT(1)));
  plhs[2] = mxCreateScalarDouble( (FIT(2)));

  mexPrintf("Errors: %g %g %g\n", ERR(0), ERR(1), ERR(2));
  

  gsl_multifit_fdfsolver_free (s);
}
