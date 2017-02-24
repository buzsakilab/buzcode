/* cross_decode.c derived from inform.c and MEX-ified for use with MATLAB
   batta 1999 
   usage
   [DEC_STIM, PROB_MAX, DEC_HIST] = cross_decode(S1, R1, n_stim1, S2, R2, n_stim2, cross_decode_opt)
*/









/* ntn/src/inform.c (derived from masters inf_bin.c and infom_new.c)
   used to extract the amount of information, in bits, stimuli decoded from 
   population responses carry about actual stimuli, and the percent correct  
   developed by Alessandro Treves, ale@limbo.sissa.it, with
   Edmund Rolls, Stefano Panzeri, Bill Skaggs, during 1994-95  */

/* 10/04/95: version for use with Tucson style data
    with cross validation, frequencies P^Q(s,s_p) & probabilities P^S(s,s') 
   01/08/95: added -dec switch for decoding and probability est. algorithms:  
    1 - PE (Bayesian) probability estimation with Gaussian fits (default)
    5 - Bayesian probability estimation with Poisson fits,
    6 - DP uses normalized dot thresholded products directly as probabilities, 
    7 - uses correlations (Pearson product moments) raised to a high power,
    8 - uses distance between current and average vectors, exponentiated. 
   30/08/95: reads ascii infom.dat if run with the -ascii option 
   08/09/95: efficient cross validation with only one trial used for testing,
    and the probability table averaged over all different test trials 
    (older cross validation remains in old programs only) 
   02/02/96: add up a scrambling routine (SP)    
   16/08/96: small tuning of the finite size corrections (gg parameter) */

/* usage         OPTION               MEANING                         DEFAULT 
        inform -seed 989765     : change random seed                   243342
	       -r filename.r    : name rate (input) file               stdin
               -o xxx.out       : name (main) output file              stdout
               -p xxx.res       : name results file for gnuplot        res
               -q xxx.pcc       : name auxiliary output file           pcc
               -max_t 10        : max # trials per stimulus allowed    2*max_s
               -min_t 10        : min #   "     "     "        "       2*max_s
	       -nc_min 4        : min population size                  1
	       -nc_max 26       : max population size                  max_c
	       -niter_fact 3    : this times (num_c-nc+1), i.e. niter,
	                          is the # of population samples of 
	                          given size averaged together         10
	       -dec 7           : decoding algorithm                   1
	       -s0err           : throw trials with stimulus = 0       False
	       -s_throw 4       : throw eg 4 stimuli at random, to
	                          try with a reduced subset            0 
	       -ascii           : input file is ascii                  False 
	       -scramble        : data scrambling                      False */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"

#define True    1
#define False   0
double  drand48();
void    Read_Input(int, int);

#undef DEBUG1
#undef  DEBUG2
#undef DEBUG3
#define max_c   77		/* max no of cells */
#define max_s   50		/* max no of stimuli */

void    Decode(void);
void    Read_Rfile(char *, double ****);
void    Read_Rfile_seq(char *, double ***, double **);
double  phims(double val, double mean, double sigma);

double  *rc;		/* current response */
int     *ic;		/* cell code */
int     nc;                     /* cells in sample being analysed */
double  **rcs;	/* prob distr mean rates */
double  **acs;	/* prob distr nonzero rates */
double  **ra;	/* (training set) averages */
double  **sd;	/* (training set) variances */
double  *ra0,*sd0;	/* mean avg. and var. across all stim. */
double  st_dev;                 /* global mean variance across cells */ 
double  **f0;	/* (tr. set) fract. of zero rates */
int     num_s, num_s1, num_s2;                  /* # stimuli found by Read_file */
int     best_s;			/* predicted (most probable) stimulus */
double  *pq_r;		/* probabilities of each stimulus */

double ***r1;			/* responses (Hz) (training set) */
double ***r2;                   /* responses (Hz) (decoding set) */
double ***r;
double **rdec;                  /* responses in the decoding set,
** listed by trial */ 
double *stim;                   /* stimuli by trial */

int num_trials, num_trials_dec;


int     max_t;
int     min_t;
int     num_c, num_c1, num_c2;                  /* actual # cells */
int    *nr1;		        /* number of trials (training set) */
int    *nr2;                    /* number fo trials (decoding set) */ 
int    *nr;

int     n_tot, n_tot1, n_tot2;	/* total of above */
double  *x_s;		/* cumul. fract. of data per stimul. */
double  *Ps;		/* Prob(s)    */
int     s0err = 0;       
int     read_ASCII = 0;
int     s_throw = 0;            /* stimuli to be discarded */
int     niter,niter_fact;       /* # samples of size nc */
int     decoding_method = 1;	/* default decoding method is PE */
int     scramble = False ;      /* scrambling of data            */
int     nonorm=0;               /* if non-zero, the histogrma is not
                                ** normalized */


int     n_trials;



double  p_max;
double  **dec_hist;

double p_thresh;
double p_thresh_fact = 2.;
int    *good_train_stim;

double x_prev, pcc0, rt, dp;
double **Pmix;
double *DecStim;
double  *ProbMax;

int s1;
double *stimuli, **responses, *stimuli2, **responses2;   /* the input array copied from the */
				/*   mxArrays */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
		 *prhs[])
{
   int     nc_min, nc_max;
   int     seed;
   int     cur_arg;
   double  r_max;		        /* max response of a cell */
   int     c, s, t;		        /* cell,stimulus,trial,bin */
   int     ic1;	        /* auxiliary cell indexes */
   int i,j;
   
   int     nt;                          /* number of trials */
   double  eps = 0.000001;              /* small number */
   int mrows, ncols;
  double *ptr;
  mxArray *field_value;
  
   

  if(nrhs != 6 && nrhs != 7)
    mexErrMsgTxt("Six (or seven) inputs required");

   if(nlhs != 3)
    mexErrMsgTxt("Requires three outputs.");
  

#ifndef DEBUG1

  /* input 0 is the training set stimuli vector */
   
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  

  if(!(mrows == 1 || ncols == 1))
    mexErrMsgTxt("Input 1 must be a row or column vector");
  if(mrows == 1)
    num_trials = ncols;
  else
    num_trials = mrows;
  mexPrintf("num_trials = %d\n", num_trials);
  
  stimuli = mxGetPr(prhs[0]);

   /* input 1 is the training set responses matrix */
   
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  
  if(mrows != num_trials)
    mexErrMsgTxt("Input 1 and 2 must have the same number of rows.");
  
  num_c1 = ncols;
  num_c = num_c1;
  
  mexPrintf("Number of cells: %d, number of trials: %d\n", num_c, num_trials);
  

  ptr = mxGetPr(prhs[1]);

  if(!(responses = (double **)mxCalloc(num_trials, sizeof(double *))))
    mexErrMsgTxt("Mallocing responses\n");
  
  for(i = 0; i < num_trials; i++)
    if(!(responses[i] = mxCalloc(num_c1, sizeof(double))))
      mexErrMsgTxt("mallocing responses 2\n");
  

  for(i = 0; i < num_trials; i++)
    for(j = 0; j < num_c1; j++)
      {
	responses[i][j] = ptr[i + num_trials * j];
	if(responses[i][j] > 70) 
	  mexPrintf("responses[%d][%d] = %g", i, j, responses[i][j]);
      }



/* input 2 is the number of stimuli (scalar) */
 
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  
  if(mrows * ncols != 1)
    mexErrMsgTxt("input 3 must be a scalar.");
  
  ptr = mxGetPr(prhs[2]);
  
  num_s1 = ptr[0];

  max_t = 2*max_s;
  min_t = 2*max_s;
  nc_min = 1;
  nc_max = max_c;
  seed = 243342;
  niter_fact = 1; 

 

  /* input 3 is the decoding set stimuli vector */
  
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  

  if(!(mrows == 1 || ncols == 1))
    mexErrMsgTxt("Input 1 must be a row or column vector");
  if(mrows == 1)
    num_trials_dec = ncols;
  else
    num_trials_dec = mrows;
  mexPrintf("I: num_trials_dec = %d\n", num_trials_dec);
  

  stimuli2 = mxGetPr(prhs[3]);


   /* input 4 is the decoding set responses matrix */

  mrows = mxGetM(prhs[4]);
  ncols = mxGetN(prhs[4]);
  
  if(mrows != num_trials_dec)
    mexErrMsgTxt("Input 1 and 2 must have the same number of rows.");
  
  num_c2 = ncols;
  mexPrintf("num_c2 = %d\n", num_c2);
  

  ptr = mxGetPr(prhs[4]);
#else 
  num_trials = 1496;
  num_trials_dec = 1496;
  
#endif 

  if(!(responses2 = (double **)mxCalloc(num_c2, sizeof(double *))))
     mexErrMsgTxt("Mallocing responses2 step 1");
  for(i = 0; i < num_c2; i++)
    if(!(responses2[i] = mxCalloc(num_trials_dec, sizeof(double))))
       mexErrMsgTxt("Mallocing responses2 step 2");

  for(i = 0; i < num_trials_dec; i++)
    for(j = 0; j < num_c2; j++)
      {
	responses2[j][i] = ptr[i + num_trials_dec * j];
	if(responses2[j][i] > 70) 
	  mexPrintf("responses2[%d][%d] = %g", i, j, responses2[j][i]);
      }


   

/* input 5 is the number of stimuli (scalar) */

  mrows = mxGetM(prhs[5]);
  ncols = mxGetN(prhs[5]);
  
  if(mrows * ncols != 1)
    mexErrMsgTxt("input 3 must be a scalar.");
  
  ptr = mxGetPr(prhs[5]);
  
  num_s2 = ptr[0];

  max_t = 2*max_s;
  min_t = 2*max_s;
  nc_min = 1;
  nc_max = max_c;
  seed = 243342;
  niter_fact = 1; 
       



  /* input 6 is a structure containing the options passed to infors,
     if a field contains the empty array the value is not changed */

  if(nrhs > 6)
    {
      if(!mxIsStruct(prhs[6]))
	mexErrMsgTxt("Input 7 must be a structure containing the options passed to informs");
      
      field_value = mxGetField(prhs[6], 0, "seed");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  seed = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "max_t");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  max_t = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "min_t");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  min_t = (int)ptr[0];
	  
	}
      field_value = mxGetField(prhs[6], 0, "nc_min");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  nc_min = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "nc_max");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  nc_max = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "s0err");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  s0err = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "s_throw");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  s_throw = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "niter_fact");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  niter_fact = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "scramble");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  scramble = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "decoding_method");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  decoding_method = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "nonorm");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  nonorm = (int)ptr[0];
	}
      field_value = mxGetField(prhs[6], 0, "thresh");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  p_thresh_fact = ptr[0];
	}
      
    }
  


      /* allocate the memory for all the variables */
      
      {
	int i;

#ifndef DEBUG1
	int num_c_max = (num_c1 > num_c2 ? num_c1 : num_c2);
	int num_s_max = (num_s1 > num_s2 ? num_s1 : num_s2);
#else
	int num_c_max = 36;
	int num_s_max = 4;
	num_trials_dec = 1496;
	
#endif
	if(!(rc = (double *)mxCalloc(num_c_max, sizeof(double))))
	  mexErrMsgTxt("Malloc error 1");

	if(!(ic = (int *)mxCalloc(num_c_max, sizeof(int))))
	  mexErrMsgTxt("Malloc error 2");


	if(!(ra0 = (double *)mxCalloc(num_c_max, sizeof(double))))
	  mexErrMsgTxt("Malloc error 3");

	if(!(sd0 = (double *)mxCalloc(num_c_max, sizeof(double))))
	  mexErrMsgTxt("Malloc error 4");


	if(!(rcs = (double **)mxCalloc(num_c_max, sizeof(double *))))
	  mexErrMsgTxt("Malloc error 5");
	for(i = 0; i < num_c_max; i++)
	  if(!(rcs[i] = (double *)mxCalloc(num_s_max, sizeof(double))))
	    mexErrMsgTxt("Malloc error 6");

	if(!(acs = (double **)mxCalloc(num_c_max, sizeof(double *))))
	  mexErrMsgTxt("Malloc error 7");
	for(i = 0; i < num_c_max; i++)
	  if(!(acs[i] = (double *)mxCalloc(num_s_max, sizeof(double))))
	    mexErrMsgTxt("Malloc error 8");

	if(!(ra = (double **)mxCalloc(num_c_max, sizeof(double *))))
	  mexErrMsgTxt("Malloc error 9");
	for(i = 0; i < num_c_max; i++)
	  if(!(ra[i] = (double *)mxCalloc(num_s_max, sizeof(double))))
	    mexErrMsgTxt("Malloc error 10");

	if(!(sd = (double **)mxCalloc(num_c_max, sizeof(double *))))
	  mexErrMsgTxt("Malloc error 11");
	for(i = 0; i < num_c_max; i++)
	  if(!(sd[i] = (double *)mxCalloc(num_s_max, sizeof(double))))
	    mexErrMsgTxt("Malloc error 12");  

	if(!(f0 = (double **)mxCalloc(num_c_max, sizeof(double *))))
	  mexErrMsgTxt("Malloc error 13");
	for(i = 0; i < num_c_max; i++)
	  if(!(f0[i] = (double *)mxCalloc(num_s_max, sizeof(double))))
	    mexErrMsgTxt("Malloc error 14");  

	if(!(pq_r = (double *)mxCalloc(num_s_max, sizeof(double))))
	  mexErrMsgTxt("Malloc error 15");
	if(!(x_s = (double *)mxCalloc(num_s_max, sizeof(double))))
	  mexErrMsgTxt("Malloc error 16");
	if(!(Ps = (double *)mxCalloc(num_s_max, sizeof(double))))
	  mexErrMsgTxt("Malloc error 17");

	if(!(DecStim = (double *)mxCalloc(num_trials_dec, sizeof(double))))
	  mexErrMsgTxt("Malloc error 18");
	if(!(ProbMax = (double *)mxCalloc(num_trials_dec, sizeof(double))))
	  mexErrMsgTxt("Malloc error 19");

	if(!(Pmix = (double **)mxCalloc(num_s_max, sizeof(double *))))
	  mexErrMsgTxt("Malloc error 20");
	for(i = 0; i < num_s_max; i++)
	  if(!(Pmix[i] = (double *)mxCalloc(num_s_max, sizeof(double))))
	    mexErrMsgTxt("Malloc error 21");  


      }
      
      




/*                  *************************************
                    read data from the training set and set constant parameters
                    *************************************    */ 
                    
   srand48(seed);
/*    Read_Rfile(rfile1_name, &r); */
#ifndef DEBUG2
    Read_Input(num_s1, num_c1); 
      mexPrintf("Done Read_Input\n");
#endif   
#ifdef DEBUG3
   {
     plhs[0] = mxCreateDoubleMatrix(1, num_trials_dec, mxREAL);
     ptr = mxGetPr(plhs[0]);
     memcpy(ptr, DecStim, num_trials_dec * sizeof(double));
     mexPrintf("Done DecStim\n");
     

     plhs[1] = mxCreateDoubleMatrix(1, num_trials_dec, mxREAL);
     ptr = mxGetPr(plhs[1]);
     memcpy(ptr, ProbMax, num_trials_dec * sizeof(double));
     mexPrintf("Done ProbMax\n");
     mexPrintf("num_s1 = %d num_s2 = %d\n", num_s1, num_s2);
     plhs[2] = mxCreateDoubleMatrix(num_s1, num_s2, mxREAL);
     mexPrintf("Created Matrix\n");
     
     ptr = mxGetPr(plhs[2]);
     mexPrintf("Got Pointer\n");
     
     {
       int i, j;
       
       dec_hist = (double **)mxCalloc(num_s2, sizeof(double *));
   for(s=0;s<num_s2;s++)
     dec_hist[s] = (double *)mxCalloc(num_s1, sizeof(double));
       for(i = 0;i < num_s2; i++)
	 for(j = 0; j < num_s1; j++)
	   ptr[i+num_s1*j] = dec_hist[i][j];
     }
     mexPrintf("Done dec_hist\n");
     
   }
   

   return;
#endif /* DEBUG3 */   
   

   

   r1 = r;
   num_c = num_c1;
   num_s = num_s1;
   
   

   nr1 = nr;
   n_tot1 = n_tot;
   


   for (c = 0; c < num_c; c++)
   {
      r_max = 0.0;
      for (s = 0; s < num_s; s++)
	if(good_train_stim[s])
      {
         nt = nr[s];
         for (t = 0; t < nt; t++)
            if (r[c][s][t] > r_max)
               r_max = r[c][s][t];
      }
      mexPrintf("max %5.1f spikes for cell %d\n",r_max,c);
      ra0[c] = 0.0;
      sd0[c] = 0.0;
      for (s = 0; s < num_s; s++)
	if(good_train_stim[s])
      {
	 nt = nr[s];
	 ra[c][s] = 0.;
	 sd[c][s] = 0.;
	 f0[c][s] = 0.;
	 for (t = 0; t < nt; t++)
	 {
	    ra[c][s] += r[c][s][t] / nt;
	    sd[c][s] += r[c][s][t] * r[c][s][t] / (nt - 1.);
	    if(r[c][s][t] < eps)
	       f0[c][s] += 1. / nt;
	 }
	 sd[c][s] -= ra[c][s] * ra[c][s] * nt / (nt - 1.);
	 ra0[c] += ra[c][s] * nt / n_tot;
	 sd0[c] += sd[c][s] * nt / n_tot;	
         acs[c][s] = 1.;
         if (ra[c][s] <= 0.2 / nt)
	    rcs[c][s] = 0.2 / nt;
         else if (sd[c][s] <= ra[c][s])
	    rcs[c][s] = ra[c][s];
	 else
         {
	    rcs[c][s] = sd[c][s] / ra[c][s] - 1. + ra[c][s];
	    acs[c][s] = ra[c][s] / rcs[c][s];
         }
      }
   }

 x_prev = 0.;
   pcc0 = 0.;
   for (s = 0; s < num_s; s++)
     if(good_train_stim[s])
   {
      rt = nr[s] - 1;                 /* number of training trials */
      for (s1 = 0; s1 < num_s; s1++)
	 Pmix[s][s1] = 0.0;           /* this is for the metric analysis */
      Ps[s] = (float) nr[s] / (float) n_tot;
      pcc0 += Ps[s] * Ps[s];
      mexPrintf("Ps[s]=%g \n", Ps[s]);
      x_s[s] = x_prev + Ps[s];
      x_prev = x_s[s];

   }
   mexPrintf("Assigned probs for each stimulus, pcc0=%4.2f\n", pcc0 * 100.0);

   dp = 1. / (float)(n_tot);
 





/*    Read_Rfile_seq(rfile2_name, &rdec, &stim); */
   stim = stimuli2;
   rdec = responses2;
   

   num_s2 = num_s;
   num_c2 = num_c;
   

   dec_hist = (double **)mxCalloc(num_s2, sizeof(double *));
   for(s=0;s<num_s2;s++)
     dec_hist[s] = (double *)mxCalloc(num_s1, sizeof(double));
   



   /* for the moment only allow nc to be num_c and check that training
   ** and decoding have the same number of cells */

   if(num_c1 != num_c2)
     {
       mexErrMsgTxt("I have a problem: training and decoding set must have the same number of cells.\nExiting.\n");
     
     }
   
   nc_min = nc_max = num_c1;
   


 
   for (nc = nc_min; nc <= nc_max; nc++)
   {

/*    pick niter different nc-plets */
/*       niter = niter_fact * (num_c - nc + 1);  */
/*       for (iter = 0; iter < niter; iter++) */
/*       { */
/*       fprintf(outfile, "%d th n-plet.\n", iter+1); */
/* 	 icc = iter; */
/* 	 while (icc >= num_c) */
/* 	    icc -= num_c; */
/* 	 ic[0] = icc; */
/* 	 for (ic1 = 1; ic1 < nc; ic1++) */
/* 	 { */
/* 	    ok = False; */
/* 	    while (!ok) */
/* 	    { */
/* 	       c = drand48() * num_c; */
/* 	       used = False; */
/* 	       ic2 = 0; */
/* 	       while (!used && ic2 < ic1) */
/* 	       { */
/* 		  if (c == ic[ic2]) */
/* 		  used = True; */
/* 		  ic2++; */
/* 	       } */
/* 	       if (!used) */
/* 	       { */
/* 		  ok = True; */
/* 		  for (s = 0; s < num_s; s++) */
/* 		     if (nr[s] <= 1) */
/* 			ok = False; */
/* 	       } */
/* 	    } */
/* 	    ic[ic1] = c; */
/* 	 } */

     for (ic1 = 1; ic1 < nc; ic1++)
       ic[ic1] = ic1;
     

     p_thresh = p_thresh_fact / num_s2;
     
     mexPrintf("fact = %g, thresh = %g, num_trials_dec = %d\n",p_thresh_fact,
	     p_thresh, num_trials_dec);
     


     for(t=0; t < num_trials_dec;t++)
       {
	 int dec_stim;
	 
	 for(c=0;c<num_c;c++)
	   rc[c] = rdec[c][t];
	 
	 num_s = num_s1;
	 Decode();

	 /*	 fprintf(stderr, "p_max:   %g\n",p_max);    DEBUG */
	 
	 if(p_max < p_thresh)
	   dec_stim = -1;
	 else
	   {
	     dec_stim = best_s;
	     dec_hist[dec_stim][(int)stim[t]]++;
	   }
	 
	 DecStim[t] = dec_stim;
	 ProbMax[t] = p_max;
	 
	 
	/*  fprintf(outfile,"%d\t%d\t%d\t%g\n",t, (int)stim[t], dec_stim, */
/* 		 p_max); */
       }
     
     /* normalize the histogram and write it into the file */
     {
       int i,j;
       double tot1;
       for(i=0;i<num_s2;i++)
	 {
	   tot1 = 0.;
	   for(j=0;j<num_s1;j++)
	     tot1 += dec_hist[i][j];
	   if(tot1 == 0.)
	     for(j=0;j<num_s1;j++)
	       dec_hist[i][j] = -1;
	   else if(!nonorm)
	     for(j=0;j<num_s1;j++) 
	       dec_hist[i][j] /= tot1;
	 }

 /*       for(i=0;i<num_s2;i++) */
/* 	 { */
/* 	   for(j=0;j<num_s1;j++) */
/* 	     fprintf(histfile,"%g\t",dec_hist[i][j]); */
/* 	   fprintf(histfile,"\n"); */
/* 	 } */
/*        fclose(histfile); */
       


     }
     
   }
   

   plhs[0] = mxCreateDoubleMatrix(1, num_trials_dec, mxREAL);
   ptr = mxGetPr(plhs[0]);
   memcpy(ptr, DecStim, num_trials_dec * sizeof(double));

   plhs[1] = mxCreateDoubleMatrix(1, num_trials_dec, mxREAL);
   ptr = mxGetPr(plhs[1]);
   memcpy(ptr, ProbMax, num_trials_dec * sizeof(double));

   plhs[2] = mxCreateDoubleMatrix(num_s1, num_s2, mxREAL);
   ptr = mxGetPr(plhs[2]);
   {
     int i, j;
     mexPrintf("num_s1 = %d num_s2 = %d\n", num_s1, num_s2);
     
     for(i = 0;i < num_s1; i++)
       for(j = 0; j < num_s2; j++)
	 ptr[j+num_s1*i] = dec_hist[j][i];
   }
     
       

   mexPrintf("\nnum_s1 = %d num_s2 = %d\n\n", num_s1, num_s2);
     
     /*    fclose(outfile); */
   
   
   

   return;
   
   
    
  
}



void    Decode(void)
{
   double  p_tot,  ee,  eps, leps, noise_fact, uppe, fact;
   int     ic1, c, s, b;
   double  x, y, xx, yy, xy, sd2x, sd2y, sdxsdy;
   double  xy_corr, xy_corr_sq, xy_corr_ave, xy_corr_max, root;

   eps = 0.000001;
   leps = 2.*log(1./eps);
   best_s = 0;
   p_tot = 0.;
   p_max = 0.;

   if (decoding_method == 1)	        /* dec 1, i.e. approximate Bayesian */
   {                                    /* likelihoods with Gaussian fits   */
      for (s = 0; s < num_s; s++)       /* for each s, estimate
      ** P(r|s)P(s)  */
	if(good_train_stim[s])
      {
	 pq_r[s] = Ps[s];	        /* start by writing P(s)            */
	 for (ic1 = 0; ic1 < nc; ic1++)	/* and multiply by the product over */
	 {			        /* cells of P(rc[c]|s)              */
	    c = ic[ic1];
	    if (rc[c] < eps)	        /* if cell c does not fire now then */
	    {
	       if (f0[c][s] > eps)	/* if there were trials in which it */
		  pq_r[s] *= f0[c][s];	/* also did not fire, then multiply */
	       else		        /* by that fraction; else, by the   */
	       {                        /* tail of the probab. integral     */
                  fact = phims(0., ra[c][s], sd[c][s]);
                  if(fact < eps) fact = eps;
                  pq_r[s] *= fact;           
               }
	    }
	    else
	    {               /* if instead c does fire in the current vector */
	       ee = rc[c] - ra[c][s];
	       ee = ee*ee;	/* if it fires at exactly the mean rate, then
				   just multiply by one (i.e., do nothing)    */
	       if (ee > eps)	/* if not at exactly the mean rate, then..... */
	       {
		  if (sd[c][s] < ee/leps) /* if there was no variance in */
		     pq_r[s] *= eps;	/* train trials, then kill the test */
		  else		   /* while if there was, do the full thing */
		  {
		     fact = exp(-ee / (sd[c][s] * 2.0));
		     pq_r[s] *= fact;
		  }
	       }
	    }
	 }
	 p_tot += pq_r[s];
	 noise_fact = 1. + eps * (drand48() - 0.5);
	 if (p_max < pq_r[s] * noise_fact)
	 {
	    p_max = pq_r[s] * noise_fact;
	    best_s = s;
	 }
      }
   }

   else if (decoding_method == 5)	/* approximate Bayesian likelihoods */
   {                                    /* with Poisson distributions       */
      for (s = 0; s < num_s; s++)	/* for each s, estimate P(r|s)P(s)  */
	if(good_train_stim[s])
	  {
	 pq_r[s] = Ps[s];	        /* start by writing P(s)            */
	 for (ic1 = 0; ic1 < nc; ic1++) /* and multiply by the product over */
	 {			        /* cells of P(rc[c]|s)              */
	    c = ic[ic1];
	    b = 0;
	    uppe = eps;
	    fact = acs[c][s] * exp(-rcs[c][s]);
	    if (rc[c] < uppe)	        /* if cell c does not fire now then.. */
	    {
	       fact += 1. - acs[c][s];
	    }
	    else
	    {
	       while (rc[c] > uppe)	/* if instead c does fire in */
	       {                        /* the current vector ...... */
		  b++;
		  uppe = eps + b;
		  fact *= rcs[c][s] / (double) b;
	       }
	    }
	    pq_r[s] *= fact;
	 }
	 p_tot += pq_r[s];
	 noise_fact = 1. + eps * (drand48() - 0.5);
	 if (p_max < pq_r[s] * noise_fact)
	 {
	    p_max = pq_r[s] * noise_fact;
	    best_s = s;
	 }
      }
   }

   else if (decoding_method == 6)	/* method 6, (cos) angle between */
   {                                    /* the vectors */
      xy_corr_max = 0.0;
      xy_corr_ave = 0.0;
      xy_corr_sq = 0.0;
      xx = 0.0;
      for (ic1 = 0; ic1 < nc; ic1++)
      {
	 c = ic[ic1];
	 xx += (rc[c] * rc[c]);              /* Sigma x_squared */
      }
      for (s = 0; s < num_s; s++)	/* for each s, estimate
      ** P(r|s)P(s)  */
	if(good_train_stim[s])
	  {
	 xy_corr = 0.0;
  	 yy = 0.0;
	 xy = 0.0;
	 for (ic1 = 0; ic1 < nc; ic1++)
	 {
	    c = ic[ic1];
	    yy += (ra[c][s] * ra[c][s]);	/* Sigma y_squared */
	    xy += (rc[c] * ra[c][s]);		/* Sigma xy */
	 }
         if (xx > 0.0)
         {
	    root = sqrt(xx * yy); /* length of the resultant */
	    if (root == 0.0)
	       xy_corr = 0.0;
	    else
	       xy_corr = (xy / root);	/* (cos) angle between the vectors */
         }
         else
	    xy_corr = -yy;                      /* length of y vectors */
	 pq_r[s] = xy_corr;
	 xy_corr_ave += xy_corr / num_s;
	 xy_corr_sq += xy_corr * xy_corr / num_s;
      }
      xy_corr_sq -= xy_corr_ave * xy_corr_ave;
      xy_corr_sq = sqrt(xy_corr_sq + eps*eps);
      xy_corr_ave += 1.0 * xy_corr_sq;
      for (s = 0; s < num_s; s++)
	if(good_train_stim[s])
	  {
	 if (pq_r[s] > xy_corr_ave)
	    pq_r[s] -= xy_corr_ave;
	 else
	    pq_r[s] = 0.0;
	 p_tot += pq_r[s];
	 noise_fact = 1. + eps * (drand48() - 0.5); 
	 if (p_max < pq_r[s] * noise_fact)
	 {
	    p_max = pq_r[s] * noise_fact;
	    best_s = s;
	 }
      }
   }

   else if (decoding_method == 7)	/* method 7, correlation between */ 
   {                                    /* the vectors */
      xy_corr_max = 0.0;
      xy_corr_ave = 0.0;
      xy_corr_sq = 0.0;
      x = 0.0;
      xx = 0.0;
      for (ic1 = 0; ic1 < nc; ic1++)
      {
	 c = ic[ic1];
	 x += rc[c];
	 xx += (rc[c] * rc[c]);	/* Sigma x_squared */
      }
      sd2x = xx - (x * x) / nc;
      if(sd2x < 0.0)sd2x = 0.0;

      for (s = 0; s < num_s; s++)
	if(good_train_stim[s])
      {
	 xy_corr = 0.0;

	 if (nc > 1 && sd2x > 0.0)     /* can only do corrln if >1 cells */
	 {   	                               /* and sd2x > 0.0 */
	    y = 0.0;
	    yy = 0.0;
	    xy = 0.0;
                                
	    for (ic1 = 0; ic1 < nc; ic1++)
	    {
	       c = ic[ic1];
	       y += ra[c][s];
	       yy += (ra[c][s] * ra[c][s]);	/* Sigma y_squared */
	       xy += (rc[c] * ra[c][s]);	/* Sigma xy */
	    }
	    sd2y = yy - (y * y) / nc;
	    if(sd2y < 0.0)sd2y = 0.0;
	    sdxsdy = xy - (x * y) / nc;
	    root = sd2x * sd2y;
	    if (root < eps*eps)
	       xy_corr = sdxsdy / eps;
	    else
	       xy_corr = sdxsdy / sqrt(root);	/* corrln between the vectors */
	 }
	 else               /* special for 1 cell or uniform or zero response*/
	 {			/* get min square of diff of firing rate */
	    for (ic1 = 0; ic1 < nc; ic1++)
	    {
	       c = ic[ic1];
	       xy_corr -= (rc[c] - ra[c][s])*(rc[c] - ra[c][s]);  /* diff xy */
	    }
	 }			/* end of special case for 1 cell */

	 pq_r[s] = xy_corr;
	 xy_corr_ave += xy_corr / num_s;
	 xy_corr_sq += xy_corr * xy_corr / num_s;
      }

      xy_corr_sq -= xy_corr_ave * xy_corr_ave;
      xy_corr_sq = sqrt(xy_corr_sq + eps*eps);
      xy_corr_ave += 1.0 * xy_corr_sq;
      for (s = 0; s < num_s; s++)
	if(good_train_stim[s])
      {
	 if (pq_r[s] > xy_corr_ave)
	    pq_r[s] -= xy_corr_ave;
	 else
	    pq_r[s] = 0.0;
	 p_tot += pq_r[s];
	 noise_fact = 1. + eps * (drand48() - 0.5); 
	 if (p_max < pq_r[s] * noise_fact)
	 {
	    p_max = pq_r[s] * noise_fact;
	    best_s = s;
	 }
      }
   }

   else if (decoding_method == 8)	/* method 8, exp - dist(rc,ra[s]) */
   {                                    
      xx = 0.0;
      st_dev = 0.0;
      for (ic1 = 0; ic1 < nc; ic1++)
      {
	 c = ic[ic1];
	 xx += (rc[c] * rc[c]);              /* Sigma x_squared */
	 st_dev += sd0[c];
      }
      st_dev = 0.1*sqrt(st_dev/(double)nc);
      for (s = 0; s < num_s; s++)	/* for each s, estimate P(r|s)P(s)  */
	if(good_train_stim[s])
      {
  	 yy = 0.0;
	 xy = 0.0;
	 for (ic1 = 0; ic1 < nc; ic1++)
	 {
	    c = ic[ic1];
	    yy += (ra[c][s] * ra[c][s]);	/* Sigma y_squared */
	    xy += (rc[c] * ra[c][s]);		/* Sigma xy */
	 }
	 root = xx + yy - 2.*xy;                /* square distance x to y */
         if (root <= 0.0)root = 0.0;
         root = sqrt(root)/st_dev;              /* normalized distance */
	 xy_corr = exp(-root);                  /* length of y vectors */
	 pq_r[s] = xy_corr;
      }
      for (s = 0; s < num_s; s++)
	if(good_train_stim[s])
      {
	 p_tot += pq_r[s];
	 noise_fact = 1. + eps * (drand48() - 0.5); 
	 if (p_max < pq_r[s] * noise_fact)
	 {
	    p_max = pq_r[s] * noise_fact;
	    best_s = s;
	 }
      }
   }

   else
   {
     mexErrMsgTxt("error - decoding method not found\n");
   }

   for (s = 0; s < num_s; s++)
     if(good_train_stim[s])
   {
      if (p_tot == 0.0)
      {
	 pq_r[s] = Ps[s];	/* if they all died, equal honor to all */
	 best_s = num_s * drand48();
      }
      else
	 pq_r[s] /= p_tot;	/* finally normalize to get P(s|r)  */
   }
   if(p_tot != 0.)
     p_max /= p_tot;
   
  
}




void    Read_Input(int num_s, int num_c)
{
   int     ix, c, s, ss, n, s_start;
   float   rate;


   
   

   mexPrintf("number of cells=%4d   number of stimuli=%4d\n", num_c, num_s);

   nr = (int *) mxCalloc(num_s, sizeof(int));

   r = (double ***) mxCalloc(num_c, sizeof(double **));
   good_train_stim = (int *)mxCalloc(num_s, sizeof(int)); 
   for (c = 0; c < num_c; c++)
   {
      r[c] = (double **) mxCalloc(num_s, sizeof(double *));

      for (s = 0; s < num_s; s++)
      {
	 r[c][s] = (double *) mxCalloc(max_t, sizeof(double ));

	 if (!r[c][s])
	 {
	   char ErrMsg[80];
	   
	   sprintf(ErrMsg, "Unable to malloc at c = %d, s = %d.\n", c, s);
	   mexErrMsgTxt(ErrMsg);
	   

	 }
	
      }

   }

   for (s = 0; s < num_s; s++)
      nr[s] = 0;


     for(ix = 0; ix < num_trials; ix++)
      {
	s = stimuli[ix];
	
	 for (c = 0; c < num_c; c++)
	 {
	    if (nr[s] < max_t)
	      {
		
		r[c][s][nr[s]] = responses[ix][c];
		if(r[c][s][nr[s]] > 70 )
		  mexPrintf("ix = %d r[%d][%d][%d] = %g\n", ix, c, s, nr[s], r[c][s][nr[s]]);
	      }
	    
		 
	 }
	 ++nr[s];
      }

      for(ix = 0; ix < num_trials; ix++) 
        mxFree(responses[ix]); 
      mxFree(responses); 
     
     



     




   for (s = 0; s < num_s; s++)
      if (nr[s] >= max_t)
	 nr[s] = max_t;

   s_start = 0;
   if (s0err)
      s_start = 1;
   n_tot = 0;
   ss = 0;

   for (s = s_start; s < num_s; s++)
   {
      if (nr[s] >= min_t)
      {
	good_train_stim[s] = 1;
	
	 for (c = 0; c < num_c; c++)
	    r[c][ss] = r[c][s];
	 nr[ss] = nr[s];
	 n_tot += nr[s];
	 ++ss;
      }
      else
      {
	char WarnMsg[80];
	sprintf(WarnMsg, "not enough samples for s = %d.\n", s);
	mexWarnMsgTxt(WarnMsg);
	
      }
   }
/*    num_s = ss;*/

   while (s_throw > 0 && num_s > 0)
   {
     mexPrintf("Within s_throw\n");
     
      s_start = drand48() * num_s;
      mexPrintf( "Reducing stimulus set by throwing s = %d.\n", s_start);
      n_tot -= nr[s_start];
      ss = s_start;
      for (s = s_start + 1; s < num_s; s++)
      {
	 for (c = 0; c < num_c; c++)
	    r[c][ss] = r[c][s];
	 nr[ss] = nr[s];
	 ++ss;
      }
      num_s -= 1;
      s_throw -= 1;
   }

   mexPrintf(" Final num_s=%d num_c=%d trials=%d\n", num_s, num_c, n_tot);

}






/* reads the rfile and generates a amtrix containing sequentially the
** firing of cells and a vectror with the stimuli */






double  phims(val, mean, sigma)
double  val, mean, sigma;
{
   double  h, z, t, erfcc, phim;

   h = val - mean;
   if (h == 0.)
      return (0.5);
   else if (h * h > 25.0 * sigma)
   {
      if (h > 0.)
      {
	 return (1.0);
      }
      else
	 return (0.0);
   }
   else
   {
      h /= sqrt(sigma);
      z = sqrt(0.5 * h * h);
      t = 1.0 / (1.0 + 0.5 * z);
      erfcc = t * exp(-z * z - 1.26551223 +
	      t * (1.00002368 + t * (0.37409196 +
	      t * (0.09678418 + t * (-0.18628806 +
	      t * (0.27886807 + t * (-1.13520398 +
	      t * (1.48851587 + t * (-0.82215223 +
	      t * 0.17087277)))))))));
      if (h < 0.)
	 phim = 0.5 * erfcc;
      else
	 phim = 1.0 - 0.5 * erfcc;
      return (phim);
   }
}




