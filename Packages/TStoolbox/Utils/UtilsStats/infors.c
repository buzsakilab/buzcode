/* infors_mex.c */
/* version modified to a me file for use in MATLAB */
/* batta 1999 */


/* infors.c
   used to extract the amount of information, in bits, rate responses of
   individual cells carry about a set of stimuli, and sparseness measures  
   developed by Alessandro Treves, ale@limbo.sissa.it, with
   Edmund Rolls, Stefano Panzeri, Bill Skaggs, during 1994-95  */

/* 16 Aug. 96: small tuning of the finite size corrections (gg parameter)
   18 Apr. 96: calculate also the stimulus specific information (S.Panzeri)
   17 Apr. 96: prints out the mean rates to each stimulus (S.P.)
   4 Dec 1995: revision in the f.s. subtraction (S.Panzeri)
   10/04/95: version for use with Tucson style data
   30/08/95: reads ascii infom.dat if run with the -ascii option  */

/* usage         OPTION               MEANING              DEFAULT 


   infors


   /* these are the option you have to take care of */
/*    -seed 989765     : change random seed                   243342 */

/*    -max_t 10        : max # trials per stimulus allowed    2*max_s */
/*    -min_t 10        : min #   "     "     "        "       2*max_s */
/*    -s0err           : throw trials with stimulus = 0       False */
/*    -s_throw 4       : throw eg 4 stimuli at random, to */
/*                       try with a reduced subset            0  */

/*    -max_b           : max number of bins                   25  */
/*    -timewin 500     : time window size in msec.            1000 */
 
/******* warning: if you do not choose the command line argument -timewin
the mean rates are printed in spikes, otherwise in spikes/sec     **********/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mex.h"



#define True    1
#define False   0
double  drand48();

#define max_c   200		/* max no of cells */
#define max_s   200		/* max no of stimuli */
#define MAX_B   100

int cmdlineOutput = 0;

void    Read_Input(void);
double  phims(double val, double mean, double sigma);

double  **rcs;	/* prob distr mean rates */
double  **acs;	/* prob distr nonzero rates */
double  **ra;	/* rate averages */
double  **ra_true;	/* rate averages before rescaling */
double  **sd;	/* rate variances */
double  **f0;	/* fraction of zero rates */

double ***r;			/* responses (total number of spikes) */

double *stimuli, **responses;   /* the input array copied from the */
				/*   mxArrays */
double *II_raw, *II_C1, *II_corr;


int     num_s;
int     num_c;
int     num_trials;
int     max_t;
int     min_t;
int     max_b = 25;
int    *nr;			/* number of trials */
int     n_tot;			/* total of above */
double  *Ps;	        	/* Prob(s)    */
int     s0err = 0;
int     read_ASCII = 0;
int     s_throw = 0;		/* stimuli to be discarded */
double  rmax_b;
int     FREQ = False;
double  timewin = 1000.;
double  *sp, *lsp;	/* sparsity and its log */
double  *It, It_av;	/* info derivative */
int     seed;
int     infors()
{

   int     cur_arg;
   double  r_max;		/* max response of a cell */
   double  *r_a1, *r_a2;	/* averages for sparseness */
   double  *r_a1_true;	/* before rescaling  */
   double  *r_a1_true_max;
   double  *r_a1_true_min;
   double  **Qsq, q1;	/* quantized probability table */
   double  *Qq;		/* quantized Prob(s_q) */
   double  **Psq, p1;	/* probability table */
   double  *Pq;		/* Prob(s_q) */
   int     ***nrb;	/* number of rates per bin */
   int     **nrb0;	/* total number of rates per bin */
   double  I_raw;		/* information values */
   double  *I_raws;
   double  I_C1;		/* first correction term from data */
   double  I_avg_aft;
   double  *I_C1s;
   double  I_mean_raw;		/* average info across cells */
   double  I_mean_C1;		/* average 1st correction term */
   double  *I_means;
   double  dp;			/* probability of one trial */
   int     c, s, t, b;		/* cell,stimulus,trial,bin */

   int     nt;
   double  rt;
   int     nb, xtr, xxx;	/* number of relevant bins */
   double  nb_x;		/* expected number of occupied bins */
   double  qc_x;		/* expected conditional probability */
   double  delta_N, delta_N_prev;	/* (bins expected - bins found) */
   double  gg;			/* gamma factor for Stefano's routine */
   int     q;
   double  eps = 0.000001;	/* small number */
   double  uppe, fact;

#define Match(arg)	(strcmp(argv[cur_arg], (arg)) == 0)

/*                    ******************************
                      process command line arguments
                      ******************************         */




   seed = 243342;

 

/*                  *************************************
                    read data and set constant parameters
                    *************************************    */

   rmax_b = ((double) max_b - 1.) - eps;

   
   srand48(seed);
   Read_Input();

   
   /* dynamically allocate memory for everything */
   {
     int i, j, k;
     
     if(!(r_a1 = (double *)mxCalloc(num_c, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     
     if(!(r_a2 = (double *)mxCalloc(num_c, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(r_a1_true = (double *)mxCalloc(num_c, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(r_a1_true_max = (double *)mxCalloc(num_c, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(r_a1_true_min = (double *)mxCalloc(num_c, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(Qsq = (double **)mxCalloc(num_s, sizeof(double *))))
       mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_s; i++)
       if(!(Qsq[i] = (double *)mxCalloc(max_b, sizeof(double))))
	  mexErrMsgTxt("Malloc error 1");
     if(!(Qq = (double *)mxCalloc(max_b, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(Psq = (double **)mxCalloc(num_s, sizeof(double *))))
       mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_s; i++)
       if(!(Psq[i] = (double *)mxCalloc(max_b, sizeof(double))))
	 mexErrMsgTxt("Malloc error 1");
     if(!(Pq = (double *)mxCalloc(max_b, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");

     if(!(nrb = (int ***)mxCalloc(num_c, sizeof(int **))))
	mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_c; i++)
       {
	 if(!(nrb[i] = (int **)mxCalloc(num_s, sizeof(int *))))
	   mexErrMsgTxt("Malloc error 1");
	 for(j = 0; j < num_s; j++)
	   if(!(nrb[i][j] = (int *)mxCalloc(max_b, sizeof(int))))
	     mexErrMsgTxt("Malloc error 1");
       }
     
     if(!(nrb0 = (int **)mxCalloc(num_c, sizeof(int *))))
       mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_c; i++)
       if(!(nrb0[i] = (int *)mxCalloc(max_b, sizeof(int))))
	 mexErrMsgTxt("Malloc error 1");
     

     if(!(I_raws = (double *)mxCalloc(num_s, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     
     if(!(I_C1s = (double *)mxCalloc(num_s, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(I_means = (double *)mxCalloc(num_s, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(rcs = (double **)mxCalloc(num_c, sizeof(double *))))
       mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_c; i++)
       if(!(rcs[i] = (double*)mxCalloc(num_s, sizeof(double))))
	 mexErrMsgTxt("Malloc error 1");

     if(!(acs = (double **)mxCalloc(num_c, sizeof(double *))))
       mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_c; i++)
       if(!(acs[i] = (double*)mxCalloc(num_s, sizeof(double))))
	 mexErrMsgTxt("Malloc error 1");

     if(!(ra = (double **)mxCalloc(num_c, sizeof(double *))))
       mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_c; i++)
       if(!(ra[i] = (double*)mxCalloc(num_s, sizeof(double))))
	 mexErrMsgTxt("Malloc error 1");

     if(!(ra_true = (double **)mxCalloc(num_c, sizeof(double *))))
       mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_c; i++)
       if(!(ra_true[i] = (double*)mxCalloc(num_s, sizeof(double))))
	 mexErrMsgTxt("Malloc error 1");

     if(!(sd = (double **)mxCalloc(num_c, sizeof(double *))))
       mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_c; i++)
       if(!(sd[i] = (double*)mxCalloc(num_s, sizeof(double))))
	 mexErrMsgTxt("Malloc error 1");

     if(!(f0 = (double **)mxCalloc(num_c, sizeof(double *))))
       mexErrMsgTxt("Malloc error 1");
     for(i = 0; i < num_c; i++)
       if(!(f0[i] = (double*)mxCalloc(num_s, sizeof(double))))
	 mexErrMsgTxt("Malloc error 1");

     if(!(Ps = (double *)mxCalloc(num_s, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(sp = (double *)mxCalloc(num_c, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(lsp = (double *)mxCalloc(num_c, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");
     if(!(It = (double *)mxCalloc(num_c, sizeof(double))))
       mexErrMsgTxt("Malloc error 1");

     
   }
   
/******  calculate true firing rates  (before rescaling)    *******/
   for (c = 0; c < num_c; c++)
   {
      r_a1_true[c] = 0.0;
      r_a1_true_max[c] = 0.0;
      r_a1_true_min[c] = 0.0;
      for (s = 0; s < num_s; s++)
      {
	 nt = nr[s];
	 ra_true[c][s] = 0.;
	 for (t = 0; t < nt; t++)
	 {
	    ra_true[c][s] += r[c][s][t] / nt;
	 }
	 if(ra_true[c][s] > r_a1_true_max[c] || r_a1_true_max[c] == 0)
	   r_a1_true_max[c] = ra_true[c][s];
	 if(ra_true[c][s] < r_a1_true_min[c] || r_a1_true_min[c] == 0)
	   r_a1_true_min[c] = ra_true[c][s];
	 
	 r_a1_true[c] += ra_true[c][s] * (double) nt / (double) n_tot;
      }
   }

   for (c = 0; c < num_c; c++)
   {
      for (s = 0; s < num_s; s++)
      {
	 ra_true[c][s] /= timewin / 1000.;	/* now in spikes/sec. */
      }
      r_a1_true[c] /= timewin / 1000.;	/* now in spikes/sec   */

      
   }



   It_av = 0.0;
   for (c = 0; c < num_c; c++)
   {

     
      r_max = 0.0;
      It[c] = 0.0;
      r_a1[c] = 0.0;
      r_a2[c] = 0.0;
      for (s = 0; s < num_s; s++)
      {
	 nt = nr[s];
	 for (t = 0; t < nt; t++)
	    if (r[c][s][t] > r_max)
	       r_max = r[c][s][t];
      }
      if (r_max > rmax_b)
	{
	  char warn_msg[80];

	  sprintf(warn_msg, "%d bins exceeded for cell %d\n", max_b,
		  c);
	  mexWarnMsgTxt(warn_msg);
	}
      

      

      for (s = 0; s < num_s; s++)
      {

	
	 nt = nr[s];
	 ra[c][s] = 0.;
	 sd[c][s] = 0.;
	 f0[c][s] = 0.;
	 for (b = 0; b < max_b; b++)
	    nrb[c][s][b] = 0;

	 
	 for (t = 0; t < nt; t++)
	 {
	   
	    if (r_max > rmax_b)
	       r[c][s][t] *= rmax_b / r_max;

	    
	    b = 0;
	    uppe = eps;
	    while (r[c][s][t] > uppe && b < max_b-2)
	    {
	       b++;
	       uppe = (double) b + eps;
	    }

	    
	    nrb[c][s][b]++;
	    nrb0[c][b]++;

	    
	    ra[c][s] += r[c][s][t] / nt;
	    sd[c][s] += r[c][s][t] * r[c][s][t] / (nt - 1.);
	    if (r[c][s][t] < eps)
	       f0[c][s] += 1. / nt;
	 }

	 sd[c][s] -= ra[c][s] * ra[c][s] * nt / (nt - 1.);	/* back to sd (at) */
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
	 r_a1[c] += ra[c][s] * (double) nt / (double) n_tot;
	 r_a2[c] += ra[c][s] * ra[c][s] * (double) nt / (double) n_tot;
	 if (ra[c][s] > eps)
	    It[c] += log(ra_true[c][s]) * ra_true[c][s] * (double) nt / (double) n_tot;

	 fact = acs[c][s] * exp(-rcs[c][s]);


	 if (f0[c][s] > eps)
	    fact = f0[c][s];
	 else
	    fact = phims(0., ra[c][s], sd[c][s]);

	 uppe = 0.0;
	 for (b = 1; b < max_b; b++)
	 {
	    fact = -phims(uppe, ra[c][s], sd[c][s]);
	    uppe += 1.;
	    fact += phims(uppe, ra[c][s], sd[c][s]);

	 }


      }
      if (r_a1[c] > eps)
	 It[c] -= log(r_a1_true[c]) * r_a1_true[c];
      It[c] /= log(2.);
      It_av += It[c] / num_c;
   }

   for (s = 0; s < num_s; s++)
   {
      rt = nr[s];
      Ps[s] = rt / (double) n_tot;

      
   }


   I_mean_raw = 0.;
   I_mean_C1 = 0.;
   dp = 1. / (double) n_tot;
   for (s = 0; s < num_s; s++)
      I_means[s] = 0.0;

   


/*                  ***********************************
                    loop calculating info for each cell
                    ***********************************    */
   for (c = 0; c < num_c; c++)
   {
      I_raw = 0.;
      I_C1 = 0.;
      for (s = 0; s < num_s; s++)
      {
	 I_raws[s] = 0.;
	 I_C1s[s] = 0.;
	 fact = acs[c][s] * exp(-rcs[c][s]);
	 Psq[s][0] = Ps[s] * (fact + 1. - acs[c][s]);
	 Qsq[s][0] = (double) nrb[c][s][0] / (double) n_tot;
	 for (b = 1; b < max_b; b++)
	 {
	    fact *= rcs[c][s] / b;
	    Psq[s][b] = Ps[s] * fact;
	    Qsq[s][b] = (double) nrb[c][s][b] / (double) n_tot;
	 }
      }
      for (b = 0; b < max_b; b++)
      {
	 Pq[b] = 0.;
	 Qq[b] = 0.;
	 for (s = 0; s < num_s; s++)
	 {
	    Pq[b] += Psq[s][b];
	    Qq[b] += Qsq[s][b];
	 }

	 
      }

      
      sp[c] = 0.0;
      lsp[c] = 0.0;
      if (r_a2[c] > 0.0)
      {
	 sp[c] = r_a1[c] * r_a1[c] / r_a2[c];
	 lsp[c] = -sp[c] * log(sp[c]);
      }

      for (s = 0; s < num_s; s++)
      {
	 nt = nr[s];
	 nb = 0;
	 for (b = 0; b < max_b; b++)
	 {
	    q1 = Qsq[s][b];
	    if (q1 > eps)
	    {
	       I_raw += q1 * log(q1 / (Qq[b] * Ps[s]));
	       I_raws[s] += (q1 / Ps[s]) * log(q1 / (Qq[b] * Ps[s]));
	       nb += 1.;
	    }
	 }
	 if (nb < max_b)
	 {
	    nb_x = 0.0;
	    for (b = 0; b < max_b; b++)
	    {
	       qc_x = ((Qsq[s][b] / Ps[s] - eps) * nt + 1.) / (nt + nb);
	       if (Qsq[s][b] > eps)
		  nb_x += 1. - exp(log(1. - qc_x) * nt);
	    }
	    delta_N_prev = max_b * max_b;
	    delta_N = (nb - nb_x) * (nb - nb_x);
	    xtr = 0;
	    while (delta_N < delta_N_prev && (nb + xtr) < max_b)
	    {
	       xtr++;
	       nb_x = 0.0;
	       gg = log(1. + .8*(double) nb / (double) nt) * xtr / nt;
	       xxx = 0;
	       for (b = 0; b < max_b; b++)
	       {
		  if (Qsq[s][b] > eps)
		  {
		     qc_x = (1. - gg) * ((Qsq[s][b] / Ps[s]) * nt + 1.) / (nt + nb);
		     nb_x += 1. - exp(log(1. - qc_x) * nt);
		  }
		  else if (xxx < xtr)
		  {
		     qc_x = gg / xtr;
		     nb_x += 1. - exp(log(1. - qc_x) * nt);
		     xxx++;
		  }
	       }
	       delta_N_prev = delta_N;
	       delta_N = (nb - nb_x) * (nb - nb_x);
	    }
	    nb += xtr - 1;
	    if (delta_N < delta_N_prev)
	       nb++;
	 }
	 if (Ps[s] > eps)
	    I_C1s[s] = (1. / Ps[s]) * (nb - 1.) - 1.;
	 I_C1 += nb;
	 for (b = 0; b < max_b; b++)
	 {
	    if ((Ps[s] > eps) && (Qsq[s][b] > eps))
	    {
	       q1 = Qsq[s][b] / Ps[s];
	       I_C1s[s] += (-q1 + 2 * q1 * q1) / (Qq[b]);
	    }
	 }
	 I_C1s[s] /= 2. * (double) n_tot *log(2.);

	 I_raws[s] /= log(2.);
      }
      for (b = 0; b < max_b; b++)
	 if (Qq[b] > eps)
	    I_C1 -= 1.;

      I_raw /= log(2.);
      I_mean_raw += I_raw / num_c;
      I_C1 -= num_s - 1.;
      I_C1 *= dp / (2. * log(2.));
      I_mean_C1 += I_C1 / num_c;

      for (s = 0; s < num_s; s++)
	 I_means[s] += (I_raws[s] - I_C1s[s]) / num_c;

      II_raw[c] = I_raw;
      II_C1[c] = I_C1;
      II_corr[c] = I_raw - I_C1;

      
      

      
      
   }
   
   
/*    { */
/*      double lambda = 0., den = 0.; */
/*      FILE *lambdafile; */
     
/*      for(c=0; c < num_c; c++) */
/*        { */
/* 	  if(sp[c] > 1.e-10) */ 
/* 	   lambda -=  r_a1_true[c] * log(sp[c]); */ 
/* 	 if(I_raw - I_C1 > 0.) */
/* 	   lambda += It[c]; */
	 
/* 	 den += (r_a1_true_max[c] - r_a1_true_min[c]); */
/*        } */
/*      lambda /= den; */
     
/*      lambdafile = fopen("lambda.dat", "a"); */
/*      fprintf(lambdafile, "%s\t%g\n", rfile_name, lambda); */
     
/*      fclose(lambdafile); */
     
  
/*    } */
   
	   

 

}

void    Read_Input(void)
{
   int     ix, c, s, ss, n, s_start;
   float   rate;

   if(cmdlineOutput)
     mexPrintf("number of cells=%4d   number of stimuli=%4d\n", num_c, num_s);

   nr = (int *) mxCalloc(num_s, sizeof(int));

   r = (double ***) mxCalloc(num_c, sizeof(double **));

   for (c = 0; c < num_c; c++)
   {
      r[c] = (double **) mxCalloc(num_s, sizeof(double *));

      for (s = 0; s < num_s; s++)
      {
	 r[c][s] = (double *) mxCalloc(max_t, sizeof(double));

	 if (!r[c][s])
	 {
	   char err_msg[80];
	   sprintf(err_msg, "Unable to malloc at c = %d, s = %d.\n", c, s);
	   mexErrMsgTxt(err_msg);
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

		if(cmdlineOutput)
		  if(r[c][s][nr[s]] > 70 || r[c][s][nr[s]] != floor(r[c][s][nr[s]]))
		    mexPrintf("ix = %d r[%d][%d][%d] = %g\n", ix, c, s, nr[s], r[c][s][nr[s]]);
	      }
	    
		 
	 }
	 ++nr[s];
      }

/*      for(ix = 0; ix < num_trials; ix++) */
/*        mxFree(responses[ix]); */
/*      mxFree(responses); */
     
     



     




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
	 for (c = 0; c < num_c; c++)
	    r[c][ss] = r[c][s];
	 nr[ss] = nr[s];
	 n_tot += nr[s];
	 ++ss;
      }
      else
      {
	 fprintf(stderr, "Warning: not enough samples for s = %d.\n", s);
      }
   }
   num_s = ss;

   while (s_throw > 0 && num_s > 0)
   {
      s_start = drand48() * num_s;
      fprintf(stderr, "Reducing stimulus set by throwing s = %d.\n", s_start);
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

   fprintf(stderr, " Final num_s=%d num_c=%d trials=%d\n", num_s, num_c, n_tot);

}

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
/*      fprintf(stderr,"h %3.1f z %3.1f t %6.3f erfcc %6.3f\n",h,z,t,erfcc); */
      if (h < 0.)
	 phim = 0.5 * erfcc;
      else
	 phim = 1.0 - 0.5 * erfcc;
/*      fprintf(stderr,"phim %6.3f\n",phim); */
      return (phim);
   }
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
		 *prhs[])
{
  int i, j;
  int mrows, ncols;
  double *ptr;
  mxArray *field_value;
  
  min_t = 5;
  max_t = 2 * max_s;
  

  




  if(nrhs != 3 && nrhs != 4)
    mexErrMsgTxt("Three (or four) inputs required");
  
  if(nlhs < 1)
    mexErrMsgTxt("At least one output required");
  
  if(nlhs > 3) 
    mexErrMsgTxt("Too many output arguments");
  
  /* input 0 is the stimuli vector  */
  
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  if(!(mrows == 1 || ncols == 1))
    mexErrMsgTxt("Input 1 must be a row or column vector");
  if(mrows == 1)
    num_trials = ncols;
  else
    num_trials = mrows;
  
  stimuli = mxGetPr(prhs[0]);
  


  

  /* input 1 is the responses matrix */

  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  
  if(mrows != num_trials)
    mexErrMsgTxt("Input 1 and 2 must have the same number of rows.");
  
  num_c = ncols;
  

  ptr = mxGetPr(prhs[1]);

  
  

  responses = (double **)mxCalloc(num_trials, sizeof(double *));
  for(i = 0; i < num_trials; i++)
    responses[i] = mxCalloc(num_c, sizeof(double));
  
  
  for(i = 0; i < num_trials; i++)
    for(j = 0; j < num_c; j++)
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
  
  num_s = ptr[0];
  
  II_raw = (double *)mxCalloc(num_c, sizeof(double));
  II_C1   = (double *)mxCalloc(num_c, sizeof(double));
  II_corr = (double *)mxCalloc(num_c, sizeof(double));


  /* input 3 is a structure containing the options passed to infors,
     if a field contains the empty array the value is not changed */

  if(nrhs > 3)
    {
      if(!mxIsStruct(prhs[3]))
	mexErrMsgTxt("Input 4 must be a structure containing the options passed to informs");
      
      field_value = mxGetField(prhs[3], 0, "seed");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  seed = (int)ptr[0];
	}
      field_value = mxGetField(prhs[3], 0, "max_t");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  max_t = (int)ptr[0];
	}
      field_value = mxGetField(prhs[3], 0, "min_t");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  min_t = (int)ptr[0];
	  
	}
      field_value = mxGetField(prhs[3], 0, "s0err");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  s0err = (int)ptr[0];
	}
      field_value = mxGetField(prhs[3], 0, "s_throw");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  s_throw = (int)ptr[0];
	}
      field_value = mxGetField(prhs[3], 0, "timewin");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  timewin = ptr[0];
	}
      field_value = mxGetField(prhs[3], 0, "max_b");
      if(field_value != NULL)
	{
	  ptr = mxGetPr(field_value);
	  max_b = ptr[0];
	}



    }
  
  

  infors();
  
  plhs[0] = mxCreateDoubleMatrix(1, num_c, mxREAL);
  
  ptr = mxGetPr(plhs[0]);
  memcpy(ptr, II_corr, num_c * sizeof(double));
  
  if(nlhs > 1)
    {
      plhs[1]  = mxCreateDoubleMatrix(1, num_c, mxREAL);
      ptr = mxGetPr(plhs[1]);
      memcpy(ptr, It, num_c * sizeof(double));
    }
  
  if(nlhs > 2)
    {
      plhs[2]  = mxCreateDoubleMatrix(1, num_c, mxREAL);
      ptr = mxGetPr(plhs[2]);
      memcpy(ptr, sp, num_c * sizeof(double));
    }
}



    
  
