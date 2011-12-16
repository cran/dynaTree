/****************************************************************************
 *
 * Dynamic Trees for Learning and Design
 * Copyright (C) 2010, Universities of Cambridge and Chicago
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
 *
 ****************************************************************************/


extern "C" {
#include "rhelp.h"
#include "R.h"
#include "matrix.h"
}
#include "cloud.h"
#include "assert.h"

extern "C" {

unsigned int NC = 0;
Cloud **clouds = NULL;

/* useful prototypes */
unsigned int get_cloud(void);
  int ** alloc_XNA(unsigned int T, unsigned int m, double **X, 
		   int *Xna_in, int *XNA_in, unsigned int *nna);

/*
 * dynaTree_R:
 *
 * main C routine for creating a particle cloud of
 * dynamic trees based on T (x,y) pairs and N 
 * particles;  An integer reference to the established
 * cloud and the T-2*minp marginal likelihood calculations
 * are passed back to R
 */

void dynaTree_R(/* inputs */
		int* m_in,
		int* T_in,
		int* N_in,
		double* X_in,
		int *Xna_in,
		int *XNA_in,
		double* y_in,
		int* model_in,
		double* params_in,
		int *nstart_in,
		int *verb_in,
		
		/* outputs */
		double *lpred_out,
		int *c_out) 
{ 
  GetRNGstate();

  /* training data */
  unsigned int m = (unsigned int) *m_in;
  unsigned int T = (unsigned int) *T_in;
  unsigned int N = (unsigned int) *N_in;
  double **X = new_matrix_bones(X_in, T, m);
  double *y = new_dup_vector(y_in, T);
  
  /* any missing data (?) */
  unsigned int nna = 0;
  int ** XNA = alloc_XNA(T, m, X, Xna_in, XNA_in, &nna);

  /* make new Pall */
  Pall *pall = new_pall(X, T, m, Xna_in, XNA, nna, 
			y, params_in, *model_in);

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* choose starting indices and bounding rectangle */
  unsigned int nstart = *nstart_in; 
  assert(nstart >= (unsigned int) params_in[4]);
  if(nstart >= T) nstart = T-1;
  int *pstart = iseq(0,nstart-1);
  
  /* get a new cloud index */
  unsigned int c = get_cloud();
  Cloud *cloud;

  /* allocate a new particle cloud */
  cloud = clouds[c] = new Cloud(N, pall, pstart, nstart);
  free(pstart);
  
  /* initialize marginal likelihood */
  zerov(lpred_out, T);
  
  /* Particle Learning steps for each new data point */
  for(unsigned int t=nstart; t<T; t++){

    /* PL Resample step, and gathering of cloud summaries */
    double logpost = cloud->Resample(t, verb);

    /* trace marginal likelihood */
    lpred_out[t] = logpost;
  
    /* PL propagate step */
    cloud->Propagate(t);
  } 

  /* free data */
  free(X);
  if(XNA) free(XNA);
  free(y);

  /* return to R */
  *c_out = (int) c;
  PutRNGstate();
}


/*
 * rejuvenate_R:
 *
 * create a new particle cloud from an old pall object,
 * obtained from a different cloud -- and then combine
 * the particle clouds
 */

void rejuvenate_R(/* inputs */
		int* c_in,
		int *o_in,
		int *n_in,
		int *verb_in,
		
		/* outputs */
		double *lpred_out)
{ 
  GetRNGstate();

  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];

  /* sanity check */
  unsigned int T = cloud->pall->n;
  if(*n_in > 0) assert(*n_in == (int) T);
  else assert(!o_in);

  /* reorder the indices in the cloud */
  if(o_in) cloud->Reorder(o_in);

  /* choose starting indices and bounding rectangle */
  unsigned int nstart = cloud->pall->minp; // no possible splits until 2*minp
  if(nstart >= T) nstart = T-1;
  int *pstart = iseq(0,nstart-1);
  
  /* print something? */
  if(*verb_in > 0) myprintf(stdout, "rejuvenating\n");

  /* allocate a new particle cloud */
  Cloud *newcloud = new Cloud(cloud->Nrevert, cloud->pall, pstart, nstart);
  free(pstart);
  
  /* initialize marginal likelihood */
  zerov(lpred_out, T);
  
  /* Particle Learning steps for each new data point */
  for(unsigned int t=nstart; t<T; t++){

    /* PL Resample step, and gathering of cloud summaries */
    double logpost = newcloud->Resample(t, *verb_in);

    /* trace marginal likelihood */
    lpred_out[t] = logpost;
  
    /* PL propagate step */
    newcloud->Propagate(t);
  } 

  /* combine particle clouds; also deletes newcloud */
  cloud->Combine(newcloud);
  newcloud = NULL;

  /* print something */
  if(*verb_in > 0) myprintf(stdout, "now %d particles\n", cloud->N);

  /* return to R */
  PutRNGstate();
}


/*
 * update_R:
 *
 * add (x,y) pairs to the data and update the particle
 * approximation (in a particular cloud) via the 
 * particle learning algorithm 
 */

void update_R(/* inputs */
	      int* c_in,
	      int* m_in,
	      int* T_in,
	      double* X_in,
	      int *Xna_in,
	      int *XNA_in,
	      double* y_in,
	      int *verb_in,
	      
	      /* outputs */
	      double *lpred_out)
{ 
  GetRNGstate();

  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;
  assert((int) m == *m_in);

  /* dimensions of new data */
  unsigned int T = (unsigned int) *T_in;
  double **X = new_matrix_bones(X_in, T, m);
  double *y = new_dup_vector(y_in, T);

  /* any missing data (?) */
  unsigned int nna = 0;
  int ** XNA = alloc_XNA(T, m, X, Xna_in, XNA_in, &nna);

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* add the data and get start/end times */
  unsigned int nstart = cloud->pall->n;
  add_data(cloud->pall, X, T, Xna_in, XNA, nna, y);
  T = cloud->pall->n;
  
  /* Particle Learning steps for each new data point */
  for(unsigned int t=nstart; t<T; t++){

    /* PL Resample step, and gathering of cloud summaries */
    double logpost = cloud->Resample(t, verb);
  
    /* trace marginal likelihood */
    lpred_out[t-nstart] = logpost;

    /* PL propagate step */
    cloud->Propagate(t);
  } 

  /* free data */
  free(X);
  if(XNA) free(XNA);
  free(y);

  /* return to R */
  PutRNGstate();
}


/*
 * retire_R:
 *
 * R-interface to a routine that removes the input/output
 * pairs corresponding to a certain indices, moving the info
 * to the prior to the extent possible
 */

void retire_R(/* inputs */
	     int *c_in,
	     int *pretire_in,
	     int *nretire_in,
	     double *lambda_in,
	     int *verb_in,
	     double *X_out,
	     double *y_out)
{
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];

  /* remove from pall and each cloud; causes changes in
     pall->X and pall->y */
  cloud->Retire(pretire_in, *nretire_in, *lambda_in, *verb_in);

  /* write out the new X matrix and y vector */
  dupv(X_out, cloud->pall->X[0], cloud->pall->n * cloud->pall->m);
  dupv(y_out, cloud->pall->y, cloud->pall->n);
}


/*
 * intervals_R:
 *
 * R-interface to a routine that extracts the rectangle bounds
 * for X[index,var] in the particle cloud
 */

void intervals_R(/* inputs */
		 int *c_in,
		 int *index_in,
		 int *var_in,
		 double *a_out,
		 double *b_out)
{
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];

  cloud->Intervals((*index_in)-1, (*var_in)-1, a_out, b_out);
}


/*
 * treestats_R:
 *
 * R-interface to a routine that extracts average tree
 * information like height and number of observations
 * in each leaf
 */

void treestats_R(/* inputs */
		 int *c_in,
		 double *height_out,
		 double *avgsize_out,
		 double *avgretire_out)
{
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];

  cloud->TreeStats(height_out, avgsize_out, avgretire_out);
}


/*
 * sameleaf_R:
 *
 * R-interface to a routine that extracts counts of the
 * number of time each X location (provided) is in the 
 * same leaf as each other x
 */

void sameleaf_R(/* inputs */
		int *c_in,
		double *X_in,
		int *n_in,
		int *counts_out)
{
  /* IF MISSING DATA THEN NEED TO GET RANDOM SEED */

  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;

  /* matrixify X */
  unsigned int n = (unsigned int) *n_in;
  double **X = new_matrix_bones(X_in, n, m);

  /* call the cloud function */
  cloud->SameLeaf(X, n, counts_out);

  /* clean up */
  free(X);
}


/*
 * predict_R:
 *
 * function to predict at new XX locations based on 
 * the c-th particle cloud; returns the particle average
 * moments of the (student-t) predictive distribution 
 * and its quantiles, and expected improvement (EI) stats
 * if desired
 */

void predict_R(/* inputs */
	       int *c_in,
	       double *XX_in,
	       double *yy_in,
	       int *nn_in,
	       int *verb_in,
	       
 	       /* outputs */
	       double *mean_out,
	       double *var_out,
	       double *q1_out,
	       double *q2_out,
	       double *yypred_out,
	       double *ei_out)
{
  /* IF MISSING DATA THEN NEED TO GET RANDOM SEED */

  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* data and storage for samples from
   the posterior predictive distribution */
  unsigned int nn = (unsigned int) *nn_in;
  double **XX = new_matrix_bones(XX_in, nn, m);

  /* predict at the XX locations for each particle */
  cloud->Predict(XX, yy_in, nn, mean_out, var_out, q1_out, 
		 q2_out, yypred_out, ei_out, verb);

  /* free predictive data */
  free(XX);
}


/*
 * alc_R:
 *
 * function to calculate ALC at new XX locations, over 
 * reference locations Xref, or a rectangle in Xref,
 *  based on the c-th particle cloud
 */

void alc_R(/* inputs */
	   int *c_in,
	   double *XX_in,
	   int *nn_in,
	   double *Xref_in,
	   int *nref_in,
	   int *cat_in,
	   int *approx_in,
	   double *probs_in,
	   int *verb_in,
	       
	   /* outputs */
	   double *alc_out)
{
  /* IF MISSING DATA THEN NEED RANDOM SEED */

  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* data and storage for samples from
   the posterior predictive distribution */
  unsigned int nn = (unsigned int) *nn_in;
  double **XX = new_matrix_bones(XX_in, nn, m);
  double **probs = NULL;
  if(probs_in) probs = new_matrix_bones(probs_in, cloud->N, *nref_in);

  /* reference locations or rectangle for ALC */
  double **Xref = NULL;
  double **rect = NULL;
  if(*nref_in > 0) { /* ALC by new reference locations */
    assert(Xref_in);
    Xref = new_matrix_bones(Xref_in, *nref_in, m);
  } else if(*nref_in == -1) /* ALC by rectangle */
    rect = new_matrix_bones(Xref_in, 2, m);
  else assert(Xref_in == NULL);

  /* cat sanity check */
  if(!rect) assert(!cat_in);
  else assert(cat_in);

  /* deal with ALC */
  assert(alc_out);
  if(Xref) cloud->ALC(XX, nn, Xref, *nref_in, probs, alc_out, verb);
  else if(rect) cloud->ALC(XX, nn, rect, cat_in, (bool) *approx_in, alc_out, verb);
  else cloud->ALC(XX, nn, XX, nn, probs, alc_out, verb);

  /* clean up ALC predictive data */
  if(Xref) free(Xref);
  if(rect) free(rect);
  if(probs) free(probs);
  free(XX);
}


/*
 * ieci_R:
 *
 * function to calculate IECI at new XX locations, over 
 * reference locations Xref (or XX if not specified)
 * based on the c-th particle cloud
 */

void ieci_R(/* inputs */
	   int *c_in,
	   double *XX_in,
	   int *nn_in,
	   double *Xref_in,
	   int *nref_in,
	   double *probs_in,
	   int *verb_in,
	       
	   /* outputs */
	   double *ieci_out)
{
  /* IF MISSING DATA THEN NEED RANDOM SEED */

  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* data and storage for samples from
   the posterior predictive distribution */
  unsigned int nn = (unsigned int) *nn_in;
  double **XX = new_matrix_bones(XX_in, nn, m);
  double **probs = NULL;
  if(probs_in) probs = new_matrix_bones(probs_in, cloud->N, *nref_in);

  /* explicit reference locations or use XX? */
  double **Xref = NULL;
  if(*nref_in > 0) { /* ALC by new reference locations */
    assert(Xref_in);
    Xref = new_matrix_bones(Xref_in, *nref_in, m);
  } else assert(Xref_in == NULL);

  /* deal with IECI */
  assert(ieci_out);
  if(Xref) cloud->IECI(XX, nn, Xref, *nref_in, probs, ieci_out, verb);
  else cloud->IECI(XX, nn, XX, nn, probs, ieci_out, verb);

  /* clean up ALC predictive data */
  if(Xref) free(Xref);
  if(probs) free(probs);
  free(XX);
}


/*
 * alcX_R:
 *
 * function to calculate the ALC statistic at the
 * data locations by integrating over a bounding rectangle
 * 
 */

void alcX_R(/* inputs */
	       int *c_in,
	       double *rect_in,
	       int *cat_in,
	       int *approx_in,
	       int *verb_in,
	       
 	       /* outputs */
	       double *alc_out)
{
  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;

  /* verbosity argument */
  unsigned int verb = *verb_in;

  assert(rect_in);
  double ** rect = new_matrix_bones(rect_in, 2, m);
  assert(alc_out != NULL);
  assert(cat_in);

  /* calculate ALC */
  cloud->ALC(rect, cat_in, (bool) *approx_in, alc_out, verb);

  /* free data */
  free(rect);
}


/*
 * relevance_R:
 *
 * function to calculate the average variance 
 * by integrating over a bounding rectangle, i.e., 
 * calculate the partial dependencies for each 
 * in put direction
 * 
 */

void relevance_R(/* inputs */
		   int *c_in,
		   double *rect_in,
		   int *cat_in,
		   int *approx_in,
		   int *verb_in,
		   
		   /* outputs */
		   double *delta_out)
{
  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;

  /* verbosity argument */
  unsigned int verb = *verb_in;
  
  /* matrix pointers */
  assert(delta_out);
  double **delta = new_matrix_bones(delta_out, cloud->N, cloud->pall->m);
  assert(rect_in);
  double **rect = new_matrix_bones(rect_in, 2, m);

  /* calculate ALC */
  cloud->Relevance(rect, cat_in, (bool) *approx_in, delta, verb);

  /* free data */
  free(rect);
  free(delta);
}

/*
 * entropyX_R:
 *
 * function to calculate the Entropy statistics at the
 * data locations 
 * 
 */

void entropyX_R(/* inputs */
	       int *c_in,
	       int *verb_in,
	       
 	       /* outputs */
	       double *entropy_out)
{
  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* sanity check */
  assert(entropy_out != NULL);

  /* calculate ALC */
  cloud->Entropy(entropy_out, verb);
}


/*
 * predclass_R:
 *
 * function to predict at new XX locations based on 
 * the c-th particle cloud -- special for the classification
 * case
 */

void predclass_R(/* inputs */
		 int *c_in,
		 double *XX_in,
		 int *yy_in,
		 int *nn_in,
		 int *verb_in,
		 
		 /* outputs */
		 double *p_out,
		 double *yypred_out,
		 double *entropy_out)
{
  /* IF MISSING DATA THEN NEED RANDOM SEED */

  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;
  unsigned int nc = cloud->pall->nc;
  assert(nc > 0);

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* matrix pointer for XX */
  unsigned int nn = (unsigned int) *nn_in;
  double **XX = new_matrix_bones(XX_in, nn, m);
  double **p = new_matrix_bones(p_out, nc, nn);

  /* predict at the XX locations for each particle */
  cloud->Predict(XX, yy_in, nn, p, yypred_out, entropy_out, verb);

  /* clean up p */
  free(p);
  free(XX);
}


/*
 * classprobs_R:
 *
 * return the particle probabilities of a particular
 * a particle class -- a sample from the posterior 
 * distribution
 */

void classprobs_R(/* inputs */
		  int *c_in,
		  int *class_in,
		  double *XX_in,
		  int *nn_in,
		  int *verb_in,
		  
		  /* outputs */
		  double *p_out,
		  double *cs_out)
{
  /* IF MISSING DATA THEN NEED RANDOM SEED */

  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;
  unsigned int cl = *class_in;
  assert(cl < (unsigned) cloud->pall->nc);

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* matrix pointer for XX */
  unsigned int nn = (unsigned int) *nn_in;
  double **XX = new_matrix_bones(XX_in, nn, m);

  /* matrix pointers for outputs */
  double **p,  **cs;
  p = cs = NULL;
  if(p_out) p = new_matrix_bones(p_out, cloud->N, nn);
  if(cs_out) cs = new_matrix_bones(cs_out, cloud->N, nn);

  /* predict at the XX locations for each particle */
  cloud->Predict(cl, XX, nn, p, cs, verb);

  /* clean up p */
  if(p) free(p);
  if(cs) free(cs);
  free(XX);
}


/*
 * sens_R:
 *
 * function to to calculate sensitivity indices for
 * the c-th particle cloud; returns main mean effects and
 * quantiles, along with S and calculations
 */

void sens_R(/* inputs */
	    int *c_in,
	    int *class_in,
	    int *nns_in,
	    int *aug_in,
	    double *rect_in,
	    double *shape_in,
	    double *mode_in,
	    int *cat_in,
	    int *ngrid_in,
	    double *span_in,
	    double *Xgrid_t_in,
	    int *verb_in,
	    
	    /* outputs */
	    double *mean_out,
	    double *q1_out,
	    double *q2_out,
	    double *S_out,
	    double *T_out)
{
  GetRNGstate();

  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;
  unsigned int aug = (unsigned) *aug_in;
  if(cloud->pall->model != LINEAR) assert(aug == 0);
  unsigned int N = cloud->N;

  /* verbosity argument */
  unsigned int verb = *verb_in;
  
  /* extract stuff to do with the latin hypercube */
  unsigned int nns = (unsigned int) *nns_in;
  double **rect = NULL;  /* indicates boot method instead */
  if(rect_in) rect = new_matrix_bones(rect_in, 2, m-aug);
  else assert(shape_in == NULL && mode_in == NULL);

  /* storage for samples from the mean effects distribution */
  unsigned int ngrid = (unsigned int) *ngrid_in;
  double **Xgrid_t = new_matrix_bones(Xgrid_t_in, m, ngrid);
  double **mean = new_matrix_bones(mean_out, m, ngrid);
  double **q1 = new_matrix_bones(q1_out, m-aug, ngrid);
  double **q2 = new_matrix_bones(q2_out, m-aug, ngrid);

  /* storage for the samples of the sensitivity indices */
  double **S = new_matrix_bones(S_out, N, m-aug);
  double **T = new_matrix_bones(T_out, N, m-aug);

  /* sensitivity calculations for each particle */
  cloud->Sens(class_in, nns, aug, rect, shape_in, mode_in, cat_in,
	      Xgrid_t, ngrid, *span_in, mean, q1, q2, S, T, verb);

  /* free bones */
  free(rect);
  free(Xgrid_t); free(mean); free(q1); free(q2);
  free(S); free(T);

  PutRNGstate();
}


/*
 * varpropuse_R:
 *
 * function to tally the number of particles that 
 * use a variable as a split location
 */

void varpropuse_R(int *c_in, double* props_out)
{
  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];

  /* call the varcount function for clounds */
  cloud->VarPropUse(props_out);
}


/*
 * varproptotal_R:
 *
 * function to tally the number of particles that 
 * use a variable as a split location
 */

void varproptotal_R(int *c_in, double* props_out)
{
  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];

  /* call the varcount function for clounds */
  cloud->VarPropTotal(props_out);
}


/*
 * copy_cloud_R:
 *
 * function to copy a cloud in memory into a 
 * newly allocated cloud reference
 */

void copy_cloud_R(int *c_in, int *c_out)
{
  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];

  /* get a new cloud index */
  unsigned int cnew = get_cloud();
  
  /* allocate a new particle cloud */
  clouds[cnew] = new Cloud(cloud);

  /* return to R */
  *c_out = (int) cnew;
}


/*
 * get_cloud:
 *
 * returns an integer reference to a free particle 
 * cloud that can be allcocated to dynaTree inference
 * by Particle Learning
 */

unsigned int get_cloud(void)
{
  if(NC == 0) {
    assert(clouds == NULL);
    clouds = (Cloud**) malloc(sizeof(Cloud*));
    clouds[0] = NULL;
    NC = 1;
    return 0;
  } else {
    for(unsigned int i=0; i<NC; i++) {
      if(clouds[i] == NULL) return i;
    }
    clouds = (Cloud**) realloc(clouds, sizeof(Cloud*) * (2*NC));
    for(unsigned int i=NC; i<2*NC; i++) clouds[i] = NULL;
    NC *= 2;
    return NC/2;
  }
}


/* 
 * delete_cloud:
 *
 * delete the i-th cloud
 */

void delete_cloud(unsigned int i)
{
  if(clouds == NULL || clouds[i] != NULL) { 
    delete clouds[i];
    clouds[i] = NULL;
  } else error("cloud %d is not allocated\n", i);
}


/*
 * delete_cloud_R:
 *
 * R-interface to delete_cloud
 */

void delete_cloud_R(int *cloud)
{
  delete_cloud(*cloud);
}


/*
 * delete_clouds:
 *
 * delete all of the clouds stored in
 * the array, and destroy the array
 */

void delete_clouds(void)
{
  for(unsigned int i=0; i<NC; i++) {
    if(clouds[i]) {
      myprintf(stdout, "removing cloud %d\n", i);
      delete clouds[i];
    }
  }
  free(clouds);
  clouds = NULL;
  NC = 0;
}


/*
 * delete_clouds_R:
 *
 * R interface to delete_clouds
 */

void delete_clouds_R(void)
{ 
  if(clouds)
    delete_clouds();
}


/*
 * alloc_XNA:
 *
 * read Xna_in and XNA_in to create to count the number
 * of missing, NA, entries and allocate the nna * m indicator
 * matrix of missingness
 */

int ** alloc_XNA(unsigned int T, unsigned int m, double **X, 
		 int *Xna_in, int *XNA_in, unsigned int *nna)
{
  int **XNA = NULL;
  *nna = 0;

  /* calculating the dimensions of XNA */
  if(Xna_in) {
    for(unsigned int i=0; i<T; i++) {
      if(Xna_in[i] == 0) Xna_in[i] = -1;
      else { Xna_in[i] = *nna; (*nna)++; }
    }

    /* allocating full matrix XNA */
    XNA = new_imatrix_bones(XNA_in, *nna, m);

    /* put Infs in missing values of X */
    /* while sanity checking XNA versus X entries */
    for(unsigned int i=0; i<T; i++) {
      if(Xna_in[i] >= 0) {
	for(unsigned int j=0; j<m; j++) {
	  if(XNA[Xna_in[i]][j]) {
	    assert(X[i][j] == -12345.0);
	    X[i][j] = -INF;
	  }
	}
      }
    }
  } else assert(XNA_in == NULL);
  return(XNA);
}
  
} // extern brace

