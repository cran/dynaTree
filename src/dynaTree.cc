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
 * Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
 *
 ****************************************************************************/


extern "C" {
#include "rhelp.h"
#include "R.h"
#include "Rmath.h"
#include "matrix.h"
}
#include "cloud.h"
#include "assert.h"

extern "C" {

unsigned int NC = 0;
Cloud **clouds = NULL;

unsigned int get_cloud(void);
double EI(const double m, const double s2, const double df,
	const double fmin);


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
		double* y_in,
		int* model_in,
		double* params_in,
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
  Pall *pall = new_pall(X, T, m, y, params_in, *model_in);

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* minimum parition size */
  unsigned int minp = (int) params_in[2];

  /* choose starting indices and bounding rectangle */
  unsigned int nstart = minp; // no possible splits until 2*minp
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
  free(y);

  /* return to R */
  *c_out = (int) c;
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

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* add the data and get start/end times */
  unsigned int nstart = cloud->pall->n;
  add_data(cloud->pall, X, T, y);
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
  free(y);

  /* return to R */
  PutRNGstate();
}


/*
 * predict_R:
 *
 * function to predict at new XX locations based on 
 * the c-th particle cloud; returns the particle average
 * moments of the (student-t) predictive distribution 
 * and its quantiles
 */

void predict_R(/* inputs */
	       int *c_in,
	       double *XX_in,
	       int *nn_in,
	       int *verb_in,
	       
 	       /* outputs */
	       double *mean_out,
	       double *var_out,
	       double *q1_out,
	       double *q2_out,
	       double *alc_out,
	       double *ei_out)
{
  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds[c] == NULL) error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;
  unsigned int N = cloud->N;

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* data and storage for samples from
   the posterior predictive distribution */
  unsigned int nn = (unsigned int) *nn_in;
  double **XX = new_matrix_bones(XX_in, nn, m);
  double **mean = new_zero_matrix(N, nn);
  double **var = new_zero_matrix(N, nn);

  /* quantiles might not be required */
  double **q1, **q2;
  if(q1_out) {
    q1 = new_zero_matrix(N, nn);
    q2 = new_zero_matrix(N, nn);
  } else q1 = q2 = NULL;

  /* EI calculations might not be required */
  double **sd, **df;
  if(ei_out) {
    sd = new_zero_matrix(N, nn);
    df = new_zero_matrix(N, nn);
  } else sd = df = NULL;

  /* predict at the XX locations for each particle */
  cloud->Predict(XX, nn, mean, sd, df, var, q1, q2, verb);

  /* average means of particles */
  wmean_of_columns(mean_out, mean, N, nn, NULL);
  
  /* average variance plus variances of means */
  wmean_of_columns(var_out, var, N, nn, NULL);
  double *temp = new_vector(nn);
  wvar_of_columns(temp, mean, N, nn, NULL);
  for(unsigned int i=0; i<nn; i++) var_out[i] += temp[i];
  free(temp);
  delete_matrix(var);
  
  /* average quantiles of particles */
  if(q1_out) {
    wmean_of_columns(q1_out, q1, N, nn, NULL);
    wmean_of_columns(q2_out, q2, N, nn, NULL);
    delete_matrix(q1);
    delete_matrix(q2);
  }
    
  /* calculate EI if required */
  if(ei_out) {
    zerov(ei_out, nn);
    unsigned int which; double fmin;
    // fmin = min(mean_out, nn, &which);
    for(unsigned int i=0; i<N; i++) {
      fmin = min(mean[i], nn, &which);
      for(unsigned int j=0; j<nn; j++) {
	ei_out[j] += EI(mean[i][j], sd[i][j], df[i][j], fmin);
      }
    }
    scalev(ei_out, nn, 1.0/N);

    /* clean up */
    delete_matrix(sd);
    delete_matrix(df);
  }
  delete_matrix(mean);

  /* deal with ALC */
  if(alc_out) cloud->ALC(XX, nn, alc_out, verb);

  /* free predictive data */
  free(XX);
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
		 int *nn_in,
		 int *verb_in,
		 
		 /* outputs */
		 double *p_out,
		 double *entropy_out)
{
  /* get the cloud */
  unsigned int c = *c_in;
  if(clouds == NULL || clouds[c] == NULL) 
    error("cloud %d is not allocated\n", c);
  Cloud *cloud = clouds[c];
  unsigned int m = cloud->pall->m;
  unsigned int N = cloud->N;
  unsigned int nc = cloud->pall->nc;
  assert(nc > 0);

  /* verbosity argument */
  unsigned int verb = *verb_in;

  /* matrix pointer for XX */
  unsigned int nn = (unsigned int) *nn_in;
  double **XX = new_matrix_bones(XX_in, nn, m);

  /* data and storage for samples from
     the posterior predictive distribution */
  double ***p = (double***) malloc(sizeof(double **) *nc);
  for(unsigned int i=0; i<nc; i++) {
    p[i] = new_zero_matrix(N, nn);
  }
  double **entropy = new_zero_matrix(N, nn);

  /* predict at the XX locations for each particle */
  cloud->Predict(XX, nn, p, entropy, verb);

  /* average, write out, and delete */
  for(unsigned int i=0; i<nc; i++) {
    wmean_of_columns(p_out + nn*i, p[i], N, nn, NULL);
    delete_matrix(p[i]);
  }
  wmean_of_columns(entropy_out, entropy, N, nn, NULL);
  delete_matrix(entropy);

  /* clean up p */
  free(p);
  free(XX);
}


/*
 * EI:
 *
 * calculates the expected improvement following
 * Williams et al by integrating over the parameters
 * to the GP predictive
 */

double EI(const double m, const double s2, const double df,
          const double fmin)
{
  double diff, sd, diffs, scale, ei;

  diff = fmin - m;
  sd = sqrt(s2);
  diffs = diff/sd;
  scale = (df*sd + sq(diff)/sd)/(df-1.0);
  ei = diff*pt(diffs, (double) df, 1, 0);
  ei += scale*dt(diffs, (double) df, 0);

  return(ei);
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
  delete_clouds();
}


} // extern brace

