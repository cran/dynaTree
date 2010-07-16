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
#include "assert.h"
#include <stdlib.h>
#include <math.h>
}
#include "cloud.h"

/*
 * Cloud:
 *
 * particle cloud constructor function 
 */

Cloud::Cloud(unsigned int N, Pall *pall, int *pstart, 
	     unsigned int nstart)
{
  /* set pall */
  this->pall = pall;

  /* allocate new particles */
  this->N = N;
  particle = (Particle **) malloc(sizeof(Particle*) * N);
  for(unsigned int i=0; i<N; i++) 
    particle[i] = new Particle(pall, pstart, nstart);

  /* allocate indices and other helper variables */
  index = new_ivector(N);
  prob = new_vector(N);
  rsmin = N;
}


/*
 * ~Cloud:
 *
 * particle cloud destructor function 
 */

Cloud::~Cloud(void)
{
  /* destroy particles */
  for(unsigned int i=0; i<N; i++) delete particle[i];
  free(particle);
  
  /* clean-up */
  delete_pall(pall);
  free(prob);
  free(index);
}


/*
 * Resample:
 *
 * Particle Learning resample function applied on the
 * t-th observation contained in the data structure 
 * using the posterior predictive calculations -- this
 * function also gathers summaries of the particle 
 * cloud trees for printing to the screen if desired
 */

double Cloud::Resample(unsigned int t, unsigned int verb)
{ 
  double prob_var, pred;
  pred = 0.0;

  /* sanity check */
  assert(t < pall->n);

  /* calculate predictive probabilites and gather statistics */
  for(unsigned int i=0; i<N; i++) {
    prob[i] = particle[i]->PostPred(pall->X[t], pall->y[t]);
    pred += prob[i];
  }
   
  /* normalize the weights and calculate their variance */
  prob_var = norm_weights(prob, N);
  
  /* resample particles */
  unsigned int np = Resample();
  
  /* print progress meter */
  if(verb > 0 && (t+1) % verb == 0){

    /* gather tree stats */
    double height, avgsize;
    avgsize = height = 0.0;
    for(unsigned int i=0; i<N; i++){
      height += particle[i]->getHeight();
      avgsize += particle[i]->AvgSize();
    }
    height /= (double) N;
    avgsize /= (double) N;

    /* print tree stats */
    myprintf(stdout, "t=%d, np=%d, v(w)=%g, avg: depth=%g, size=%g\n",
	     t+1, np, prob_var, height, avgsize);
    myflush(stdout);
  }

  /* return marginal likelihood contribution */
  return log(pred) - log((double) N);
}


/*
 * Resample:
 *
 * the actual random re-sampling of particles.  
 * Must be called *after* the probabilities have
 * been calculated
 */

unsigned int Cloud::Resample()
{
  /* first re-sample indices */
  unsigned int np = indexsample(index, N, N, prob);
  if(np < rsmin) rsmin = np;

  /* resample and keep track of which particles have been
     resampled (multiple) times */
  Particle ** resampled = (Particle **) malloc(sizeof(Particle*) * N);
  bool *in = (bool*) malloc(sizeof(bool) * N);
  for(unsigned int i=0; i<N; i++) in[i] = false;
  for(unsigned int i=0; i<N; i++) { 
    if(in[index[i]]) resampled[i] = new Particle(particle[index[i]]);
    else {	
      resampled[i] = particle[index[i]];
      in[index[i]] = true;
    }
  }

  /* swap and delete old */
  for(unsigned int i=0; i<N; i++) { 
    if(!in[i]) delete particle[i];
    particle[i] = resampled[i];
  }
  
  /* clean up */
  free(in);
  free(resampled);

  /* return number of new particles */
  return np;
}


/*
 * Propagate:
 *
 * Particle Learning propagate step applied to the
 * t-th observation in the data
 */

void Cloud::Propagate(unsigned int t)
{
  /* sanity check */
  assert(t < pall->n);

  /* propagate particles */
  for(unsigned int i=0; i<N; i++) particle[i]->Propagate(t);
}


/*
 * Posterior:
 *
 * calculate the average marginal posterior (over all
 * particles)
 */

double Cloud::Posterior(void)
{
  double post = 0.0;
  for(unsigned int i=0; i<N; i++) 
    post += exp(particle[i]->Posterior());
  return(log(post) - log((double) N));
}


/*
 * Predict:
 *
 * collect the moments and quantiles of the particle 
 * predictive distribution at XX predictive locations
 */

void Cloud::Predict(double **XX, unsigned int nn, double **mean,
		    double **sd, double **df, double **var, double **q1, 
		    double **q2, unsigned int verb)
{
  /* predict at the XX locations for each particle */
  for(unsigned int i=0; i<N; i++){

    /* print progress if required */
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }

    /* quantiles pointers if required */
    double *q1i, *q2i;
    if(q1) { q1i = q1[i]; q2i = q2[i]; }
    else q1i = q2i = NULL;

    /* t-variance and DoF pointers if required */
    double *sdi, *dfi;
    if(sd) {
      assert(df);
      sdi = sd[i]; dfi = df[i];
    } else { sdi = dfi = NULL; }

    /* predict at XX */
    particle[i]->Predict(XX, nn, mean[i], sdi, dfi, var[i], q1i, q2i);
  }
}


/*
 * Predict:
 *
 * collect the moments and quantiles of the particle 
 * predictive distribution at XX predictive locations
 */

void Cloud::Predict(double **XX, unsigned int nn, double ***p,
		    double **entropy, unsigned int verb)
{
  /* temporary storage for per-particle predictions */
  double **ptemp = new_matrix(pall->nc, nn);

  /* predict at the XX locations for each particle */
  for(unsigned int i=0; i<N; i++){

    /* print progress if required */
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }

    /* predict at XX */
    particle[i]->Predict(XX, nn, ptemp, entropy[i]);

    /* copy ptemp to the right plot in p */
    for(unsigned int j=0; j<pall->nc; j++) {
      dupv(p[j][i], ptemp[j], nn);
    }
  }

  /* clean up */
  delete_matrix(ptemp);
}


/*
 * ALC:
 *
 * compute the Active Learning Cohn statistic for sequential
 * design, and store in the n-vector alc
 */

void Cloud::ALC(double **XX, unsigned int nn, double *alc_out, 
		unsigned int verb)
{
  /* allocate for all xx pairs */
  double **alc = new_zero_matrix(nn, nn);

  /* accumulate ALC at XX locations for each particle */
  for(unsigned int i=0; i<N; i++){
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }
    particle[i]->ALC(XX, nn, alc);
  }

  /* now average over the number of particles */
  scalev(*alc, nn*nn, 1.0/((double) N));

  /* now average over the rows */
  for(unsigned int i=0; i<nn; i++) {
    alc_out[i] = meanv(alc[i], nn);
  }

  /* clean up */
  delete_matrix(alc);
}


/*
 * GetRsmin:
 *
 * return the minimum number of re-sampled particles
 * encountered so far 
 */

unsigned int Cloud::GetRsmin(void)
{
  return rsmin;
}


/* 
 * indexsample:
 *
 * Just simple "with replacement" index sampling
 * Returns the number of unique components.
 * 
 */

int indexsample(int *ind, int n, int num_probs, double *probs)
{
  double pick;
  int i, counter;
  double *cumprob = new_vector(num_probs);
  double *selected = new_zero_vector(num_probs);

  assert(num_probs > 0);
  assert(n > 0);

  assert(probs[0] >= 0);
  cumprob[0] = probs[0];
  for(i=1; i<num_probs; i++) {
    assert(probs[i] >= 0);
    cumprob[i] = cumprob[i-1] + probs[i];
  }
  if(cumprob[num_probs-1] < 1.0) cumprob[num_probs-1] = 1.0;

  for(i=0; i<n; i++) {
    counter = 0;
    pick = unif_rand();
    while(cumprob[counter] < pick) counter++;
    ind[i] = counter;
    selected[counter]++;
  }

  counter = 0;
  for(i=0; i<num_probs; i++) if(selected[i] > 0) counter++;
  free(cumprob);
  free(selected);
  return counter;
}


/*
 * norm_weights:
 *
 * normalize a weight vector and return the variance
 * of the weight vector as a summary statistic that
 * can be used to track the performance of the SMC
 */

double norm_weights(double *v, int n)
{
  int i;
  double vsum = 0.0;
  for(i=0; i<n; i++) vsum += v[i];
  if(vsum == 0.0 || isnan(vsum)) {
    assert(0);
    printf("zero/nan vector sum in normalize, replacing with unif\n");
    for(i=0; i<n; i++) v[i] = 1.0/n;
    vsum = 1.0;
  }
  double vavg = 0.0;
  double var = 0.0;
  for(i=0; i<n; i++){
    v[i] = v[i]/vsum;
    vavg += v[i];
    var += v[i]*v[i];
  }
  vavg = vavg/((double) n);
  var = var/((double) n) - (vavg*vavg);
  if(var < 0.0) var = 0.0;
  return var;
}
