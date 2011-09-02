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
#include "assert.h"
#include "linalg.h"
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
  this->N = this->Nrevert = N;
  particle = (Particle **) malloc(sizeof(Particle*) * N);
  for(unsigned int i=0; i<N; i++) 
    particle[i] = new Particle(pall, pstart, nstart);

  /* allocate indices and other helper variables */
  index = new_ivector(N);
  prob = new_vector(N);
}


/*
 * Cloud:
 *
 * particle cloud new copy constructor function 
 */

Cloud::Cloud(Cloud *cold) 
{
  /* set pall */
  this->pall = copy_pall(cold->pall);

  /* allocate new particles */
  this->N = cold->N;
  this->Nrevert = cold->Nrevert;
  particle = (Particle **) malloc(sizeof(Particle*) * N);
  for(unsigned int i=0; i<N; i++) 
    particle[i] = new Particle(cold->particle[i], pall);

  /* allocate indices and other helper variables */
  index = new_dup_ivector(cold->index, N);
  prob = new_dup_vector(cold->prob, N);
}


/*
 * Cobine:
 *
 * combine two particle clouds; must have identical palls
 */


void Cloud::Combine(Cloud *c2) 
{
  /* sanity check */
  assert(this->pall == c2->pall);

  /* allocate new particles */
  particle = (Particle **) realloc(particle, sizeof(Particle*) * (N + c2->N));
  for(unsigned int i=0; i<c2->N; i++) {
    particle[N+i] = c2->particle[i];
    c2->particle[i] = NULL;
  }
  c2->pall = NULL;

  /* combine indices and other helper variables */
  prob = (double*) realloc(prob, sizeof(double) * (N+c2->N));
  dupv(prob + N, c2->prob, c2->N);

  /* update N */
  this->N += c2->N;
}


/*
 * ~Cloud:
 *
 * particle cloud destructor function 
 */

Cloud::~Cloud(void)
{
  /* destroy particles */
  for(unsigned int i=0; i<N; i++) {
    if(particle[i]) delete particle[i];
  }
  if(particle) free(particle);
  
  /* clean-up */
  if(pall) delete_pall(pall);
  if(prob) free(prob);
  if(index) free(index);
}


/*
 * Reorder:
 *
 * re-order the X-y paris in pall, and the pointers to
 * those entries recorded in the particles/trees 
 */

void Cloud::Reorder(int *o)
{
  /* sanity check */
  assert(o);

  /* reorder the trees in each particle */
  for(unsigned int i=0; i<N; i++) particle[i]->Reorder(o);
  
  /* re-rder the X-y entries in pall */
  reorder(pall, o);
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

  /* special handling for missing data */
  int *xna = NULL;
  if(pall->Xna && pall->Xna[t] >= 0) xna = pall->XNA[pall->Xna[t]];

  /* calculate predictive probabilites and gather statistics */
  for(unsigned int i=0; i<N; i++) {
    prob[i] = particle[i]->PostPred(pall->X[t], pall->y[t], xna);
    pred += prob[i];
  }

  /* normalize the weights and calculate their variance */
  prob_var = norm_weights(prob, N);
  
  /* resample particles */
  unsigned int np = N; /* no need to resample if all particles
			  are the same for new data */
  /* only PROBLEM with this is Nrevert via rejuvinate */
  if(prob_var > 0) np = Resample();

  /* print progress meter */
  if(verb > 0 && (t+1+pall->g) % verb == 0){

    /* gather tree stats */
    double height, avgsize, avgretire;
    avgretire = avgsize = height = 0.0;
    for(unsigned int i=0; i<N; i++){
      height += particle[i]->getHeight();
      avgsize += particle[i]->AvgSize();
      if(pall->g > 0) avgretire += particle[i]->AvgRetired();
    }
    height /= (double) N;
    avgsize /= (double) N;
    avgretire /= (double) N;

    /* print tree stats */
    if(pall->g > 0) myprintf(stdout, "t=%d[%d]", t+1+pall->g, t+1);
    else myprintf(stdout, "t=%d", t+1+pall->g);
    myprintf(stdout, ", np=%d, v(w)=%g, avg: depth=%g, size=%g",
	     np, prob_var, height, avgsize);
    if(pall->g > 0) myprintf(stdout, "(%g)", avgretire);
    myprintf(stdout, "\n");
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
  // indexsample(index, N, N, prob);
  ressample(index, Nrevert, N, prob);
  
  /* resample and keep track of which particles have been
     resampled (multiple) times */
  Particle ** resampled = (Particle **) malloc(sizeof(Particle*) * Nrevert);
  bool *in = (bool*) malloc(sizeof(bool) * N);
  for(unsigned int i=0; i<N; i++) in[i] = false;
  
  /* copy/move particles */
  unsigned int np = 0;
  for(unsigned int i=0; i<Nrevert; i++) { 
    if(in[index[i]]) resampled[i] = new Particle(particle[index[i]]);
    else {	
      np++;
      resampled[i] = particle[index[i]];
      in[index[i]] = true;
    }
  }

  /* delete particles that are not being carried over */
  for(unsigned int i=0; i<N; i++) if(!in[i]) delete particle[i];
    
  /* deal with possible N & Nrevert discrepency */
  if(N != Nrevert) {
    particle = (Particle **) realloc(particle, sizeof(Particle*) * Nrevert);
    prob = (double*) realloc(prob, sizeof(double) * Nrevert);
    N = Nrevert;
  }
  
  /* copy resampled particles */
  for(unsigned int i=0; i<N; i++) particle[i] = resampled[i];

  /* clean up */
  free(in);
  free(resampled);

  /* return unique number of new particles */
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
 * Retire:
 *
 * function to remove the input/output pair for a certain index
 * and moving the info to the prior(s) to the extent possible;
 * the indices must be sorted
 */

void Cloud::Retire(int *pretire, unsigned int nretire, double lambda, 
		   unsigned int verb)
{
  /* ensure that we're not removing too many entries */
  assert(nretire < pall->n);

  /* print how many will be removed */
  if(verb > 0) {
    myprintf(stdout, "Retiring %d observations: ");
    printIVector(pretire, nretire, stdout);
  }

  /* for each index to be chosted */
  for(unsigned int i=0; i<nretire; i++) {

    /* ensure that we're removing a valid entry */
    assert(pretire[i] < (int) pall->n);

    /* print the index and X and y values that are being removed */
    if(verb > 0) {
      myprintf(stdout, "removing y=%g and X=", pall->y[pretire[i]]);
      printVector(pall->X[pretire[i]], pall->m, stdout, HUMAN);
    }

    /* remove from each particle */
    for(unsigned int p=0; p<N; p++)
      particle[p]->Retire(pretire[i], lambda);

    /* remove from pall->X and pall->Y */
    retire(pall, pretire[i]);

    /* decrement pretire */
    for(unsigned int j=i+1; j<nretire; j++) 
      if(pretire[j] == (int) pall->n) { pretire[j] = pretire[i]; break; }
  }
}


/*
 * Intervals:
 *
 * extracts the intervals (rectangle bounds) for X[index,var]
 * in the particle cloud
 */


void Cloud::Intervals(unsigned int index, unsigned int var, double *a, double *b)
{
  /* sanity check */
  assert(index < pall->n);
  assert(var < pall->m);

  for(unsigned int p=0; p<N; p++) 
    particle[p]->Interval(index, var, a+p, b+p);

}


/*
 * Predict:
 *
 * collect the moments and quantiles of the particle 
 * predictive distribution at XX predictive locations
 */

void Cloud::Predict(double **XX, double *yy, unsigned int nn, double *mean,
		    double *var, double *q1, double *q2, double *yypred, 
		    double *ei, unsigned int verb)
{

  /* quantiles pointers if required */
  double *q1i, *q2i;
  if(q1) { 
    assert(q2);
    q1i = new_vector(nn); zerov(q1, nn);
    q2i = new_vector(nn); zerov(q2, nn); 
  }
  else q1i = q2i = NULL;

  /* var and mean should be allocated */
  assert(mean); zerov(mean, nn);
  assert(var); zerov(var, nn);
  double *meani = new_vector(nn);
  double *vari = new_vector(nn);
  double *m2 = new_zero_vector(nn);

  /* predictive probability if required */
  double *yypredi;
  if(yy) {
    assert(yypred);
    yypredi = new_vector(nn); zerov(yypred, nn);
  } else yypredi = NULL;

  /* t-variance and DoF pointers if required */
  double *sdi, *dfi;
  if(ei) {
    sdi = new_vector(nn);
    zerov(ei, nn);
    dfi = new_vector(nn);
  } else { sdi = dfi = NULL; }
  unsigned int which;
  double fmin;

  /* predict at the XX locations for each particle */
  for(unsigned int i=0; i<N; i++){

    /* print progress if required */
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }

    /* predict at XX */
    particle[i]->Predict(XX, yy, nn, meani, sdi, dfi, vari, 
			 q1i, q2i, yypredi, NULL);

    /* accumulate for mean */
    linalg_daxpy(nn, 1.0, meani, 1, mean, 1);

    /* accumulate yypred */
    if(yy) linalg_daxpy(nn, 1.0, yypredi, 1, yypred, 1);

    /* accumulate for EI */
    /* must be done outside of particle predict due to fmin */
    if(ei) {
      fmin = min(meani, nn, &which);
      for(unsigned int j=0; j<nn; j++)
	ei[j] += EI(meani[j], sdi[j], dfi[j], fmin);
    }

    /* accumulate for variance */
    linalg_daxpy(nn, 1.0, vari, 1, var, 1); /* mean of var */
    for(unsigned int j=0; j<nn; j++) meani[j] *= meani[j];
    linalg_daxpy(nn, 1.0, meani, 1, m2, 1); /* second moment */

    /* accumulate for quantiles */
    if(q1) linalg_daxpy(nn, 1.0, q1i, 1, q1, 1);
    if(q2) linalg_daxpy(nn, 1.0, q2i, 1, q2, 1);
  }

  /* finish variance calculations */
  scalev(mean, nn, 1.0/N);
  if(yy) scalev(yypred, nn, 1.0/N);
  if(ei) scalev(ei, nn, 1.0/N);
  scalev(m2, nn, 1.0/N);
  scalev(var, nn, 1.0/N);
  for(unsigned int j=0; j<nn; j++) var[j] += m2[j] - sq(mean[j]);
  if(q1) scalev(q1, nn, 1.0/N);
  if(q2) scalev(q2, nn, 1.0/N);

  /* clean up */
  if(yypredi) free(yypredi);
  if(q1i) free(q1i);
  if(q2i) free(q2i);
  free(meani);
  free(vari);
  free(m2);
  if(sdi) free(sdi);
  if(dfi) free(dfi);
}


/*
 * Predict:
 *
 * collect the predictive distribution of class 
 * labels at XX predictive locations, as well as 
 * entropy, etc.
 */

void Cloud::Predict(double **XX, int *yy, unsigned int nn, double **p,
		    double *yypred, double *entropy, unsigned int verb)
{
  /* temporary storage for per-particle predictions */
  double **pi = new_matrix(pall->nc, nn);
  zerov(*p, pall->nc*nn);
  double *ei = new_vector(nn);
  double *yypredi = NULL;
  if(yy) {
    yypredi = new_vector(nn);
    zerov(yypred, nn);
  }

  /* predict at the XX locations for each particle */
  for(unsigned int i=0; i<N; i++){

    /* print progress if required */
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }

    /* predict at XX */
    particle[i]->Predict(XX, yy, nn, pi, yypredi, ei);

    /* accumulate */
    linalg_daxpy(nn*pall->nc, 1.0, *pi, 1, *p, 1);
    linalg_daxpy(nn, 1.0, ei, 1, entropy, 1);
    if(yy) linalg_daxpy(nn, 1.0, yypredi, 1, yypred, 1);
  }

  /* normalize */
  scalev(*p, nn*pall->nc, 1.0/N);
  scalev(entropy, nn, 1.0/N);
  if(yy) scalev(yypred, nn, 1.0/N); 

  /* clean up */
  if(yypredi) free(yypredi);
  delete_matrix(pi);
  free(ei);
}


/*
 * Predict:
 *
 * collect the full sample of predictive distribution of a 
 * particular class label at XX predictive locations
 */

void Cloud::Predict(unsigned int cl, double **XX, unsigned int nn, double **p,
		    double **c, unsigned int verb)
{
  double *ci, *pi;

  /* predict at the XX locations for each particle */
  for(unsigned int i=0; i<N; i++){

    /* print progress if required */
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }

    /* predict at XX */
    
    if(p) pi = p[i]; else pi = NULL;
    if(c) ci = c[i]; else ci = NULL;
    particle[i]->Predict(cl, XX, nn, pi, ci);
  }
}


/*
 * Sens:
 *
 * sensitivity analysis of each particle, accumulating the main
 * effect results, and retaining samples of the S and T indices
 */

void Cloud::Sens(int *cls, unsigned int nns, unsigned int aug, 
		 double **rect, double *shape, double *mode, int *cat,
		 double **Xgrid_t, unsigned int ngrid, double span, 
		 double **mean, double **q1, double **q2, 
		 double **S, double **T, unsigned int verb)
{
  /* allocate temporary space for accumulation of main effects */
  double ***M = (double***) malloc(sizeof(double **) * (pall->m - aug));
  for(unsigned int j=0; j<pall->m-aug; j++) M[j] = new_matrix(N, ngrid);
  double **mtemp = new_matrix(pall->m - aug, ngrid);
  
  /* for each particle */
  for(unsigned int i=0; i<N; i++) {

    /* print progress if required */
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }

    /* sensitivity within the particle */
    if(!cls) particle[i]->Sens(nns, aug, rect, shape, mode, cat, Xgrid_t,
			       ngrid, span, mtemp, S[i], T[i]); /* regression */
    else particle[i]->Sens(*cls, nns, aug, rect, shape, mode, cat, 
			   Xgrid_t, ngrid, span, mtemp, S[i], T[i]); /* class */

    /* accumulate */
    for(unsigned int j=0; j<pall->m - aug; j++) 
      dupv(M[j][i], mtemp[j], ngrid);
  }
  delete_matrix(mtemp);

  /* now summary statistics */
  double q[2] = {0.05, 0.95};
  for(unsigned int j=0; j<pall->m-aug; j++) {
    wmean_of_columns(mean[j], M[j], N, ngrid, NULL);
    double **Q = (double**) malloc(sizeof(double*) * 2);
    Q[0] = q1[j];  Q[1] = q2[j];
    quantiles_of_columns(Q, q, 2, M[j], N, ngrid, NULL);
    delete_matrix(M[j]);
    free(Q);
  }
  free(M);
}


/*
 * VarPropUse:
 *
 * count the proportion of particles which use var
 * somewhere in a tree partition
 */

void Cloud::VarPropUse(double *prop)
{
  /* make a new count vector to be used by each parcile */
  int* c = new_ivector(pall->m);
  int *counts = new_zero_ivector(pall->m);

  /* for each particle */
  for(unsigned int i=0; i<N; i++) {
    particle[i]->VarCountUse(c);
    add_ivector(counts, c, pall->m);
  }

  /* normalize */
  for(unsigned int j=pall->smin; j<pall->m; j++)
    prop[j] = ((double) counts[j]) / ((double) N);

  /* invalid values */
  for(unsigned int j=0; j<pall->smin; j++)
    prop[j] = -1.0;

  /* clean up */
  free(c); free(counts);
}


/*
 * VarPropTotal:
 *
 * count the average proportion of the tree in each particle
 * which ivolves splits on each predictor
 */

void Cloud::VarPropTotal(double *prop)
{
  /* make a new count vector to be used by each parcile */
  double* c = new_vector(pall->m);

  /* initialize the totals vector */
  double *totals = new_zero_vector(pall->m);

  /* for each particle */
  for(unsigned int i=0; i<N; i++) {
    particle[i]->VarCountTotal(c);
    add_vector(1.0, totals, 1.0, c, pall->m);
  }

  /* normalize */
  for(unsigned int j=0; j<pall->m; j++)
    prop[j] = totals[j] / ((double) N);

  /* invalid values */
  for(unsigned int j=0; j<pall->smin; j++)
    prop[j] = -1.0;

  /* clean up */
  free(c); free(totals);
}


/*
 * IECI:
 *
 * compute the Integrated Expected Conditional Improvement 
 * statistic for sequential design, and store in the n-vector 
 * ieci
 */

void Cloud::IECI(double **XX, unsigned int nn, double **Xref, 
		 unsigned int nref, double **probs, 
		 double *ieci_out, unsigned int verb)
{
  /* allocate for all xx pairs */
  double **ecis = new_zero_matrix(nn, nref);

  /* accumulate EImECI at XX locations for each particle */
  double *pi;
  for(unsigned int i=0; i<N; i++){
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }
    if(probs) pi = probs[i];
    else pi = NULL;
    particle[i]->EImECI(XX, nn, Xref, nref, pi, ecis);
  }

  /* now average over the number of particles */
  scalev(*ecis, nn*nref, 1.0/((double) N));

  /* now average over the rows */
  for(unsigned int i=0; i<nn; i++) 
    ieci_out[i] = meanv(ecis[i], nref);

  /* clean up */
  delete_matrix(ecis);
}


/*
 * ALC:
 *
 * compute the Active Learning Cohn statistic for sequential
 * design, and store in the n-vector alc
 */

void Cloud::ALC(double **XX, unsigned int nn, double **Xref, 
		unsigned int nref, double **probs, 
		double *alc_out, unsigned int verb)
{
  /* allocate for all xx pairs */
  double **alc = new_zero_matrix(nn, nref);

  /* accumulate ALC at XX locations for each particle */
  double *pi;
  for(unsigned int i=0; i<N; i++){
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }
    if(probs) pi = probs[i];
    else pi = NULL;
    particle[i]->ALC(XX, nn, Xref, nref, pi, alc);
  }

  /* now average over the number of particles */
  scalev(*alc, nn*nref, 1.0/((double) N));

  /* now average over the rows */
  for(unsigned int i=0; i<nn; i++) alc_out[i] = meanv(alc[i], nref);

  /* clean up */
  delete_matrix(alc);
}


/*
 * ALC:
 *
 * compute the Active Learning Cohn statistic for sequential
 * design by integrating over the rectangle of reference locations
 * and store in the n-vector alc_out
 */

void Cloud::ALC(double **XX, unsigned int nn, double **rect,
		int *cat, bool approx, double *alc_out, unsigned int verb)
{
  /* initialize */
  zerov(alc_out, nn);

  /* accumulate ALC at XX locations for each particle */
  for(unsigned int i=0; i<N; i++){
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }
    particle[i]->ALC(XX, nn, rect, cat, approx, alc_out);
  }

  /* calculate the area of the rectangle for normalization */
  double area = 1.0;
  if(approx) area = (double) (pall->n + pall->g); 
  else {
    for(unsigned int j=0; j<pall->bmax; j++) {
      if(cat[j] || rect[1][j] - rect[0][j] < DOUBLE_EPS) continue;
      area *= rect[1][j] - rect[0][j];
    }
  }

  /* now average over the number of particles, and normalize */
  scalev(alc_out, nn, 1.0/(((double) N)*area));
}


/*
 * ALC:
 *
 * compute the Active Learning Cohn statistic for sequential
 * design at the input locations pall->X by integrating over 
 * the rectangle of reference locations and store in the 
 * n-vector alc_out
 */

void Cloud::ALC(double **rect, int *cat, bool approx, double *alc_out,
		unsigned int verb)
{
  /* initialize */
  zerov(alc_out, pall->n);

  /* accumulate ALC at XX locations for each particle */
  for(unsigned int i=0; i<N; i++) {
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }
    particle[i]->ALC(rect, cat, approx, alc_out);
  }

  /* calculate the area of the rectangle for normalization */
  double area = 1.0; 
  if(approx) area = (double) (pall->n + pall->g);
  else {
    for(unsigned int j=0; j<pall->bmax; j++) {
      if(cat[j] || rect[1][j] - rect[0][j] < DOUBLE_EPS) continue;
      area *= rect[1][j] - rect[0][j];
    } 
  }

  /* now average over the number of particles, and normalize */
  scalev(alc_out, pall->n, 1.0/(((double) N)*area));
}


/*
 * Relevance:
 *
 * sample the reduction in Average Variance statistic by integrating 
 * over the rectangle of reference locations, and return; i.e., calculate
 * the partial dependencies
 */

void Cloud::Relevance(double **rect, int *cat, bool approx, double **delta, 
		   unsigned int verb)
{
  /* get ALC at XX locations for each particle */
  for(unsigned int i=0; i<N; i++) {
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "relevance %d of %d done\n", i+1, N);
      myflush(stdout);
    }
    particle[i]->Relevance(rect, cat, approx, delta[i]);
  }

  /* calculate the area of the rectangle for normalization */
  double area = 1.0;
  if(approx) area = (double) (pall->n + pall->g);
  else {
    for(unsigned int j=0; j<pall->bmax; j++) {
      if(cat[j] || rect[1][j] - rect[0][j] < DOUBLE_EPS) continue;
      area *= rect[1][j] - rect[0][j];
    }
  }

  /* now average over the number of particles, and normalize */
  scalev(*delta, N*(pall->m), 1.0/area);
}

/*
 * Entropy:
 *
 * compute Entropy statistics for sequential
 * design at the input locations pall->X and store in the 
 * n-vector alc_out
 */

void Cloud::Entropy(double *entropy_out, unsigned int verb)
{
  /* initialize */
  zerov(entropy_out, pall->n);

  /* accumulate ALC at XX locations for each particle */
  for(unsigned int i=0; i<N; i++) {
    if(verb > 0 && (i+1) % verb == 0) {
      myprintf(stdout, "prediction %d of %d done\n", i+1, N);
      myflush(stdout);
    }
    particle[i]->Entropy(entropy_out);
  }

  /* now average over the number of particles, and normalize */
  scalev(entropy_out, pall->n, 1.0/((double) N));
}


/* 
 * indexsample:
 *
 * Just simple "with replacement" index sampling
 * Returns the number of unique components.
 * 
 */

void indexsample(int *ind, int n, int num_probs, double *probs)
{
  double pick;
  int i, counter;
  double *cumprob = new_vector(num_probs);

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
  }

  free(cumprob);
}


/* 
 * ressample:
 *
 * residual resampling to keep variance and duplication
 * of particles to a minimum
 */

void ressample(int *ind, int n, int num_probs, double *probs)
{
  /* calculate expected number for each particle */
  double *expect = new_dup_vector(probs, num_probs);
  scalev(expect, num_probs, ((double) n));

  /* calculate deterministic number */
  int *deter = new_ivector(num_probs);
  int sum = 0;
  for(int i=0; i< num_probs; i++) {
    deter[i] = (int) expect[i];
    sum += deter[i];
  }
  
  /* deterministic indices */
  unsigned int k=0;
  for(int i=0; i<num_probs; i++)
    for(int j=0; j<deter[i]; j++, k++)
      ind[k] = i;
  assert((int) k <= n);

  /* calculate residual size */
  int residual = n - sum;
  
  /* residual indices */
  if(residual > 0) {
    /* adjusted probabilities */
    double *reprobs = new_vector(num_probs);
    for(int i=0; i<num_probs; i++)
      reprobs[i] = (expect[i] - deter[i])/((double) residual);

    /* index sample for the residuals */
    indexsample(ind + k, residual, num_probs, reprobs);
    free(reprobs);
  }

  /* finish clean up */
  free(expect);
  free(deter);
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
    // assert(0);
    myprintf(stderr, "zero/nan vector or sum in normalize, replacing with unif\n");
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
