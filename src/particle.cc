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


#include <Rmath.h>
#include <R.h>
extern "C" {
#include "matrix.h"
#include "rhelp.h"
#include "lh.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
}
#include "particle.h"


/*
 * Particle:
 * 
 * constructor for new data
 */

Particle::Particle(Pall *pall_in, int *pstart, unsigned int nstart)
{
  assert(pall_in->minp > 0);
  int *p = new_dup_ivector(pstart, nstart);
  this->pall = pall_in;
  tree = new Tree(this, p, nstart, NULL);
}


/*
 * Particle:
 *
 * pointer-based copy constructor for existing
 * particle
 */

Particle::Particle(Particle *p) 
{
  pall = p->pall;
  tree = new Tree(p->tree, this, NULL);
}


/*
 * Particle:
 *
 * pointer-based copy constructor for existing
 * particle, also swapping t->pall
 */

Particle::Particle(Particle *p, Pall *pall_in) 
{
  this->pall = pall_in;
  tree = new Tree(p->tree, this, NULL);
}


/*
 * ~Particle:
 *
 * particle de-constructor
 */

Particle::~Particle(void)
{
  delete tree;
}


/*
 * Reorder:
 *
 * the pall X-y pairs recorded in the particles/trees 
 */

void Particle::Reorder(int *o)
{
  tree->ReorderP(o);
}
 

/*
 * Propagate:
 *
 * stochastic propagate rule for dynaTree based on
 * stay, grow, and prune moves, incorporating the
 * data at pall->X[index]
 */

void Particle::Propagate(unsigned int index)
{
  /* add the new data point to the tree */
  Tree* leaf = tree->AddDatum(index);

  /* calculate the probability of staying at the same tree */
  double pstay = leaf->stayProb();
  // MYprintf(stdout, "pstay=%g\n", pstay);

  /* calculate the probability of growing to a new tree by 
     splitting at the leaf node that pall->X[index] belongs to */
  int var; double val;
  double pgrow = leaf->growProb(&var, &val);
  // MYprintf(stdout, "pgrow=%g\n", pgrow);
 
  /* calculate the probability of pruning back from the
     parent of the leaf node that pall->X[index] belongs to */
  double pprune = leaf->pruneProb();
  // MYprintf(stdout, "pprune=%g\n\n", pprune);

  /* clever normalization to reduce numerical error */
  double lnorm;
  if(R_FINITE(pprune)) {
    if(R_FINITE(pgrow)) {
      lnorm = fmin2(pstay, fmin2(pprune, pgrow));
      pstay -= lnorm; pgrow -=lnorm; pprune -= lnorm;
    } else {
      lnorm = fmin2(pstay, pprune);
      pstay -= lnorm; pprune -= lnorm;
    }
  } else {
    if(R_FINITE(pgrow)) {
      lnorm = fmin2(pstay, pgrow);
      pstay -= lnorm; pgrow -= lnorm;
    } else {
      pstay = 0.0;
    }
  }
  pstay = exp(pstay); pgrow = exp(pgrow); pprune = exp(pprune);

  /* normalize the probabilities of the three (two) moves */
  double norm = pstay + pprune + pgrow;
  if(!R_FINITE(norm) || ISNAN(norm) || norm==0.0) return; 
  pstay = pstay/norm;
  pprune = pprune/norm;

  /* choose randomly between the three moves */
  double u = unif_rand();

  /* do nothing if staying */
  if(u < pstay) return;

  /* if not prune then automatically grow */
  if(u < pstay + pprune) leaf->Parent()->prune();
  else leaf->grow(var, val);
  /* grow assumes missing pall->X values un-(randomly)-changed
     since growProb calculation */
}


/*
 * PostPred:
 *
 * calculate the posterior probability of a particular
 * (x.y) pair under the tree model in the particle --
 * used in the PL resample ste
 */

double Particle::PostPred(double *xx, double yy, int *xna)
{
  assert(tree != NULL);
  // if(xna) printIVector(xna, pall->m, stdout);
  Tree* leaf = tree->GetLeaf(xx, xna);
  double post = leaf->PostPred(xx, yy);
  return post;
}


/*
 * Posterior:
 *
 * calculate the (log) posterior probability of the
 * entire tree
 */

double Particle::Posterior(void)
{
  return tree->FullPosterior();
}


/*
 * Predict:
 *
 * return the predictive mean, variance, and quantiles for
 * each of the XX locations, under the tree model in the
 * particle
 */

void Particle::Predict(double **XX, double *yy, unsigned int nn, 
		       double *mean, double *sd, double *df, 
		       double *var, double *q1, double *q2,
		       double *yypred, double *ZZ)
{
  /* dumMY pointers memory and pointers */
  double meani, sdi, dfi;

  for(unsigned int i=0; i<nn; i++) {

    /* actually predict */
    tree->Predict(XX[i], &meani, &sdi, &dfi);

    /* maybe save the predictive variances */
    if(var) var[i] = dfi*sq(sdi)/(dfi - 2.0);

    /* maybe save the quantiles of the t */
    if(q1) q1[i] = qt(0.05, dfi, 1, 0)*sdi + meani;
    if(q2) q2[i] = qt(0.95, dfi, 1, 0)*sdi + meani;
    
    /* posterior predictive probability */
    if(yy) yypred[i] = dt((yy[i] - meani)/sdi, dfi, 0)/sdi;
    /* sdi needed for Jacobian */
    
    /* samples from the predictive distribution */
    if(ZZ) ZZ[i] = rt(dfi)*sdi + meani;
    
    /* save mean, s2 and df if allocated */
    if(mean) mean[i] = meani;
    if(sd) sd[i] = sdi;
    if(df) df[i] = dfi; 
  }
}


/*
 * Coef:
 *
 * return the linear regression coefficients for each of the XX locations, 
* under the tree model in the particle
 */

void Particle::Coef(double **XX, unsigned int nn, double **beta)
{
  /* sanity check */
  assert(pall->model == LINEAR);

  for(unsigned int i=0; i<nn; i++) {

    /* actually predict */
    tree->Coef(XX[i], beta[i]);
  }
}

/*
 * Predict:
 *
 * classification based prediction and entropy calculation
 */

void Particle::Predict(double **XX, int *yy, unsigned int nn, double **p,
		       double *yypred, double *entropy)
{
  double *ptemp = new_vector(pall->nc);
  for(unsigned int i=0; i<nn; i++) {
    tree->Predict(XX[i], ptemp);
    entropy[i] = 0.0;
    for(unsigned int j=0; j<pall->nc; j++) {
      p[j][i] = ptemp[j];
      entropy[i] += 0.0 - ptemp[j] * log(ptemp[j]);
    }
    if(yy) yypred[i] = p[yy[i]][i];
  }
  free(ptemp);
}


/*
 * Predict:
 *
 * classification based prediction
 */

void Particle::Predict(unsigned int cls, double **XX, unsigned int nn, 
		       double *probs, double *ZZ)
{
  double p, u;
  assert(cls < pall->nc);
  for(unsigned int i=0; i<nn; i++) {
    p = tree->Predict(XX[i], cls);
    if(probs) probs[i] = p;
    if(ZZ) {
      u = unif_rand();
      ZZ[i] = (double) (u < p);
    }
  }
}


/*
 * Retire:
 *
 * remove the corresponding x-y pair stored at X[index] y[index] 
 * and update the sufficientinformation (moving to the prior) 
 * using as few operations as possible
 */

void Particle::Retire(unsigned int index, double lambda)
{
  /* retire */
  Tree* collapse = tree->RetireDatum(index, lambda);

  /* check to see if we need to remove/prune a leaf node */
  if(collapse) {
    collapse->Collapse();
    delete collapse->Parent(); /* should also delete collapse */
  }

  /* adjust the data pointers */
  tree->DecrementP(pall->n - 1, index);
}


/*
 * Interval:
 *
 * extracts the interval (rectangle bounds) for X[index,var]
 * in the particle cloud
 */


void Particle::Interval(unsigned int index, unsigned int var, double *a, double *b)
{
  Tree *leaf = tree->GetLeaf(index);
  assert(leaf);
  *a = leaf->Min(var);
  *b = leaf->Max(var);
}



/*
 * Sens: (for regression)
 *
 * within-particle calculation of main effects and S and T 
 * sensitivity indices
 */

void Particle::Sens(unsigned int nns, unsigned int aug, double **rect, 
		    double *shape, double *mode, int *cat, double **Xgrid_t, 
		    unsigned int ngrid, double span, double **main, 
		    double *S, double *T)
{
  /* sanity check */
  assert(pall->model != CLASS);

  /* input column dimension */
  unsigned int m = pall->m;

  /* generate XX locations from two LHS samples */
  int nn;
  double **XX;
  if(rect == NULL) 
    XX = sens_boot(nns, m, aug, &nn, pall->X, pall->n);
  else XX = sens_lhs(nns, m, aug, rect, shape, mode, &nn);
  /* first nns rows is M1, second is M2, and each nns group
     thereafter is one of the N_js */

  /* now do the prediction at all nn locations */
  double *ZZ = new_vector(nn);
  double *pmean = new_vector(nn);
  Predict(XX, NULL, nn, pmean, NULL, NULL, NULL, NULL, NULL, NULL, ZZ);

  /* mean main effects (mean, and quantiles) */
  main_effects(XX, 2*nns, m, aug, cat, pmean, Xgrid_t, ngrid, span, main);

  /* sobol indices */
  sobol_indices(ZZ, nns, m-aug, S, T);

  /* clean up */
  delete_matrix(XX);
  free(ZZ); free(pmean);
}


/*
 * Sens: (for classification; could be combined with the above function)
 *
 * within-particle calculation of main effects and S and T 
 * sensitivity indices
 */

void Particle::Sens(unsigned int cls, unsigned int nns, unsigned int aug,
		    double **rect, double *shape, double *mode, int *cat,
		    double **Xgrid_t, unsigned int ngrid, double span, 
		    double **main, double *S, double *T)
{
  /* sanity check */
  assert(pall->model == CLASS && cls < pall->nc);

  unsigned int m = pall->m;

  /* generate XX locations from two LHS samples */
  int nn;
  double **XX;
  if(rect == NULL) XX = sens_boot(nns, m, aug, &nn, pall->X, pall->n);
  else XX = sens_lhs(nns, m, aug, rect, shape, mode, &nn);
  /* first nns rows is M1, second is M2, and each nns group
     thereafter is one of the N_js */

  /* now do the prediction at all nn locations */
  double *ZZ = new_vector(nn);
  double *pcls = new_vector(nn);
  Predict(cls, XX, nn, pcls, ZZ);

  /* mean main effects (mean, and quantiles) */
  main_effects(XX, 2*nns, m, aug, cat, pcls, Xgrid_t, ngrid, 
	       span, main);

  /* sobol indices */
  sobol_indices(ZZ, nns, m-aug, S, T);

  /* clean up */
  delete_matrix(XX);
  free(ZZ); free(pcls);
}


/*
 * VarCountUse:
 *
 * fill c with a Booleans indicating if the tree in this
 * particle uses the corresponding variable as a split 
 * location or not
 */

void Particle::VarCountUse(int *c)
{
  int len;
  Tree **leaves =  tree->internalsList(&len);

  unsigned int var;
  unsigned int used = 0;
  zeroiv(c, pall->m);
  for(int i=0; i<len; i++) {
    var = leaves[i]->GetVar();
    if(c[var] == 0) {
      c[var] = 1;
      used++;
      if(used == pall->m - pall->smin) break;
    }
  }

  free(leaves);
}


/*
 * VarCountTotal:
 *
 * fill c with a the cout of the number of times the tree
 * in this particle uses the corresponding variable as a split 
 * location
 */

void Particle::VarCountTotal(double *c)
{
  int len;
  Tree **leaves =  tree->internalsList(&len);

  zerov(c, pall->m);
  for(int i=0; i<len; i++)
    c[leaves[i]->GetVar()] += 1.0;

  /* normalize */
  if(len > 0) {
    for(unsigned int j=pall->smin; j<pall->m; j++)
      c[j] /= ((double) len);
  }
    
  free(leaves);
}


/*
 * SameLeaf:
 *
 * return a count of the number of other X values that
 * are in the same leaf node as each individual X
 */

void Particle::SameLeaf(double **X, unsigned int n, int *counts)
{
  /* indicators */
  int *p = iseq(0,n-1);

  /* recursive calculation */
  tree->SameLeaf(X, p, n, counts);

  /* clean up */
  free(p);
}


/*
 * ALC:
 *
 * accumulate the (un-normalized) Active Learning Cohn
 * statistic for sequential design
 */


void Particle::ALC(double **XX, unsigned int nn, double **Xref, 
		   unsigned int nref, double *probs, double **alc) 
{
  double ealc;
  for(unsigned int j=0; j<nref; j++) {
    if(probs && probs[j] <= DOUBLE_EPS) continue;
    for(unsigned int i=0; i<nn; i++) {
      ealc = tree->ALC(XX[i], Xref[j]);
      if(probs) alc[i][j] += probs[j] * ealc;
      else alc[i][j] += ealc;
    }
  }
}


/*
 * ECImECI:
 *
 * accumulate the reduction in Expected Conditional Improvement 
 * (from the expected improvement), used for optimization under
 * unknown constraints
 */

void Particle::EImECI(double **XX, unsigned int nn, double **Xref, 
		    unsigned int nref, double *probs, double **eimeci) 
{
  /* first get the predictive moments at the Xref locations */
  double *mean = new_vector(nref);
  double *sd = new_vector(nref);
  double *df = new_vector(nref);
  for(unsigned int j=0; j<nref; j++)
    tree->Predict(Xref[j], mean+j, sd+j, df+j);
  
  /* get fmin */
  unsigned int which;
  double fmin = min(mean, nref, &which);

  /* then calculate EI and ECIs */
  double eci; 
  for(unsigned int j=0; j<nref; j++) {
    if(probs && probs[j] <= DOUBLE_EPS) continue;
    double ei = EI(mean[j], sd[j], df[j], fmin);
    for(unsigned int i=0; i<nn; i++) {
      eci = tree->ECI(XX[i], Xref[j], mean[j], sd[j], fmin, ei);
      if(probs) eimeci[i][j] += probs[j] * (ei - eci);
      else eimeci[i][j] += ei - eci;
    }
  }

  /* clean up */
  free(mean);
  free(sd);
  free(df);
}


/*
 * qEI:
 *
 * quantile expected improvement calculation
 */

void Particle::qEI(double q, double alpha, double **XX, unsigned int nn, 
		      double *qei) 
{
  double mean, sd, df, epsilon, u1, u2;

  for(unsigned int i=0; i<nn; i++) {
    tree->Predict(XX[i], &mean, &sd, &df);
    
    epsilon = alpha * sd * sqrt(df/(df-2.0));
    u1 = (q - mean - epsilon)/sd;
    u2 = (q - mean + epsilon)/sd;

    /* assume qei initialized at zeros */
    qei[i] += (sq(epsilon) - sq(mean - q) - sq(sd))*(pt(u2,df,1,0) - pt(u1,df,1,0));
    qei[i] += sq(sd)*(u2*dt(u2,df,0) - u1*dt(u1,df,0));
    qei[i] += 2.0*(mean - q)*sd*(dt(u2,df,0) - dt(u1,df,0));
  }
}


/*
 * ALC:
 *
 * accumulate the (un-normalized) Active Learning Cohn
 * statistic for sequential design based on integrating over
 * the rectangle of reference locations
 */

void Particle::ALC(double **XX, unsigned int nn, double **rect, 
		   int *cat, bool approx, double *alc)
{
  for(unsigned int i=0; i<nn; i++)
    alc[i] += tree->ALC(XX[i], rect, cat, approx);
}


/*
 * ALC:
 *
 * accumulate the (un-normalized) Active Learning Cohn
 * statistic for sequential design at the input locations
 * pall->X based on integrating over the rectangle of 
 * reference locations
 */

void Particle::ALC(double **rect, int *cat, bool approx, double *alc)
{
  tree->ALC(rect, cat, approx, alc);
}


/*
 * Relevance:
 *
 * accumulate the (un-normalized) reduction in average variance
 * statistic for the each split in tree based on integrating 
 * over the rectangle of reference locations described
 * by the partition, i.e., calculate the partial dependencies
 * for each input direction
 */

double Particle::Relevance(double **rect, int *cat, bool approx, double *delta)
{
  zerov(delta, pall->m);
  return tree->Relevance(rect, cat, approx, delta);
}


/*
 * Entropy:
 *
 * accumulate the (un-normalized) Entropy statistics
 * for sequential design at the input locations
 * pall->X 
 */

void Particle::Entropy(double *entropy)
{
    tree->Entropy(entropy);
}


/*
 * Simple accessor functions follow 
 */

/* return the height of the tree in the particle */
int Particle::getHeight(void){ return tree->Height(); }

/* return the height of the tree in the particle */
int Particle::numLeaves(void){ return tree->numLeaves(); }

/* calculate the average size of the leaves of the tree */
double Particle::AvgSize(void){ return tree->leavesAvgSize(); }

/* calculate the average size of the leaves of the tree */
double Particle::AvgRetired(void){ return tree->leavesAvgRetired(); }

/* calculate the total number of data points under the 
   treed partition */
double Particle::getT(void){ return tree->leavesN(); }

/* print the contents of the particle (i.e., of the tree) */
void Particle::Print(int l){ tree->Print(); }


/*
 * main_effects:
 *
 * use nonparametric (move_avg) smoothing to estimate
 * main effects on a pre-defined Xgrid using the 
 * samples at nn LHS locations
 */

void main_effects(double **XX, unsigned int nn, unsigned int m, 
		  unsigned int aug, int *cat, double *ZZm, 
		  double **Xgrid_t, unsigned int ngrid, double span, 
		  double **main)
{
  /* for extracting each column */
  double *XXj = new_vector(nn);

  /* loop over columns */
  for(unsigned int j=0; j<m-aug; j++) {

    /* real-valued predictor */
    if(cat[j] == 0) {

      /* extract the j-th column of XX */
      for(unsigned int k=0; k<nn; k++) XXj[k] = XX[k][j+aug];

      /* nonparametric regression */
      move_avg(ngrid, Xgrid_t[j], main[j], nn, XXj, ZZm, span);

    } else { /* categorical variables */

      /* accumulate for each value of the binary predictor */
      unsigned int n0 = 0;
      double b0 = 0.0, b1 = 0.0;
      for(unsigned int k=0; k<nn; k++){
	if(XX[k][j+aug] == 0) { n0++; b0 += ZZm[k]; }
	else b1 += ZZm[k];
      }
      
      /* normalize; and fill for all grid locations */
      for(unsigned int i=0; i<ngrid; i++) {
	if(Xgrid_t[j][i] < 0.5) main[j][i] = b0/((double) n0);
	else main[j][i] = b1/((double) (nn-n0));
      }
    }
  }

  /* clean up */
  free(XXj);
}


/* 
 * move_avg:
 * 
 * simple moving average smoothing.  
 * Assumes that XX is already in ascending order!
 * Uses a squared difference weight function.
 */
 
void move_avg(int nn, double* XX, double *YY, int n, double* X, 
              double *Y, double frac)
{
  int q, i, j, l, u, search;
  double dist, range, sumW;
  double *Xo, *Yo, *w;
  int *o;
	
  /* frac is the portion of the data in the moving average
     window and q is the number of points in this window */
  assert( 0.0 < frac && frac < 1.0);
  q = (int) floor( frac*((double) n));
  if(q < 2) q=2;
  if(n < q) q=n;

  /* assume that XX is already in ascending order. 
   * put X in ascending order as well (and match Y) */
  Xo = new_vector(n);
  Yo = new_vector(n);
  o = order(X, n);
  for(i=0; i<n; i++) { Xo[i] = X[o[i]-1]; Yo[i] = Y[o[i]-1]; }
  free(o);

  /* window paramters */
  w = new_vector(n);  /* window weights */
  l = 0;              /* lower index of window */
  u = q-1;            /* upper index of the window */
  
  /* now slide the window along */
  for(i=0; i<nn; i++){

    /* find the next window */
    search=1;
    while(search){
      if(u==(n-1)) search = 0;
      else if( MYfmax(fabs(XX[i]-Xo[l+1]), fabs(XX[i]-Xo[u+1])) > 
               MYfmax(fabs(XX[i]-Xo[l]), fabs(XX[i]-Xo[u]))) search = 0;
      else{ l++; u++; }
    }

    /* width of the window in X-space */
    range = MYfmax(fabs(XX[i]-Xo[l]), fabs(XX[i]-Xo[u]));
    
    /* calculate the weights in the window; 
     * every weight outside the window will be zero */
    zerov(w,n);
    for(j=l; j<=u; j++){
      dist = fabs(XX[i]-Xo[j])/range;
      w[j] = (1.0-dist)*(1.0-dist);
    }

    /* record the (normalized) weighted average in the window */
    sumW = sumv(&(w[l]), q);
    YY[i] = vmult(&(w[l]), &(Yo[l]), q)/sumW;
  }
  
  /* clean up */
  free(w);
  free(Xo);
  free(Yo);
}


/*
 * sobol_indices:
 *
 * calculate the Sobol S and T indices using samples of the
 * posterior predictive distribution (ZZm and ZZvar) at 
 * nn*(d+2) locations
 */
 
void sobol_indices(double *ZZ, unsigned int nn, unsigned int m, 
		   double *S, double *T)
{
  /* pointers to responses for the original two LHSs */
  double *fM1 = ZZ; 
  double *fM2 = ZZ + nn;

  /* accumilate means and variances */
  double EZ, EZ2, Evar;
  Evar = EZ = EZ2 = 0.0;
  for(unsigned int j=0; j<nn; j++){
    EZ += fM1[j] + fM2[j];
    EZ2 += sq(fM1[j]) + sq(fM2[j]);
  }

  /* normalization for means and variances */
  double dnn = (double) nn;
  EZ = EZ/(dnn*2.0);
  EZ2 = EZ2/(dnn*2.0);
  Evar = Evar/(dnn*2.0);
  double sqEZ = sq(EZ);
  double lVZ = log(EZ2 - sqEZ); 
  
  /* fill S and T matrices */
  double ponent;
  double *fN;
  double U, Uminus;
  for(unsigned int k=0; k<m; k++) { /* for each column */
    
    /* accumulate U and Uminus for each k: the S and T dot products */
    fN = ZZ + (k+2)*nn;
    U = Uminus = 0.0;
    for(unsigned int j=0; j<nn; j++) {
      U += fM1[j]*fN[j];
      Uminus += fM2[j]*fN[j];
    }

    /* normalization for U and Uminus */
    U = U/(dnn - 1.0);
    Uminus = Uminus/(dnn - 1.0);
    
    /* now calculate S and T */
    ponent = U - sqEZ;
    if(ponent < 0.0) ponent = 0;
    S[k] = exp(log(ponent) - lVZ); 
    ponent = Uminus - sqEZ;
    if(ponent < 0.0) ponent = 0;
    T[k] = 1 - exp(log(ponent) - lVZ);
  }
}
