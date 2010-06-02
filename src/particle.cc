extern "C" {
#include "matrix.h"
#include "rhelp.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <Rmath.h>
}
#include "particle.h"


/*
 * Particle:
 * 
 * constructor for new data
 */

Particle::Particle(Pall *pall, int *pstart, unsigned int nstart)
{
  assert(pall->minp > 0);
  int *p = new_dup_ivector(pstart, nstart);
  double **rect = get_data_rect(pall->X, nstart, pall->m);
  tree = new Tree(pall, p, nstart, rect, NULL);
}


/*
 * Particle:
 *
 * copy constructor for existing particle
 */

Particle::Particle(const Particle &p) 
{
  tree = new Tree(p.tree, NULL);
}


/*
 * Particle:
 *
 * pointer-based copy constructor for existing
 * particle#
 */

Particle::Particle(Particle *p) 
{
  tree = new Tree(p->tree, NULL);
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
 * Particle =:
 *
 * copy constructor overloading the equals operator 
 */

Particle& Particle::operator=(const Particle &p)
{
  tree = new Tree(p.tree, NULL);
  return *this;
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

  /* calculate the probability of growing to a new tree
     by splitting at the leaf node that pall->X[index] 
     belongs to */
  int var;
  double val;
  double pgrow = leaf->growProb(&var, &val);
 
  /* calculate the probability of pruning back from the
     parent of the leaf node that pall->X[index] belongs
     to */
  double pprune = leaf->pruneProb();

  /* normalize the probabilities of the three (two) moves */
  double norm = pstay + pprune + pgrow;
  if(isnan(norm) || norm==0.0) return; 
  pstay = pstay/norm;
  pprune = pprune/norm;

  /* choose randomly between the three moves */
  double u = unif_rand();

  /* do nothing if staying */
  if(u < pstay) return;

  /* if not prune then automatically grow */
  if(u < pstay + pprune) leaf->Parent()->prune();
  else leaf->grow(var, val);
}


/*
 * PostPred:
 *
 * calculate the posterior probability of a particular
 * (x.y) pair under the tree model in the particle --
 * used in the PL resample ste
 */

double Particle::PostPred(double *xx, double yy)
{
  assert(tree != NULL);
  Tree* leaf = tree->GetLeaf(xx);
  double post = leaf->PostPred(xx,yy);
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

void Particle::Predict(double **XX, unsigned int nn, 
		       double *mean, double *sd, double *df, 
		       double *var, double *q1, double *q2)
{
  /* dummy pointers for s2 and df */
  double sdi, dfi;

  for(unsigned int i=0; i<nn; i++) {
    double *q1i, *q2i;

    /* pointers for quantiles if needed */
    if(q1) { q1i = &(q1[i]); q2i = &(q2[i]); }
    else q1i = q2i = NULL;

    /* actually predict */
    tree->Predict(XX[i], &mean[i], &sdi, &dfi, &var[i], q1i, q2i);

    /* save s2 and df if allocated */
    if(sd) {
      assert(df);
      sd[i] = sdi; df[i] = dfi;
    }
    
  }
}


/*
 * Predict:
 *
 * classification based prediction
 */

void Particle::Predict(double **XX, unsigned int nn, double **p,
		       double *entropy)
{
  double *ptemp = new_vector(tree->pall->nc);
  for(unsigned int i=0; i<nn; i++) {
    tree->Predict(XX[i], ptemp);
    entropy[i] = 0.0;
    for(unsigned int j=0; j<tree->pall->nc; j++) {
      p[j][i] = ptemp[j];
      entropy[i] += 0.0 - ptemp[j] * log(ptemp[j]);
    }
  }
  free(ptemp);
}


/*
 * ALC:
 *
 * compute the (un-normalized) Active Learning Cohn
 * statistic for sequential design
 */

void Particle::ALC(double **XX, unsigned int nn, double **alc) 
{
  for(unsigned int i=0; i<nn; i++) {
    for(unsigned int j=0; j<nn; j++) {
      alc[i][j] += tree->ALC(XX[i], XX[j]);
    }
  }
}



/*
 * Simple accessor functions follow 
 */

/* return the height of the tree in the particle */
int Particle::getHeight(void){ return tree->Height(); }

/* calculate the average size of the leaves of the tree */
double Particle::AvgSize(void){ return tree->leavesAvgSize(); }

/* calculate the total number of data points under the 
   treed partition */
double Particle::getT(void){ return tree->leavesN(); }

/* print the contents of the particle (i.e., of the tree) */
void Particle::Print(int l){ tree->Print(); }
