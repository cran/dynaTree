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


#ifndef __CLOUD_H__
#define __CLOUD_H__

#include "particle.h"
#include "matrix.h"

using namespace std;

class Cloud
{
 private:
  
  Particle** particle;         /* the particle cloud */

  int *index;                  /* indices for re-sampling particles */
  double *prob;                /* probabilities for re-sampling */

 protected:

  unsigned int Resample(void);

 public:

  Pall *pall;                  /* data common to each particle */
  unsigned int N;              /* number of particles */
  unsigned int Nrevert;        /* remember number of particles */
  
  /* constructor and destructor */
  Cloud(unsigned int N, Pall *pall, int *pstart, unsigned int nstart);
  Cloud(Cloud *cold);
  ~Cloud(void);

  /* rejuvination functions */
  void Reorder(int *o);
  void Combine(Cloud *c2);

  /* Particle Learning steps */
  double Resample(unsigned int t, unsigned int verb);
  void Propagate(unsigned int t);
  double Posterior(void);

  /* prediction and sequential design */
  void Predict(double **XX, double *yy, unsigned int nn, double *mean,
	       double *vmean, double *var, double *df, double *q1, double *q2, 
         double *yypred, double *ei, unsigned int verb);
  void IECI(double **XX, unsigned int nn, double **Xref, 
	    unsigned int nref, double **probs, double *ieci_out, 
	    unsigned int verb);
  void ALC(double **XX, unsigned int nn, double **Xref, 
	   unsigned int nref, double **probs, double *alc_out, 
	   unsigned int verb);
  void ALC(double **XX, unsigned int nn, double **rect, int *cat,
	   bool approx, double *alc_out, unsigned int verb);
  void ALC(double **rect, int *cat, bool approx, double *alc_out, 
	   unsigned int verb);
  void Sens(int *cls, unsigned int nns, unsigned int aug, 
	    double **rect, double *shape, double *mode, int *cat,
	    double **Xgrid, unsigned int ngrid, double span, 
	    double **mean, double **q1, double **q2, 
	    double **S, double **T, unsigned int verb);
  void Relevance(double **rect, int *cat, bool approx, double **delta, 
	      unsigned int verb);
  void qEntropy(double q, double **XX, unsigned int nn, double *qentropy, 
		unsigned int verb);
  void qEI(double q, double alpha, double **XX, unsigned int nn, 
	   double *qei, unsigned int verb);

  /* retiring particles for online learning */
  void Retire(int *pretire, unsigned int nretire, double lambda, 
	      unsigned int verb);

  /* variable selection */
  void VarPropUse(double *counts);
  void VarPropTotal(double *counts);
  void Intervals(unsigned int index, unsigned int var, double *a, double *b);

  /* prediction for classification */
  void Predict(double **XX, int *yy, unsigned int nn, double **p,
         double *yypred, double *entropy, unsigned int verb);
  void Coef(double **XX, unsigned int nn, double **beta, unsigned int verb);
  void Predict(unsigned int cl, double **XX, unsigned int nn, double **p,
	       double **c, unsigned int verb);
  void Entropy(double *entropy_out, unsigned int verb);

  /* information extracting about particles */
  void TreeStats(double *height, double *leaves, double *avgsize, double *avgretire);
  void SameLeaf(double **X, unsigned int n, int *counts);
};


/* utility functions for the ::Resample step */
double norm_weights(double *v, int n);
void indexsample(int *ind, int n, int num_probs, double *probs);
void ressample(int *ind, int n, int num_probs, double *probs);

#endif
