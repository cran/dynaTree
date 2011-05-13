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


#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "tree.h"
#include "matrix.h"

class Particle
{
 private:

  Tree *tree;             /* pointer to the tree */
  
 public:
  
  Pall *pall;		  /* holding stuff common to all particles */

  /* constructors, destructors and copying */
  Particle(Pall *pall, int *pstart, unsigned int nstart);
  Particle(Particle *p);
  Particle(Particle *p, Pall *pall);
  ~Particle(void);

  /* rejuvination */
  void Reorder(int *o);

  /* Particle Learning */
  void Propagate(unsigned int index);
  double Posterior(void);
 
  /* access tree characteristics */
  int getHeight();
  double AvgSize();
  double AvgRetired();
  double getT();
  void Print(int l);

  /* prediction */
  double PostPred(double *xx, double yy);
  void Predict(double **XX, double *yy, unsigned int nn, double *mean, 
	       double *sd, double *df, double *var, double *q1, 
	       double *q2, double *yypred, double *ZZ);
  void EImECI(double **XX, unsigned int nn, double **Xref, 
	      unsigned int nref, double *probs, double **eimeci);
  void ALC(double **XX, unsigned int nn, double **Xref, unsigned int nref, 
	   double *probs, double **alc);
  void ALC(double **XX, unsigned int nn, double **rect, double *alc);
  void ALC(double **rect, double *alc);
  void Sens(unsigned int nns, unsigned int aug, double **rect, double *shape, 
	    double *mode, int *cat, double **Xgrid, unsigned int ngrid, 
	    double span, double **main, double *S, double *T);

  /* retire particles for online learning */
  void Retire(unsigned int index, double lambda);

  /* prediction for classification */
  void Predict(double **XX, int *yy, unsigned int nn, double **p, 
	       double *yypred, double *entropy);
  void Predict(unsigned int cls, double **XX, unsigned int nn,  
	       double *p, double *ZZ);
  void Sens(unsigned int cls, unsigned int nns, unsigned int aug, 
	    double **rect, double *shape, double *mode, int *cat, 
	    double **Xgrid, unsigned int ngrid, double span, 
	    double **main, double *S, double *T);
  void Entropy(double *entropy);

  /* variable selection */
  void VarCountUse(int *c);
  void VarCountTotal(double *c);
};

void main_effects(double **XX, unsigned int nn, unsigned int m, 
		  unsigned int adj, int *cat, double *ZZm, 
		  double **Xgrid_t, unsigned int ngrid, double span, 
		  double **main);
void move_avg(int nn, double* XX, double *YY, int n, double* X, 
              double *Y, double frac);
void sobol_indices(double *ZZ, unsigned int nn, unsigned int m, 
		   double *S, double *T);

#endif
