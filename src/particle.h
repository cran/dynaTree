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

  /* constructors, destructors and copying */
  Particle(Pall *pall, int *pstart, unsigned int nstart);
  Particle(const Particle &p);
  Particle(Particle *p);
  ~Particle(void);
  Particle& operator=(const Particle &p);

  /* Particle Learning */
  void Propagate(unsigned int index);
  double Posterior(void);
 
  /* access tree characteristics */
  int getHeight();
  double AvgSize();
  double getT();
  void Print(int l);

  /* prediction */
  double PostPred(double *xx, double yy);
  void Predict(double **XX, unsigned int nn, double *mean, 
	       double *sd, double *df, double *var, double *q1, 
	       double *q2);
  void ALC(double **XX, unsigned int nn, double **alc);

  /* prediction for classification */
  void Predict(double **XX, unsigned int nn, double **p, 
	       double *entropy);
};


#endif
  

 
