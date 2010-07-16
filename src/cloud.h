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
  unsigned int rsmin;          /* min number of re-sampled particles */

 protected:

  unsigned int Resample(void);

 public:

  Pall *pall;                  /* data common to each particle */
  unsigned int N;              /* number of particles */
  
  /* constructor and destructor */
  Cloud(unsigned int N, Pall *pall, int *pstart, unsigned int nstart);
  ~Cloud(void);

  /* Particle Learning steps */
  double Resample(unsigned int t, unsigned int verb);
  void Propagate(unsigned int t);
  double Posterior(void);

  /* prediction and sequential design */
  void Predict(double **XX, unsigned int nn, double **mean,
	       double **sd, double **df, double **var, double **q1, 
	       double **q2, unsigned int verb);
  void ALC(double **XX, unsigned int nn, double *alc_out, 
		unsigned int verb);

  /* prediction for classification */
  void Predict(double **XX, unsigned int nn, double ***p,
	       double **entropy, unsigned int verb);
  /* access */
  unsigned int GetRsmin(void);

};


/* utility functions for the ::Resample step */
double norm_weights(double *v, int n);
int indexsample(int *ind, int n, int num_probs, double *probs);

#endif
