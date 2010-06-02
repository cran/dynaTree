/*  Written by Matt Taddy & Bobby Gramacy */

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
  

 
