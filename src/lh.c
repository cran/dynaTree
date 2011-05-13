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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <Rmath.h>
#include "lh.h"
#include "matrix.h"
#include "rhelp.h"

int compareRank(const void* a, const void* b);
int compareDouble(const void* a, const void* b);


/*
 * structure for ranking
 */

typedef struct rank
{
  double s;
  int r;
} Rank;


/*
 * rect_sample_lh:
 *
 * returns a unidorm sample of (n) points
 * within a regular (dim)-dimensional cube.
 * (n*dim matrix returned)
 */

double** rect_sample(int dim, int n)
{
  int i,j;
  double **s = new_matrix(dim, n);
  for(i=0; i<dim; i++) {
    for(j=0; j<n; j++) {
      s[i][j] = unif_rand();
    }
  }

  return s;
}


/*
 * sens_boot:
 *
 * use a Bootstrap to create a sample of M1 and M2 matrices
 * from the rows of X for a sensitivity analysis
 */

double ** sens_boot(int nn_boot, int d, int aug, int *nn, 
		    double **X, int n)
{
  double **M1, **M2;

  /* get the two latin hypercube samples */
  M1 = boot_sample(d-aug, aug, nn_boot, X, n);
  assert(M1);
  M2 = boot_sample(d-aug, aug, nn_boot, X, n);
  assert(M2);
  
  /* create one big XX matrix from M1 and M2 */
  return Ms_to_XX(nn_boot, d, aug, M1, M2, nn);
}


/*
 * Ms_to_XX:
 *
 * create one big XX matrix from M1 and M2 as created
 * by either of sens_lhs or sens_boot
 */

double **Ms_to_XX(unsigned int nns, int d, int aug, 
		  double **M1, double **M2, int *nn)
{
  double **XX, **XX2;
  int i, j;

  /* allocate space for the new matrix */
  /* full size of XX matrix */
  *nn = nns * (d-aug + 2);
  XX = new_matrix(*nn, d-aug);
  assert(XX);

  /* copy the first M1 and M2 into XX */
  dup_matrix(XX, M1, nns, d-aug);
  dupv(XX[nns], M2[0], nns*(d-aug));

  /* now do d copies of M2 into XX with the appropriate Nj
     replacement */
  for(j=0; j<d-aug; j++) {
    dup_matrix(&XX[nns*(2+j)], M2, nns, d-aug);
    for(i=0; i<nns; i++) XX[nns*(2+j)+i][j] = M1[i][j];
  }

  /* clean up M1 & M2 matrices */
  delete_matrix(M1);
  delete_matrix(M2);

  /* now augment with a column of ones if needed */
  if(aug > 0) {
    assert(aug == 1);
    XX2 = new_matrix(*nn, d);
    for(i=0; i<*nn; i++) {
      XX2[i][0] = 1.0;
      for(j=0; j<d-aug; j++) XX2[i][j+1] = XX[i][j];
    }
    delete_matrix(XX);
    XX = XX2;
  }

  /* return new XX matrix */
  return(XX);
}


/*
 * boot_sample:
 *
 * create a new matrix, M, filled with a bootstrap sample
 * of size nn from the n x d matrix X;
 */

double ** boot_sample(int d, int aug, int nn, double **X, int n)
{
  int i, indx;
  double **M;

  /* allocate new matrix for bootstrap */
  M = new_matrix(nn, d);

  /* each row of nre matrix */
  for(i=0; i<nn; i++) {
    /* sample a new index  uniformly */
    indx = (int) (unif_rand()*((double) n));
    dupv(M[i], X[indx]+aug, d);
  }
  
  return(M);
}


/*
 * sens_lhs:
 *
 * use a Latin Hypercube to create a sample of M1 and M2 matrices
 * from the rows of X for a sensitivity analysis
 */

double ** sens_lhs(int nn_lhs, int d, int aug, double **bnds, 
		   double *shape, double *mode, int *nn)
{
  double **M1, **M2;

  /* get the two latin hypercube samples */
  M1 = beta_sample_lh(d-aug, nn_lhs, bnds, shape, mode);
  assert(M1);
  M2 = beta_sample_lh(d-aug, nn_lhs, bnds, shape, mode);
  assert(M2);

  /* create one big XX matrix from M1 and M2 */
  return Ms_to_XX(nn_lhs, d, aug, M1, M2, nn);
}


/*
 * beta_sample_lh:
 *
 * returns a latin hypercube sample of (n) points
 * within a regular (dim)-dimensional cube, proportional
 * to independant scaled beta distributions over the cube,
 * with specified modes and shape parameters.
 */

double** beta_sample_lh(int dim, int n, double** rect, double* shape, 
			double* mode)
{
  int i,j;
  double **z, **s, **zout;
  double** e;
  int **r;
  Rank ** sr;

  double alpha, mscaled;
  
  assert(n >= 0);
  if(n == 0) return NULL;
  z = e = s = NULL;
  
  /* We could just draw random permutations of (1..n) here, 
     which is effectively what we are doing. 
     This ranking scheme could be valuable, though, 
     in drawing lhs for correlated variables.
     In that case, s would instead be a sample from the correct
     joint distribution, and the quantile functions at the end
     would have to correspond to the marginal distributions
     for each variable.  See Stein, 1987 (Technometrics).
     This would have to be coded on a case to case basis though. */

  /* get initial sample */
  s = rect_sample(dim, n);
  
  /* get ranks */
  r = (int**) malloc(sizeof(int*) * dim);
  for(i=0; i<dim; i++) {
    sr = (Rank**) malloc(sizeof(Rank*) * n);
    r[i] = new_ivector(n);
    for(j=0; j<n; j++) {
      sr[j] = (Rank*) malloc(sizeof(Rank));
      sr[j]->s = s[i][j];
      sr[j]->r = j;
    }
    
    qsort((void*)sr, n, sizeof(Rank*), compareRank);
    
    /* assign ranks	*/
    for(j=0; j<n; j++) {
      r[i][sr[j]->r] = j+1;
      free(sr[j]);
    }
    free(sr);
  }
  
  /* Draw random variates */
  e = rect_sample(dim, n);
  /* Obtain latin hypercube sample on the unit cube:
   The alpha parameters for each beta quantile function are calculated
   from the (re-scaled) mode and the shape parameter.  */
  z = new_matrix(dim,n);
  for(i=0; i<dim; i++) {

    if(shape[i]==0){ /* for binary variables, draw 0-1. */
      if(mode==NULL || mode[i] > 1.0 || mode[i] < 0) mscaled=0.5;
      else mscaled = mode[i];
      for(j=0; j<n; j++){
	z[i][j] = 0.0;
	if(unif_rand() < mscaled) z[i][j] = 1.0; 
      }
      free(r[i]);
      continue;
    }

    if(mode==NULL) mscaled = 0.5;
    else mscaled = (mode[i]-rect[0][i])/(rect[1][i] - rect[0][i]);
    if( 0 > mscaled || 1 < mscaled ) mscaled=0.5;
    if(shape[i] < 1) shape[i] = 1; /* only concave betas, else uniform */
    alpha = (1 + mscaled*(shape[i]-2))/(1-mscaled);
    assert( alpha > 0 );
    for(j=0; j<n; j++) {
      z[i][j] = qbeta( ( ((double)r[i][j]) - e[i][j]) / n, alpha, shape[i], 1, 0);
    }
    free(r[i]);
  }
    /* Shift and scale from the unit cube to rect */
  rect_scale(z, dim, n, rect);

  /* Wrap up */
  free(r);
  delete_matrix(s);
  delete_matrix(e);
  zout = new_t_matrix(z, dim, n);
  delete_matrix(z);
  return zout;
}

/*
 * rect_sample_lh:
 *
 * returns a (uniform) latin hypercube sample of (n) points
 * within a regular (dim)-dimensional cube.
 */

double** rect_sample_lh(int dim, int n, double** rect, int er)
{
  int i,j;
  double **z, **s, **zout;
  double** e;
  int **r;
  Rank ** sr;
  
  assert(n >= 0);
  if(n == 0) return NULL;
  z = e = s = NULL;
  
  /* get initial sample */
  s = rect_sample(dim, n);
  
  /* get ranks */
  r = (int**) malloc(sizeof(int*) * dim);
  for(i=0; i<dim; i++) {
    sr = (Rank**) malloc(sizeof(Rank*) * n);
    r[i] = new_ivector(n);
    for(j=0; j<n; j++) {
      sr[j] = (Rank*) malloc(sizeof(Rank));
      sr[j]->s = s[i][j];
      sr[j]->r = j;
    }
    
    qsort((void*)sr, n, sizeof(Rank*), compareRank);
    
    /* assign ranks	*/
    for(j=0; j<n; j++) {
      r[i][sr[j]->r] = j+1;
      free(sr[j]);
    }
    free(sr);
  }
  
  /* Draw random variates */
  if(er) e = rect_sample(dim, n);
  
  /* Obtain latin hypercube sample */
  z = new_matrix(dim,n);
  for(i=0; i<dim; i++) {
    for(j=0; j<n; j++) {
      if(er) z[i][j] = (r[i][j] - e[i][j]) / n;
      else z[i][j] = (double)r[i][j] / n;
    }
    free(r[i]);
  }
  
  /* Wrap up */
  free(r);
  delete_matrix(s);
  if(er) delete_matrix(e);
  
  rect_scale(z, dim, n, rect);
  
  zout = new_t_matrix(z, dim, n);
  delete_matrix(z);
  
  return zout;
}


/*
 * compareRank:
 *
 * comparison function for ranking
 */

int compareRank(const void* a, const void* b)
{
  Rank* aa = (Rank*)(*(Rank **)a); 
  Rank* bb = (Rank*)(*(Rank **)b); 
  if(aa->s < bb->s) return -1;
  else return 1;
}


/*
 * compareDouble:
 *
 * comparison function double sorting ranking
 */

int compareDouble(const void* a, const void* b)
{
  double aa = (double)(*(double *)a); 
  double bb = (double)(*(double *)b); 
  if(aa < bb) return -1;
  else return 1;
}


/*
 * rect_scale:
 *
 * shift/scale a draws from a unit cube into 
 * the specified rectangle
 */ 

void rect_scale(double** z, int d, int n, double** rect)
{
  int i,j;
  double scale, shift;
  for(i=0; i<d; i++) {
    scale = rect[1][i] - rect[0][i];
    shift = rect[0][i];
    for(j=0; j<n; j++) {
      z[i][j] = z[i][j]*scale + shift;
    }
  }
}


/*
 * printRect:
 *
 * print 2xd double rectangle to (FILE* outfile)
 */

void printRect(FILE* outfile, int d, double** rect)
{
  int j,i;
  for(j=0; j<2; j++) {
    for(i=0; i<d; i++) {
      myprintf(outfile, " %5.4g", rect[j][i]);
    }
    myprintf(outfile, "\n");
  }
}


/*
 * errorBadRect:
 *
 * Bad rectangle (argv[2]) error message
 * uses printUsage();;
 */

void errorBadRect(void)
{
  error("bad rectangle format"); 
}


/*
 * sortDouble:
 *
 * sort an array of doubles
 */

void sortDouble(double *s, unsigned int n)
{
  qsort((void*)s, n, sizeof(double), compareDouble);
}


/*
 * order:
 *
 * obtain the integer order of the indices of s
 * from least to greatest.  the returned indices o
 * applied to s, (e.g. s[o]) would resort in a sorted list
 */

int* order(double *s, unsigned int n)
{
  int j;
  int *r;
  Rank ** sr;
  
  r = new_ivector(n);
  sr = (Rank**) malloc(sizeof(Rank*) * n);
  for(j=0; j<n; j++) {
    sr[j] = (Rank*) malloc(sizeof(Rank));
    sr[j]->s = s[j];
    sr[j]->r = j;
  }
  
  qsort((void*)sr, n, sizeof(Rank*), compareRank);
  
  /* assign ranks */
  for(j=0; j<n; j++) {
    r[j] = sr[j]->r +1;
    free(sr[j]);
  }
  free(sr);
  
  return r;
}


/*
 * rank:
 *
 * obtain the integer rank of the elemts of s
 */

int* rank(double *s, unsigned int n)
{
  int j;
  int *r;
  Rank ** sr;
  
  r = new_ivector(n);
  sr = (Rank**) malloc(sizeof(Rank*) * n);
  for(j=0; j<n; j++) {
    sr[j] = (Rank*) malloc(sizeof(Rank));
    sr[j]->s = s[j];
    sr[j]->r = j;
  }
  
  qsort((void*)sr, n, sizeof(Rank*), compareRank);
  
  /* assign ranks */
  for(j=0; j<n; j++) {
    r[sr[j]->r] = j+1;
    free(sr[j]);
  }
  free(sr);
  
  return r;
}
