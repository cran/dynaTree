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
#include <assert.h>
#include "pall.h"
#include "rhelp.h"
#include "Rmath.h"
#include "matrix.h"

/*
 * new_pall:
 *
 * create a new pall structure with (x,y) pairs and
 * bounding rectangle; the pointers to X
 */

Pall *new_pall(double **X, unsigned int n, unsigned int m, 
	       int *Xna, int **XNA, unsigned int nna, 
	       double *y, double *params, int model_in)
{
  /* copy parameters common to all particles */
  unsigned int w;
  Pall *pall;
  pall = (Pall*) malloc(sizeof(struct pall));
  pall->X = new_dup_matrix(X, n, m);
  pall->n = n;
  pall->g = 0;
  pall->m = m;
  pall->y = new_dup_vector(y, n);
  pall->nu0 = params[0];
  pall->s20 = params[1];
  pall->a = params[2];
  pall->b = params[3];
  pall->minp = (unsigned int) params[4];
  pall->smin = (unsigned int) params[5] - 1;
  pall->bmax = (unsigned int) params[6];

  /* deal with NAs */
  pall->nna = nna;
  if(nna > 0) {
    pall->Xna = new_dup_ivector(Xna, n);
    pall->XNA = new_dup_imatrix(XNA, nna, m);
  } else { pall->Xna = NULL; pall->XNA = NULL; }

  /* determine the growing rectangle proposal typr */
  switch((int) params[8]) {
  case 1: pall->rprop = LUALL; break;
  case 2: pall->rprop = LUVAR; break;
  case 3: pall->rprop = REJECT; break;
  default: error("rprop %d not defined\n", (int) params[8]);
  }
  pall->icept = (unsigned int) params[7];
  pall->nc = 0;

  /* determine the model */
  if(model_in == 1) pall->model = CONSTANT;
  else if(model_in == 2) pall->model = LINEAR;
  else if(model_in == 3) {
    pall->model = CLASS;
    pall->nc = (unsigned int) max(y, n, &w) + 1;
  } else if(model_in == 4) {
    pall->model = PRIOR;
  }
  else {
    warning("model %d not defined, using constant\n", model_in);
    pall->model = CONSTANT;
  }

  /* allocate temporary vector of length bmax in LINEAR model */
  if(pall->model == LINEAR) pall->bmaxv = new_vector(pall->bmax);
  else pall->bmaxv = NULL;

  return pall;
}


/*
 * copy_pall:
 *
 * copy an existing pall
 */

Pall *copy_pall(Pall *pold)
{
  Pall *pall;
  assert(pold);
  pall = (Pall*) malloc(sizeof(struct pall));
  pall->X = new_dup_matrix(pold->X, pold->n, pold->m);
  pall->n = pold->n;
  pall->g = pold->g;
  pall->m = pold->m;
  pall->y = new_dup_vector(pold->y, pold->n);
  pall->nna = pold->nna;
  if(pold->Xna) pall->Xna = new_dup_ivector(pold->Xna, pold->n);
  else pall->Xna = NULL;
  if(pold->XNA) pall->XNA = new_dup_imatrix(pold->XNA, pold->nna, pold->m);
  else pall->XNA = NULL;
  pall->nu0 = pold->nu0;
  pall->s20 = pold->s20;
  pall->a = pold->a;
  pall->b = pold->b;
  pall->minp = pold->minp;
  pall->smin = pold->smin;
  pall->bmax = pold->bmax;
  pall->icept = pold->icept;
  pall->rprop = pold->rprop;
  pall->icept = pold->icept;
  pall->nc = pold->nc;
  pall->model = pold->model;
  if(pold->bmaxv) pall->bmaxv = new_dup_vector(pold->bmaxv, pold->bmax);
  else pall->bmaxv = NULL;
  return pall;
}


/*
 * delete_pall:
 *
 * delete the pall structure and its components 
 */

void delete_pall(Pall *pall)
{
  delete_matrix(pall->X);
  if(pall->Xna) free(pall->Xna);
  if(pall->XNA) delete_imatrix(pall->XNA);
  free(pall->y);
  if(pall->bmaxv) free(pall->bmaxv);
  free(pall);
}


/*
 * add_data:
 *
 * add new rows to X and components to y
 */

void add_data(Pall *pall, double **X, unsigned int n, 
	      int *Xna, int **XNA, unsigned int nna, double *y)
{
  unsigned int i, nnew;

  /* sanity check for classification */
  if(pall->model == CLASS)
    for(i=0; i<n; i++)
      assert(y[i] < pall->nc);

  /* new data size */
  nnew = pall->n + n;

  /* make new bigger matrices and vectors */
  pall->X = new_bigger_matrix(pall->X, pall->n, pall->m, nnew, pall->m);
  pall->y = realloc(pall->y, sizeof(double)*nnew);

  /* fillin the new X bits */
  for(i=pall->n; i<nnew; i++) 
    dupv(pall->X[i], X[i-pall->n], pall->m);

  /* and now the new y bits */
  dupv(pall->y+pall->n, y, n);

  /* deal with NAs */
  if(pall->Xna || Xna) {
    int newXna = pall->Xna == NULL;
    pall->Xna = realloc(pall->Xna, sizeof(int)*nnew);
    if(newXna) for(unsigned int i=0; i<pall->n; i++) pall->Xna[i] = -1;
    if(!Xna) for(unsigned int i=pall->n; i<nnew; i++) pall->Xna[i] = -1;
    else { /* augment new Xnas >= 0 by previous nna count */
      for(unsigned int i=0; i<n; i++) {
	pall->Xna[i+pall->n] = Xna[i];
	if(Xna[i] != -1) pall->Xna[i+pall->n] += pall->nna;
      }
    }
    /* copy new NA indicator rows if there are new NAs */
    if(XNA) {
      pall->XNA = new_bigger_imatrix(pall->XNA, pall->nna + nna, 
					   pall->m, nna, pall->m);
      pall->nna += nna;
    } assert(nna == 0);

    /* BOBBY: delete this later */
    printIVector(pall->Xna, nnew, MYstdout);
    printIMatrix(pall->XNA, pall->nna, pall->m, MYstdout);
  }

  /*  update n */
  pall->n = nnew;
}


/*
 * retire:
 *
 * retire/remove the index(ed) entry of the X and y
 * elements of the pall structure 
 */

void retire(Pall *pall, unsigned int index)
{
  double **Xnew;
  int **XNAnew;
  unsigned int k;
  int xna;

  /* sanity check */
  assert(index < pall->n);

  /* remove from X and Y and al (if used) */
  (pall->n)--;
  if(index < (int) pall->n) {
    dupv(pall->X[index], pall->X[pall->n], pall->m);
    pall->y[index] = pall->y[pall->n];
  }

  /* reduce the size of y and X */
  pall->y = (double*) realloc(pall->y, sizeof(double)*pall->n);
  Xnew = (double**) malloc(sizeof(double*) * pall->n);
  Xnew[0] = (double*) realloc(pall->X[0], sizeof(double) * pall->n * pall->m);
  free(pall->X);
  for(k=1; k<pall->n; k++) Xnew[k] = Xnew[k-1] + pall->m;
  pall->X = Xnew;

  /* deal with NA entries */
  if(pall->Xna) {
    xna = pall->Xna[index];
    /* maybe re-size XNA */
    if(xna >= 0) {
      (pall->nna)--;
      if(pall->nna == 0) { /* no more NAs left */
	delete_imatrix(pall->XNA); pall->XNA = NULL; 
      } else {
	if(xna < (int) pall->nna) /* move last into xna position */
	  dupiv(pall->XNA[xna], pall->XNA[pall->nna], pall->m);
	/* reduce the size of XNA */
	/* BOBBY: NEED FUNCTION FOR THIS SINCE THIS IS DUPLICATED 
	   FROM A FEW LINES ABOVE */
	XNAnew = (int**) malloc(sizeof(int*) * pall->nna);
	XNAnew[0] = (int*) realloc(pall->XNA[0], sizeof(int) * pall->nna * pall->m);
	free(pall->XNA);
	for(k=1; k<pall->nna; k++) XNAnew[k] = XNAnew[k-1] + pall->m;
	pall->XNA = XNAnew;
      }
    }
    
    /* resize and reset the highest number in pall->Xna */
    if(pall->nna > 0) {
      if(index < (int) pall->n) pall->Xna[index] = pall->Xna[pall->n];
      pall->Xna = (int*) realloc(pall->Xna, sizeof(int)*pall->n);
      for(k=0; k<pall->n; k++) 
	if(pall->Xna[k] == pall->nna) pall->Xna[k] = xna;
    }
  }

  /* keep track of the fact that we've removed an entry */
  (pall->g)++;
}


/*
 * reorder:
 *
 * reorder the X-y pairs 
 */

void reorder(Pall *pall, int *o)
{
  unsigned int i;
  double **Xnew;
  double *ynew;
  int *Xnanew;

  /* sanity check */
  assert(o);

  /* allocate new space */
  Xnew = new_matrix(pall->n, pall->m);
  ynew = new_vector(pall->n);
  if(pall->Xna) Xnanew = new_ivector(pall->n);
  else Xnanew = NULL;

  /* fill that space */
  for(i=0; i<pall->n; i++) {
    dupv(Xnew[o[i]], pall->X[i], pall->m);
    ynew[o[i]] = pall->y[i];
    if(Xnanew) Xnanew[o[i]] = pall->Xna[i];
  }

  /* free old space and copy pointer */
  delete_matrix(pall->X); pall->X = Xnew;
  free(pall->y); pall->y = ynew;
  if(Xnanew) { free(pall->Xna); pall->Xna = Xnanew; }
}


/*
 * EI:
 *
 * calculates the expected improvement following
 * Williams et al by integrating over the parameters
 * to the GP predictive
 */

double EI(const double m, const double sd, const double df,
          const double fmin)
{
  double diff, diffs, scale, ei;

  diff = fmin - m;
  diffs = diff/sd;
  scale = (df*sd + sq(diff)/sd)/(df-1.0);
  ei = diff*pt(diffs, (double) df, 1, 0);
  ei += scale*dt(diffs, (double) df, 0);

  return(ei);
}
