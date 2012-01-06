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


#ifndef __PALL_H__
#define __PALL_H__ 

#define INF 1e300*1e300

/* types of models supported */
typedef enum Model {CONSTANT=1001, LINEAR=1002, CLASS=1003} Model;
typedef enum RProp {LUALL=1001, LUVAR=1002, REJECT=1003} RProp;

typedef struct pall {
  double **X;              /* input matrix n x m */
  double *y;               /* n-vector of outputs or class labels */
  unsigned int n;          /* number of (x,y) pairs */
  unsigned int g;          /* number of retired inputs */
  unsigned int m;          /* number of predictors; cols of X */
  unsigned int nna;        /* number of NA rows in X */
  int *Xna;                /* NA index matrix; (-1) for row not NA */
  int **XNA;               /* col NA indicator; nna x m matrix */
  unsigned int nc;         /* number of class labels */
  double nu0;              /* variance prior sample size */
  double s20;              /* prior variance estimate */
  double a;                /* tree prior alpha parameter */
  double b;                /* tree prior beta parameter */
  unsigned int smin;       /* tree prior splitmin index for linear models */
  unsigned int bmax;       /* tree prior basemax index for linear models */
  double *bmaxv;           /* temporary vector of length bmax */
  unsigned int icept;      /* Boolean intercept indicator for linear models */
  unsigned int minp;       /* tree prior min-part parameter */
  Model model;             /* model indicator */
  RProp rprop;             /* grow rectangle proposal inticator */
} Pall;

Pall *new_pall(double **X, unsigned int n, unsigned int m, 
	       int *Xna, int **XNA, unsigned int nna,
	       double *y, double *params, int model);
Pall *copy_pall(Pall *pold);
void add_data(Pall *data, double **X, unsigned int n, 
	      int *Xna, int **XNA, unsigned int nna, double *y);
void retire(Pall *data, unsigned int index);
void delete_pall(Pall *data);
void reorder(Pall *pall, int *o);
double EI(const double m, const double sd, const double df,
	const double fmin);

#endif
