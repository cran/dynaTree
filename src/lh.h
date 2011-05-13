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


#ifndef __LH_H__
#define __LH_H__

#include <stdio.h>

double** sens_lhs(int nn_lhs, int d, int aug, double **bnds, 
		     double *shape, double *mode, int *nn);
double ** sens_boot(int nn_lhs, int d, int aug, int *nn, double **X, 
		    int n);
double **Ms_to_XX(unsigned int nns, int d, int aug, 
		  double **M1, double **M2, int *nn);
double** rect_sample(int dim, int n);
double** rect_sample_lh(int dim, int n, double** rect, int er);
double** beta_sample_lh(int dim, int n, double** rect, double* shape, 
			double* mode);
double ** boot_sample(int d, int aug, int nn, double **X, int n);
void rect_scale(double** z, int n, int d, double **rect);
void errorBadRect(void);
void printRect(FILE* outfile, int d, double** rect);
void errorBadRect(void);
int* order(double *s, unsigned int n);
void sortDouble(double *s, unsigned int n);

#endif
