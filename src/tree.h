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


#ifndef __TREE_H__
#define __TREE_H__ 

extern "C" {
#include "pall.h"
#include "matrix.h"
}
#include <cstdio>

class Particle;

class Tree
{
 private: 
  
  /* pointer to particle this tree lives in */
  Particle *particle;

  /* data */
  unsigned int n;               /* size of data in current parition */
  int *p;			/* n-vector of indices into original data */
  double *al;                   /* n-vector of stored active learning stats */
  double ng;                    /* (drifting) number of retired observations */

  /* sufficient stats for the classification model */
  unsigned int *counts;         /* counts in each class */

  /* retireed stats for the classification model */
  double *gcounts;              /* (drifting) retired counts in each class */
  /* sufficient stats for constant model */
  double sy;                    /* data sum; n*ybar */
  double syy;                   /* sum of y-squared  */

  /* retired stats for constant model */
  double syg;                   /* retireed ng*ybar */
  double syyg;                  /* retireed sum of y-squared */

  /* extra sufficient stats for linear model */
  double **XtX;                 /* t(X) %*% X, NULL when icept = FALSE */
  double *Xty;                  /* t(X) %*% y, NULL when icept = FALSE */
  double **XtXi;                /* inv(XtX) */
  double ldet_XtXi;             /* log determiniant of XtXi */
  double *bmu;                  /* inv(t(X) %*% X) %*% t(X) %*% y */
  double bb;                    /* t(bmu) %*% XtX %*% bmu */
  double *xmean;                /* mean of the cols of X */

  /* extra retired stats for linear model; must have icept = FALSE */
  double **XtXg;                 /* retireed t(X) %*% X */
  double *Xtyg;                  /* retireed t(X) %*% y */
  
  /* splits */
  int var;	                /* split point dimension */
  double val;		        /* split point value */
  
  /* structural */
  Tree* parent;		        /* parent partition */
  Tree* leftChild;	        /* partition LEQ (<=) split point */
  Tree* rightChild;	        /* partition GT (>) split point */
  Tree* next;		        /* used for making lists of tree nodes */
  int depth;	                /* depth of partition in tree */
   
  /* functions */
  
  /* create lists of tree nodes, 
   * and traverse them from first to next ... to last */
  int leaves(Tree** first, Tree** last);
  int internals(Tree **first, Tree **last);
  
  /* creating new leaves, and removing them */
  unsigned int grow_child(Tree** child, FIND_OP op);
  int part_child(FIND_OP op, int **pnew, unsigned int *plen);
  bool grow_children(bool missrand);
  void Missing(void);
  bool Missing(unsigned int index, unsigned int var);

  /* calculating split locations */
  bool ChooseVarVal(void);
  unsigned int GetXcol(unsigned int var, double *x);

 public:

  /* constructor, destructor and misc partition initialization */
  Tree(Particle *particle, int *p, unsigned int n, Tree* parent_in);
  Tree(const Tree *oldt, Particle *particle, Tree *parent);
  ~Tree(void);
  void IEconomy(void);
  
  /* propose tree operations */
  void grow(int var, double val);
  double growProb(int *var, double *val);
  double stayProb(void);
  double pruneProb(void);
  void prune(void);

  /* access functions:
   * return current values of the parameters */
  int getDepth(void) const;
  int getN(void) const;
  Tree* GetLeaf(double *x, int *xna);
  Tree* GetLeaf(unsigned int index);
  double Min(unsigned int var);
  double Max(unsigned int var);

  /* access function: info about nodes */
  bool goLeft(unsigned int index, bool always);
  bool goLeft(double *x, int *xna);
  double LeftBal(void);
  bool isLeaf(void) const;
  bool isRoot(void) const;
  char* State(int which);
  void Clear(void);

  /* create an array of typed tree nodes,
   * passing back the length of the array */
  Tree** leavesList(int* len);
  Tree** internalsList(int* len);
  Tree** buildTreeList(int len);
  int numLeaves(void);
  Tree* Parent(void) const; 
  Tree* Sibling(void) const; 

  /* size checks */
  bool wellSized(void) const;
  int Height(void) const;
   
  /* computing the full posterior or likelihood of the tree */
  double Prior(void);
  double FullPosterior(void);
  double Posterior(void);
  double Regression(double *mean, double *s2numer, double *df, double *s2p_out);
  double PostPred(double *x, double y);

  /* calculating sufficient statistics for leaf nodes */
  void Calc(void);
  void CalcClass(void);
  void CalcConst(void);
  void CalcLinear(void);
  void ReCalcLinear(void);

  /* adjustments for linear prediction */
  double LinearAdjust(double *x, double *bmean, double *xtXtXix, 
		      double *XtXix, double *y);

  /* adding and retiring (moving to prior) data to a leaf node */
  Tree* AddDatum(unsigned int index);
  Tree* RetireDatum(unsigned int index, double lambda);
  void DecrementP(unsigned int oldi, unsigned int newi);
  void Collapse(void);

  /* for rejuvination */
  void ReorderP(int *o);

  /* collecting info from leaves to construct parents */
  int* GetP(unsigned int *n);
  void AccumClass(unsigned int *counts, double *gcounts);
  void AccumNg(double *ng);
  void AccumConst(double *sy, double *syy, double *ng, double *syg, 
		  double *syyg);
  void AccumLinear(double **XtX, double *Xty, double **XtXg, double *Xtyg);
  void AccumCalc(void);

  /* summary statistics */
  int leavesN(void);
  double leavesAvgSize(void);
  double leavesAvgRetired(void);
  void Print(void);
  unsigned int GetVar(void);
  void SameLeaf(double **X, int *p, unsigned int n, int *count);
  
  /* prior capping/forgetting factor for online learning */
  void CapRetired(void);

  /* prediction and adaptive sampling */
  void Predict(double *xx, double *mean, double *sd, double *df);
  double ALC(double *x, double *y);
  double ALC(double *x, double **rect, int *cat, bool approx);
  void ALC(double **rect, int *cat, bool approx, double *alc_out);
  double ECI(double *x, double *y, double ymean, double ysd, 
	     double fmin, double ei);
  double AvgVar(double **rect, int *cat, bool approx);
  double AvgEntropy(double **rect, int *cat, bool approx);
  double Relevance(double **rect, int *cat, bool approx, double *delta);

  /* prediction for classification */
  void Predict(double *x, double *pred);
  void Coef(double *x, double *beta);
  void Predict(double *pred);
  double Predict(double *x, unsigned int cls);
  void Entropy(double *entropy_out);
};

/* calculating subroutines */
void Collapse(Tree *collapse);
double calculate_linear(unsigned int m, double **XtX, double*Xty, 
			double **XtXi, double *ldet_XtXi, double *bmu);
double log_determinant_chol(double **M, const unsigned int n);
double intdot2(unsigned int m, double c, double *y, double *a, double *b,
	       int *cat, double approx);
double intdot(unsigned int m, double c, double **g, double *a, double *b,
	      int *cat, double approx);

#endif
