#ifndef __TREE_H__
#define __TREE_H__ 

extern "C" {
#include <stdio.h>
#include "pall.h"
#include "matrix.h"
}

class Tree
{
 private: 
  
  /* prior parameter settings */
  double as2;                  
  double bs2;

  /* data */
  unsigned int n;               /* size of data in current parition */
  int *p;			/* n, indices into original data */
  double **rect;                /* bounding rectangle of data */

  /* sufficient stats for the classification model */
  unsigned int *counts;         /* counts in each class */

  /* sufficient stats for constant model */
  double mu;                    /* data mean */
  double syy;                   /* sum of squared residuals (from mu)  */

  /* extra sufficient stats for linear model */
  double bb;                    /* t(bmu) %*% XtX bmu */
  double mm;                    /* for mu adjustment */
  double **XtXi;                /* solve(XtX) */
  double ldet_XtXi;             /* log determiniant of XtXi */
  double *bmu;                  /* mean of the beta posterior */
  double *xmean;                /* mean of the cols of X */

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
  int part_child(FIND_OP op, int **pnew, unsigned int *plen, 
		 double ***newRect);
  bool grow_children(void);

  /* compute lost of the posterior
   * (likelihood + plus some prior stuff) 
   * of a particular lef, or all leaves */
    
 public:

  Pall *pall;		        /* holding stuff common to all particles */
  
  /* constructor, destructor and misc partition initialization */
  Tree(Pall *pall, int *p, unsigned int n, double **rect, 
       Tree* parent_in);
  Tree(const Tree *oldt, Tree *parent);
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
  Tree* GetLeaf(double *x);

  /* access function: info about nodes */
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

  /* size checks */
  bool wellSized(void) const;
  int Height(void) const;
   
  /* computing the full posterior or likelihood of the tree */
  double Prior(void);
  double FullPosterior(void);
  double Posterior(void);
  double PostPred(double *x, double y);
  void Update(void);
  void UpdateLinear(void);
  double LinearAdjust(double *x, double *bmean, double *xtXtXix, double *y);

  /* adding data to a leaf node */
  Tree* AddDatum(unsigned int index);
  void ExpandRect(double *x);

  /* collecting info from leaves to construct parents */
  int* GetP(unsigned int *n);
  double** GetRect(void);

  /* summary statistics */
  int leavesN(void);
  double leavesAvgSize(void);
  void Print(void);
  
  /* prediction and adaptive sampling */
  void Predict(double *xx, double *mean, double *sd, double *df, 
	       double *var, double *q1, double *q2);
  double ALC(double *x, double *y);

  /* prediction for classification */
  void Predict(double *x, double *pred);
};

/* calculating determinants */
double log_determinant_chol(double **M, const unsigned int n);

#endif
