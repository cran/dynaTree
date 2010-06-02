#ifndef __PALL_H__
#define __PALL_H__ 


/* types of models supported */
typedef enum Model {CONSTANT=1001, LINEAR=1002, CLASS=1003} Model;

typedef struct pall {
  double **X;              /* input matrix */
  double *y;               /* outputs or class labels */
  unsigned int n;          /* number of (x,y) pairs */
  unsigned int m;          /* number of predictors; cols of X */
  unsigned int nc;         /* number of class labels */
  double a;                /* tree prior alpha parameter */
  double b;                /* tree prior beta parameter */
  double minp;             /* tree prior min-part parameter */
  Model model;             /* model indicator */
} Pall;

Pall *new_pall(double **X, unsigned int n, unsigned int m, 
	       double *y, double *params, int model);
void add_data(Pall *data, double **X, unsigned int n, double *y);
void delete_pall(Pall *data);

#endif
