#include <stdlib.h>
#include <assert.h>
#include "pall.h"
#include "rhelp.h"
#include "matrix.h"

/*
 * new_pall:
 *
 * create a new pall structure with (x,y) pairs and
 * bounding rectangle; the pointers to X
 */

Pall *new_pall(double **X, unsigned int n, unsigned int m, 
	       double *y, double *params, int model_in)
{
  unsigned int w;
  Pall *pall;
  pall = (Pall*) malloc(sizeof(struct pall));
  pall->X = new_dup_matrix(X, n, m);
  pall->n = n;
  pall->m = m;
  pall->y = new_dup_vector(y, n);
  pall->a = params[0];
  pall->b = params[1];
  pall->minp = (unsigned int) params[2];
  pall->nc = 0;

  /* get the model */
  if(model_in == 1) pall->model = CONSTANT;
  else if(model_in == 2) pall->model = LINEAR;
  else if(model_in == 3) {
    pall->model = CLASS;
    pall->nc = (unsigned int) max(y, n, &w) + 1;
  }
  else {
    warning("model %d not defined, using constant\n", model_in);
    pall->model = CONSTANT;
  }

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
  free(pall->y);
  free(pall);
}


/*
 * add_data:
 *
 * add new rows to X and components to y
 */

void add_data(Pall *pall, double **X, unsigned int n, double *y)
{
  unsigned int i, nnew;

  /* sanity check for classification */
  if(pall->model == CLASS)
    for(unsigned int i=0; i<n; i++)
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

  /*  update n */
  pall->n = nnew;
}
