extern "C" 
{
#include "matrix.h"
#include "linalg.h"
#include "pall.h"
}

#include <iomanip>
#include "tree.h"
#include <stdlib.h>
#include <assert.h>
#include <Rmath.h>
#include <math.h>

/*
 * Tree:
 * 
 * tree constructor based on external parameters indicating the
 * partition 
 */

Tree::Tree(Pall *pall_in, int *p, unsigned int n,
	   double **rect, Tree* parent_in)
{
  /* data storage */
  this->rect = rect;
  this->pall = pall_in;
  this->n = n;
  this->p = p;

  /* prior parameters */
  this->as2 = 0;//10;
  this->bs2 = 0;//1;

  /* dummy initial values for const sufficient stats */
  syy = mu = 0.0;

  /* dummy initial values for linear extended stats */
  bb = mm = 0.0;
  XtXi = NULL;
  ldet_XtXi = 0.0;
  bmu = NULL;
  xmean = NULL;

  /* dummy initial values for classifications stats */
  counts = NULL;

  /* tree pointers */
  leftChild = NULL;
  rightChild = NULL;
  if(parent_in != NULL) depth = parent_in->depth+1;
  else depth = 0;
  this->parent = parent_in;

  /* changepoint (split) variables */
  var = 0; val = 0;

  /* update the sufficient parameters */
  Update();
}


/*
 * Tree:
 * 
 * pointer-copy constructor
 */
Tree::Tree(const Tree *told, Tree *parentold)
{
  unsigned int m = told->pall->m;

  /* tree parameters */
  var = told->var; 	
  val = told->val;
  depth = told->depth; 	
  leftChild = rightChild = next = NULL;
  this->parent = parentold;

  /* data */
  pall = told->pall;
  n = told->n;
  if(told->p) p = new_dup_ivector(told->p, n); 
  else p = NULL;
  if(told->rect) rect = new_dup_matrix(told->rect, 2, m);
  else rect = NULL;

  /* priors parameters */
  as2 = told->as2;
  bs2 = told->bs2;
  
  /* sufficient statistics for constant model */
  mu = told->mu;
  syy = told->syy;

  /* extended sufficient stats for linear model */
  bb = told->bb;
  mm = told->mm;
  if(told->XtXi) XtXi = new_dup_matrix(told->XtXi, m, m);
  else XtXi = NULL;
  ldet_XtXi = told->ldet_XtXi;
  if(told->bmu) bmu = new_dup_vector(told->bmu, m);
  else bmu = NULL;
  if(told->xmean) xmean = new_dup_vector(told->xmean, m);
  else xmean = NULL;

  /* sufficient stats for classification */
  if(told->counts) counts = new_dup_uivector(told->counts, pall->nc);
  else counts = NULL;

  /* recurse down the leaves */
  if(! told->isLeaf()) {
    leftChild = new Tree(told->leftChild, this);
    rightChild = new Tree(told->rightChild, this);
  } 
}


/* 
 * ~Tree:
 * 
 * destructor
 */

Tree::~Tree(void)
{
  IEconomy();
  if(leftChild) delete leftChild;
  if(rightChild) delete rightChild;
}


/*
 * IEconomy:
 *
 * free quantities that would not be used by
 * internal nodes 
 */

void Tree::IEconomy(void)
{
  if(p) { free(p); p = NULL; }
  if(rect) { delete_matrix(rect); rect = NULL; }
  if(XtXi) { delete_matrix(XtXi); XtXi = NULL; }
  if(bmu) { free(bmu); bmu = NULL; }
  if(xmean) { free(xmean); xmean = NULL; }
  if(counts) { free(counts); counts = NULL; }
}



/* 
 * getDepth:
 * 
 * return the node's depth in the tree
 */

int Tree::getDepth(void) const
{
  return depth;
}


/*
 * isLeaf:
 * 
 * TRUE if the node is a leaf,
 * FALSE otherwise
 */

bool Tree::isLeaf(void) const
{
  /* assert(!(leftChild != NULL && rightChild == NULL));
     assert(!(leftChild == NULL && rightChild != NULL)); */
  if(leftChild == NULL && rightChild == NULL) return true;
  else return false;
}


/*
 * isRoot:
 * 
 * TRUE if the node is the root (parent == NULL),
 * FALSE otherwise
 */

bool Tree::isRoot(void) const
{
  if(parent == NULL) return true;
  else return false;
}


/*
 * internals:
 * 
 * get a list of internal (non-leaf) nodes, where the first in
 * list is pointed to by the first pointer, and the last by the 
 * last pointer.  The length of the list is returned.
 */

int Tree::internals(Tree **first, Tree **last)
{
  if(isLeaf()) {
    *first = *last = NULL;
    return 0;
  }

  Tree *leftFirst, *leftLast, *rightFirst, *rightLast;
  leftFirst = leftLast = rightFirst = rightLast = NULL;
  
  int left_len = leftChild->internals(&leftFirst, &leftLast);
  int right_len = rightChild->internals(&rightFirst, &rightLast);
  
  if(left_len == 0) {
    this->next = rightFirst;
    *first = this;
    if(right_len > 0) {
      *last = rightLast;
      (*last)->next = NULL;
    } else {
      *last = this;
      (*last)->next = NULL;
    }
    return right_len + 1;
  } else {
    leftLast->next = rightFirst;
    this->next = leftFirst;
    *first = this;
    if(right_len == 0) *last = leftLast;
    else *last = rightLast;
    (*last)->next = NULL;
    return left_len + right_len + 1;
  }
}



/*
 * leaves:
 * 
 * get a list of leaf nodes, where the first in list is 
 * pointed to by the first pointer, and the last by the 
 * last pointer.  The length of the list is returned.
 */

int Tree::leaves(Tree **first, Tree **last)
{
  if(isLeaf()) {
    *first = this;
    *last = this;
    (*last)->next = NULL;
    return 1;
  }
  
  Tree *leftFirst, *leftLast, *rightFirst, *rightLast;
  leftFirst = leftLast = rightFirst = rightLast = NULL;
  
  int left_len = leftChild->leaves(&leftFirst, &leftLast);
  int right_len = rightChild->leaves(&rightFirst, &rightLast);
  
  leftLast->next = rightFirst;
  *first = leftFirst;
  *last = rightLast;
  return left_len + right_len;
}


/*
 * leavesN:
 * 
 * get the partition sizes (n) at all
 * leaf children of this node
 */

int Tree::leavesN(void)
{
  Tree *first, *last;
  int numLeaves = leaves(&first, &last);
  assert(numLeaves > 0);
  int N = 0;
  while(first) {
    N += first->n;
    first = first->next;
  }
  return N;
}


/* 
 * grow:
 * 
 * attempt to add two children to this LEAF 
 * by splitting at a particular var and val
 */

void Tree::grow(int var, double val)
{
  assert(isLeaf());
    
  /* assign the split */
  this->var = var;
  this->val = val;

  /* grow the children; stop if partition too small */
  assert(grow_children());
  assert(leftChild->n + rightChild->n == n);

  /* clear p and rect */
  IEconomy();
}


/*
 * growProb:
 *
 * randomly try growing, and return the probability,
 * but clean up afterwards
 */

double Tree::growProb(int *gvar, double *gval)
{
  assert(isLeaf());

  /* check if we're allowing grows */
  if(pall->a <= 0 || pall->b <= 0) 
    return 0.0;

  /* pick a var and val to grow at */
  this->var = (int) floor(((double)(pall->m)) * runif(0,1));
  this->val = runif(rect[0][var], rect[1][var]);

  /* try growing */
  bool success = grow_children();
  if(!success) return 0.0;

  /* already done in grow_child in new Tree(...) call */

  /* calculate full posterior */
  double prob;
  if(parent != NULL) prob = exp(parent->FullPosterior());
  else prob = exp(FullPosterior());

  /* clean up */
  delete leftChild; leftChild = NULL;
  delete rightChild; rightChild = NULL;

  /* pass back */
  *gvar = var; *gval = val;
  return(prob);
}


/* 
 * pruneProb:
 *
 * try pruning, and return the probability, 
 * but clean up afterwards -- called from a child
 */

double Tree::pruneProb(void)
{
  /* cannot prune if no parent, i.e., if root */
  if(parent == NULL) return 0.0;

  /* get the tree process prior parameters */
  assert(parent->p == NULL); 
  parent->p = parent->GetP(&(parent->n));
  parent->Update();

  /* save old L and R children */
  Tree *oldLC = parent->leftChild;
  Tree *oldRC = parent->rightChild;

  /* calculate the probability of the prune */
  parent->leftChild = NULL;
  parent->rightChild = NULL;
  double prob = exp(parent->FullPosterior());

  /* undo */
  parent->leftChild = oldLC;
  parent->rightChild = oldRC;
  // free(parent->p); parent->p = NULL;
  parent->IEconomy();

  return prob;
}


/*
 * prune:
 *
 * actually perform a prune for keeps
 */

void Tree::prune(void)
{
  assert(!isLeaf());

  /* collect the partition descriptors from 
     the children nodes */
  assert(p == NULL);
  p = GetP(&n);
  rect = GetRect();

  /* update the sufficient information */
  Update();

  /* perform the prune by cutting off the children */
  delete leftChild; leftChild = NULL;
  delete rightChild; rightChild = NULL;
}


/*
 * stayProb:
 *
 * calculate the probability of the current tree, as is,
 * by calculating the posteriors at the leaves of the parent
 */

double Tree::stayProb(void)
{
  assert(isLeaf());
  double pstay;
  if(parent != NULL) pstay = exp(parent->FullPosterior());
  else pstay = exp(FullPosterior());
  assert(!isinf(pstay));
  return pstay;
 }


/*
 * grow_children:
 * 
 * grow both left and right children based on splitpoint --
 * essentially copied from tgp
 */

bool Tree::grow_children(void)
{
  /* can't grow in this case */
  if(n < 2*(pall->minp)) return false;

  unsigned int suc1 = grow_child(&leftChild, LEQ);
  if(!suc1 || !(leftChild->wellSized())) {
    if(leftChild) delete leftChild;
    leftChild = NULL; 
    assert(rightChild == NULL);
    return false;
  }
  unsigned int suc2 = grow_child(&rightChild, GT);
  if(!suc2 || !(rightChild->wellSized())) {
    delete leftChild;
    if(rightChild) delete rightChild;
    leftChild = rightChild = NULL;
    return false;
  }
  assert(suc1 + suc2 == n);
  return true;
}


/*
 * part_child:
 * 
 * creates the data according to the current partition
 * the current var and val parameters, and the operation "op"
 */

int Tree::part_child(FIND_OP op, int **pnew, unsigned int *plen,  
		     double ***newRect)
{
  int *pchild = find_col(pall->X, p, n, var, op, val, plen);
  if(*plen == 0) return 0;

  /* check for big enough partition */
  if(*plen < pall->minp) { free(pchild); return 0; }
  
  /* partition the data and predictive locations */
  *pnew = new_ivector(*plen);
  for(unsigned int j=0; j<*plen; j++) (*pnew)[j] = p[pchild[j]];
  if(pchild) free(pchild); 
  
  /* record the boundary of this partition */
  *newRect = new_matrix(2, pall->m);
  for(unsigned int i=0; i<pall->m; i++) {
    (*newRect)[0][i] = rect[0][i];
    (*newRect)[1][i] = rect[1][i];
  }
  if(op == LEQ) (*newRect)[1][var] = val; 
  else { assert(op == GT); (*newRect)[0][var] = val; }
  
  return (*plen);
}

/*
 * grow_child:
 * 
 * based on current val and var variables, create the corresponding 
 * leftChild partition returns the number of points in the grown region
 */

unsigned int Tree::grow_child(Tree** child, FIND_OP op)
{
  assert(!(*child));
	
  /* find partition indices */
  unsigned int plen; 
  double **newRect = NULL;
  int *pnew = NULL; 
  
  /* construct the partition, assuming it can be made */
  int success = part_child(op, &pnew, &plen, &newRect);
  if(success == 0) return success; /* bad partition */
  
  /* grow the Child */
  (*child) = new Tree(pall, pnew, plen, newRect, this);
  return plen;
}



/* 
 * getN:
 * 
 * return the number of input locations, N
 */

int Tree::getN(void) const
{
  return n;
}


/* 
 * internalsList:
 * 
 * get an array containing the internal nodes of the tree --
 * essentially copied from tgp
 */

Tree** Tree::internalsList(int* len)
{
  Tree *first, *last;
  first = last = NULL;
  *len = internals(&first, &last);
  if(*len == 0) return NULL;
  return first->buildTreeList(*len);
}


/* 
 * leavesList:
 * 
 * get an array containing the leaves of the tree --
 * essentially copied from tgp
 */

Tree** Tree::leavesList(int* len)
{
  Tree *first, *last;
  first = last = NULL;
  *len = leaves(&first, &last);
  if(*len == 0) return NULL;
  return first->buildTreeList(*len);
}


/* 
 * numLeaves:
 * 
 * get a count of the number of leaves in the tree --
 * essentially copied from tgp
 */

int Tree::numLeaves(void)
{
  Tree *first, *last;
  first = last = NULL;
  int len = leaves(&first, &last);
  return len;
}


/*
 * buildTreeList:
 * 
 * takes a pointer to the first element of a Tree list and a 
 * length parameter and builds an array style list --
 * essentially copied from tgp
 */

Tree** Tree::buildTreeList(int len)
{
  int i;
  Tree* first = this;
  Tree** list = (Tree**) malloc(sizeof(Tree*) * (len));
  for(i=0; i<len; i++) {
    assert(first);
    list[i] = first;
    first = first->next;
  }
  return list;
}



/*
 * wellSized:
 * 
 * return true if this node (leaf) is well sized (nonzero 
 * area and > t_minp points in the partition)
 */

bool Tree::wellSized(void) const
{
  if(n <= pall->minp) return false;
  else return true;
}


/* 
 * Height:
 *
 * compute the height of the the tree -- essentially 
 * copied from tgp
 */

int Tree::Height(void) const
{
  if(isLeaf()) return 1;
  
  int lh = leftChild->Height();
  int rh = rightChild->Height();
  if(lh > rh) return 1 + lh;
  else return 1 + rh;
}


/*
 * Prior:
 *
 * calculate the prior probabilities from this node
 * downwards in the tree
 */

double Tree::Prior(void)
{
  double prior;

  /* get the tree process prior parameters */
  double a = pall->a;
  double b = pall->b;

  if(isLeaf()) {

    /* probability of not growing this branch */
    prior = log(1.0 - a*pow(1.0+depth,0.0-b));

  } else {
    
    /* probability of growing here */
    prior = log(a) - b*log(1.0 + depth);

    /* probability of the children */
    prior += leftChild->Prior();
    prior += rightChild->Prior();
  }

  return prior;
}



/*
 * FullPosterior:
 *
 * Calculate the full posterior of (the leaves of) 
 * the tree using the base particles and the probability
 * of growing (or not) at internal (leaf) nodes with
 * process prior determined by a and b
 *
 * returns a log posterior probability
 */

double Tree::FullPosterior( void )
{
  double post;

  /* get the tree process prior parameters */
  double a = pall->a;
  double b = pall->b;

  if(isLeaf()) {

    /* probability of not growing this branch */
    post = log(1.0 - a*pow(1.0+depth,0.0-b));

    post += Posterior();

  } else {
    
    /* probability of growing here */
    post = log(a) - b*log(1.0 + depth);

    /* probability of the children */
    post += leftChild->FullPosterior();
    post += rightChild->FullPosterior();
  }

  return post;
}

/* 
 * Posterior:
 *
 * compute the LOG posterior of the leaf node under
 * whatever model is being used */

double Tree::Posterior(void) /* log post! */
{
  /* sanity checks */
  assert(isLeaf());
  assert(n >= pall->minp);

  /* extract double versions of m and n */
  double dm = (double) 1.0; 
  double dn = (double) n;

  /* special classification treatment */
  if(pall->model == CLASS) {
    double post = 0.0 - lgamma(dn + 1.0);
    dm = pall->nc;
    for(unsigned int i=0; i<pall->nc; i++) {
      post += lgamma(counts[i] + 1.0/dm) - lgamma(1.0/dm);
    }
    return post;
  }
  /* otherwise continue with constan and linear treatment */

  /* linear model adjustment of dm */
  if(pall->model == LINEAR && bb > 0) dm += (double) pall->m;

  /* invalid leaf node */
  if(syy-bb <= 0.0) return -1e300*1e300;

  /* integrated likelihood for this branch */
  double post = 0.0 - 0.5*((as2 + dn-dm) * log(M_2_PI) + log(dn));
  post -= 0.5*(as2 + dn-dm) * log(0.5*(syy-bb + bs2));
  post += lgamma(0.5*(as2 + dn-dm));

  /* include prior if proper */
  if(bs2 > 0 && as2 > 0)
    post += 0.5*as2*log(0.5*bs2) - lgamma(0.5*as2);

  /* further linear augmentations */
  if(pall->model == LINEAR) post += 0.5*ldet_XtXi;

  /* sanity check */
  assert(!isnan(post));

  return post;
}


/*
 * AddDatum:
 *
 * add a new x-y pair to the tree and update the sufficient
 * information using as few operations as possible
 */

Tree* Tree::AddDatum(unsigned int index)
{
  if(isLeaf()) {
    
    /* augment the p-vector */
    p = (int*) realloc(p, sizeof(int)*(n+1));
    p[n] = index;
    n++;

    /* check to see if the rectangle needs expanding */
    ExpandRect(pall->X[index]);

    /* update the classification model */
    if(pall->model == CLASS) counts[(int) pall->y[index]]++;
    else { /* constant and linear models */

      /* update the constant model sufficient statistics */
      double y = pall->y[index];
      double nold = (double) (n-1);
      syy += sq(y) + nold*sq(mu);
      mu = (nold*mu + y)/(nold + 1.0);
      syy -= (nold + 1.0)*sq(mu);
      
      /* update linear part of model */
      if(pall->model == LINEAR) UpdateLinear();
    }

    /* done */
    return this;

  } else {
    assert(p == NULL);
    if(pall->X[index][var] <= val) return leftChild->AddDatum(index);
    else return rightChild->AddDatum(index);
  }
}


/*
 * ExpandRect:
 *
 * possibly expand the rectangle to contain this new point
 */

void Tree::ExpandRect(double *x)
{
  for(unsigned int i=0; i<pall->m; i++) {
    if(x[i] > rect[1][i]) rect[1][i] = x[i];
    else if(x[i] < rect[0][i]) rect[0][i] = x[i];
  }
}


/*
 * GetP:
 *
 * combine the p-vectors from the left and right children --
 * mainly used in the ::prune operation
 */

int* Tree::GetP(unsigned int *n)
{
  if(isLeaf()) {
    *n = this->n;
    assert(p);
    return new_dup_ivector(p, this->n); 
  } else {
    unsigned int nl, nr;
    int *pl = leftChild->GetP(&nl);
    int *pr = rightChild->GetP(&nr);
    *n = nl + nr;
    assert(*n > 0);
    pl = (int*) realloc(pl, sizeof(int) * (*n));
    dupiv(pl+nl, pr, nr); free(pr);
    return pl;
  }
}


/*
 * GetP:
 *
 * combine the rectangles from the left and right children --
 * mainly used in the ::prune operation
 */

double** Tree::GetRect(void)
{
  if(isLeaf()) {
    return new_dup_matrix(this->rect, 2, pall->m);
  } else {
    double **rl = leftChild->GetRect();
    double **rr = rightChild->GetRect();
    for(unsigned int i=0; i<pall->m; i++) {
      rl[0][i] = myfmin(rl[0][i], rr[0][i]);
      rl[1][i] = myfmax(rl[1][i], rr[1][i]);
    }
    delete_matrix(rr);
    return rl;
  }
}


/*
 * GetLeaf:
 *
 * return the leaf node that contains x
 */

Tree* Tree::GetLeaf(double *x)
{
  if(isLeaf()) return this;
  else {
    if(x[var] <= val) return leftChild->GetLeaf(x);
    else return rightChild->GetLeaf(x);
  }
}


/*
 * PostPred:
 *
 * return the posterior predictive probability of the
 * (x,y) pair under whatever model is being used at this
 * leaf 
 */

double Tree::PostPred(double *x, double y)
{
  /* sanity checks */
  assert(isLeaf());
  assert(n >= pall->minp);

  /* double versions of data sizes */
  double dn = (double) n;
  double dm = (double) 1.0;

  /* special handling of classification */
  if(pall->model == CLASS) {
    dm = (double) pall->nc;
    return (((double) counts[(int) y]) + 1.0/dm)/(1.0 + dn);
  } /* otherwise constant and linear model */

  /* adjustment for linear model */
  if(pall->model == LINEAR && bb > 0) dm += (double) pall->m;

  /* check for proper within-partition model */
  if(syy-bb <= 0.0) return 0.0;

  /* dummy variables that will get more use in linear model, later */
  double xtXtXix = 1.0/dn;
  double bmean = 0.0;

  /* adjustments under the linear model */
  if(pall->model == LINEAR) LinearAdjust(x, &bmean, &xtXtXix, NULL);

  /* student t moments */
  double nt = (as2 + dn-dm)/(1.0 + xtXtXix);
  double prec = sqrt(nt/(syy-bb + bs2));
  double yt = (y-mu-bmean)*prec;

  /* calculate the predictve probability */
  double prob = dt(yt, dn-dm, 0)*prec;

  return prob;
}


/*
 * LinearAdjust:
 *
 * adjustments to the latter two arguments, used for prediction,
 * that must be made under the linear model
 */

double Tree::LinearAdjust(double *x, double *bmean, double *xtXtXix, 
			  double *y)
{
  assert(pall->model == LINEAR);

  unsigned int m = pall->m;
  
  /* adjustments under the linear model */
  if(bmean) {
    *bmean -= mm;
    for(unsigned int i=0; i<m; i++) *bmean += x[i]*bmu[i];
  }

  /* subtract off mean of X */
  double *xdup = new_dup_vector(x, m);
  //add_vector(1.0, xdup, 0.0-1.0, xmean, m);
  linalg_daxpy(m,0.0-1.0,xmean,1,xdup,1);
  
  /* t(x) %*% XtXi %*% x */
  double *XtXix = new_zero_vector(m);
  linalg_dsymv(m, 1.0, XtXi, m, xdup, 1, 0.0, XtXix, 1);
  *xtXtXix += linalg_ddot(m, xdup, 1, XtXix, 1);
  free(xdup);

  /* deal with y (for ALC) */
  double yadj;
  if(y) {
    double *ydup = new_dup_vector(y, m);
    // add_vector(1.0, ydup, 0.0-1.0, xmean, m);
    linalg_daxpy(m,0.0-1.0,xmean,1,ydup,1);
    yadj = linalg_ddot(m, ydup, 1, XtXix, 1);
    free(ydup);
  } else yadj = 0.0;

  /* final cleanup */
  free(XtXix);

  return yadj;
 }


/*
 * Parent:
 *
 * return the parent of this node
 */

Tree* Tree::Parent(void) const
{
  return parent;
}


/*
 * Update:
 *
 * update the sufficient paramters for this tree node
 * based on the data under management in the corresponding
 * partition
 */

void Tree::Update(void)
{
  /* sanity check */
  assert(p);

  /* special handling of classification */
  if(pall->model == CLASS) {
    if(counts == NULL) counts = new_zero_uivector(pall->nc);
    else zerouiv(counts, pall->nc);
    for(unsigned int i=0; i<n; i++) 
      counts[(int) pall->y[p[i]]]++;
    return;
  }

  /* calculate intercept */
  mu = 0.0;
  for(unsigned int i=0; i<n; i++) mu += pall->y[p[i]];
  mu /= n;

  /* calculate un-scaled standard error */
  syy = 0.0;
  for(unsigned int i=0; i<n; i++) syy += sq(pall->y[p[i]] - mu);

  /* possibly update for the linear model too */
  if(pall->model == LINEAR) UpdateLinear();
}


/* 
 * UpdateLinear:
 *
 * finish off mu (b) calculation -- CAN MAKE THIS MORE
 * EFFICIENT BY RE-USING VECTORS THAT HAVE BEEN FREED TO
 * CUT DOWN ON MALLOCS
 */

void Tree::UpdateLinear(void)
{
  /* sanity checks */
  assert(pall->model == LINEAR);
  
  /* build y and subtract the (constant) mean */
  double *y = new_sub_vector(p, pall->y, n);
  for(unsigned int i=0; i<n; i++) y[i] -= mu;

  /* build X and subtract off mean */
  unsigned int m = pall->m;
  double **Xorig = new_matrix(n, m);
  for(unsigned int i=0; i<n; i++) dupv(Xorig[i], pall->X[p[i]], m);
  double **X = new_dup_matrix(Xorig, n, m);
  if(!xmean) xmean = new_vector(m);
  wmean_of_columns(xmean, X, n, m, NULL);
  for(unsigned int i=0; i<n; i++) {
    // add_vector(1.0, X[i], 0.0-1.0, xmean, m);
    linalg_daxpy(m,0.0-1.0,xmean,1,X[i],1);
  }
  
  /* calculate XtX */
  double **XtX = new_zero_matrix(m, m);
  linalg_dgemm(CblasNoTrans,CblasTrans,m,m,n,1.0,
		 X,m,X,m,0.0,XtX,m);

  /* calculate XtXi = inv(XtX) and determinant */
  if(!isZero(XtX, m, 1)) {
    double ** XtX_chol = new_dup_matrix(XtX, m, m);
    if(!XtXi) XtXi = new_id_matrix(m);
    else id(XtXi, m);
    int info = linalg_dposv(m, XtX_chol, XtXi);
    if(info != 0) { /* zero values for bad Gram matrix */
      zero(XtXi, m, m);
      zero(XtX, m, m);
      ldet_XtXi = 0.0;
    } else ldet_XtXi = 0.0-log_determinant_chol(XtX_chol, m);
    delete_matrix(XtX_chol);
  } else {
    if(!XtXi) XtXi = new_zero_matrix(m, m);
    else zero(XtXi, m, m);
    ldet_XtXi = 0.0; 
  }

  /* calculate Xty = t(X) %*% y */
  double *Xty = new_zero_vector(m);
  linalg_dgemv(CblasNoTrans,m,n,1.0,X,m,y,1,0.0,Xty,1);
  delete_matrix(X);
  /* free(y); use later for Xbeta */

  /* now calculate beta-hat (bmu) = inv(t(X) %*% X) %*% t(X) %*% y */
  if(!bmu) bmu = new_zero_vector(m);
  else zerov(bmu, m);
  linalg_dsymv(m, 1.0, XtXi, m, Xty, 1, 0.0, bmu, 1);
  /* free(Xty); use later for XtXbmu */

  /* compute: bb = t(bmu) %*% (XtX) %*% bmu */
  double *XtXbmu = Xty; zerov(XtXbmu, m); //new_zero_vector(m); 
  linalg_dsymv(m, 1.0, XtX, m, bmu, 1, 0.0, XtXbmu, 1);
  bb = linalg_ddot(m, bmu, 1, XtXbmu, 1);
  free(XtXbmu);
  delete_matrix(XtX);

  /* augment mu */
  double *Xbeta = y; zerov(Xbeta, n); //new_zero_vector(n);
  linalg_dgemv(CblasTrans,m,n,1.0,Xorig,m,bmu,1,0.0,Xbeta,1);
  delete_matrix(Xorig);
  mm = meanv(Xbeta, n);
  free(Xbeta);
}


/*
 * Predict:
 *
 * extract the predictive distribution at x in terms of its
 * moments and, perhaps, quantiles
 */

void Tree::Predict(double *x, double *mean, double *sd, double *df, 
		   double *var, double *q1, double *q2)
{
   if(isLeaf()) {

     /* this is the wrong predict function for classification */
     assert(pall->model != CLASS);

     double dm = (double) 1.0;
     double dn = (double) n;
     if(pall->model == LINEAR && bb > 0) dm += (double) pall->m;

     double xtXtXix = 1.0/dn;
     double bmean = 0.0;

     /* adjustments under the linear model */
     if(pall->model == LINEAR) LinearAdjust(x, &bmean, &xtXtXix, NULL);

     /* student t parameters */
     double nt = (as2 + dn - dm)/(1.0 + xtXtXix);
     *sd = sqrt((syy-bb + bs2)/nt);
     *df = dn - dm;
     
     /* mean and variance of the t */
     *var = (*df)*(*sd)/(*df - 2.0);
     *mean = mu + bmean;

     /* quantiles of the t */
     if(q1) *q1 = qt(0.05, *df, 1, 0)*(*sd) + *mean;
     if(q2) *q2 = qt(0.95, *df, 1, 0)*(*sd) + *mean;

   } else {
     if(x[this->var] <= val) 
       return leftChild->Predict(x, mean, sd, df, var, q1, q2);
     else 
       return rightChild->Predict(x, mean, sd, df, var, q1, q2);
   }
}


/*
 * Predict:
 *
 * prediction for classification
 */

void Tree::Predict(double *x, double *pred)
{
  if(isLeaf()) {
    assert(pall->nc > 0);

    /* double versions of integers */
    double dm = (double) pall->nc;
    double dn = (double) n;

    /* probability for each class */
    for(unsigned int i=0; i<pall->nc; i++)
      pred[i] = (((double) counts[i]) + 1.0/dm)/(1.0 + dn);

    return;
  } else {
     if(x[this->var] <= val) return leftChild->Predict(x, pred);
     else return rightChild->Predict(x, pred);
   }
}


/* 
 * ALC:
 *
 * calculate the Active Learning Cohn statistic for 
 * sequential design
 */

double Tree::ALC(double *x, double *y)
{
   if(isLeaf()) {

     /* ALC not supported for classification */
     assert(pall->model != CLASS);

     /* get double versions of m and n */
     double dn = (double) n;
     double dm = (double) 1.0;
     if(pall->model == LINEAR && bb > 0) dm += (double) pall->m;

     /* important prediction quanties */
     double xtXtXix = 1.0/n;
     double ytXtXix = 1.0/n;

     /* adjustments to above under the linear model */
     if(pall->model == LINEAR) 
       ytXtXix += LinearAdjust(x, NULL, &xtXtXix, y);

     /* wrapping up */
     assert(dn-dm-2.0 > 0);
     double s2 = (syy-bb + bs2)/(dn-dm+as2);
     return ((dn-dm)*s2/(dn-dm-2.0)) * sq(ytXtXix)/(1.0 + xtXtXix);

   } else { /* recurse into leaves */
     if(x[var] <= val && y[var] <= val) return leftChild->ALC(x, y);
     else if(x[var] > val && y[var] > val) return rightChild->ALC(x, y);
     else return 0.0; /* x and y aren't in same region */
   }
}


/*
 * Print:
 *
 * printing debugging information about the leaf nodes --
 * currently out of use
 */

void Tree::Print(void)
{
  if(isLeaf()) {
    assert(0);
  } else {
    leftChild->Print();
    rightChild->Print();
  }
}


/*
 * leavesAvgSize:
 *
 * collecting the average size of leaf nodes 
 */

double Tree::leavesAvgSize(void)
{
  Tree *first, *last;
  int numLeaves = leaves(&first, &last);
  assert(numLeaves > 0);
  double size = 0.0;
  while(first) {
    size += first->n;
    first = first->next;
  }
  return size/((double) numLeaves);
}


/*
 * log_determinant_chol:
 *
 * returns the log determinant of the n x n
 * choleski decomposition of a matrix M
 */

double log_determinant_chol(double **M, const unsigned int n)
{
  double log_det;
  unsigned int i;

  /* det = prod(diag(R)) .^ 2 */
  log_det = 0;
  for(i=0; i<n; i++) log_det += log(M[i][i]);
  log_det = 2*log_det;

  return log_det;
}
