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
 * License alog1 with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
 *
 ****************************************************************************/


extern "C" 
{
#include "matrix.h"
#include "rhelp.h"
#include "linalg.h"
#include "pall.h"
}

#include <iomanip>
#include "tree.h"
#include "particle.h"
#include <stdlib.h>
#include <assert.h>
#include <Rmath.h>
#include <math.h>

// #define DEBUG

/*
 * Tree:
 * 
 * tree constructor based on external parameters indicating the
 * partition 
 */

Tree::Tree(Particle *particle_in, int *p, unsigned int n, Tree* parent_in)
{
  /* data storage */
  this->particle = particle_in;
  Pall *pall = particle->pall;
  this->n = n;
  this->p = p;
  this->al = NULL;

  /* retired model sufficient stats start NULL/zero */
  this->Xtyg = NULL;
  this->XtXg = NULL;
  this->gcounts = NULL;
  this->syg = this->syyg = 0.0;
  this->ng = 0.0;

  /* retired/prior parameters */
  /* SHOULD MAKE INTO OWN FUNCTION */
  if(parent_in != NULL && parent_in->ng != 0) {

    /* size of number retired in this leaf */
    double mult = ((double) n)/((double) parent_in->n);
    this->ng = mult * parent_in->ng;

    /* different retiring for each model type */
    if(this->ng > 0) {

      /* re-calculated weight of prior from parent */
      mult = ((double) this->ng)/((double) parent_in->ng);

      /* retired classification model sufficient stats */
      if(pall->model == CLASS) {	
	gcounts = new_vector(pall->nc);
	for(unsigned int i=0; i<pall->nc; i++)
	  gcounts[i] = mult * parent_in->gcounts[i];
      } else {
	/* retired constant model sufficient stats */
	this->syg = mult * parent_in->syg;
	this->syyg = mult * parent_in->syyg;
	
	/* and then retired linear model sufficient stats */
	if(pall->model == LINEAR) {
	  assert(pall->icept == FALSE);
	  unsigned int m = pall->bmax;
	  this->XtXg = new_dup_matrix(parent_in->XtXg, m, m);
	  scalev(*XtXg, m*m, mult);
	  this->Xtyg = new_dup_vector(parent_in->Xtyg, m);
	  scalev(Xtyg, m, mult);
	}
      }
    }
  }

  /* cap contribution of retired stats */
  // CapRetired();

  /* dummy initial values for const sufficient stats */
  syy = sy = 0.0;

  /* dummy initial values for linear extended stats */
  bb = ldet_XtXi = 0.0;
  XtX = XtXi = NULL;
  Xty = bmu = xmean = NULL;

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
  Calc();
}


/*
 * Tree:
 * 
 * pointer-copy constructor
 */

Tree::Tree(const Tree *told, Particle *particle, Tree *parentold)
{
  /* tree parameters */
  var = told->var; 	
  val = told->val;
  depth = told->depth; 	
  leftChild = rightChild = next = NULL;
  this->parent = parentold;

  /* data */
  this->particle = particle;
  n = told->n;
  if(told->p) p = new_dup_ivector(told->p, n); 
  else p = NULL;

  /* active learning stats */
  if(told->al) al = new_dup_vector(told->al, n);
  else al = NULL;

  /* sufficient statistics for constant model */
  sy = told->sy;
  syy = told->syy;

  /* retired constant model sufficient stats */
  ng = told->ng;
  syg = told->syg;
  syyg = told->syyg;
  
  /* extended sufficient stats for linear model */
  Pall *pall = particle->pall;
  bb = told->bb;
  if(told->XtXi) XtXi = new_dup_matrix(told->XtXi, pall->bmax, pall->bmax);
  else XtXi = NULL;
  if(told->XtX) XtX = new_dup_matrix(told->XtX, pall->bmax, pall->bmax);
  else XtX = NULL;
  ldet_XtXi = told->ldet_XtXi;
  if(told->Xty) Xty = new_dup_vector(told->Xty, pall->bmax);
  else Xty = NULL;
  if(told->bmu) bmu = new_dup_vector(told->bmu, pall->bmax);
  else bmu = NULL;
  if(told->xmean) xmean = new_dup_vector(told->xmean, pall->bmax);
  else xmean = NULL;

  /* extended retired sufficient stats for linear model */
  if(told->XtXg) XtXg = new_dup_matrix(told->XtXg, pall->bmax, pall->bmax);
  else XtXg = NULL;
  if(told->Xtyg) Xtyg = new_dup_vector(told->Xtyg, pall->bmax);
  else Xtyg = NULL;

  /* sufficient stats for classification */
  if(told->counts) counts = new_dup_uivector(told->counts, pall->nc);
  else counts = NULL;

  /* retired sufficient stats for classification */
  if(told->gcounts) gcounts = new_dup_vector(told->gcounts, pall->nc);
  else gcounts = NULL;

  /* recurse down the leaves */
  if(! told->isLeaf()) {
    leftChild = new Tree(told->leftChild, particle, this);
    rightChild = new Tree(told->rightChild, particle, this);
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
  if(XtXi) { delete_matrix(XtXi); XtXi = NULL; }
  if(XtX) { delete_matrix(XtX); XtX = NULL; }
  if(Xty) { free(Xty); Xty = NULL; }
  if(bmu) { free(bmu); bmu = NULL; }
  if(xmean) { free(xmean); xmean = NULL; }
  if(XtXg) { delete_matrix(XtXg); XtXg = NULL; }
  if(Xtyg) { free(Xtyg); Xtyg = NULL; }
  if(counts) { free(counts); counts = NULL; }
  if(gcounts) { free(gcounts); gcounts = NULL; }
  if(al) { free(al); al = NULL; }
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
#ifdef DEBUG
  assert(isLeaf());
#endif
    
  /* assign the split */
  assert(var >= (int) particle->pall->smin);
  this->var = var;
  this->val = val;

  /* grow the children; stop if partition too small */
  assert(grow_children());
  assert(leftChild->n + rightChild->n == n);

  /* clear p and and other data */
  IEconomy();
}


/*
 * ChooseVarVal:
 *
 * pick a var and val to grow at: either propose
 * only from valid split locations (LUALL) or 
 * blindly from the bounding rectangle of the data 
 * in REJECT case
 *
 */

bool Tree::ChooseVarVal(void)
{
  /* get pall */
  Pall *pall = particle->pall;

  /* do nothing if we can't possibly split */
  if(n < 2*pall->minp) return false;

  /* only use ones that make valid partitions */
  if(pall->rprop == LUALL) {

    /* allocations */
    double **rect = new_matrix(2, pall->m - pall->smin);
    double *x = new_vector(n);
    int *vars = iseq(pall->smin, pall->m - 1);
    
    /* calculate the growable dimensions, simultaneously,
       the edges of the growable rectangle */
    int growable = 0;
    for(unsigned int j=0; j<pall->m - pall->smin; j++) {
      for(unsigned int i=0; i<n; i++) x[i] = pall->X[p[i]][vars[j]];
      rect[0][j] = quick_select(x, n, pall->minp - 1);
      rect[1][j] = quick_select(x, n, n - pall->minp);
      if(rect[0][j] < rect[1][j]) vars[growable++] = vars[j];
    }
    free(x);
    
    /* return failure if there are no growable dims, which could 
       happen even if n >= 2*minp if there are repeated x-values */
    if(growable == 0) { 
      free(vars); delete_matrix(rect); return false; 
    } 
    
    /* choose var and val from growable */
    var = vars[(int) floor(((double)(growable)) * runif(0,1))];
    val = runif(rect[0][var-pall->smin], rect[1][var-pall->smin]);
    
    /* clean up */
    free(vars);
    delete_matrix(rect);
    
  } else {  /* blindly choose from rectangle */
  
    /* choose dimension uniformly */
    var = (int) floor(((double)(pall->m - pall->smin)) * runif(0,1));
    var += pall->smin;
    double mn, mx;
    mn = mx = pall->X[0][var];
    
    if(pall->rprop == LUVAR) { /* use LU on the var */
      double *x = new_vector(n);
      for(unsigned int i=0; i<n; i++) x[i] = pall->X[p[i]][var];
      mn = quick_select(x, n, pall->minp - 1);
      mx = quick_select(x, n, n - pall->minp);
      free(x);
      if(mn >= mx) return false;
    } else { /* find side of the rectangle */
      for(unsigned int i=1; i<n; i++) {
	if(pall->X[i][var] < mn) mn = pall->X[i][var];
	else if(pall->X[i][var] > mx) mx = pall->X[i][var];
      }
    }
    
    /* select val */
    val = runif(mn, mx);
  }

  return true;
}

/*
 * growProb:
 *
 * randomly try growing, and return the probability,
 * but clean up afterwards
 */

double Tree::growProb(int *gvar, double *gval)
{
#ifdef DEBUG
  assert(isLeaf());
#endif

  /* pall for theparticle this tree belongs to */
  Pall *pall = particle->pall;

  /* check if we're allowing grows */
  if(pall->a <= 0 || pall->b <= 0) 
    return 0.0;

  /* choose the split var and val */
  if(!ChooseVarVal()) return -1e300*1e300;

  /* try growing */
  bool success = grow_children();
  /* should have returned above if !success */
  if(pall->rprop == LUALL) assert(success);  
  else if(!success) return -1e300*1e300;

  /* calculate full posterior */
  double prob;
  if(parent != NULL) prob = parent->FullPosterior();
  else prob = FullPosterior();

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
  if(parent == NULL) return -1e300*1e300;

  /* get the node data indices */
  assert(parent->p == NULL); 
  parent->p = parent->GetP(&(parent->n));

  /* check if the prune is allowed */
  if(parent->n < 2*particle->pall->minp) {
    parent->IEconomy();
    return -1e300*1e300;
  }
  
  /* update new parameters */
  parent->AccumCalc();

  /* save old L and R children */
  Tree *oldLC = parent->leftChild;
  Tree *oldRC = parent->rightChild;

  /* calculate the probability of the prune */
  parent->leftChild = NULL;
  parent->rightChild = NULL;
  double prob = parent->FullPosterior();

  /* undo */
  parent->leftChild = oldLC;
  parent->rightChild = oldRC;
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

  /* update the sufficient information */
  AccumCalc();

  /* perform the prune by cutting off the children */
  delete leftChild; leftChild = NULL;
  delete rightChild; rightChild = NULL;
}


/*
 * stayProb:
 *
 * calculate the (log) robability of the current tree, as is,
 * by calculating the posteriors at the leaves of the parent
 */

double Tree::stayProb(void)
{
#ifdef DEBUG
  assert(isLeaf());
#endif
  double pstay;
  if(parent != NULL) pstay = parent->FullPosterior();
  else pstay = FullPosterior();
  assert(isfinite(pstay));
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
  if(n < 2*(particle->pall->minp)) return false;

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

int Tree::part_child(FIND_OP op, int **pnew, unsigned int *plen)  
{
  int *pchild = find_col(particle->pall->X, p, n, var, op, val, plen);
  if(*plen == 0) return 0;

  /* check for big enough partition */
  if(*plen < particle->pall->minp) { free(pchild); return 0; }
  
  /* partition the data and predictive locations */
  *pnew = new_ivector(*plen);
  for(unsigned int j=0; j<*plen; j++) (*pnew)[j] = p[pchild[j]];
  if(pchild) free(pchild); 
  
  return (*plen);
}

/*
 * grow_child:
 * 
 * based on current val and var variables, create the corresponding 
 * leftChild partition returns the number of points in the grown 
 * region
 */

unsigned int Tree::grow_child(Tree** child, FIND_OP op)
{
  assert(!(*child));
	
  /* find partition indices */
  unsigned int plen; 
  int *pnew = NULL; 
  
  /* construct the partition, assuming it can be made */
  int success = part_child(op, &pnew, &plen);
  if(success == 0) return success; /* bad partition */
  
  /* grow the Child */
  (*child) = new Tree(particle, pnew, plen, this);
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
  if(n < particle->pall->minp) return false;
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
  double a = particle->pall->a;
  double b = particle->pall->b;

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
  double a = particle->pall->a;
  double b = particle->pall->b;

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
#ifdef DEBUG
  assert(isLeaf());
#endif
  Pall *pall = particle->pall;
  assert(n + ng >= pall->minp);

  /* special classification treatment */
  if(pall->model == CLASS) {
    double post = 0.0 - lgamma(((double) n) + 1.0);
    double dm = (double) pall->nc;
    if(ng > 0) { /* add in retires */
      assert(gcounts);
      double asum = 0; double zsum = 0;
      for(unsigned int i=0; i<pall->nc; i++) {
	double dpi = gcounts[i] + 1.0/dm; asum += dpi;
	double dci = (double) counts[i]; zsum += dci;
	post += lgamma(dci + dpi) - lgamma(dpi) - lgamma(dci + 1.0);
      }
      post += lgamma(((double) n) + 1.0) + lgamma(asum) - lgamma(asum + zsum);
    } else { /* no retires */
      for(unsigned int i=0; i<pall->nc; i++) {
	double dci = (double) counts[i];
	post += lgamma(dci + 1.0/dm) - lgamma(dci + 1.0);
      }
      post -= dm * lgamma(1.0/dm);
    }
    return post;
  }

  /* initial regression calculations */
  double df, s2numer, s2p;
  double dm = Regression(NULL, &s2numer, &df, &s2p);
  if(s2numer <= 0.0) return -1e300*1e300;

  /* integrated likelihood for this branch */
  double dn = (double) n;
  double post = 0.0 - 0.5*((dn-dm) * log(M_2_PI));
  if(pall->icept) {
    if(n > 0) post -= 0.5*log(dn);
    if(ng > 0) post += 0.5*log(((double)ng));
  }
  post -= 0.5*df*log(0.5*(s2numer));
  post += lgamma(0.5*df);

  /* include sigma2 prior if proper */
  double nup = pall->nu0 + ng;
  if(nup > 0 && s2p > 0)
    post += 0.5*nup*log(0.5*s2p) - lgamma(0.5*nup);

  /* further linear augmentations */
  if(pall->model == LINEAR) post += 0.5*ldet_XtXi;

  /* sanity check */
  assert(!isnan(post));

  return post;
}


/*
 * Regression:
 *
 * initial calculation of regression parameters shared
 * by most leaf-operating (regression, i.e., not classification) 
 * functions in the Tree class; the mean is allowed to be null
 * since it is not always needed
 */

double Tree::Regression(double *mean, double *s2numer, double *df, 
			double *s2p_out)
{
  /* pall for particle this tree belongs to */
  Pall *pall = particle->pall;

  /* get double versions of m and n */
  double dn = (double) n;
  double dm = (double) pall->icept;
  
  /* linear model adjustment of dm */
  if(pall->model == LINEAR && bb > 0) dm += (double) pall->bmax;
  
  /* calculate prior part */
  double s2p = pall->nu0*pall->s20 + syyg;
  if(ng > 0) s2p -= sq(syg)/ng;
  
  /* calculate the data part */
  double ybar = 0.0;
  if(n > 0) ybar = sy/dn;
  double s2y = syy - dn*sq(ybar);
  
  /* combined prior and data part; (valid node check removed) */
  if(ng > 0 && n > 0) 
    s2p += (ng * dn) * sq(ybar - syg/ng) / (ng + dn);
  
  /* parts of student-t */
  *df = pall->nu0 + ng + dn - dm;
  *s2numer = s2y - bb + s2p;
  if(mean) {
    if(ng <= 0.0) *mean = ybar;
    else *mean = (syg + sy)/(ng + dn);
  }
  
  /* return the adjusted m calculation */
  if(s2p_out) *s2p_out = s2p;
  return dm;
}


/*
 * AddDatum:
 *
 * add a new x-y pair to the tree and update the sufficient
 * information using as few operations as possible
 */

Tree* Tree::AddDatum(unsigned int index)
{
  /* pall of the particle this tree belongs to */
  Pall *pall = particle->pall;

  if(isLeaf()) {
    
    /* augment the p-vector */
    p = (int*) realloc(p, sizeof(int)*(n+1));
    p[n] = index;
    n++;

    /* update the classification model */
    if(pall->model == CLASS) counts[(int) pall->y[index]]++;
    else { /* constant and linear models */

      /* update the constant model sufficient statistics */
      double y = pall->y[index];
      if(n == 0) {
	syy = sq(y);
	if(pall->icept) sy = y;
	else assert(pall->model == LINEAR && sy == 0.0);
      } else {
	syy += sq(y);
	if(pall->icept) sy += y;
	else assert(pall->model == LINEAR && sy == 0.0);
      }
	
      /* update linear part of model */
      if(pall->model == LINEAR) {
	if(pall->icept)	CalcLinear();
	else {
	  /* sanity check */
	  assert(XtX && Xty);

	  /* accumulate XtX and Xty */
	  unsigned int m = pall->bmax;
	  double **X = &(pall->X[index]);
	  linalg_dgemm(CblasNoTrans,CblasTrans,m,m,1,1.0,
		       X,m,X,m,1.0,XtX,m);
	  linalg_dgemv(CblasNoTrans,m,1,1.0,X,m,&y,1,1.0,Xty,1);

	  /* calculate XtXi, bmu and bb */
	  assert(XtXi && bmu); 
	  bb = calculate_linear(m, XtX, Xty, XtXi, &ldet_XtXi, bmu);
	}
      }
    }

    /* clear the al calculations */
    if(al) { free(al); al = NULL; }

    /* done */
    return this;

  } else {
    assert(p == NULL);
    if(pall->X[index][var] <= val) return leftChild->AddDatum(index);
    else return rightChild->AddDatum(index);
  }
}


/*
 * RetireDatum:
 *
 * remove the corresponding x-y pair and update the sufficient
 * information using as few operations as possible
 */

Tree* Tree::RetireDatum(unsigned int index, double lambda)
{
  /* pall of the particle this tree belongs to */
  Pall *pall = particle->pall;

  if(isLeaf()) {
    
    /* save old information */
    double y = pall->y[index];
    /* Need to save X-row in linear case */
    unsigned int pi=0;
    for(pi=0; pi<n; pi++) if(p[pi] == (int) index) break;
    assert(pi < n);
    
    /* remove p[index] */
    n--;
    p[pi] = p[n];
    p = (int*) realloc(p, sizeof(int) * n);
    /* indices get adjusted later, after Linear updates */

    /* maybe remove al[index] */
    if(al) {
      al[pi] = al[n];
      if(n == 0) { free(al); al = NULL; } 
      else al = (double*) realloc(al, sizeof(double) * n);
    }

    /* increment the count of retires */
    ng = lambda*ng + 1.0;
    
    /* update the classification model */
    if(pall->model == CLASS) { 
      counts[(int) y]--; /* take from likelihood */
      if(ng == 1.0) { /* special initalization case */
	assert(!gcounts); gcounts = new_zero_vector(pall->nc); 
      }
      scalev(gcounts, pall->nc, lambda);
      gcounts[(int) y] += 1.0; /* add to prior */
    } else { /* constant and linear models */

      /* update the constant model sufficient statistics */
      if(pall->icept) { /* implicit intercept */
	assert(pall->model != LINEAR); 	/* can't do linear in this case */

	/* remove from likelihood */
	if(n == 0) { sy = syy = 0.0; }
	else { syy -= sq(y); sy -= y; }

	/* add to prior */
	syg = lambda*syg + y; 
	syyg = lambda*syyg + sq(y);

      } else { /* explicit intercept */
	if(n == 0) { syy = 0; }
	else syy -= sq(y);
	syyg = lambda*syyg + sq(y);
	assert(pall->model == LINEAR);
      }

      /* update linear part of model */
      if(pall->model == LINEAR) {
	
	unsigned int m = pall->bmax;

	/* special initialization case */
	if(ng == 1.0) {
	  assert(!XtXg); XtXg = new_zero_matrix(m, m);
	  assert(!Xtyg); Xtyg = new_zero_vector(m);
	}

	/* accumulate retire by adding to prior */
	double **X = &(pall->X[index]);
	linalg_dgemm(CblasNoTrans,CblasTrans,m,m,1,1.0,
		     X,m,X,m,lambda,XtXg,m);
	linalg_dgemv(CblasNoTrans,m,1,1.0,X,m,&y,1,lambda,Xtyg,1);

	/* no explicit removal from likelihood since that
	   can be fascilitated in ::grow via ::Calc; therefore bb, 
	   bu and ldet_XtXi don't change */
      }
    }

    /* done */
    return this;

  } else {
    assert(p == NULL);
    if(pall->X[index][var] <= val) return leftChild->RetireDatum(index, lambda);
    else return rightChild->RetireDatum(index, lambda);
  }
}

/* 
 * DecrementP:
 *
 * make sure we move any p[i] pointing to the last position
 * to point to the index position, which is where the 
 * corresponging pall->X[n] will be moved to outside this
 * function 
 */

void Tree::DecrementP(unsigned int oldi, unsigned int newi)
{
  if(isLeaf()) {
    
    /* sanity checks */
    assert(oldi == particle->pall->n - 1);
    assert(newi <= oldi);
    if(newi == oldi) return;
    
    unsigned int i;
    for(i=0; i<n; i++) {
      if(p[i] == (int) oldi) { p[i] = newi; break; }
    }
    assert(i != n);

  } else {
    assert(p == NULL);
    if(particle->pall->X[oldi][var] <= val) 
      return leftChild->DecrementP(oldi, newi);
    else return rightChild->DecrementP(oldi, newi);
  }
}


/* 
 * ReorderP:
 *
 * reorder the pointers p to new entries in pall->X and pall->y;
 * used in the rejuvination steps from Cloud via Particle
 */

void Tree::ReorderP(int *o)
{
  if(isLeaf()) for(unsigned int i=0; i<n; i++) p[i] = o[p[i]]; 
  else {
    assert(p == NULL);
    leftChild->ReorderP(o);
    rightChild->ReorderP(o);
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
    if(!p) return NULL;
    else return new_dup_ivector(p, this->n); 
  } else {
    unsigned int nl, nr;
    int *pl = leftChild->GetP(&nl);
    int *pr = rightChild->GetP(&nr);
    *n = nl + nr;
    if(*n == 0) return NULL;
    else {
      pl = (int*) realloc(pl, sizeof(int) * (*n));
      dupiv(pl+nl, pr, nr); free(pr);
      return pl;
    }
  }
}


/*
 * AccumClass:
 *
 * Accumulate the sufficient statistics for the classification
 * model their corresponding retire stats too
 */

void Tree::AccumClass(unsigned int *counts, double *gcounts)
{
  if(isLeaf()) {
    unsigned int nc = particle->pall->nc;
    for(unsigned int i=0; i<nc; i++)
      counts[i] += this->counts[i];
    if(gcounts && this->gcounts) { /* and for retires */
      for(unsigned int i=0; i<nc; i++)
	gcounts[i] += this->gcounts[i];
    }
  } else {
    /* get stats from both children */
    leftChild->AccumClass(counts, gcounts);
    rightChild->AccumClass(counts, gcounts);
  }
}


/*
 * AccumNg:
 *
 * Accumulate the number of retires from leaves
 */

void Tree::AccumNg(double *ng)
{
  if(isLeaf()) *ng += this->ng;
  else {
    /* get stats from both children */
    leftChild->AccumNg(ng);
    rightChild->AccumNg(ng);
  }
}


/*
 * AccumConst:
 *
 * Accumulate the sufficient statistics for the constant model and
 * their corresponding retire stats too
 */

void Tree::AccumConst(double *sy, double *syy, double *ng, 
		      double *syg, double *syyg)
{
  if(isLeaf()) {
    *sy += this->sy;
    *syy += this->syy;
    *ng += this->ng;
    *syg += this->syg;
    *syyg += this->syyg;
  } else {
    /* get stats from both children */
    leftChild->AccumConst(sy, syy, ng, syg, syyg);
    rightChild->AccumConst(sy, syy, ng, syg, syyg);
  }
}


/*
 * AccumLinear:
 *
 * Accumulate the sufficient statistics for the linear model and
 * their corresponding retire stats too -- only applies when
 * icelt = FALSE -- automatically accumulates to save on
 * memory allocations
 */

void Tree::AccumLinear(double **XtX, double *Xty, double **XtXg, double *Xtyg)
{
  if(isLeaf()) {
    assert(particle->pall->icept == FALSE);
    unsigned int m = particle->pall->bmax;
    linalg_daxpy(m*m,1.0,*(this->XtX),1,*XtX,1);
    linalg_daxpy(m,1.0,this->Xty,1,Xty,1);
    if(XtXg && this->XtXg) linalg_daxpy(m*m,1.0,*(this->XtXg),1,*XtXg,1);
    if(Xtyg && this->Xtyg) linalg_daxpy(m,1.0,this->Xtyg,1,Xtyg,1);
  } else {
    leftChild->AccumLinear(XtX, Xty, XtXg, Xtyg);
    rightChild->AccumLinear(XtX, Xty, XtXg, Xtyg);
  }
}


/*
 * CapRetired:
 *
 * cap prior ESS parameters to forget
 */

void Tree::CapRetired(void)
{
  Pall *pall = particle->pall;
  unsigned int cap = 1;//pall->minp;
  if(pall->model == LINEAR) cap = pall->bmax;
  if(ng > cap) { 
    double dcap = (double) cap;
    double scale = dcap/ng;
    syg = scale * syg;
    syyg = scale * syyg;
    ng = dcap;
    if(pall->model == LINEAR) {
      unsigned int m = pall->bmax;
      scalev(*XtXg, m*m, scale);
      scalev(Xtyg, m, scale);
    }
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
 * GetVar:
 *
 * return the splitting variable used -- must be an internal
 * node
 */

unsigned int Tree::GetVar(void)
{
  assert(!isLeaf());
  return var;
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
#ifdef DEBUG
  assert(isLeaf());
#endif

  /* pall of the particle this tree belongs to */
  Pall *pall = particle->pall;

  /* sanity check and extract double n*/
  assert(n + ng >= pall->minp);
  double dn = (double) n;

  /* special handling of classification */
  if(pall->model == CLASS) {
    double dm = (double) pall->nc;
    if(ng > 0) {
      assert(gcounts);
      double dc = (double) (counts[(int) y] + gcounts[(int) y]);
      return (dc + 1.0/dm)/(1.0 + dn + ng);
    } else return (((double) counts[(int) y]) + 1.0/dm)/(1.0 + dn);
  } /* otherwise constant and linear model */

  /* initial regression calculations */
  double mean, df, s2numer;
  Regression(&mean, &s2numer, &df, NULL);
  if(s2numer <= 0.0) return 0.0;

  /* adjustments under the linear model */
  double xtXtXix = (1.0/(((double) ng) + dn)) * ((double) pall->icept);
  double bmean = 0.0;
  if(pall->model == LINEAR) 
    LinearAdjust(x, &bmean, &xtXtXix, pall->bmaxv, NULL);
  mean += bmean;

  /* student t moments */
  double nt = df/(1.0 + xtXtXix);
  double prec = sqrt(nt/(s2numer));
  double yt = (y-mean)*prec;

  /* calculate the predictve probability */
  double prob = dt(yt, df, 0)*prec; /* prec needed for Jacobian */

  return prob;
}


/*
 * LinearAdjust:
 *
 * adjustments to the latter two arguments, used for prediction,
 * that must be made under the linear model; if XtXix is allocated,
 * then that space is used to hold that value (for use outside), 
 * if not then new space is allocated; if y is allocated then 
 * yXtXix is returned, otherwize zero is returned
 */

double Tree::LinearAdjust(double *x, double *bmean, double *xtXtXix, 
			  double *XtXix, double *y)
{
  /* pall of the particle this tree belongs to */
  Pall *pall = particle->pall;
  
  /* sanity check */
  assert(pall->model == LINEAR);
  assert(XtXix);
  
  /* subtract off mean of X */
  unsigned int m = pall->bmax;
  if(xmean) linalg_daxpy(m,0.0-1.0,xmean,1,x,1);
  else assert(pall->icept == 0);
  
  /* the part of the linear predictor to do with the slope */
  if(bmean) *bmean = linalg_ddot(m, x, 1, bmu, 1);

  /* now the variance of the linear predictor */
  /* t(x) %*% XtXi %*% x */
  zerov(XtXix, m);
  linalg_dsymv(m, 1.0, XtXi, m, x, 1, 0.0, XtXix, 1);
  *xtXtXix += linalg_ddot(m, x, 1, XtXix, 1);
  assert(*xtXtXix > 0);
  if(xmean) linalg_daxpy(m,1.0,xmean,1,x,1); /* mean back in x */

  /* deal with y (for ALC) */
  double yadj;
  if(y) {
    if(xmean) linalg_daxpy(m,0.0-1.0,xmean,1,y,1);
    yadj = linalg_ddot(m, y, 1, XtXix, 1);
    if(xmean) linalg_daxpy(m,1.0,xmean,1,y,1); /* mean back in y */
  } else yadj = 0.0;

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
 * Calc:
 *
 * calculate the sufficient paramters *from scratch*
 * for this tree node based on the data under management 
 * in the corresponding partition
 */

void Tree::Calc(void)
{
  /* pall of the particle this tree belongs to */
  Pall *pall = particle->pall;

  /* special handling of classification */
  if(pall->model == CLASS) { CalcClass(); return; }

  /* calculate intercept and un-scaled error */
  CalcConst();
  
  /* possibly update for the linear model too */
  if(pall->model == LINEAR) {
    if(pall->icept) CalcLinear();
    else {
      /* build X and y */
      unsigned int m = pall->bmax;
      double **X = new_matrix(n, m);
      double *y = new_sub_vector(p, pall->y, n);
      for(unsigned int i=0; i<n; i++) dupv(X[i], pall->X[p[i]], m);

      /* calculate XtX */
      assert(!XtX); XtX = new_zero_matrix(m, m);
      linalg_dgemm(CblasNoTrans,CblasTrans,m,m,n,1.0,
		   X,m,X,m,0.0,XtX,m);
      /* possibly add retire */
      if(ng > 0) { assert(XtXg); linalg_daxpy(m*m, 1.0, *XtXg, 1, *XtX, 1); }

      /* calculate Xty = t(X) %*% y */
      assert(!Xty); Xty = new_zero_vector(m);
      linalg_dgemv(CblasNoTrans,m,n,1.0,X,m,y,1,0.0,Xty,1);
      delete_matrix(X);
      free(y);
      /* possibly add retire */
      if(ng > 0) { assert(Xtyg); linalg_daxpy(m, 1.0, Xtyg, 1, Xty, 1); }

      /* calculate XtXi, bmu and bb */
      assert(!XtXi); XtXi = new_matrix(m, m);
      assert(!bmu); bmu = new_vector(m);
      bb = calculate_linear(m, XtX, Xty, XtXi, &ldet_XtXi, bmu);
    }
  }
}


/*
 * AccumCalc:
 *
 * for prune: gather the sufficient stats from the leaves
 * and then do the relevant calculations if needed
 *
 * NEED TO CODE UP FOR CLASSIFICATION
 */

void Tree::AccumCalc(void)
{
  /* pall of the particle this tree belongs to */
  Pall *pall = particle->pall;

  /* special handling of classification */
  if(pall->model == CLASS) { 
    ng = 0.0; AccumNg(&ng); /* initialize retire count */
    assert(!counts); 
    counts = new_zero_uivector(pall->nc);
    if(ng > 0) { /* maybe allocate retires */
      assert(!gcounts); 
      gcounts = new_zero_vector(pall->nc);
    }
    AccumClass(counts, gcounts);
    return; 
  }

  /* calculate intercept and un-scaled error */
  sy = syy = syg = syyg = 0.0; ng = 0.0; 
  AccumConst(&sy, &syy, &ng, &syg, &syyg);
  
  /* possibly update for the linear model too */
  if(pall->model == LINEAR) {
    if(pall->icept) CalcLinear();
    else {
      /* build XtX and Xty from children */
      unsigned int m = pall->bmax;
      assert(!XtX); XtX = new_zero_matrix(m, m);
      assert(!Xty); Xty = new_zero_vector(m);
      if(ng > 0) {
	assert(!XtXg); XtXg = new_zero_matrix(m, m);
	assert(!Xtyg); Xtyg = new_zero_vector(m);
      }
      AccumLinear(XtX, Xty, XtXg, Xtyg);

      /* calculate XtXi, bmu and bb */
      assert(!XtXi); XtXi = new_matrix(m, m);
      assert(!bmu); bmu = new_vector(m);
      bb = calculate_linear(m, XtX, Xty, XtXi, &ldet_XtXi, bmu);
    }
  }
}


/*
 * CalcClass:
 *
 * do the necessary calculations of sufficient statistics
 * for the classification leaf model
 */

void Tree::CalcClass(void)
{
  Pall *pall = particle->pall;
  assert(pall->model == CLASS);
  if(counts == NULL) counts = new_zero_uivector(pall->nc);
  else zerouiv(counts, pall->nc);
  for(unsigned int i=0; i<n; i++) 
    counts[(int) pall->y[p[i]]]++;
  return;
}


/*
 * CalcConst:
 *
 * do the necessary calculations of the sufficient statistics 
 * for the regression leaf model
 */

void Tree::CalcConst(void)
{
  Pall *pall = particle->pall;
  assert(p);

  /* calculate intercept and un-scaled standard error */
  sy = syy = 0.0;
  if(pall->icept) for(unsigned int i=0; i<n; i++) sy += pall->y[p[i]];
  else assert(pall->model == LINEAR);
  for(unsigned int i=0; i<n; i++) syy += sq(pall->y[p[i]]);
  
  /* check for degenerate but non-zero syy due to numerical imprecision */
  if(syy < DOUBLE_EPS) syy = 0.0;
}


/* 
 * CalcLinear:
 *
 * finish off bmu calculation and XtXi, etc.
 */

void Tree::CalcLinear(void)
{
  /* pall of the particle this tree belongs to */
  Pall *pall = particle->pall;

  /* sanity checks */
  assert(pall->model == LINEAR);
  
  /* build y and subtract the (constant) mean */
  double ybar = sy/((double) n);
  double *y = new_sub_vector(p, pall->y, n);
  for(unsigned int i=0; i<n; i++) y[i] -= ybar;

  /* build X */
  unsigned int m = pall->bmax;
  double **X = new_matrix(n, m);
  for(unsigned int i=0; i<n; i++) dupv(X[i], pall->X[p[i]], m);

  /* maybe subtract off column means */
  if(pall->icept) {
    if(!xmean) xmean = new_vector(m);
    wmean_of_columns(xmean, X, n, m, NULL);
    for(unsigned int i=0; i<n; i++)
      linalg_daxpy(m,0.0-1.0,xmean,1,X[i],1);
  }
  
  /* calculate XtX */
  double **XtX = new_zero_matrix(m, m);
  linalg_dgemm(CblasNoTrans,CblasTrans,m,m,n,1.0,
		 X,m,X,m,0.0,XtX,m);

  /* calculate Xty = t(X) %*% y */
  double *Xty = new_zero_vector(m);
  linalg_dgemv(CblasNoTrans,m,n,1.0,X,m,y,1,0.0,Xty,1);
  delete_matrix(X);
  free(y);

  /* calculate XtXi, bmu and bb */
  if(!XtXi) XtXi = new_matrix(m, m);
  if(!bmu) bmu = new_vector(m);
  bb = calculate_linear(m, XtX, Xty, XtXi, &ldet_XtXi, bmu);

  /* clean up */
  free(Xty);
  delete_matrix(XtX);
}


/*
 * calculate_linear:
 *
 * inverts XtX, calculates its determinant, the MLE esitimator
 * bmu, and returns the component bb of the posterior variance
 * parameter 
 */

double calculate_linear(unsigned int m, double **XtX, double*Xty, 
			double **XtXi, double *ldet_XtXi, double *bmu)
{
  /* calculate inv Gram matrix XtXi = inv(XtX) and determinant */
  assert(XtXi);
  if(!isZero(XtX, m, 1)) { /* non-trivial Gram */
    double ** XtX_chol = new_dup_matrix(XtX, m, m);
    id(XtXi, m);
    int info = linalg_dposv(m, XtX_chol, XtXi);
    if(info != 0) { /* zero values for bad Gram matrix */
      zero(XtXi, m, m);
      zero(XtX, m, m);
      *ldet_XtXi = 0.0;
    } else *ldet_XtXi = 0.0 - log_determinant_chol(XtX_chol, m);
    delete_matrix(XtX_chol);
  } else { /* trivial Gram matrix */
    zero(XtXi, m, m);
    *ldet_XtXi = 0.0; 
  }

  /* now calculate beta-hat (bmu) = inv(XtX) %*% t(X) %*% y */
  assert(bmu);
  zerov(bmu, m);
  linalg_dsymv(m, 1.0, XtXi, m, Xty, 1, 0.0, bmu, 1);

  /* compute: bb = t(bmu) %*% (XtX) %*% bmu */
  double *XtXbmu = new_zero_vector(m); 
  linalg_dsymv(m, 1.0, XtX, m, bmu, 1, 0.0, XtXbmu, 1);
  double bb = linalg_ddot(m, bmu, 1, XtXbmu, 1);
  free(XtXbmu);

  /* return bb */
  return bb;
}


/*
 * Predict:
 *
 * extract the predictive distribution at x in terms of its
 * moments and, perhaps, quantiles; mean, sd and df must be 
 * allocated; The others can be NULL
 */

void Tree::Predict(double *x, double *mean_out, double *sd_out, double *df_out)
{
   if(isLeaf()) {

     /* pall of the particle this tree belongs to */
     Pall *pall = particle->pall;

     /* this is the wrong predict function for classification */
     assert(pall->model != CLASS);

     /* initial regression calculations */
     double mean, df, s2numer;
     Regression(&mean, &s2numer, &df, NULL);

     /* adjustments under the linear model */
     double dn = (double) n;
     double xtXtXix = (1.0/(((double) ng) + dn)) * ((double) pall->icept);
     double bmean = 0.0;
     if(pall->model == LINEAR) 
       LinearAdjust(x, &bmean, &xtXtXix, pall->bmaxv, NULL);
     mean += bmean;

     /* student t parameters */
     *df_out = df;
     double nt = (*df_out)/(1.0 + xtXtXix);
     *sd_out = sqrt((s2numer)/nt);
     assert(*sd_out > 0);
     *mean_out = mean;

     /* myprintf(stdout, "mean=%g, s2numer=%g, df=%g, xtXtXix=%g, sd=%g\n", 
	mean, s2numer, df, xtXtXix, *sd_out); */

   } else {
     if(x[this->var] <= val) 
       return leftChild->Predict(x, mean_out, sd_out, df_out);
     else 
       return rightChild->Predict(x, mean_out, sd_out, df_out);
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
    assert(particle->pall->nc > 0);

    /* double versions of integers */
    unsigned int nc = particle->pall->nc;
    double dm = (double) nc;
    double dn = (double) n;

    /* probability for each class */
    if(ng > 0) { /* using retires */
      assert(gcounts);
      for(unsigned int i=0; i<nc; i++) 
	pred[i] = (((double) (counts[i] + gcounts[i])) + 1.0/dm)/(1.0 + dn + ng);
    } else /* no retires */
      for(unsigned int i=0; i<nc; i++)
	pred[i] = (((double) counts[i]) + 1.0/dm)/(1.0 + dn);

    return;
  } else {
     if(x[this->var] <= val) return leftChild->Predict(x, pred);
     else return rightChild->Predict(x, pred);
   }
}


/*
 * Predict:
 *
 * prediction for classification, for leaves only
 */

void Tree::Predict(double *pred)
{
  assert(isLeaf());
  unsigned int nc = particle->pall->nc;
  assert(nc > 0);

  /* double versions of integers */
  double dm = (double) nc;
  double dn = (double) n;
  
  /* probability for each class */
  if(ng > 0) { /* using retires */
    assert(gcounts);
    for(unsigned int i=0; i<nc; i++) 
      pred[i] = (((double) (counts[i] + gcounts[i])) + 1.0/dm)/(1.0 + dn + ng);
  } else /* no retires */
    for(unsigned int i=0; i<nc; i++)
      pred[i] = (((double) counts[i]) + 1.0/dm)/(1.0 + dn);
}


/*
 * Predict:
 *
 * prediction for classification (for a single class)
 */

double Tree::Predict(double *x, unsigned int cls)
{
  if(isLeaf()) {
    unsigned int nc = particle->pall->nc;
    assert(nc > 0 && cls < nc);

    /* double versions of integers */
    double dm = (double) nc;
    double dn = (double) n;

    /* probability for each class */
    if(ng > 0) { /* using retires */
      assert(gcounts);
      return (((double) (counts[cls] + gcounts[cls])) + 1.0/dm)/(1.0 + dn + ng);
    } else /* no retires */
      return (((double) counts[cls]) + 1.0/dm)/(1.0 + dn);

  } else {
     if(x[this->var] <= val) return leftChild->Predict(x, cls);
     else return rightChild->Predict(x, cls);
   }
}


/* 
 * ECI:
 *
 * calculate the Expected Conditional Improvement for 
 * sequential design at x relative to reference location y
 * y
 */

double Tree::ECI(double *x, double *y, double ymean, double ysd, 
		 double fmin, double ei)
{
   if(isLeaf()) {

     /* pall of the particle this tree belongs to */
     Pall *pall = particle->pall;

     /* ALC not supported for classification */
     assert(pall->model != CLASS);

     /* initial regression calculations */
     double df, s2numer;
     Regression(NULL, &s2numer, &df, NULL);

     /* adjustments to above under the linear model */
     double xtXtXix, ytXtXix;
     double dn = (double) n;
     xtXtXix = (1.0/(ng + dn)) * ((double) pall->icept);
     ytXtXix = xtXtXix;
     if(pall->model == LINEAR) 
       ytXtXix += LinearAdjust(x, NULL, &xtXtXix, pall->bmaxv, y);
     
     /* student t parameters for the predictive at y */
     double sd = sqrt(sq(ysd) - (s2numer/df) * sq(ytXtXix)/(1.0 + xtXtXix));
     assert(sd > 0);
     
     /* do EI calculation based on adjusted sd */
     return EI(ymean, sd, df, fmin);

   } else { /* recurse into leaves */
     if(x[var] <= val && y[var] <= val) 
       return leftChild->ECI(x, y, ymean, ysd, fmin, ei);
     else if(x[var] > val && y[var] > val) 
       return rightChild->ECI(x, y, ymean, ysd, fmin, ei);
     else return ei; /* x and y aren't in same region */
   }
}


/* 
 * ALC:
 *
 * calculate the Active Learning Cohn statistic for 
 * sequential design at x relative to reference location y
 * y
 */

double Tree::ALC(double *x, double *y)
{
   if(isLeaf()) {

     /* pall of the particle this tree belongs to */
     Pall *pall = particle->pall;

     /* ALC not supported for classification */
     assert(pall->model != CLASS);

     /* initial regression calculations */
     double df, s2numer;
     Regression(NULL, &s2numer, &df, NULL);

     /* adjustments to above under the linear model */
     double xtXtXix, ytXtXix;
     double dn = (double) n;
     xtXtXix = (1.0/(ng + dn)) * ((double) pall->icept);
     ytXtXix = xtXtXix;
     if(pall->model == LINEAR) 
       ytXtXix += LinearAdjust(x, NULL, &xtXtXix, pall->bmaxv, y);

     /* wrapping up */
     assert(df-2.0 > 0);
     return (s2numer/(df-2.0)) * sq(ytXtXix)/(1.0 + xtXtXix);

   } else { /* recurse into leaves */
     if(x[var] <= val && y[var] <= val) 
       return leftChild->ALC(x, y);
     else if(x[var] > val && y[var] > val) 
       return rightChild->ALC(x, y);
     else return 0.0; /* x and y aren't in same region */
   }
}


/* 
 * ALC:
 *
 * calculate the Active Learning Cohn statistic for 
 * sequential design by integrating over a rectangle
 * of reference locations
 */

double Tree::ALC(double *x, double **rect)
{
   if(isLeaf()) {

     /* pall of the particle this tree belongs to */
     Pall *pall = particle->pall;

     /* ALC not supported for classification */
     assert(pall->model != CLASS);

     /* initial regression calculations */
     double df, s2numer;
     Regression(NULL, &s2numer, &df, NULL);

     /* constant model ALC calculations */
     double xtXtXix;
     double dn = (double) n;
     xtXtXix = (1.0/(ng + dn)) * ((double) pall->icept);
     double c = xtXtXix;

     /* adjustments to above under the linear model */
     double *y = pall->bmaxv;
     if(pall->model == LINEAR) {
       LinearAdjust(x, NULL, &xtXtXix, y, NULL);
       /* now y is XtXix */
       if(xmean) { /* adjust rectangle by xmean */
	 linalg_daxpy(pall->bmax,0.0-1.0,xmean,1,rect[0],1);
	 linalg_daxpy(pall->bmax,0.0-1.0,xmean,1,rect[1],1);
       }
     }

     /* do the ALC integration analytically */
     double gral = intdot2(pall->bmax, c, y, rect[0], rect[1]);
     
     /* undo any rect adjustments */
     if(pall->model == LINEAR && xmean) {
       linalg_daxpy(pall->bmax,1.0,xmean,1,rect[0],1);
       linalg_daxpy(pall->bmax,1.0,xmean,1,rect[1],1);
     }

     /* factor in the constant part */
     assert(df-2.0 > 0);
     return s2numer*gral/((df-2.0)*(1.0 + xtXtXix));

   } else { /* recurse into leaves */
     if(x[var] <= val) { /* adjust right edge of rect */
       double save = rect[1][var]; rect[1][var] = val;
       double alc = leftChild->ALC(x, rect);
       rect[1][var] = save; /* then put back */
       return alc;
     } else { /* adjust left edge of rect */
       double save = rect[0][var]; rect[0][var] = val;
       double alc = rightChild->ALC(x, rect);
       rect[0][var] = save; /* then put back */
       return alc;
     }
   }
}


/* 
 * ALC:
 *
 * calculate the Active Learning Cohn statistic for 
 * sequential design at the input locations pall->X 
 * by integrating over a rectangle of reference locations
 */

void Tree::ALC(double **rect, double *alc_out)
{
  if(isLeaf()) {

    /* pall of the particle this tree belongs to */
    Pall *pall = particle->pall;

    /* do nothing if no data */
    if(n == 0) return;

    /* maybe ALC is already stored */
    if(al) { add_p_vector(1.0, alc_out, p, 1.0, al, n); return; }
    else al = new_vector(n);
    
    /* ALC not supported for classification */
    assert(pall->model != CLASS);
    
    /* initial regression calculations */
    double df, s2numer;
    Regression(NULL, &s2numer, &df, NULL);
    assert(df-2.0 > 0);
    
    /* constant model ALC calculations */
    double xtXtXix;
    double dn = (double) n;
    xtXtXix = (1.0/(ng + dn)) * ((double) pall->icept);
    double c = xtXtXix;
    
    /* adjust rectangle mean */
    if(pall->model == LINEAR && xmean) {
      linalg_daxpy(pall->bmax,0.0-1.0,xmean,1,rect[0],1);
      linalg_daxpy(pall->bmax,0.0-1.0,xmean,1,rect[1],1);
    }
    
    double *y = pall->bmaxv;
    /* for each X location in this region */
    for(unsigned int i=0; i<n; i++) {
      
      /* adjustments to above under the linear model */
      if(pall->model == LINEAR) { 
	xtXtXix = c;
	LinearAdjust(pall->X[p[i]], NULL, &xtXtXix, y, NULL);
	/* now y is XtXix */
      }
      
      /* do the ALC integration analytically */
      double gral = intdot2(pall->bmax, c, y, rect[0], rect[1]);
      
      /* factor in the constant part */
      alc_out[p[i]] += al[i] = s2numer*gral/((df-2.0)*(1.0 + xtXtXix));
    }
    
    /* undo any rect adjustments */
    if(pall->model == LINEAR && xmean) {
      linalg_daxpy(pall->bmax,1.0,xmean,1,rect[0],1);
      linalg_daxpy(pall->bmax,1.0,xmean,1,rect[1],1);
    }
    
  } else { /* recurse into leaves */
    
    /* left leaf: adjust right edge of rect */
    double save = rect[1][var]; rect[1][var] = val;
    leftChild->ALC(rect, alc_out);
    rect[1][var] = save; /* then put back */
    
    /* right leaf: adjust left edge of rect */
    save = rect[0][var]; rect[0][var] = val;
    rightChild->ALC(rect, alc_out);
    rect[0][var] = save; /* then put back */
  }
}


/* 
 * Entropy:
 *
 * calculate the Entropy statistic for 
 * sequential design at the input locations pall->X 
 */

void Tree::Entropy(double *entropy_out)
{
  if(isLeaf()) {

    /* do nothing if no data */
    if(n == 0) return;

    /* maybe ALC is already stored */
    if(al) { add_p_vector(1.0, entropy_out, p, 1.0, al, n); return; }
    else al = new_vector(n);
    
    /* obtain prediction labels */
    unsigned int nc = particle->pall->nc;
    double *pred = new_vector(nc);
    Predict(pred);

    /* do entropy calculation */
    double entropy = 0.0;
    for(unsigned int j=0; j<nc; j++)
      entropy += 0.0 - pred[j] * log(pred[j]);
    free(pred);

    /* for each X location in this region */
    for(unsigned int i=0; i<n; i++) {

      /* factor in the constant part */
      entropy_out[p[i]] += al[i] = entropy;
    }
    
  } else { /* recurse into leaves */   
    leftChild->Entropy(entropy_out);
    rightChild->Entropy(entropy_out);
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
 * leavesAvgRetired:
 *
 * collecting the average number of retired data points
 * in leaf nodes 
 */

double Tree::leavesAvgRetired(void)
{
  Tree *first, *last;
  int numLeaves = leaves(&first, &last);
  assert(numLeaves > 0);
  double size = 0.0;
  while(first) {
    size += first->ng;
    first = first->next;
  }
  // myprintf(stderr, "%g ", size/((double)numLeaves));
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


/*
 * intdot2:
 *
 * integral (dx, m-dimensional) of the squared dot product 
 * (a + t(x) %*% y)^2 over the m-rectangle given by cbind(a,b)
 */

double intdot2(unsigned int m, double c, double *y, double *a, double *b)
{
  unsigned int i,j;
  double gral, pinv, area;

  /* calculate the area of the rectangle */
  area = 1.0;
  for(i=0; i<m; i++) {
    pinv = b[i] - a[i];
    if(pinv <= 0) { c += b[i]*y[i]; continue; }
    area *= pinv;
  }
  assert(area > 0);

  /* start calculation of the intrgral with a^2 part */
  gral = c*c;

  /* skips the hard integration part if not needed */
  if(y) { 

    for(i=0; i<m; i++) {

      /* then axy part */

      /* pinv <- prod(b[-i] - a[-i]) */
      pinv = (b[i] - a[i]);
      if(pinv <= DOUBLE_EPS) continue;

      /* gral <- gral + c*sum(x^2 * y)|_a^b */
      gral += c*y[i]*(b[i]*b[i] - a[i]*a[i])/pinv;
    
      /* then the (xy)^2 part */
      
      /* gral <- gral + p*(y[i]^2)*(b[i]^3 - a[i]^3)/3 */
      gral += y[i]*y[i]*(b[i]*b[i]*b[i]-a[i]*a[i]*a[i])/(pinv*3.0);
      
      for(j=0; j<i; j++) {
	
	/* pinv <- prod(b[-c(i,j)] - a[-c(i,j)]) */
	pinv = ((b[i] - a[i])*(b[j] - a[j]));
	if(pinv <= DOUBLE_EPS) continue;

	/* gral <- gral + p*y[i]*y[j]*(b[i]^2-a[i]^2)*(b[j]^2-a[j]^2)/2 */
	gral += y[i]*y[j]*(b[i]*b[i]-a[i]*a[i])*(b[j]*b[j]-a[j]*a[j])/(pinv*2.0);
      }
      
    }
  }

  return area*gral;
}
