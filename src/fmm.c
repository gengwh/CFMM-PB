/*
 *  fmm.c
 *    routines related to fmm-style calculations
 *
 *    for piecewise constant elements
 *
 *  Author: Johannes Tausch
 *  Modified by J. Chen
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"

#define STOREM2L 0
#define SETUPONLY 0

/* blas: matrix times vector */
void dgemv_(char *tr, int *m, int *n, double *alpha, double *A, int *lda,
          double *x, int *incx, double *beta, double *y, int *incy);

void kernelKER4( double *x, double *y);
double *panelIA0(panel *pnlX, panel *pnlY );

double *calcMoments0(ssystem *sys, int order, cube *cb, int job);
void transM2M(ssystem *sys, cube *cbIn, cube *cbOut);
void transM2L(ssystem *sys, double *G0, double *Gk, cube *cbIn, cube *cbOut );
void transL2L(ssystem *sys, cube *cbIn, cube *cbOut );
void kernelDC0( double r, int p, double *G );
void kernelDS0( double r, int p, double *G );
void kernelDG0DGk( double r, int p, double *G0, double *Gk );
void setupDerivs(int order, double *x );

double curvature(double *x, double *h20, double *h11, double *h02);
double paramEllip( panel *pnl, double x, double y, double *r, double *nrm );
void dumpStats(ssystem *sys);
int nrCommonVtx( panel *p, panel *q, int *idxX, int *idxY );

extern double **dG0;     /* workspace for setupDerivs */
extern double **dGk;     /* workspace for setupDerivs */
extern int normErr;
extern void (*kernel)(), (*kernelD)(), (*kernelDC)(), (*kernelDS)();
extern double kappa;
extern double epsilon;

double **Gp0, **Gpk;            /* translation matrices */
double **Q2PK1, **Q2PK2, **Q2PK3, **Q2PK4;
double **Q2M0, **Q2M1;
double **L2P0, **L2P1;

/* preconditioning variables */
extern double **matrixA, *rhs;
extern int *ipiv;


/*
 * setup everything related to FMM-style matrix-vector multiply
 * note that the orders used in the FMM are given by sys->ordM2L[]
 */
void setupFMM(ssystem *sys) {
  cube *cb, *cb1, *kid;
  int depth=sys->depth, height=sys->height;
  int lev, idx, i, j, k;
  int inbr, nNbrs;
  int order, nMoments, nCubesL;
  double r[3];


  for ( nCubesL=0, cb=sys->cubeList[depth]; cb != NULL; cb=cb->next ) {
    nCubesL++;
  }
  CALLOC(Q2M0, nCubesL, double*);
  CALLOC(Q2M1, nCubesL, double*);
  CALLOC(L2P0, nCubesL, double*);
  CALLOC(L2P1, nCubesL, double*);

  order=sys->ordMom[depth];

  /* moments depend on which layer operator we have */
  for ( idx=0, cb=sys->cubeList[depth]; cb!=NULL; cb=cb->next, idx++ ) {
    L2P0[idx] = Q2M0[idx] = calcMoments0(sys, order, cb, 0);
    L2P1[idx] = Q2M1[idx] = calcMoments0(sys, order, cb, 1);
  }


  for ( lev=sys->depth; lev>=sys->height; lev-- ) {
    order = sys->ordMom[lev];
    nMoments  = sys->nMom[order];
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      CALLOC(cb->mom_pot, nMoments, double);
      CALLOC(cb->mom_dpdn, nMoments, double);
      CALLOC(cb->lec_k1, nMoments, double);
      CALLOC(cb->lec_k2, nMoments, double);
      CALLOC(cb->lec_k3, nMoments, double);
      CALLOC(cb->lec_k4, nMoments, double);
    }
  }

  kernelDC = kernelDC0;
  kernelDS = kernelDS0;
  kernelD  = kernelDG0DGk;


} /* setupFMM */


/*
 * Add the nearfield.
 * Pcw constant case. The nearfield coefficients are computed
 * by at Here.
 * Same parameters as applyFMM().
 */
void applyNearfield1(ssystem *sys, double *alpha, double *sgm, double *beta, double *pot) {
  cube *cb, *cb1;
  double *x_pot, *y_pot, *x_dpdn, *y_dpdn, *KER;
  int inbr, nNbrs, idx, nPnls=sys->nPnls;
  int i, j, k, n, n1, inc=1;
  panel *pnlX, *pnlY;

  /* set up kernel */
  kernel = kernelKER4;

  for ( idx=0, cb=sys->cubeList[sys->depth]; cb!=NULL; cb=cb->next ) {
    y_pot = &(pot[cb->pnls->idx]);
    y_dpdn = &(pot[cb->pnls->idx+nPnls]);
    n = cb->nPnls;
    nNbrs = cb->nNbrs;
    for ( inbr=0; inbr<nNbrs; inbr++, idx++ ) {
      cb1 = cb->nbrs[inbr];
      n1 = cb1->nPnls;
      x_pot = &(sgm[cb1->pnls->idx]);
      x_dpdn = &(sgm[cb1->pnls->idx+nPnls]);
      for ( i=0, pnlX=cb->pnls; i<n; i++, pnlX=pnlX->nextC ) {
        for ( j=0, pnlY=cb1->pnls; j<n1; j++, pnlY=pnlY->nextC ) {
          KER = panelIA0(pnlX, pnlY);
          y_pot[i] += (KER[0]*x_dpdn[j] + KER[1]*x_pot[j])*(*alpha);
          y_dpdn[i] += (KER[2]*x_dpdn[j] + KER[3]*x_pot[j])*(*alpha);
        }
      }
    }
  }
} /* applyNearfield1 */


/*
 * FMM-style matrix-vector multiply for panels
 * Parameters
 *    sgm      input density (function on panels)
 *    pot      output potential (function on panels)
 */
void applyFMM(ssystem *sys, double *alpha, double *sgm, double *beta, double *pot) {
  cube *cb, *cb1;
  double *x, *y, *lec;
  double r[3], *self;
  int depth=sys->depth, height=sys->height, nPnls=sys->nPnls;
  int nKid, nKid1, iNbr, iPnl, nNbrs, idx, nMom, order;
  int i, k, lev, n, n1, inc = 1;
  time_t time1, time2;

  /* zero out mom's and lec's */
  for ( lev=depth; lev>=height; lev-- ) {
    nMom  = sys->nMom[sys->ordMom[lev]];
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      for ( k=0; k<nMom;  k++ ) {
        cb->mom_pot[k] = 0.0;
        cb->mom_dpdn[k] = 0.0;
        cb->lec_k1[k] = 0.0;
        cb->lec_k2[k] = 0.0;
        cb->lec_k3[k] = 0.0;
        cb->lec_k4[k] = 0.0;
      }
    }
  }

  /* Q2M transformations */
  time1 = time(&time1);
  nMom = sys->nMom[sys->ordMom[depth]];
  for ( idx=0, cb=sys->cubeList[depth]; cb != NULL; cb=cb->next, idx++ ) {
    x = &(sgm[cb->pnls->idx]);
    n = cb->nPnls;
    y = cb->mom_pot;
    dgemv_(&nChr, &nMom, &n, &one, Q2M1[idx], &nMom, x, &inc, &one, y, &inc);
    x = &(sgm[cb->pnls->idx+nPnls]);
    y = cb->mom_dpdn;
    dgemv_(&nChr, &nMom, &n, &one, Q2M0[idx], &nMom, x, &inc, &one, y, &inc);
  }
  time2 = time(&time2);
  fmmQ2MTime += difftime(time2, time1);

  /* upward pass */
  time1 = time(&time1);
  for ( lev=depth-1; lev>=height; lev-- ) {
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      for ( nKid=0; nKid<cb->nKids; nKid++ ) {
        transM2M(sys, cb->kids[nKid], cb);
      }
    }
  }
  time2 = time(&time2);
  fmmM2MTime += difftime(time2, time1);

  /* Interaction phase */
  time1 = time(&time1);
  for ( idx=0, lev=depth; lev>=height; lev-- ) {
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      nNbrs = cb->nNbrs;
#if !STOREM2L
      order=sys->ordM2L[lev];
#endif
      for ( iNbr=cb->n2Nbrs-1; iNbr>=nNbrs; iNbr--, idx++ ) {
#if STOREM2L
        transM2L(sys, Gp0[idx], Gpk[idx], cb->nbrs[iNbr], cb);
#else
        cb1 = cb->nbrs[iNbr];
        for ( k=0; k<3; k++ ) r[k] = cb->x[k] - cb1->x[k];
        setupDerivs(order, r);
        transM2L(sys, dG0[0], dGk[0], cb1, cb);
#endif
      }
    }
  }
  time2 = time(&time2);
  fmmM2LTime += difftime(time2, time1);

  /* downward pass */
  time1 = time(&time1);
  for ( lev=height; lev<depth; lev++ ) {
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      for ( nKid=0; nKid<cb->nKids; nKid++ ) {
        transL2L(sys, cb, cb->kids[nKid]);
      }
    }
  }
  time2 = time(&time2);
  fmmL2LTime += difftime(time2, time1);

  /* L2P transformations */
  time1 = time(&time1);
  nMom = sys->nMom[sys->ordMom[depth]];
  for ( idx=0, cb=sys->cubeList[depth]; cb != NULL; cb=cb->next, idx++ ) {
    y = &(pot[cb->pnls->idx]);
    n = cb->nPnls;
    x = cb->lec_k1;
    dgemv_(&hChr, &nMom, &n, alpha, L2P0[idx], &nMom, x, &inc, beta, y, &inc);
    x = cb->lec_k2;
    dgemv_(&hChr, &nMom, &n, alpha, L2P0[idx], &nMom, x, &inc, &one, y, &inc);

    y = &(pot[cb->pnls->idx+nPnls]);
    x = cb->lec_k3;
    dgemv_(&hChr, &nMom, &n, alpha, L2P1[idx], &nMom, x, &inc, beta, y, &inc);
    x = cb->lec_k4;
    dgemv_(&hChr, &nMom, &n, alpha, L2P1[idx], &nMom, x, &inc, &one, y, &inc);
  }
  time2 = time(&time2);
  fmmL2PTime += difftime(time2, time1);

  time1 = time(&time1);
  applyNearfield1(sys, alpha, sgm, beta, pot);
  time2 = time(&time2);
  fmmNearTime += difftime(time2, time1);

} /* applyFMM */
