#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"

/* blas: matrix times vector */
void dgemv_(char *tr, int *m, int *n, double *alpha, double *A, int *lda,
          double *x, int *incx, double *beta, double *y, int *incy);

void setupDerivs(int order, double *x );
double *panelPotential(int qOrder, double *r, panel *pnl);
void transM2M(ssystem *sys, cube *cbIn, cube *cbOut);

extern double **dG0;     /* workspace for setupDerivs */
extern double **dGk;     /* workspace for setupDerivs */
extern double kappa;
extern double epsilon;
extern void (*kernel)();
extern double **Q2M0, **Q2M1;   /* moments */
extern double *ifact3;
extern int *sgn3;

double partcluster( ssystem *sys, double *G0, double *Gk, cube *cb ) {
  int i, k1, k2, k3, k, j;
  int order = sys->ordM2L[cb->level];
  int nMom  = sys->nMom[order];
  double tmp, tmp1, tmp2, pot;
  double *mom_pot, *mom_dpdn, *lec_k1, *lec_k2;

  mom_pot = cb->mom_pot; mom_dpdn = cb->mom_dpdn;
  for ( tmp1=tmp2=0.,i=0; i<nMom; i++ ) {
    tmp = sgn3[i];
    tmp1 += tmp*mom_dpdn[i]*(G0[i]-Gk[i]);
    tmp2 += tmp*mom_pot[i] *(epsilon*Gk[i]-G0[i]);
    //printf("%f %f %f %e %e\n",ifact3[i],mom_dpdn[i],mom_pot[i],G0[i],Gk[i]);
  }
  pot = tmp1+tmp2;

  return pot;
}

double Treecode( ssystem *sys, cube *cb, double *pos, double *sgm ) {
  int lev = cb->level, order=sys->ordM2L[lev];
  int depth=sys->depth, height=sys->height;
  int nPnls=sys->nPnls, qOrder=sys->maxQuadOrder;
  int i, k;
  double theta = sys->maxSepRatio;
  double dist, pot=0., r[3], *intgr;
  cube *cbKid;
  panel *pnl;
  theta = 0.2;

  for (k=0; k<3; k++) r[k] = pos[k]-cb->x[k];
  dist = sqrt(SQR(r[0])+SQR(r[1])+SQR(r[2]));

  if ( cb->eRad < theta*dist & cb->level>=height ) {
    //printf("target: %f %f %f\n", pos[0], pos[1], pos[2]);
    //printf("cubect: %f %f %f\n", cb->x[0], cb->x[1], cb->x[2]);
    //printf("dist: %f, radius: %f\n", dist, cb->eRad);
    //printf("lev=%d (i,j,k)=(%d,%d,%d)\n", lev, cb->i, cb->j, cb->k);
    //printf("#pnls: %d\n", cb->nPnls);

    //double pot1 = 0.;
    //for ( i=0,pnl=cb->pnls; i<cb->nPnls; i++,pnl=pnl->nextC ) {
    //  intgr = panelPotential(qOrder, pos, pnl);
    //  pot1 += intgr[0]*sgm[pnl->idx+sys->nPnls]+intgr[1]*sgm[pnl->idx];
    //}

    // particle-cluster interaction
    setupDerivs(order, r);
    pot = partcluster(sys, dG0[0], dGk[0], cb);
    //printf("%e %e\n\n",pot1, pot);
    //exit(0);
    //if ( fabs(pot1-pot) > 1e-4 ) exit(0);
    return pot;
  } else {
    if ( cb->level == sys->depth ) {
      //printf("cube lvl: %d, pnl #:%d\n", cb->level, cb->nPnls);
      // direct summation
      for ( i=0, pnl=cb->pnls; i<cb->nPnls; i++,pnl=pnl->nextC ) {
        intgr = panelPotential(qOrder, pos, pnl);
        pot += intgr[0]*sgm[pnl->idx+sys->nPnls]+intgr[1]*sgm[pnl->idx];
      }
      return pot;
    } else {
      for ( i=0; i<cb->nKids; i++ ) {
        cbKid = cb->kids[i];
        pot += Treecode(sys, cbKid, pos, sgm);
      }
      return pot;
    }
  }
}

/*
 * this subroutine apply the particle-cluster interaction
 * work on source term and solvation energy
 * input should specify the kernels (can't apply it on RHS)
*/
//void applyTreecode( ssystem *sys, double *sgm, double *pot, void (*kernelType)() ) {
void applyTreecode( ssystem *sys, double *sgm, double *pot ) {
  cube *cb, *Topcb=sys->cubeList[0];
  int depth=sys->depth, height=sys->height, nPnls=sys->nPnls;
  int nKid, nKid1, iNbr, iPnl, nNbrs, idx, nMom, order;
  int i, k, lev, n, n1, inc = 1;
  double r[3], *self;
  double *x, *y, *k1, *k2;

  /* zero out mom's and lec's */
  for ( lev=depth; lev>=height; lev-- ) {
    nMom  = sys->nMom[sys->ordMom[lev]];
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      for ( k=0; k<nMom;  k++ ) {
        cb->mom_pot[k] = 0.;
        cb->mom_dpdn[k] = 0.;
      }
    }
  }

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
  for ( lev=depth-1; lev>=height; lev-- ) {
    for ( cb=sys->cubeList[lev]; cb != NULL; cb=cb->next ) {
      for ( nKid=0; nKid<cb->nKids; nKid++ ) {
        transM2M(sys, cb->kids[nKid], cb);
      }
    }
  }
  for ( *pot = 0.,i=0; i<sys->nChar; i++ ) {
    *pot += sys->chr[i]*Treecode(sys, Topcb, &sys->pos[3*i], sgm);
  }

}
