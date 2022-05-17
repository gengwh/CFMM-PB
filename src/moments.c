/*
 * moments.c:
 *    calculate three variable moments
 *
 * Important:  the moments of a panel are defined as
 *
 * m^\alpha =
 *   {1\over \alpha! |S_i|^{1/2}} \int_{pnl} (x - \bar x)^\alpha dS_x
 *
 *
 *  Written by Johannes Tausch
 *  Modified by J. Chen
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"

extern double *fact3, *ifact3, *fact, *ifact;  /* setup by initExpan() */
extern int ***idx3, *sgn3;
/* globals in this file */
double *convVec1, *convVec2, *convVec3; /* workspace for convolutions */
double *fcnBuf1, *fcnBuf2, *fcnBuf3;
double *MomWrk, *MomLoc;


/*
 * initialize stuff needed to evaluate the moments
 */
void initCalcMoments(ssystem *sys) {
  int k,l, nMomType;
  int order=sys->maxOrder, nMoments = sys->nMom[order];
  int ordL=sys->ordMom[sys->depth], nMomL=sys->nMom[ordL];



  CALLOC(MomWrk, nMomL, double);
  CALLOC(MomLoc, nMomL, double);
  CALLOC(convVec1, nMomL, double);
  CALLOC(convVec2, nMomL, double);
  CALLOC(convVec3, nMomL, double);
  CALLOC(fcnBuf1, order+1, double);
  CALLOC(fcnBuf2, order+1, double);
  CALLOC(fcnBuf3, order+1, double);

} /* initCalcMoments */


/*
 * computes the convolution
 *   c_\alpha = \sum_{\beta+\gamma = \alfa} a_\beta * b_\gamma
 * for \abs{\alpha} <= order
 *
 */
void convolution(int order, double *a, double *b, double *c){
  int i, j, i1, i2, i3, j1, j2, j3, k1, k2, k3, n, m;
  double tmp;

  for ( i=n=0; n<=order; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        i3=n-i1-i2;
        for ( tmp=0.0,k1=0; k1<=i1; k1++ ) {
          j1 = i1 - k1;
          for ( k2=0; k2<=i2; k2++ ) {
            j2 = i2 - k2;
            for ( k3=0; k3<=i3; k3++ ) {
              j3 = i3 - k3;
              tmp += a[idx3[j1][j2][j3]]*b[idx3[k1][k2][k3]];
            }
          }
        }
        c[i] = tmp;
      }
    }
  }
} /* convolution */

/*
 * setup the M2M translation vector.
 */
void setupConVect(int order, double *trns, double *convVec) {
  int i, i1, i2, i3, n;

  fcnBuf1[0] = fcnBuf2[0] = fcnBuf3[0] = 1.0;

  for ( i=1; i<=order; i++ ) {
    fcnBuf1[i] = fcnBuf1[i-1]*trns[0];
    fcnBuf2[i] = fcnBuf2[i-1]*trns[1];
    fcnBuf3[i] = fcnBuf3[i-1]*trns[2];
  }

  for ( i=n=0; n<=order; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        i3=n-i1-i2;
        convVec[i] = fcnBuf1[i1]*fcnBuf2[i2]*fcnBuf3[i3];
        convVec[i] *= ifact3[i];
      }
    }
  }
} /* setupConVect */




/*
 * calculate all noments of one panel up to order given by order,
 * using the previously calculated moment vector of the same panel
 * which is passed in Mom.
 * returns moment vector in parameter Nom
 * Note that the input moments need only be defined up to order-1
 */
void calcOneNoment(ssystem *sys, int order, panel *pnl, double *Mom, double *Nom) {
  int i, i1, i2, i3, n;
  int nMom = sys->nMom[order];
  double *nrm = pnl->normal;

  for ( i=0; i<nMom; i++ ) Nom[i] = 0.0;

  for ( i=n=0; n<order; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        i3=n-i1-i2;
        Nom[idx3[i1+1][i2][i3]] += nrm[0]*Mom[i];
        Nom[idx3[i1][i2+1][i3]] += nrm[1]*Mom[i];
        Nom[idx3[i1][i2][i3+1]] += nrm[2]*Mom[i];
      }
    }
  }

} /* calcOneNoment */




/*
 * calculate all moments of one panel.
 * The center is the panel center.
 * The moment vector is returned in parameter Mom
 * parameters:
 *     pnl     panel of which the moments are calculated
 *     order   order of moments
 *     Mom     returned moment
 *     job     see calcMoments()
 */
void calcOneMoment0(ssystem *sys, int order, panel *pnl, double *Mom, int job) {
  int i, i1, i2, i3, ind1, ind2, ind3, n;
  int order1;
  double area2 = 2.0*pnl->area;
  double *nrm = pnl->normal;
  double a1[3], a2[3];
  double *v0 = pnl->vtx[0], *v1 = pnl->vtx[1], *v2 = pnl->vtx[2];
  double *MomBuf, fac;

  if ( job == 0 ) {
    order1 = order;
    MomBuf = Mom;
  }
  else {
    order1 = order - 1;
    MomBuf = MomWrk;
  }

  /* jtdeb: old edge convention, could be re-written to match new one. */
  for ( i=0; i<3; i++ ) {
    a1[i] = v1[i] - v0[i];
    a2[i] = v2[i] - v0[i];
  }

  /* setup local moment vector */
  setupConVect(order1, a1, convVec1);
  setupConVect(order1, a2, convVec2);

  for ( i=n=0; n<=order1; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        convVec1[i] *= fact[n];
        convVec2[i] *= fact[n];
      }
    }
  }

  convolution(order, convVec1, convVec2, MomBuf);
  for ( i=n=0; n<=order; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        MomBuf[i] *= area2*ifact[n+2];
      }
    }
  }

  if ( job == 0 ) {
    return; /* MomBuf and Mom point to the same memory cell */
  }
  else if ( job == 1 ) {
    calcOneNoment(sys, order, pnl, MomBuf, Mom);
  }
} /* calcOneMoment0 */



/*
 * Allocate and calculate the Moments/Noments
 * for all box-functions contained in a  finest-level cube.
 * That is, this routine allocates and computes the Q2M matrix
 * Parameters
 *    cb        cube of which Q2M matrix is computed
 *    order     maximal order of moments
 *    job       if job==0 calculate moments only
 *              if job==1 calculate noments only
 *
 * Returns pointer to Q2M matrix
 */
double *calcMoments0(ssystem *sys, int order, cube *cb, int job) {
  int idx;
  int nMom = sys->nMom[order], nPnls=cb->nPnls, ldM;
  double trns[3], *ctr=cb->x, *v0, *nrm;
  double *Moments, *MomTrn1;
  panel *pnl;

  ldM = nMom;
  CALLOC(Moments, ldM*nPnls, double);

  for (idx=0, pnl=cb->pnls; idx<nPnls; pnl=pnl->nextC, idx++ ) {
    v0 = pnl->vtx[0];
    calcOneMoment0( sys, order, pnl, MomLoc, job );
    /* translate to center of cube */
    trns[0] = v0[0]-ctr[0]; trns[1] = v0[1]-ctr[1]; trns[2] = v0[2]-ctr[2];
    setupConVect(order, trns, convVec1);
    convolution(order, convVec1, MomLoc, &Moments[idx*ldM]);
  }

  return Moments;
} /* calcMoments0 */
