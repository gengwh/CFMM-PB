/*
 *  expan.c
 *  use recurrence formulas to calculate derivatives of Green's function
 *  FMM translation operators
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"



double *Gvals0;   /* workspace for setupDerivs for Coulomb */
double **dG0;     /* workspace for setupDerivs for Coulomb */

double *Gvalsk;   /* workspace for setupDerivs for screened Coulomb */
double **dGk;     /* workspace for setupDerivs for screened Coulomb */

double *fact, *ifact, *fact3, *ifact3, *cfact3;
int ***idx3, *sgn3;
double *convVec1, *convVec2;
double *convVeck1, *convVeck2, *convVeck3;
double *convVecR;
extern double epsilon;
extern double *fcnBuf1, *fcnBuf2, *fcnBuf3;  /* set up by initCalcMoments() */
extern void (*kernelD)(), (*kernelDC)(), (*kernelDS)();


/*
 * Generate a look-up table for finding the (i1,i2,i3)-th entry of an
 * array with three indices in a 1-d buffer. Vector are stored in
 * row-major format (i.e., first index varies most slowly, C-convention).
 */
void mkIndex(int order) {
  int i1, i2, i3, i, nn;

  CALLOC(idx3, order+1, int**);
  for (i1=0; i1<=order; i1++) {
    CALLOC(idx3[i1], order+1, int*);
    for (i2=0; i2<=order; i2++) {
      CALLOC(idx3[i1][i2], order+1, int);
    }
  }

  for ( i=nn=0; nn<=order; nn++ ) {
    for ( i1=0; i1<=nn; i1++ ) {
      for ( i2=0; i2<=nn-i1; i2++, i++ ) {
        i3=nn-i1-i2;
        idx3[i1][i2][i3] = i;
      }
    }
  }
} /* mkIndex */





/*
 * initialize everything related to the (Taylor) expansion.
 * Taylor coefficients using the trapezoidal rule and FFTs
 */
void initExpan(ssystem *sys) {
  int order, nMoments;
  int p, k, n, k1, k2, k3;
  int sgn=1;

  order = sys->maxOrder;
  nMoments = sys->nMom[order];
  /*
   * fact[i] = i!                    factorial
   * ifact[i] = 1/fact[i]
   * fact3[i] = i_1! * i_2! * i_3!   3D factorial in moments order
   * ifact3[i] = 1/fact3[i]
   * sgn3[i] = (-1)^i1 * (-1)^i2 * (-1)^i3
   */
  CALLOC(fact, order+4, double);
  CALLOC(ifact, order+4, double);
  CALLOC(fact3, nMoments, double);
  CALLOC(ifact3, nMoments, double);
  CALLOC(sgn3, nMoments, int);

  for ( fact[0]=ifact[0]=1.0, k=1; k<=order+3; k++ ) {
    fact[k] = fact[k-1]*k;
    ifact[k] = 1.0/fact[k];
  }

  for ( k=n=0; n<=order; n++, sgn*=-1 ) {
    for ( k1=0; k1<=n; k1++ ) {
      for ( k2=0; k2<=n-k1; k2++, k++ ) {
        k3=n-k1-k2;
        fact3[k] = fact[k1]*fact[k2]*fact[k3];
        ifact3[k] = 1.0/fact3[k];
        sgn3[k] = sgn;
      }
    }
  }

  mkIndex(order);

  CALLOC(Gvals0, order+1, double);
  CALLOC(dG0, order+1, double*);
  CALLOC(Gvalsk, order+1, double);
  CALLOC(dGk, order+1, double*);
  for ( k=0; k<=order; k++ ) {
    p = order-k;
    CALLOC(dG0[k], sys->nMom[p], double);
    CALLOC(dGk[k], sys->nMom[p], double);
  }

  CALLOC(convVecR, nMoments, double);
  CALLOC(convVec1, nMoments, double);
  CALLOC(convVec2, nMoments, double);
  CALLOC(convVeck1, nMoments, double);
  CALLOC(convVeck2, nMoments, double);
  CALLOC(convVeck3, nMoments, double);
} /* initTaylorCoeffs */



/*
 * Calculate the derivatives D^\alpha f(x),
 * using recurrence for higher derivatives of radially symmetric functions.
 * The result is stored in the global vector dG[0][*].
 */
void setupDerivs(int order, double *x ) {
  double r;
  int p, p1, iRow, iRow1, i1, i2, i3, idx, idx1, idx2, n;

  r = sqrt(SQR(x[0]) + SQR(x[1]) + SQR(x[2]));
  kernelD( r, order, Gvals0, Gvalsk );

  /* 0th order derivatives */
  for (p=0; p<=order; p++ ) {
    dG0[p][0] = Gvals0[p];
    dGk[p][0] = Gvalsk[p];
  }

  /* 1st order derivatives */
  for (p=0; p<order; p++ ) {
    p1 = p + 1;
    dG0[p][1] = dG0[p1][0]*x[2];
    dG0[p][2] = dG0[p1][0]*x[1];
    dG0[p][3] = dG0[p1][0]*x[0];
    dGk[p][1] = dGk[p1][0]*x[2];
    dGk[p][2] = dGk[p1][0]*x[1];
    dGk[p][3] = dGk[p1][0]*x[0];
  }

  /* higher order derivatives */
  for ( iRow=2; iRow<=order; iRow++ ) {
    for ( p=0; p<=order-iRow; p++ ) {
      p1 = p + 1;
      idx = idx3[0][0][iRow];  /* iRow-th order derivatives start here */
      iRow1 = iRow-1;
      /* i1 = i2 = 0, i3 = iRow */
      idx1 = idx3[0][0][iRow1];
      idx2 = idx3[0][0][iRow-2];
      dG0[p][idx] = dG0[p1][idx1]*x[2] + dG0[p1][idx2]*iRow1;
      dGk[p][idx] = dGk[p1][idx1]*x[2] + dGk[p1][idx2]*iRow1;
      idx++;
      /* i1 = 0, i2 = 1, i3 = iRow-1 */
      idx1 = idx3[0][0][iRow1];
      dG0[p][idx] = dG0[p1][idx1]*x[1];
      dGk[p][idx] = dGk[p1][idx1]*x[1];
      idx++;
      /* i1 = 0, i2 >= 2 */
      for ( i2=2; i2<=iRow; i2++, idx++ ) {
        i3 = iRow-i2;
        idx1 = idx3[0][i2-1][i3];
        idx2 = idx3[0][i2-2][i3];
        dG0[p][idx] = dG0[p1][idx1]*x[1] + dG0[p1][idx2]*(i2-1);
        dGk[p][idx] = dGk[p1][idx1]*x[1] + dGk[p1][idx2]*(i2-1);
      }
      /* i1 = 1 */
      for ( i2=0; i2<=iRow1; i2++, idx++ ) {
        i3 = iRow1-i2;
        idx1 = idx3[0][i2][i3];
        dG0[p][idx] = dG0[p1][idx1]*x[0];
        dGk[p][idx] = dGk[p1][idx1]*x[0];
      }
      /* i1 >= 2 */
      for ( i1=2; i1<=iRow; i1++ ) {
        for ( i2=0; i2<=iRow-i1; i2++, idx++ ) {
          i3 = iRow-i1-i2;
          idx1 = idx3[i1-1][i2][i3];
          idx2 = idx3[i1-2][i2][i3];
          dG0[p][idx] = dG0[p1][idx1]*x[0] + dG0[p1][idx2]*(i1-1);
          dGk[p][idx] = dGk[p1][idx1]*x[0] + dGk[p1][idx2]*(i1-1);
        }
      }
    }
  }
  //int i;
  //for (i=0;i<4;i++) printf("%f %f\n",dG0[0][i],dGk[0][i]);
  //exit(0);
} /* setupDerivs */


/*
 * computes the convolution
 *   c_\alpha = \sum_{\beta+\gamma = \alfa} (binomial coef)*a_\beta * b_\gamma
 *                      a_1!           a_2!          a_3!
 * (binomial coef)=--------------*--------------*--------------
 *                 b_1!(a_1-b_1)! b_2!(a_2-b_2)! b_3!(a_3-b_3)!
 * Here, beta is with coef before get in
 * for \abs{\alpha} <= order
 * This routine does the same as convolution() in moments.c, with
 * the only difference that convM2M is additive.
 */
void convM2M(int order, double *a, double *b, double *c, double *d, double *e){
  int i, j, i1, i2, i3, j1, j2, j3, k1, k2, k3, n, m;
  double tmp, tmp1, tmp2;

  for ( i=n=0; n<=order; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        i3=n-i1-i2;
        for ( tmp1=0.0,tmp2=0.0,k1=0; k1<=i1; k1++ ) {
          j1 = i1 - k1;
          for ( k2=0; k2<=i2; k2++ ) {
            j2 = i2 - k2;
            for ( k3=0; k3<=i3; k3++ ) {
              j3 = i3 - k3;
              tmp1 += a[idx3[j1][j2][j3]]*b[idx3[k1][k2][k3]];
              tmp2 += a[idx3[j1][j2][j3]]*d[idx3[k1][k2][k3]];
            }
          }
        }
        c[i] += tmp1;
        e[i] += tmp2;
      }
    }
  }
} /* convM2M */



/*
 * compute the convolution
 *   c_\alpha =
 *     \sum_{\beta: \abs{\beta+\alpha}\leq ordIn} a_\beta b_{\beta+\alpha}
 * for \abs{\alpha} \leq ordOut
 * This routine is additive.
 */
void convM2L(int ordIn, int ordOut, double *c1, double *c2,
             double *veck1, double *veck2, double *veck3,
             double *l1, double *l2, double *l3, double *l4) {
  int i, j, m, n, i1, i2, i3, j1, j2, j3, k1, k2, k3;
  double tmp, tmp1, tmp2, tmp3, tmp4;

  for ( i=n=0; n<=ordOut; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        i3=n-i1-i2;
        for ( tmp1=0.0,tmp2=0.0,tmp3=0.0,tmp4=0.0,j=m=0; m<=ordIn-n; m++) {
          for ( j1=0; j1<=m; j1++ ) {
            k1 = i1 + j1;
            for ( j2=0; j2<=m-j1; j2++, j++ ) {
              j3 = m-j1-j2;
              k2 = i2 + j2;
              k3 = i3 + j3;

              tmp1 += c2[j]*veck1[idx3[k1][k2][k3]];
              tmp2 += c1[j]*veck2[idx3[k1][k2][k3]];
              tmp3 += c2[j]*veck3[idx3[k1][k2][k3]];
              tmp4 += c1[j]*(-veck1[idx3[k1][k2][k3]]);
            }
          }
        }
        l1[i] += tmp1;
        l2[i] += tmp2;
        l3[i] += tmp3;
        l4[i] += tmp4;
      }
    }
  }
} /* convM2L */

/*
 * difference of L2L and M2L is the binomial coefficients for L2L
 */
void convL2L(int ordIn, int ordOut, double *a,
             double *In1, double *In2, double *In3, double *In4,
             double *Out1, double *Out2, double *Out3, double *Out4 ) {
  int i, j, m, n, i1, i2, i3, j1, j2, j3, k1, k2, k3;
  double tmp, tmp1, tmp2, tmp3, tmp4;

  for ( i=n=0; n<=ordOut; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        i3=n-i1-i2;
        for ( tmp1=0.0,tmp2=0.0,tmp3=0.0,tmp4=0.0,j=m=0; m<=ordIn-n; m++) {
          for ( j1=0; j1<=m; j1++ ) {
            k1 = i1 + j1;
            for ( j2=0; j2<=m-j1; j2++, j++ ) {
              j3 = m-j1-j2;
              k2 = i2 + j2;
              k3 = i3 + j3;
              tmp1 += a[j]*In1[idx3[k1][k2][k3]];
              tmp2 += a[j]*In2[idx3[k1][k2][k3]];
              tmp3 += a[j]*In3[idx3[k1][k2][k3]];
              tmp4 += a[j]*In4[idx3[k1][k2][k3]];
            }
          }
        }
        Out1[i] += tmp1;
        Out2[i] += tmp2;
        Out3[i] += tmp3;
        Out4[i] += tmp4;
      }
    }
  }
} /* convL2L */





/*
 * transM2M() translates a moment vector to a new center.
 * Parameters
 *    cbIn:   input cube
 *    cbOut:  output cube
 * The input order can be smaller than the output order. In that case,
 * the higher order input moments are assumed to be zero, i.e., the higher
 * output moments are not computed exactly.
 * Note: transM2M is additive.
 */
void transM2M(ssystem *sys, cube *cbIn, cube *cbOut) {
  int i, n, i1, i2, i3;
  int ordIn = sys->ordMom[cbIn->level], ordOut = sys->ordMom[cbOut->level];
  int nMomIn = sys->nMom[ordIn], nMomOut = sys->nMom[ordOut];
  double *momIn_pot = cbIn->mom_pot, *momOut_pot = cbOut->mom_pot;
  double *momIn_dpdn = cbIn->mom_dpdn, *momOut_dpdn = cbOut->mom_dpdn;
  double trns[3];

  for ( i=0; i<3; i++ ) trns[i] = cbIn->x[i] - cbOut->x[i];


  /* setup convolution vectors */
  fcnBuf1[0] = fcnBuf2[0] = fcnBuf3[0] = 1.0;
  for ( i=1; i<=ordOut; i++ ) {
    fcnBuf1[i] = fcnBuf1[i-1]*trns[0];
    fcnBuf2[i] = fcnBuf2[i-1]*trns[1];
    fcnBuf3[i] = fcnBuf3[i-1]*trns[2];
  }

  for ( i=n=0; n<=ordOut; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        i3=n-i1-i2;
        convVecR[i] = fcnBuf1[i1]*fcnBuf2[i2]*fcnBuf3[i3];
        convVecR[i] *= ifact3[i];
      }
    }
  }

  memcpy(convVec1, momIn_pot, nMomIn*sizeof(double));
  memcpy(convVec2, momIn_dpdn, nMomIn*sizeof(double));
  for ( i=nMomIn; i<nMomOut; i++ ) {
    convVec1[i] = 0.0;
    convVec2[i] = 0.0;
  }
  convM2M(ordOut, convVecR, convVec1, momOut_pot, convVec2, momOut_dpdn);


} /* transM2M */


/*
 * transL2L() translates a local expansion to a new center.
 * Parameters:
 *    cbIn:  input cube
 *    cbOut: output cube
 * The output order can be smaller than the input order. In that case,
 * the higher order output coefficients are simply not computed
 * Note that transL2L is additive
 */
void transL2L(ssystem *sys, cube *cbIn, cube *cbOut ) {
  int i, n, i1, i2, i3;
  int ordIn = sys->ordMom[cbIn->level], ordOut = sys->ordMom[cbOut->level];
  int nMomIn = sys->nMom[ordIn], nMomOut = sys->nMom[ordOut];
  double *lecInk1, *lecInk2, *lecInk3, *lecInk4;
  double *lecOutk1, *lecOutk2, *lecOutk3, *lecOutk4;
  double trns[3];

  lecInk1 = cbIn->lec_k1; lecInk2 = cbIn->lec_k2;
  lecInk3 = cbIn->lec_k3; lecInk4 = cbIn->lec_k4;
  lecOutk1 = cbOut->lec_k1; lecOutk2 = cbOut->lec_k2;
  lecOutk3 = cbOut->lec_k3; lecOutk4 = cbOut->lec_k4;

  for ( i=0; i<3; i++ ) trns[i] = cbOut->x[i] - cbIn->x[i];

  /* setup convolution vectors */

  fcnBuf1[0] = fcnBuf2[0] = fcnBuf3[0] = 1.0;
  for ( i=1; i<=ordOut; i++ ) {
    fcnBuf1[i] = fcnBuf1[i-1]*trns[0];
    fcnBuf2[i] = fcnBuf2[i-1]*trns[1];
    fcnBuf3[i] = fcnBuf3[i-1]*trns[2];
  }

  for ( i=n=0; n<=ordIn; n++ ) {
    for ( i1=0; i1<=n; i1++ ) {
      for ( i2=0; i2<=n-i1; i2++, i++ ) {
        i3=n-i1-i2;
        convVecR[i] = fcnBuf1[i1]*fcnBuf2[i2]*fcnBuf3[i3];
        convVecR[i] *= ifact3[i];
      }
    }
  }

  convL2L(ordIn, ordOut, convVecR, lecInk1, lecInk2, lecInk3, lecInk4,
          lecOutk1, lecOutk2, lecOutk3, lecOutk4);

} /* transL2L */



/*
 * translate a moment vector to a local expansion vector
 * the vector G[] contains the derivatives
 * Note that the M2L translation is ADDED to the vector lec.
 * Contrary to applyM2L() complex conjugation is NOT necessary.
 */
void transM2L(ssystem *sys, double *G0, double *Gk, cube *cbIn, cube *cbOut ) {
  int i;
  int order = sys->ordM2L[cbIn->level];
  int nMom  = sys->nMom[order];
  double tmp;
  double *mom_pot, *mom_dpdn;
  double *lec_k1, *lec_k2, *lec_k3, *lec_k4;

  mom_pot = cbIn->mom_pot; mom_dpdn = cbIn->mom_dpdn;
  lec_k1 = cbOut->lec_k1; lec_k2 = cbOut->lec_k2;
  lec_k3 = cbOut->lec_k3; lec_k4 = cbOut->lec_k4;

  for ( i=0; i<nMom; i++ ) {
    convVec1[i] = sgn3[i]*mom_pot[i];
    convVec2[i] = sgn3[i]*mom_dpdn[i];
    convVeck1[i] = G0[i]-Gk[i];
    convVeck2[i] = epsilon*Gk[i]-G0[i];
    convVeck3[i] = G0[i]-Gk[i]/epsilon;
  }
  convM2L(order, order, convVec1, convVec2,
          convVeck1, convVeck2, convVeck3,
          lec_k1, lec_k2, lec_k3, lec_k4);

} /* transM2L */
