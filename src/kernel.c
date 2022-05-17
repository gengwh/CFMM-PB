/*
 *  kernelHelm.c
 *  define routines related to the Helmholtz kernel
 *
 */
#include <stdio.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"

extern double kappa;
extern double *nrmX, *nrmY;
extern double epsilon;
extern int nKerl;

/*
 * Calculate the functions G^(k)(r), which are recursively defined as
 *
 *       G^(0)(r) = G(|r|)
 *     G^(k+1)(r) = ( d/dr G^(k)(r) )/r,   k=0..p-1
 *
 * where G(.) is the Green's function. For the Laplace kernel, this becomes
 *
 *       G^(0)(r) = 1/r
 *     G^(k+1)(r) = -(2k+1)*G^(k+1)(r)/r^2,   k=0..p-1
 */
void kernelDC0( double r, int p, double *G ) {
  double r2;
  int k;

  G[0] = 1/r*fourPiI;
  r2 = -1.0/(r*r);

  for ( k=0; k<p; k++ ) {
    G[k+1] = (2*k+1)*r2*G[k];
  }

#if NUMKERNEL
    numKernCplxEval += (p+1);
#endif

} /* kernelD0 */

/*
 *
 * Calculate the functions G^(k)(r), which are recursively defined as
 *
 *       G^(0)(r) = G(|r|)
 *     G^(k+1)(r) = ( d/dr G^(k)(r) )/r,   k=0..p-1
 * Here
 *  (1/r*d/dr)^n  exp( -kappa*r )/r
 * This is based on the recurrence derived in the Contemporary Mathematics
 * paper, and kappa -> i*kappa.
 * Note that the paper has a typo in the recurrence formula
 */
// Would change the subroutine name later
void kernelDS0( double r, int p, double *G ) {
  double rexp, r1, r2, kappa2;
  int k;

  kappa2 = SQR(kappa);
  rexp = exp( -kappa*r );
  r1 = 1.0/r;
  r2 = SQR(r1);

  G[0] = rexp*r1*fourPiI;

  G[1] = -G[0]*(r2 + kappa*r1);
  for ( k=1; k<p; k++ ) {
    G[k+1] = r2*(kappa2*G[k-1] - (2*k+1)*G[k]) ;
  }

#if NUMKERNEL
  numKernCplxEval += (p+1);
#endif

} /* kernelD */

void kernelDG0DGk( double r, int p, double *G0, double *Gk ) {
  double rexp, r1, r2, nr2, kappa2;
  int k;

  kappa2 = SQR(kappa);
  rexp = exp( -kappa*r );
  r1 = 1.0/r;
  r2 = SQR(r1);
  nr2 = -r2;

  G0[0] = r1*fourPiI;
  Gk[0] = rexp*G0[0];

  G0[1] = nr2*G0[0];
  Gk[1] = -Gk[0]*(r2 + kappa*r1);

  for ( k=1; k<p; k++ ) {
    G0[k+1] = (2*k+1)*nr2*G0[k];
    Gk[k+1] = r2*(kappa2*Gk[k-1] - (2*k+1)*Gk[k]);
  }
}

/*
 * kernel of the single layer operator
 */
void kernelKER4( double *x, double *y) {
  double r, r2, expKa, ip0, ipX, ipY, coef;
  double G0, Gk, dG0dx, dGkdx, dG0dy, dGkdy;
  double ddG0dxdy, ddGkdxdy;

  //printf("%f %f %f\n",x[0],x[1],x[2],x[3]);
  r2 = SQR(x[0]) + SQR(x[1]) + SQR(x[2]);
  r = sqrt(r2);
  expKa = exp(-kappa*r);
  //printf("%f %f %f\n",nrmY[0],nrmY[1],nrmY[2]);

  ip0 = nrmX[0]*nrmY[0] + nrmX[1]*nrmY[1] + nrmX[2]*nrmY[2];
  ipX = nrmX[0]*x[0] + nrmX[1]*x[1] + nrmX[2]*x[2];
  ipY = nrmY[0]*x[0] + nrmY[1]*x[1] + nrmY[2]*x[2];

  G0 = 1.0/r*fourPiI;
  Gk = expKa*G0;

  coef = (kappa*r + 1.0)*expKa;

  dG0dy = ipY*G0/r2;
  dGkdy = coef*dG0dy;

  dG0dx = ipX*G0/r2;
  dGkdx = coef*dG0dx;

  ddG0dxdy = (ip0*G0-dG0dy*ipX*3.0)/r2;
  ddGkdxdy = coef*ddG0dxdy-kappa*kappa*expKa*ipX*dG0dy;

  y[0] = G0-Gk;
  y[1] = epsilon*dGkdy - dG0dy;
  y[2] = dGkdx/epsilon - dG0dx;
  y[3] = ddGkdxdy - ddG0dxdy;
//  printf("%.12e %.12e\n",nrmY[0],nrmY[1]);

#if NUMKERNEL
  numKernRealEval++;
#endif
} /* kernelKER4 */

/*
 * kernel of the double layer operator
 */
void kernelRHS( double *x, double *y) {
  double ri, r2, r3i, ip;

  r2 = SQR(x[0]) + SQR(x[1]) + SQR(x[2]);
  ri = 1.0/sqrt(r2);
  r3i = ri/r2;

  ip = nrmX[0]*x[0] + nrmX[1]*x[1] + nrmX[2]*x[2];

  y[0] = ri;
  y[1] = -ip*r3i;

#if NUMKERNEL
  numKernRealEval++;
#endif
} /* kernelRHS */

void kernelPtl( double *x, double *y) {
  double r, r2, expKa, ip0, ipX, ipY, coef;
  double G0, Gk, dG0dy, dGkdy;

  //printf("%f %f %f\n",x[0],x[1],x[2],x[3]);
  r2 = SQR(x[0]) + SQR(x[1]) + SQR(x[2]);
  r = sqrt(r2);
  expKa = exp(-kappa*r);
  //printf("%f %f %f\n",nrmY[0],nrmY[1],nrmY[2]);

  ipY = nrmY[0]*x[0] + nrmY[1]*x[1] + nrmY[2]*x[2];
  //printf("kernel %f\n",ipY);

  G0 = 1.0/r*fourPiI;
  Gk = expKa*G0;

  coef = (kappa*r + 1.0)*expKa;

  dG0dy = ipY*G0/r2;
  dGkdy = coef*dG0dy;

  y[0] = G0-Gk;
  y[1] = epsilon*dGkdy - dG0dy;
  //y[0] = G0;
  //y[1] = -dG0dy;

#if NUMKERNEL
  numKernRealEval++;
#endif
} /* kernelPtl */

/*
 * kernel of the adjoint double layer operator
 */
void kernelC10( double *x, double *y) {
  double r, ip;

  r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  r = pow(r, -1.5);

  ip = nrmX[0]*x[0] + nrmX[1]*x[1] + nrmX[2]*x[2];
  *y = -ip*r;

#if NUMKERNEL
  numKernRealEval++;
#endif
} /* kernel10 */

/*
 * kernel of the double derivative operator
 */
void kernelC11( double *x, double *y) {
  double r, r3, r5, ip0, ipx, ipy;

  r = SQR(x[0]) + SQR(x[1]) + SQR(x[2]);
  r5 = pow(r, -2.5);

  ipx = nrmX[0]*x[0] + nrmX[1]*x[1] + nrmX[2]*x[2];
  ipy = nrmY[0]*x[0] + nrmY[1]*x[1] + nrmY[2]*x[2];

  *y = -3*ipx*ipy*r5;

#if NUMKERNEL
  numKernRealEval++;
#endif
} /* kernel11 */

/*
 * kernel of the single layer operator
 */
void kernelS00( double *x, double *y) {
  double r;

  r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  *y = exp(-kappa*r)/r;

#if NUMKERNEL
  numKernRealEval++;
#endif
} /* kernel00 */

/*
 * kernel of the double layer operator
 */
void kernelS01( double *x, double *y) {
  double r, r2, ip, expr2;

  r2 = SQR(x[0]) + SQR(x[1]) + SQR(x[2]);
  r = sqrt(r2);
  ip = nrmY[0]*x[0] + nrmY[1]*x[1] + nrmY[2]*x[2];
  expr2 = exp(-kappa*r)/r2;

  *y = ip*expr2*(kappa + 1.0/r);

#if NUMKERNEL
  numKernRealEval++;
#endif
} /* kernel01 */


/*
 * kernel of the adjoint double layer operator
 */
void kernelS10( double *x, double *y) {
  double r, r2, ip, expr2;

  r2 = SQR(x[0]) + SQR(x[1]) + SQR(x[2]);
  r = sqrt(r2);
  ip = nrmX[0]*x[0] + nrmX[1]*x[1] + nrmX[2]*x[2];
  expr2 = exp(-kappa*r)/r2;

  *y = -ip*expr2*(kappa + 1.0/r);

#if NUMKERNEL
  numKernRealEval++;
#endif
} /* kernel10 */

/*
 * kernel of the double derivative operator
 */
void kernelS11( double *x, double *y) {
  double r, r2, r3, r4, expr3, ipx, ipy;

  r2 = SQR(x[0]) + SQR(x[1]) + SQR(x[2]);
  r = sqrt(r2);
  r3 = pow(r,3);

  ipx = nrmX[0]*x[0] + nrmX[1]*x[1] + nrmX[2]*x[2];
  ipy = nrmY[0]*x[0] + nrmY[1]*x[1] + nrmY[2]*x[2];

  expr3 = exp(-kappa*r)/r3;

  *y = -ipx*ipy*expr3*(kappa*kappa+3*kappa/r+3/r2);

#if NUMKERNEL
  numKernRealEval++;
#endif
} /* kernel11 */
