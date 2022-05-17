/*
 * numQuad.c
 *    calculation of direct panel interactions in the Galerkin method
 *    via numerical quadrature.
 *
 *   Copyright  Johannes Tausch
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gkGlobal.h"
#include "gk.h"

/* function prototypes */
void Jacobi(int, double, double, double*, double*);  /* Gauss quad rule */

/* global variables */
double *nrmX, *nrmY;
double **tLegA, **wLegA;
//double intgrY[9];
int maxQuadOrder, nKerl;
double *intgrP, *intgr, *intgrY, intgrRHS[2], intgrPtl[2];
double *fcn1, *fcn2, *fcn3, *fcn4, *fcn5, *fcn6;

void kernelRHS( double *x, double *y);
void kernelPtl( double *x, double *y);

/*
 *  current kernel must be turned on before using quadrature
 */
extern void (*kernel)();


/*
 *  initialize Legendre quadrature rules of various orders
 */
#if 0  /* using the pierreQP library*/
void initQuad(ssystem *sys) {
  double *tLeg, *wLeg;
  int qOrder;

  maxQuadOrder = sys->maxQuadOrder;
  /* abscissas, weights of Gauss-Legendre, scaled for [0,1] */
  CALLOC(tLegA, maxQuadOrder+1, double*);
  CALLOC(wLegA, maxQuadOrder+1, double*);

  for ( qOrder=1; qOrder<=maxQuadOrder; qOrder++ ) {
    CALLOC(tLeg, qOrder, double);
    CALLOC(wLeg, qOrder, double);
    tLegA[qOrder] = tLeg;
    wLegA[qOrder] = wLeg;
    Jacobi(qOrder, 0.0, 0.0, tLeg, wLeg);
  }
}
#else  /* using a table */
void initQuad(ssystem *sys) {
  double *tLeg, *wLeg;
  int qOrder;

  maxQuadOrder = sys->maxQuadOrder;
  nKerl = sys->nKerl;

  CALLOC(intgrP, nKerl, double);
  CALLOC(intgr, nKerl, double);
  CALLOC(intgrY, nKerl, double);
  CALLOC(fcn1, nKerl, double);
  CALLOC(fcn2, nKerl, double);
  CALLOC(fcn3, nKerl, double);
  CALLOC(fcn4, nKerl, double);
  CALLOC(fcn5, nKerl, double);
  CALLOC(fcn6, nKerl, double);

  CALLOC(tLegA, 11, double*);
  CALLOC(wLegA, 11, double*);
  for ( qOrder=1; qOrder<=10; qOrder++ ) {
    CALLOC(tLeg, qOrder, double);
    CALLOC(wLeg, qOrder, double);
    tLegA[qOrder] = tLeg;
    wLegA[qOrder] = wLeg;
  }

  /* order=1 */
  tLegA[1][0]=0.500000000000000;  wLegA[1][0]=1.000000000000000;

  /* order=2 */
  tLegA[2][0]=0.211324865405187;  wLegA[2][0]=0.500000000000000;
  tLegA[2][1]=0.788675134594813;  wLegA[2][1]=0.500000000000000;

  /* order=3 */
  tLegA[3][0]=0.112701665379258;  wLegA[3][0]=0.277777777777778;
  tLegA[3][1]=0.500000000000000;  wLegA[3][1]=0.444444444444444;
  tLegA[3][2]=0.887298334620742;  wLegA[3][2]=0.277777777777778;

  /* order=4 */
  tLegA[4][0]=0.069431844202974;  wLegA[4][0]=0.173927422568727;
  tLegA[4][1]=0.330009478207572;  wLegA[4][1]=0.326072577431273;
  tLegA[4][2]=0.669990521792428;  wLegA[4][2]=0.326072577431273;
  tLegA[4][3]=0.930568155797026;  wLegA[4][3]=0.173927422568727;

  /* order=5 */
  tLegA[5][0]=0.046910077030668;  wLegA[5][0]=0.118463442528095;
  tLegA[5][1]=0.230765344947158;  wLegA[5][1]=0.239314335249683;
  tLegA[5][2]=0.500000000000000;  wLegA[5][2]=0.284444444444444;
  tLegA[5][3]=0.769234655052841;  wLegA[5][3]=0.239314335249683;
  tLegA[5][4]=0.953089922969332;  wLegA[5][4]=0.118463442528095;

  /* order=6 */
  tLegA[6][0]=0.033765242898424;  wLegA[6][0]=0.085662246189585;
  tLegA[6][1]=0.169395306766868;  wLegA[6][1]=0.180380786524070;
  tLegA[6][2]=0.380690406958402;  wLegA[6][2]=0.233956967286346;
  tLegA[6][3]=0.619309593041598;  wLegA[6][3]=0.233956967286346;
  tLegA[6][4]=0.830604693233132;  wLegA[6][4]=0.180380786524069;
  tLegA[6][5]=0.966234757101576;  wLegA[6][5]=0.085662246189585;

  /* order=7 */
  tLegA[7][0]=0.025446043828621;  wLegA[7][0]=0.064742483084435;
  tLegA[7][1]=0.129234407200303;  wLegA[7][1]=0.139852695744638;
  tLegA[7][2]=0.297077424311301;  wLegA[7][2]=0.190915025252560;
  tLegA[7][3]=0.500000000000000;  wLegA[7][3]=0.208979591836735;
  tLegA[7][4]=0.702922575688699;  wLegA[7][4]=0.190915025252559;
  tLegA[7][5]=0.870765592799697;  wLegA[7][5]=0.139852695744638;
  tLegA[7][6]=0.974553956171379;  wLegA[7][6]=0.064742483084435;

  /* order=8 */
  tLegA[8][0]=0.019855071751232;  wLegA[8][0]=0.050614268145188;
  tLegA[8][1]=0.101666761293187;  wLegA[8][1]=0.111190517226687;
  tLegA[8][2]=0.237233795041836;  wLegA[8][2]=0.156853322938944;
  tLegA[8][3]=0.408282678752175;  wLegA[8][3]=0.181341891689181;
  tLegA[8][4]=0.591717321247825;  wLegA[8][4]=0.181341891689181;
  tLegA[8][5]=0.762766204958164;  wLegA[8][5]=0.156853322938944;
  tLegA[8][6]=0.898333238706813;  wLegA[8][6]=0.111190517226687;
  tLegA[8][7]=0.980144928248768;  wLegA[8][7]=0.050614268145188;

  /* order=9 */
  tLegA[9][0]=0.015919880246187;  wLegA[9][0]=0.040637194180787;
  tLegA[9][1]=0.081984446336682;  wLegA[9][1]=0.090324080347429;
  tLegA[9][2]=0.193314283649705;  wLegA[9][2]=0.130305348201467;
  tLegA[9][3]=0.337873288298095;  wLegA[9][3]=0.156173538520002;
  tLegA[9][4]=0.500000000000000;  wLegA[9][4]=0.165119677500630;
  tLegA[9][5]=0.662126711701904;  wLegA[9][5]=0.156173538520002;
  tLegA[9][6]=0.806685716350295;  wLegA[9][6]=0.130305348201468;
  tLegA[9][7]=0.918015553663318;  wLegA[9][7]=0.090324080347429;
  tLegA[9][8]=0.984080119753813;  wLegA[9][8]=0.040637194180787;

 /* order=10 */
  tLegA[10][0]=0.013046735741414;  wLegA[10][0]=0.033335672154344;
  tLegA[10][1]=0.067468316655508;  wLegA[10][1]=0.074725674575290;
  tLegA[10][2]=0.160295215850488;  wLegA[10][2]=0.109543181257991;
  tLegA[10][3]=0.283302302935376;  wLegA[10][3]=0.134633359654998;
  tLegA[10][4]=0.425562830509184;  wLegA[10][4]=0.147762112357376;
  tLegA[10][5]=0.574437169490816;  wLegA[10][5]=0.147762112357376;
  tLegA[10][6]=0.716697697064624;  wLegA[10][6]=0.134633359654998;
  tLegA[10][7]=0.839704784149512;  wLegA[10][7]=0.109543181257991;
  tLegA[10][8]=0.932531683344492;  wLegA[10][8]=0.074725674575290;
  tLegA[10][9]=0.986953264258586;  wLegA[10][9]=0.033335672154344;
} /* initQuad */
#endif

/*
 * For two triangles p and q, find out whether they are the same
 * (return value = 3 ). If they are not the same, find out the number of
 * common vertices and assign vertices v0, v1, v2 for p and
 * w0, w1, w2 for q such that,
 * - if there is one common vertex, the common vertex is (v0,w0) (return
 * value=1)
 * - if there are two vertices, the common vertices
 * are in addition (v1,w2) (return value=2) or (v2,w1) (return value=-2).
 * Note that
 *  - the orientation of the triangle is preserved.
 *  - v0 = p->vtx[idxX[0]], v1 = p->vtx[idxX[1]], v2 = p->vtx[idxX[2]]
 *    w0 = p->vtx[idxY[0]], w1 = p->vtx[idxY[1]], w2 = p->vtx[idxY[2]]
 */
int nrCommonVtx( panel *p, panel *q, int *idxX, int *idxY ) {
  int n1, n2, k, flag, cnt=0;

  idxX[0] = idxY[0] = 0;
  idxX[1] = idxY[1] = 1;
  idxX[2] = idxY[2] = 2;

  if ( p == q ) return 3;

  /* find the first common vertex, if it exists */
  for (n1=0; n1<3; n1++ ) {
    flag=0;
    for (n2=0; n2<3; n2++ ) {
      if ( p->vtx[n1][0] == q->vtx[n2][0] &&
           p->vtx[n1][1] == q->vtx[n2][1] &&
           p->vtx[n1][2] == q->vtx[n2][2] ) {
        flag=1;
        cnt++;
        break;
      }
    }
    if ( flag ) break;
  }
  if ( flag==0 ) return 0;  /* no common vertex found */

  /* reorder vertices such that the first common one appears first,
   * or just assign the vertices. Orientation is preserved */

  idxX[0] = n1%3; idxX[1] = (n1+1)%3; idxX[2] = (n1+2)%3;
  idxY[0] = n2%3; idxY[1] = (n2+1)%3; idxY[2] = (n2+2)%3;

  /* find the second common vertex, note that v1==w1 and v2==w2
     are impossible since we assume the same orientation */
  if ( p->vtx[idxX[1]][0] == q->vtx[idxY[2]][0] &&
       p->vtx[idxX[1]][1] == q->vtx[idxY[2]][1] &&
       p->vtx[idxX[1]][2] == q->vtx[idxY[2]][2] ) cnt = 2;
  if ( p->vtx[idxX[2]][0] == q->vtx[idxY[1]][0] &&
       p->vtx[idxX[2]][1] == q->vtx[idxY[1]][1] &&
       p->vtx[idxX[2]][2] == q->vtx[idxY[1]][2] ) cnt = -2;

  return cnt;
} /* nrCommonVtx */


/*
 * calculate the potential due to panel pnlY at point xC
 * using the Duffy transform
 */
double potentialP0(panel *pnlY, double *xC) {
  double *a0=pnlY->a[0], *a2=pnlY->a[2], *v0=pnlY->vtx[0];
  double r[3], r1[3], *tLeg, *wLeg;
  int i, j, k, qOrder;
  double tmp, intgr;

  qOrder = maxQuadOrder-1;
  tLeg = tLegA[qOrder];
  wLeg = wLegA[qOrder];

  for ( k=0; k<3; k++ ){
    r1[k] = xC[k] - v0[k];
  }

  intgr = 0.0;
  for ( i=0; i<qOrder; i++ ) {
    for ( j=0; j<qOrder; j++ ) {
      for ( k=0; k<3; k++ ){
        r[k] = r1[k] - tLeg[i]*(a2[k] + tLeg[j]*a0[k]);
      }
      kernel( r, &tmp );
      intgr += tmp*tLeg[i]*wLeg[i]*wLeg[j];
    }
  }
  intgr *= 2.0*pnlY->area;
  return intgr;

} /* potentialP0 */



/*
 *  calculate the self interaction of a panel [v0,v1,v2]. Piecewise constants.
 *  The method is based on Stefan Sauter's transformation to six non-singular
 *  integrals, which are approximated using Gauss quadrature.
 *  Parameters:
 *     qOrder: quadrature order
 *     a2:     the edge v1-v0
 *     a0:     the edge v2-v1
 */
double *pnlThr0(int qOrder, double *a2, double *a0) {
  int i, ix, iy, jy, ky, k;
  double r1[3], r2[3], r3[3], r4[3], r5[3], r6[3];
  double omg1, eta1, fta1;
  double *tLeg, *wLeg;

  tLeg = tLegA[qOrder];
  wLeg = wLegA[qOrder];

  for ( i=0; i<nKerl; i++ ) intgr[i] = 0.0;
  for ( iy=0; iy<qOrder; iy++ ) {
    for ( jy=0; jy<qOrder; jy++ ) {
      omg1 = tLeg[iy]*(1.0-tLeg[jy]);
      for ( ky=0; ky<qOrder; ky++ ) {
        for ( i=0; i<nKerl; i++ ) intgrY[i] = 0.0;
        for ( ix=0; ix<qOrder; ix++ ) {
          eta1 = tLeg[ix];
          fta1 = 1.0 - eta1;
          for ( k=0; k<3; k++ ) {
            r1[k] = omg1*( eta1*a2[k] + a0[k] );
            r2[k] = omg1*( fta1*a0[k] - eta1*a2[k] );
            r3[k] = omg1*( a2[k] + eta1*a0[k] );
            r4[k] = -omg1*( eta1*a2[k] + a0[k] );
            r5[k] = omg1*( eta1*a2[k] - fta1*a0[k] );
            r6[k] = -omg1*( a2[k] + eta1*a0[k] );
          }
          kernel( r1, fcn1 );
          kernel( r2, fcn2 );
          kernel( r3, fcn3 );
          kernel( r4, fcn4 );
          kernel( r5, fcn5 );
          kernel( r6, fcn6 );
          //double dist = r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2];
          //double ipY = nrmY[0]*r1[0] + nrmY[1]*r1[1] + nrmY[2]*r1[2];
          //printf("%f %f\n", dist, ipY);

          for ( i=0; i<nKerl; i++ )
            intgrY[i] += (fcn1[i]+fcn3[i]+fcn4[i]
                       + fcn2[i]+fcn5[i]+fcn6[i])*wLeg[ix];
        }
        for ( i=0; i<nKerl; i++ )
          intgr[i] += intgrY[i]*omg1*SQR(tLeg[iy])*tLeg[jy]*wLeg[iy]*wLeg[jy]*wLeg[ky];
      }
    }
  }
  return intgr;
} /* pnlThr0 */




/*
 *  calculate the interaction of two panels which share exacty one common
 *  edge. Piecewise constants. The method is based on Stefan Sauter's
 *  transformation to six non-singular integrals, which are approximated
 *  using Gauss quadrature.
 *  Parameters:
 *     qOrder: quadrature order
 *     a:      the edge which is common to both panels
 *     b:      edge that points out of the second common vertex
 *     c:      edge that points into the second common vertex
 *          /\
 *         /  \b (points up)
 *        /    \
 *        --a--->
 *        \    /
 *         \  /c (points up)
 *          \/
 */
double *pnlTwo0(int qOrder, double *a, double *b, double *c) {
  int i, ix, jx, iy, jy, k;
  double r1[3], r2[3], r3[3], r4[3], r5[3], r6[3];
  double omg1, eta1, eta2, eta3, fta1, fta2, fta3;
  double *tLeg, *wLeg;
  double intgrd;

  tLeg = tLegA[qOrder];
  wLeg = wLegA[qOrder];

  for ( i=0; i<nKerl; i++ ) intgr[i] = 0.0;
  for ( ix=0; ix<qOrder; ix++ ) {
    for ( jx=0; jx<qOrder; jx++ ) {
      omg1 = tLeg[ix]*tLeg[jx];
      for ( i=0; i<nKerl; i++ ) intgrY[i] = 0.0;
      for ( iy=0; iy<qOrder; iy++ ) {
        eta1 = tLeg[iy];
        fta1 = 1.0 - eta1;
        for ( jy=0; jy<qOrder; jy++ ) {
          eta2 = tLeg[jy];
          eta3 = eta1*eta2;
          fta2 = 1.0 - eta2;
          fta3 = 1.0 - eta3;
          for ( k=0; k<3; k++ ) {
            r1[k] = omg1*( eta1*(a[k] + eta2*b[k]) + fta1*c[k] );
            r2[k] = omg1*( eta1*(b[k] + eta2*a[k]) + fta3*c[k] );
            r3[k] = omg1*( fta1*a[k] + b[k] + eta3*c[k] );
            r4[k] = omg1*( eta1*(eta2*b[k] - fta2*a[k]) + c[k] );
            r5[k] = omg1*( eta3*b[k] - fta3*a[k] + eta1*c[k] );
            r6[k] = omg1*( eta1*b[k] - fta1*a[k] + eta3*c[k] );
          }
          kernel( r1, fcn1 );
          kernel( r2, fcn2 );
          kernel( r3, fcn3 );
          kernel( r4, fcn4 );
          kernel( r5, fcn5 );
          kernel( r6, fcn6 );
          for ( i=0; i<nKerl; i++ ) {
            intgrd = (fcn1[i]+fcn3[i]+fcn4[i]+fcn2[i]+fcn5[i]+fcn6[i])*eta1;
            intgrY[i] += intgrd*wLeg[iy]*wLeg[jy];
          }
        }
      }
      for ( i=0; i<nKerl; i++ )
        intgr[i] += intgrY[i]*tLeg[jx]*SQR(omg1)*wLeg[jx]*wLeg[ix];
    }
  }
  return intgr;
} /* pnlTwo0 */



/*
 *  calculate the interaction of two panels which share exacty one common
 *  vertex. Piecewise constants, The method is based on Stefan Sauter's
 *  transformation to two non-singular integrals, which are approximated
 *  using Gauss quadrature.
 *  Parameters:
 *     qOrder: quadrature order
 *     a2:     the edge v1-v0 of the outer-integration panel
 *     a0:     the edge v2-v1 of the outer-integration panel
 *     b2:     the edge v1-v0 of the inner-integration panel
 *     b0:     the edge v2-v1 of the inner-integration panel
 *             (v0 is the common vertex)
 */
double *pnlOne0(int qOrder, double *a2, double *a0, double *b2, double *b0) {
  int i, iy, jy, ky, ix, k;
  double r1[3], r2[3];
  double omg, wij, eta1, eta2, eta3;
  double *tLeg, *wLeg;

  tLeg = tLegA[qOrder];
  wLeg = wLegA[qOrder];

  for ( i=0; i<nKerl; i++ ) intgr[i] = 0.0;
  for ( ix=0; ix<qOrder; ix++ ) {
    omg = tLeg[ix];
    for ( i=0; i<nKerl; i++ ) intgrY[i] = 0.0;
    for ( iy=0; iy<qOrder; iy++ ) {
      eta1 = tLeg[iy];
      for ( jy=0; jy<qOrder; jy++ ) {
        eta2 = tLeg[jy];
        wij = wLeg[iy]*wLeg[jy]*eta1;
        for ( ky=0; ky<qOrder; ky++ ) {
          eta3 = tLeg[ky];
          for ( k=0; k<3; k++ ) {
            r1[k] = omg*( a2[k] + eta2*a0[k] - eta1*(b2[k] + eta3*b0[k]) );
            r2[k] = omg*( eta1*(a2[k] + eta3*a0[k]) - (b2[k] + eta2*b0[k]) );
          }
          kernel( r1, fcn1 );
          kernel( r2, fcn2 );
          for ( i=0; i<nKerl; i++ )
            intgrY[i] += (fcn1[i]+fcn2[i])*wLeg[ky]*wij;
        }
      }
    }
    for ( i=0; i<nKerl; i++ )
      intgr[i] += intgrY[i]*SQR(omg)*omg*wLeg[ix];
  }
  return intgr;

}/* pnlOne0 */

/*
 * calculate the interaction of two panels which are well separated.
 * Piecewise constants. The triangles are Duffy transformed on squares,
 * where Gauss-Legendre quadrature is applied.
 *
 *  Parameters:
 *     qOrder: quadrature order
 *     ax2:    the edge v1-v0 of the outer-integration panel
 *     ax0:    the edge v2-v1 of the outer-integration panel
 *     ay1:    the edge v1-v0 of the inner-integration panel
 *     ay0:    the edge v2-v1 of the inner-integration panel
 *     r0:     the vector v0(outer) - v0(inner)
 */
double *pnlNil0(int qOrder, double *r0, double *ax2, double *ax0,
                     double *ay1, double *ay2) {
  int i, ix, jx, iy, jy, k, nVtx;
  double r[3], r1[3];
  double *tLeg, *wLeg;

  tLeg = tLegA[qOrder];
  wLeg = wLegA[qOrder];

  for ( i=0; i<nKerl; i++ ) intgr[i] = 0.0;
  for ( ix=0; ix<qOrder; ix++ ) {
    for ( jx=0; jx<qOrder; jx++ ) {
      for ( k=0; k<3; k++ ){
        r1[k] = r0[k] + tLeg[ix]*(ax2[k] + tLeg[jx]*ax0[k]);
      }
      for ( i=0; i<nKerl; i++ ) intgrY[i]=0.0;
      for ( iy=0; iy<qOrder; iy++ ) {
        for ( jy=0; jy<qOrder; jy++ ) {
          for ( k=0; k<3; k++ ){
            r[k] = r1[k] - tLeg[iy]*(ay1[k] + tLeg[jy]*ay2[k]) ;
          }
          kernel( r, fcn1 );
          for ( i=0; i<nKerl; i++ )
            intgrY[i] += fcn1[i]*tLeg[iy]*wLeg[iy]*wLeg[jy];
        }
      }
      for ( i=0; i<nKerl; i++ )
        intgr[i] += intgrY[i]*tLeg[ix]*wLeg[ix]*wLeg[jx];
    }
  }
  /* printf("---end panel---\n"); */
  return intgr;
} /* pnlNil0 */



/*
 * Calculate the interaction of two flat panels
 * Piecewise constant Galerkin discretization.
 * Parameters:
 *    pnlX   panel of the outside integral
 *    pnlY   panel of the inside integral
 */
double *panelIA0(panel *pnlX, panel *pnlY ) {
  int i, k, nVtx, qOrder, idxX[3], idxY[3];
  double a0[3], b0[3], *vx0, *vy0, r0[3];
  double dist2, diam2;
  double *intgr;
  /* normals are only used for double, adjoint and hypersing */
  nrmX = pnlX->normal;
  nrmY = pnlY->normal;
  //printf("%f %f %f\n",nrmY[0],nrmY[1],nrmY[2]);
  //printf("%f %f %f\n",pnlY->x[0],pnlY->x[1],pnlY->x[2]);
  for ( i=0; i<nKerl; i++ ) intgrP[i] = 0.0;

  nVtx = nrCommonVtx( pnlX, pnlY, idxX, idxY );
  //printf("%d\n",nVtx);

  if ( nVtx==0 ) {
//    qOrder = MAX(maxQuadOrder-2,2);
    qOrder = maxQuadOrder;
    vx0 = pnlX->vtx[0];
    vy0 = pnlY->vtx[0];
    for ( k=0; k<3; k++ ){
      r0[k] = vx0[k] - vy0[k];
    }
    intgr = pnlNil0(qOrder, r0, pnlX->a[2], pnlX->a[0], pnlY->a[2], pnlY->a[0]);
    //printf("%f\n",intgr[0]);
  }
  else if ( nVtx==1 ) {
//    qOrder = MAX(maxQuadOrder-1,2);
    qOrder = maxQuadOrder;
    intgr = pnlOne0(qOrder, pnlX->a[idxX[2]], pnlX->a[idxX[0]],
                         pnlY->a[idxY[2]], pnlY->a[idxY[0]] );
  }
  else if ( nVtx==2 ) {
    qOrder = maxQuadOrder;

    intgr = pnlTwo0(qOrder, pnlX->a[idxX[2]], pnlX->a[idxX[0]],
                          pnlY->a[idxY[0]] );
  }
  else if ( nVtx==-2 ) {
    qOrder = maxQuadOrder;
    for ( k=0; k<3; k++ ) {
      b0[k] = -pnlY->a[idxY[0]][k];
      a0[k] = -pnlX->a[idxX[0]][k];
    }
    intgr = pnlTwo0(qOrder, pnlX->a[idxX[1]], b0, a0 );
    /*    intgr = pnlTwo0(qOrder, pnlY->a[idxY[2]], pnlY->a[idxY[0]],
          pnlX->a[idxX[0]] );*/
  }
  else if ( nVtx==3 ) {
    qOrder = maxQuadOrder;
    intgr = pnlThr0(qOrder, pnlX->a[2], pnlX->a[0]);
    //for ( i=0; i<nKerl; i++ ) printf("%f ", intgr[i]);
    //printf("\n");
  }

  //printf("%.12e %.12e %.12e %.12e\n",4.0*intgr[1],4.0*intgr[0],4.0*intgr[3],4.0*intgr[2]);

  for ( i=0; i<nKerl; i++ ) {
    intgrP[i] = 4.0*pnlX->area*pnlY->area*intgr[i];
    //intgrP[i] = 4.0*pnlY->area*intgr[i];
  }

  return intgrP;
} /* panelIA0 */


/*
 * Calculate the interaction of two flat panels for preconditioning
 * Piecewise constant Galerkin discretization.
 * Parameters:
 *    pnlX   panel of the outside integral
 *    pnlY   panel of the inside integral
 */
double *panelIA1(panel *pnlX, panel *pnlY ) {
  int i, k, nVtx, qOrder, idxX[3], idxY[3];
  double a0[3], b0[3], *vx0, *vy0, r0[3];
  double dist2, diam2;
  double *intgr;
  /* normals are only used for double, adjoint and hypersing */
  nrmX = pnlX->normal;
  nrmY = pnlY->normal;
  //printf("%f %f %f\n",nrmY[0],nrmY[1],nrmY[2]);
  //printf("%f %f %f\n",pnlY->x[0],pnlY->x[1],pnlY->x[2]);
  for ( i=0; i<nKerl; i++ ) intgrP[i] = 0.0;

  nVtx = nrCommonVtx( pnlX, pnlY, idxX, idxY );
  //printf("%d\n",nVtx);

  if ( nVtx==0 ) {
//    qOrder = MAX(maxQuadOrder-2,2);
    qOrder = maxQuadOrder;
    vx0 = pnlX->vtx[0];
    vy0 = pnlY->vtx[0];
    for ( k=0; k<3; k++ ){
      r0[k] = vx0[k] - vy0[k];
    }
    intgr = pnlNil0(qOrder, r0, pnlX->a[2], pnlX->a[0], pnlY->a[2], pnlY->a[0]);
    //printf("%f\n",intgr[0]);
  }
  else if ( nVtx==1 ) {
//    qOrder = MAX(maxQuadOrder-1,2);
    qOrder = maxQuadOrder;
    intgr = pnlOne0(qOrder, pnlX->a[idxX[2]], pnlX->a[idxX[0]],
                         pnlY->a[idxY[2]], pnlY->a[idxY[0]] );
  }
  else if ( nVtx==2 ) {
    qOrder = maxQuadOrder;

    intgr = pnlTwo0(qOrder, pnlX->a[idxX[2]], pnlX->a[idxX[0]],
                          pnlY->a[idxY[0]] );
  }
  else if ( nVtx==-2 ) {
    qOrder = maxQuadOrder;
    for ( k=0; k<3; k++ ) {
      b0[k] = -pnlY->a[idxY[0]][k];
      a0[k] = -pnlX->a[idxX[0]][k];
    }
    intgr = pnlTwo0(qOrder, pnlX->a[idxX[1]], b0, a0 );
    /*    intgr = pnlTwo0(qOrder, pnlY->a[idxY[2]], pnlY->a[idxY[0]],
          pnlX->a[idxX[0]] );*/
  }
  else if ( nVtx==3 ) {
    qOrder = maxQuadOrder;
    intgr = pnlThr0(qOrder, pnlX->a[2], pnlX->a[0]);
    //for ( i=0; i<nKerl; i++ ) printf("%f ", intgr[i]);
    //printf("\n");
  }

  //printf("%.12e %.12e %.12e %.12e\n",4.0*intgr[1],4.0*intgr[0],4.0*intgr[3],4.0*intgr[2]);

  for ( i=0; i<nKerl; i++ ) {
    intgrP[i] = 4.0*pnlX->area*pnlY->area*intgr[i];
    //intgrP[i] = 4.0*pnlY->area*intgr[i];
  }

  return intgrP;
} /* panelIA0 */



/*
 * Calculate the interaction of one flat panels and charges
 * Piecewise constant Galerkin discretization.
 * Parameters:
 *    pnlX   panel of the outside integral
 *    chrY   position of the charge
 */
double *panelRHS( int qOrder, panel *pnlX, double *chrY ) {
  int ix, jx, k, nVtx;
  double fcn[2];
  double r0[3], r[3], *ax2, *ax0;
  double *tLeg, *wLeg;

  tLeg = tLegA[qOrder];
  wLeg = wLegA[qOrder];

  ax2 = pnlX->a[2];
  ax0 = pnlX->a[0];

  nrmX = pnlX->normal;

  for ( k=0; k<3; k++ ) {
    r0[k] = pnlX->vtx[0][k]-chrY[k];
  }

  intgrRHS[0] = 0.0; intgrRHS[1] = 0.0;
  for ( ix=0; ix<qOrder; ix++ ) {
    for ( jx=0; jx<qOrder; jx++ ) {
      for ( k=0; k<3; k++ )
        r[k] = r0[k] + tLeg[ix]*(ax2[k] + tLeg[jx]*ax0[k]);
      kernelRHS( r, fcn );
      intgrRHS[0] += fcn[0]*tLeg[ix]*wLeg[ix]*wLeg[jx];
      intgrRHS[1] += fcn[1]*tLeg[ix]*wLeg[ix]*wLeg[jx];
    }
  }
  intgrRHS[0] *= 2.0*pnlX->area;
  intgrRHS[1] *= 2.0*pnlX->area;

  return intgrRHS;
} /* singlepanel */

double *panelPotential( int qOrder, double *chrX, panel *pnlY ) {
  int ix, jx, k, nVtx;
  double fcn[2];
  double r0[3], r[3], *ax2, *ax0;
  double *tLeg, *wLeg;

  tLeg = tLegA[qOrder];
  wLeg = wLegA[qOrder];

  ax2 = pnlY->a[2];
  ax0 = pnlY->a[0];

  nrmY = pnlY->normal;

  for ( k=0; k<3; k++ )
    r0[k] = chrX[k]-pnlY->vtx[0][k];

  intgrPtl[0] = 0.0; intgrPtl[1] = 0.0;
  for ( ix=0; ix<qOrder; ix++ ) {
    for ( jx=0; jx<qOrder; jx++ ) {
      for ( k=0; k<3; k++ )
        r[k] = r0[k] - tLeg[ix]*(ax2[k] + tLeg[jx]*ax0[k]);
      kernelPtl( r, fcn );
      intgrPtl[0] += fcn[0]*tLeg[ix]*wLeg[ix]*wLeg[jx];
      intgrPtl[1] += fcn[1]*tLeg[ix]*wLeg[ix]*wLeg[jx];
    }
  }
  intgrPtl[0] *= 2.0*pnlY->area;
  intgrPtl[1] *= 2.0*pnlY->area;

  return intgrPtl;
}
