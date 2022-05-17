/*
 * coulomb.c: main driver
 * This program computes the boundary integral PB equation with fmm method
 * usage:
 *   coulomb [options] panelfile [options]
 *
 * Based on Tausch's Code
 * Copyright Jiahui Chen
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gkGlobal.h"
#include "gk.h"

/* global variables */
int orderMom=0;
double kappa, epsilon, epsilon1=1.0, epsilon2=80.0;

ssystem *sys;

/* function pointers to kernel routines */
void (*kernel)(), (*kernelD)(), (*kernelDC)(), (*kernelDS)();
int (*MtV)(), (*PtV)();

/* routines used by the main routine */
panel *loadPanel(char *panelfile, char *density, int *numSing, ssystem *sys);
void gkInit(ssystem *sys, panel *pnlList, int order, int orderMom);
void setupFMM(ssystem *sys);
void applyFMM( ssystem *sys, double *alpha, double *sgm, double *beta, double *pot );
void setupPreconditioning(ssystem *sys);

double *panelRHS(int qOrder, panel *pnlX, double *chrY );

int MtVmain(double *alpha, double *sgm, double *beta, double *pot);
int PtVfmm(double *pot, double *sgm);

int gmres_(long int *n, double *b, double *x, long int *arnoldiSz, double *work,
           long int *ldw, double *h, long int *ldh, long int *iter, double *resid,
           int (*matvec)(), int (*psolve)(), long int *info);

void applyTreecode( ssystem *sys, double *sgm, double *pot );
/*
 *  setup right hand side (exterior Neumann problem)
 */
void setupRHS(ssystem *sys, double *sgm) {
  int i, j;
  int nPnls = sys->nPnls, nChar = sys->nChar, qOrder=sys->maxQuadOrder;
  double *intgr, fac;
  panel *pnl;

  fac = fourPiI/epsilon1;

  /* triangles order for Direct */
  for ( i=0, pnl=sys->pnlLst; pnl!=NULL; pnl=pnl->nextC, i++ ) {
    sgm[i] = 0.0; sgm[nPnls+i] = 0.0;
    for ( j=0; j<nChar; j++ ) {
      intgr=panelRHS(qOrder, pnl, &sys->pos[3*j]);
      sgm[i] += sys->chr[j]*intgr[0];
      sgm[i+nPnls] += sys->chr[j]*intgr[1];
    }
    sgm[i] *= fac;
    sgm[nPnls+i] *= fac;
  }

} /* setupRHS */



int main(int nargs, char *argv[]){
  char panelfile[80], density[80];
  int order=-1, image=0, refineLev=0, numSurfOne=1;
  int i, j, k, n, nPnls, nChar;
  long int numItr=100, arnoldiSz=30, ldw, ldh;
  panel *inputLst, *pnl;
  cube *cb;
  double tolpar=1.0e-4, para=332.0716;
  double *sgm, *pot, *GMRES_work, *GMRES_h, ptl;
  static long int info;

  clock_t start_t, end_t;

  CALLOC(sys, 1, ssystem);
  sys->height = 2;
  sys->maxSepRatio = 0.8;
  sys->maxQuadOrder = 1;
  sys->nKerl = 4;
  sys->depth = 5;
  sys->mesh_flag = 1;
  sprintf(density,"1");
  epsilon = epsilon2/epsilon1;
  double bulk_strength = 0.15;
  //kappa = sqrt(8.430325455*bulk_strength/epsilon2);
  kappa = 0.1257;

  /* parse the command line */
  panelfile[0] = 0;
  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
        case 'S': sys->maxSepRatio = atof( argv[i]+3 );
          break;
        case 'o': tolpar = atof( argv[i]+3 );
          break;
        case 'p':
          if ( argv[i][2] == '=' ) order = atoi( argv[i]+3 );
          if ( argv[i][2] == 'm' ) orderMom = atoi( argv[i]+4 );
          break;
        case 'q':
          sys->maxQuadOrder = atoi( argv[i]+3 );
          break;
        case 't': sys->depth = atoi( argv[i]+3 );
          break;
        case 'd': strcpy(density,argv[i]+3);
          break;
        case 'e':
          if ( argv[i][3] == '1' ) epsilon1 = atof( argv[i]+4 );
          if ( argv[i][3] == '2' ) epsilon2 = atof( argv[i]+4 );
          break;
        case 'k': kappa = atof( argv[i]+3 );
          break;
        case 'm': sys->mesh_flag = atoi( argv[i]+3 );
          break;
      }
    else {
      strcpy(panelfile,argv[i]);
    }

  if ( panelfile[0] == 0 ) {
    printf("\n Name of the panel file > ");
    if ( scanf("%s",panelfile) < 1 ) {
      printf("PDB name input failed\n");
      exit(0);
    }
  }
  if ( sys->depth < 0 ) {
    printf("Select tree depth > ");
    if ( scanf("%d", &sys->depth) < 1 ) {
      printf("PDB density input failed\n");
      exit(0);
    }
    if( sys->depth < 1 ) {
      printf("Bad tree depth: %d\n", sys->depth );
      exit(0);
    }
  }
  //printf("PDB id: %s, MSMS density: %s\n", panelfile, density);
  printf("----------------------------\n");
  printf("FMM variables: nLev=%d ord=%d SepRat=%lg qOrd=%d\n",
    sys->depth, order, sys->maxSepRatio, sys->maxQuadOrder );
  printf("GMRES variables: tol=%1.e arnoldiSz=%ld maxIt=%ld\n",
    tolpar, arnoldiSz, numItr);
  printf("kappa=%f, eps1=%.0f, eps2=%.0f\n", kappa, epsilon1, epsilon2);
  //printf("----------------------------\n");


  /*
   * get panels by msms from pqr
   * or use the panel on sphere test example
   */
  start_t = clock();
  inputLst = loadPanel(panelfile, density, &nPnls, sys);
  sys->pnlOLst = inputLst;

  gkInit(sys, inputLst, order, orderMom);

  CALLOC(sgm, 2*nPnls, double);
  CALLOC(pot, 2*nPnls, double);

  setupFMM(sys);
  setupPreconditioning(sys);

  setupRHS(sys, sgm);
  for ( i=0; i<2*sys->nPnls; i++ ) pot[i] = sgm[i];

  MtV = MtVmain;
  PtV = PtVfmm;
  ldw = 2*nPnls;
  ldh = arnoldiSz+1;

  CALLOC(GMRES_work, ldw*(arnoldiSz+4), double);
  CALLOC(GMRES_h, ldh*(arnoldiSz+2), double);

  gmres_(&ldw, pot, sgm, &arnoldiSz, GMRES_work, &ldw, GMRES_h, &ldh,
         &numItr, &tolpar, MtV, PtV, &info);

  applyTreecode( sys, sgm, &ptl );
  ptl *= twoPi*para;
  end_t = clock()-start_t;
  printf("ttl time: %f, gmres-its=%ld\n", (double)end_t / CLOCKS_PER_SEC, numItr);
  printf("solvation energy: %f\n", ptl);

}



/*
 * Matrix times Vector, subroutine of the iterative solver
 * the vector sgm and the result pot are ordered contiguously within cubes
 */
int MtVmain(double *alpha, double *sgm, double *beta, double *pot) {
  int i, lev, inc=1, nPnls = sys->nPnls;
  cube *cb;
  panel *pnl;
  double scale1, scale2, inv_beta;

  scale1 = (1.0+epsilon)/2.0*(*alpha);
  scale2 = (1.0+1.0/epsilon)/2.0*(*alpha);

  inv_beta = -(*beta);
  applyFMM(sys, alpha, sgm, &inv_beta, pot);
  for (  i=0, pnl=sys->pnlLst; pnl!=NULL; pnl=pnl->nextC, i++ ) {
    pot[i] = (scale1*pnl->area*sgm[i]-pot[i]);
    pot[i+nPnls] = scale2*pnl->area*sgm[i+nPnls]-pot[i+nPnls];
  }

} /* MtVmain */
