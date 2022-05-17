#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"

extern double **Q2PK1, **Q2PK2, **Q2PK3, **Q2PK4;
extern void (*kernel)();
double *panelIA0(panel *pnlX, panel *pnlY );
void kernelKER4( double *x, double *y);

/* lapack: LU solver */
void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
void dgetrs_(char* C, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
void dgesdd_(char* C, int* N, int* M, double* A, int* LDA, double* s, double* u, int * LDU,
  double* VT, int* LDVT, double* work, int* lwork, int* iwork, int* info);

/* those variables are allocated at setupNearfield0 at fmm.c */
double *matrixA, *rhs;
int *ipiv;

extern double epsilon;
extern ssystem *sys;

int nlevel;

void setupPreconditioning(ssystem *sys) {

  int i, maxnPnls=0, idx, ttlcube=0;
  cube *cb;

  nlevel=sys->depth-1;
  //nlevel=sys->depth;
  kernel = kernelKER4;

  for ( idx=0, cb=sys->cubeList[nlevel]; cb!=NULL; cb=cb->next,idx++ ) {
    maxnPnls = cb->nPnls > maxnPnls ? cb->nPnls : maxnPnls;
    ttlcube += cb->nPnls;
  }

  maxnPnls *= 2;
  //CALLOC(matrixA, maxnPnls, double*, OFF, ASOLVER);
  //for ( i=0; i<maxnPnls; i++ ) CALLOC(matrixA[i], maxnPnls, double, OFF, ASOLVER);
  CALLOC(matrixA, maxnPnls*maxnPnls, double);
  CALLOC(ipiv, maxnPnls, int);
  CALLOC(rhs, maxnPnls, double);
  printf("Maximum number of elements in finest cluster: %d\n", maxnPnls);
  printf("----------------------------\n");
}


/*
 * Preconditioner by using the direct summation matrix
*/
int PtVfmm(double *pot, double *sgm) {
//void preconditioningFMM(ssystem *sys) {
  int i, j, idx, Msize, HMsize, inc;
  int nPnls=sys->nPnls;
  double scale1, scale2, *KER;
  cube *cb;
  panel *pnlX, *pnlY;

  scale1 = (1.0+epsilon)/2.0;
  scale2 = (1.0+1.0/epsilon)/2.0;

  for ( idx=0, cb=sys->cubeList[nlevel]; cb != NULL; cb=cb->next ) {
    Msize = 2*cb->nPnls;
    HMsize = cb->nPnls;

    for ( i=0, pnlY=cb->pnls; i<HMsize; i++, pnlY=pnlY->nextC ) {
      for ( j=0, pnlX=cb->pnls; j<HMsize; j++, pnlX=pnlX->nextC ) {
        KER = panelIA0(pnlX, pnlY);
        matrixA[i*Msize+j]                 = -KER[1];
        matrixA[i*Msize+j+HMsize]          = -KER[0];
        matrixA[(i+HMsize)*Msize+j]        = -KER[3];
        matrixA[(i+HMsize)*Msize+j+HMsize] = -KER[2];
      }
      matrixA[i*Msize+i] += scale1*pnlY->area;
      matrixA[(i+HMsize)*Msize+i+HMsize] += scale2*pnlY->area;

      rhs[i] = sgm[cb->pnls->idx+i];
      rhs[i+HMsize] = sgm[nPnls+cb->pnls->idx+i];
    }

    dgetrf_( &Msize, &Msize, matrixA, &Msize, ipiv, &inc );
    dgetrs_( &nChr, &Msize, &oneI, matrixA, &Msize, ipiv, rhs, &Msize, &inc );

    for ( i=0; i<HMsize; i++ ) {
      pot[cb->pnls->idx+i] = rhs[i];
      pot[nPnls+cb->pnls->idx+i] = rhs[i+HMsize];
    }
    idx += cb->nNbrs;

    for ( i=0; i<Msize; i++ ) {
      for ( j=0; j<Msize; j++ ) matrixA[i*Msize+j] = 0.0;
      rhs[i] = 0.0;
    }
  }

  //free(rhs);
  //free(ipiv);
  //for ( i=0; i<maxnPnls; i++ ) free(matrixA[i]);
  //free(matrixA);
  //exit(0);
}
