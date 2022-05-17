/*
 *  gkGlobal.h
 *
 *  Copyright by Johannes Tausch
 *  some stuff was inspired by K.Nabor's and J.White's fastcap package
 *
 */


/*
 * begin configuration parameters
 */
#define ON 1
#define OFF 0
#define NUMKERNEL ON        /* count number of kernel evaluations */

/*
 * end configuration parameters
 */
#include <time.h>

#define ONESIXTH 0.16666666666666666666666667
#define ONETHIRD 0.33333333333333333333333333
#define TWOTHIRD 0.66666666666666666666666667
#define ONENINTH 0.11111111111111111111111111
#define ONETWLF  0.08333333333333333333333333
#define ONESQTWO 0.70710678118654752440084436



#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

#ifndef MIN
#define MIN(A,B)  ( (A) > (B) ? (B) : (A) )
#endif

#define ABS(A) ( ( (A) > 0 ) ? (A) : (-(A)) )
#define SQR(x) ( (x)*(x) )
#define CUB(x) ( (x)*(x)*(x) )
#define DIST2(V,W) (SQR((V)[0]-(W)[0]) + SQR((V)[1]-(W)[1]) + SQR((V)[2]-(W)[2]))
#define INNR(V,W) ((V)[0]*(W)[0] + (V)[1]*(W)[1] + (V)[2]*(W)[2]);
#define SWITCH(x,y,tmp) tmp = y; y = x; x = tmp;

#if COMPLEX
#define FABS(A)  sqrt( __real__(A)*__real__(A) + __imag__(A)*__imag__(A) )
#define FABS2(A)  __real__(A)*__real__(A) + __imag__(A)*__imag__(A)
#else
#define FABS(A)  fabs( (A) )
#define FABS2(A) (A)*(A)
#endif

/* memory types */
#define APVE 1      /* panels, vertices, edges */
#define ACUBES 2    /* cube tree */
#define AQ2P 3      /* Q2P transformations */
#define AQ2M 4      /* Q2M/L2P transformations */
#define AM2L 5      /* M2L transformations */
#define ASOLVER 6   /* extra memory allocated in linear solver */
#define AMISC 7     /* misc. memory such as workspaces */


#define NOT !
#define  ABORT() {						                                           \
  (void)fflush(stdout);					                                         \
  (void)fprintf(stdout, "panic in file `%s' at line %d.\n",              \
	  __FILE__, __LINE__);				                                         \
  (void)fflush(stdout);					                                         \
  abort();							                                                 \
}

#define ASSERT(condition) if(NOT(condition)) ABORT()



#define CALLOC(PNTR, NUM, TYPE)                                          \
{                                                                        \
  if((NUM)*sizeof(TYPE)==0)                                              \
    (void)fprintf(stdout,                                                \
		  "zero element request in file `%s' at line %d\n",                  \
		  __FILE__, __LINE__);	                                             \
  else if(((PNTR)=(TYPE*)calloc(NUM, sizeof(TYPE)))==NULL) {             \
    (void)fprintf(stdout,                                                \
	    "\ngk: out of memory in file `%s' at line %d\n",                   \
	    __FILE__, __LINE__);                                               \
    (void)fflush(stdout);                                                \
    (void)fflush(stdout);                                                \
  }                                                                      \
}


#define TRUE 1
#define FALSE 0
#define QUADRILAT 4
#define TRIANGLE 3
#define PCWCONST 0
#define PCWLINEAR 1
#define SINGLE 0
#define DOUBLE 1
#define ADJOINT 2
#define HYPERSING 3
#define COMBINED 4

/* timers */
extern double setupBasTime, setupMatTime, setupPCTime, solveTimeNoPC, solveTimePC;
extern double setupQ2PTime, setupQ2MTime, setupM2LTime;
extern double gkInitTime, setupFMMTime, setupRHSTime, gmresTime;
extern double fmmQ2MTime, fmmM2MTime, fmmM2LTime, fmmL2LTime, fmmL2PTime, fmmNearTime;
/* counts of memory usage */
extern long memcount;
extern long memPVE, memCUBES, memQ2P, memQ2M, memM2L, memSOLVER, memMISC;
/* counts of kernel evaluations */
extern long numKernRealEval, numKernCplxEval;
/* misc constants */
extern double twoPi, twoPiI, fourPi, fourPiI;
extern double zero, one;
extern int oneI;
extern char hChr, nChr;
