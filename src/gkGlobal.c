#include <math.h>
#include <time.h>
#include "gkGlobal.h"
double gkInitTime, setupFMMTime, setupRHSTime, gmresTime;
double solveTimeNoPC, solveTimePC;
double setupQ2PTime, setupQ2MTime, setupM2LTime;
double fmmQ2MTime, fmmM2MTime, fmmM2LTime, fmmL2LTime, fmmL2PTime, fmmNearTime;

/* memory counters */
long memcount=0;    /* total memory */
long memPVE=0;      /* panels, vertices, edges */
long memCUBES=0;    /* cube tree */
long memQ2P=0;      /* Q2P transformations */
long memQ2M=0;      /* Q2M/L2P transformations */
long memM2L=0;      /* M2L transformations */
long memSOLVER=0;   /* allocated by solver */
long memMISC=0;     /* all other stuff */

/* matrix entry counters */
long numKernRealEval=0, numKernCplxEval=0;

/* misc constants */
double fourPi = 4*M_PI, fourPiI = 1.0/(4*M_PI);
double twoPi = 2*M_PI, twoPiI = 1.0/(2*M_PI);
double zero=0.0, one=1.0;
int oneI = 1;

char hChr='T';             /* transpose */
char nChr='N';             /* normal */
