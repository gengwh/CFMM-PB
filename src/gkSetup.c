/*
 *  gkSetup.c
 *    routines that
 *    place panels in cubes, define the cube hierarchy, cube relations
 *    etc
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"

extern int doSphere;
extern void **sortBuf1, **sortBuf2, **LList;
#define EPS 1e-12

/* globals in this file */
double *sideLen, eta0;
int maxLev, maxPnlVtx=0;
cube **firstCbLst, **lastCbLst, *firstCube=NULL, *lastCube=NULL;

double minx, miny, minz;         /* corner of level-0 cube */
double eps;                      /* tolerance for equality tests */

void initCalcMoments(ssystem *sys);
void initQuad(ssystem *sys);
void initExpan(ssystem *sys);

/*
 * determine which expansion orders to use as a function of the level
 * order: orders for M2L's
 *      <0: variable order with -order in the finest level
 *      >0: fixed order with order in all levels
 * orderMom: orders for M2M's and L2L's (only relevant with variable order M2L's)
 *      <0: use exact moments (i.e., M2M's and L2L's are highest order
 *          in all levels).
 *     >=0: variable order M2M's and L2L's, the order in the finest level
 *          is given by order + orderMom.
 */

void getOrders( ssystem *sys, int order, int orderMom ) {
  int depth=sys->depth, height=sys->height;
  int lev, k, cnt;

  CALLOC(sys->ordMom, depth+1, int);
  CALLOC(sys->ordM2L, depth+1, int);

  if ( order > 0 ) { /* fixed order */
    for ( lev=0; lev<=depth; lev++ ) {
      sys->ordM2L[lev] = order;
      sys->ordMom[lev] = order;
    }
    sys->maxOrder = order;
  }
  else {            /* variable order */
    sys->ordM2L[depth] = -order;
    for ( lev=depth-1; lev>=height; lev-- ) {
      sys->ordM2L[lev] = sys->ordM2L[lev+1] + 1;
    }
    sys->maxOrder = sys->ordM2L[height];
    if ( orderMom >= 0 ) {
      for ( lev=depth; lev>=height; lev-- ) {
        sys->ordMom[lev] = MIN(sys->ordM2L[lev] + orderMom, sys->maxOrder);
      }
    }
    else {
      for ( lev=depth; lev>=height; lev-- ) {
        sys->ordMom[lev] = sys->maxOrder;
      }
    }
  } /* variable order */


  /* determine the number of moments as a function of the order */
  CALLOC(sys->nMom, sys->maxOrder+3, int);
  for ( cnt=0; cnt<=sys->maxOrder+2; cnt++ ) {
    sys->nMom[cnt] = ((cnt+1)*(cnt+2)*(cnt+3))/6;
  }

} /* getOrders */

/*
 * loops over all panels:
 *   - finds the 'smallest' corner of the coarsest-level cube
 *   - calculates the sidelengths of the cubes in each level
 *   - finds the number of surfaces
 *   - returns number of panels
*/
int loopPanels(ssystem *sys, panel *panels) {
  double maxx, maxy, maxz, length;
  panel *pnl;
  cube *cb;
  int nPnls, k, lev;

  /* Figure out the length of lev 0 cube */
  pnl = panels;
  minx = maxx = pnl->vtx[0][0];
  miny = maxy = pnl->vtx[0][1];
  minz = maxz = pnl->vtx[0][2];

  for ( nPnls=0, pnl=panels; pnl != NULL; pnl = pnl->next) {
    for (k=0; k<pnl->shape; k++ ) {
      maxx = MAX(pnl->vtx[k][0], maxx);
      minx = MIN(pnl->vtx[k][0], minx);
      maxy = MAX(pnl->vtx[k][1], maxy);
      miny = MIN(pnl->vtx[k][1], miny);
      maxz = MAX(pnl->vtx[k][2], maxz);
      minz = MIN(pnl->vtx[k][2], minz);
    }
    for (k=0; k<3; k++) {
      pnl->x[k] = ONETHIRD*(pnl->vtx[0][k] + pnl->vtx[1][k] + pnl->vtx[2][k]);
    }
    nPnls++;
  }

  CALLOC(sideLen, maxLev+1, double);

  length = MAX((maxx - minx), (maxy - miny));
  length = MAX((maxz - minz), length);

  for ( lev=0; lev<=maxLev; lev++ ) {
    sideLen[lev] = length;
    length *= 0.5;
  }

  return nPnls;
} /* loopPanels */


/*
 * refPt (reference point) is a point located in cube *cb
 * The routine determines the kid of cb in which refPt is located
 * If the kid already exists the routine calls itself with kid
 * instead of cube. If the kid does not exist it allocates the kid
 * and fills in the basic members, then recurses.
 * The returned value is the pointer to the finest-level cube in which
 * refPt is located.
 * Also, generates a linked list of finest level cubes.
 * Globals referenced:   sideLen,  maxLev, firstCube, lastCube;
 */
cube* goDown( double* refPt, cube* cb ) {
  int ix, iy, iz, idx;
  cube* kid;
  /* find index of kid */
  ix = iy = iz = 0;
  if ( refPt[0] > minx + (cb->i+0.5)*sideLen[cb->level] ) ix = 1;
  if ( refPt[1] > miny + (cb->j+0.5)*sideLen[cb->level] ) iy = 1;
  if ( refPt[2] > minz + (cb->k+0.5)*sideLen[cb->level] ) iz = 1;
  idx = ix + 2*iy + 4*iz;

  if ( cb->kids[idx] == NULL ) {
    CALLOC(cb->kids[idx], 1, cube);
    kid = cb->kids[idx];
    kid->level = cb->level+1;
    kid->i = 2*cb->i + ix;
    kid->j = 2*cb->j + iy;
    kid->k = 2*cb->k + iz;
    kid->parent = cb;
    if ( kid->level == maxLev ) {
      kid->x[0] = minx + (kid->i+0.5)*sideLen[maxLev];
      kid->x[1] = miny + (kid->j+0.5)*sideLen[maxLev];
      kid->x[2] = minz + (kid->k+0.5)*sideLen[maxLev];
      if ( firstCube == NULL ) {
        firstCube = kid;
      }
      else {
        lastCube->next = kid;
      }
      lastCube = kid;
    }
  }
  else {
    kid = cb->kids[idx];
  }

  if ( kid->level == maxLev ) {
    return kid;
  }
  else {
    return goDown( refPt, kid );
  }
} /* goDown */



/*
 * This does three things:
 * 1.) Generate linked lists of cubes in the same level. The links are contigous
 *     within higher level cubes.
 * 2.) generate the links between cubes for the "nextC" linked lists of panels
 * 3.) eliminate empty kids and count non-empty kids
 */
void linkcubes(cube *cb) {
  int i1, i2;
  int lev = cb->level;
  panel *pnl;

  if ( lastCbLst[lev] == NULL ) {
    firstCbLst[lev] = lastCbLst[lev] = cb;
  }
  else {
    if ( lev == maxLev ) {
      /* connect panel lists */
      pnl = lastCbLst[lev]->pnls;
      while ( pnl->nextC!=NULL ) pnl=pnl->nextC;
      pnl->nextC = cb->pnls;

      cb->next = NULL; /* necessary because we've already generated links before */
    }
    /* append current cube to cube list */
    lastCbLst[lev]->next = cb;
    lastCbLst[lev] = cb;
  }

  if ( lev < maxLev ) {
    /* count number of kids, remove empty kids, recurse  */
    for ( i1=0; i1<8; i1++ ) {
      if ( cb->kids[i1] != NULL ) cb->nKids++;
    }
    for ( i1=i2=0; i1<8; i1++ ) {
      if ( cb->kids[i1] != NULL ) cb->kids[i2++] = cb->kids[i1];
    }
    for ( i1=cb->nKids; i1<8; i1++ ) {
      cb->kids[i1] = NULL ;
    }

    for ( i2=0; i2<cb->nKids; i2++ ) {
      linkcubes(cb->kids[i2]);
      cb->nPnls += cb->kids[i2]->nPnls;
    }
    /* let parents link their children's pnl */
    if ( cb->pnls == NULL )
      cb->pnls = cb->kids[0]->pnls;
  }
} /* linkcubes */


/*
 * adds up the number of panels in coarser levels
 * creates links to panels in coarser levels
 * by recursing through the cube-tree
 */
void countRecurse(cube *cb) {
  int i;
  cube *kid;

  if ( cb->level == maxLev ) return;
  for ( i=0; i<cb->nKids; i++ ) {
    kid=cb->kids[i];
    countRecurse(kid);
    cb->nPnls += kid->nPnls;
    if ( cb->pnls == NULL ) cb->pnls = cb->kids[i]->pnls;
  }
} /* countRecurse */



/*
 * find the center and the radius (=half-length of enclosing box's main
 * diagonal) for all cubes. Furthermore, find the enclosing boxes for cubes
 * in the  coarser levels.
 */
void getEnclBoxs(ssystem *sys){
  double *up, *lo, *x;
  cube *cb, *kid;
  int i, j, k, lev;
  panel *pnl;

  /* finest level cubes */
  for ( cb=firstCbLst[maxLev]; cb!=NULL; cb=cb->next ) {
    up = cb->eBoxUp;
    lo = cb->eBoxLo;
    pnl = cb->pnls;
    x = pnl->vtx[0];
    for ( k=0; k<3; k++ ) {
      up[k] = x[k];
      lo[k] = x[k];
    }

    for ( i=0; i<cb->nPnls; i++, pnl=pnl->nextC ) {
      for ( j=0; j<pnl->shape; j++ ) {
        x = pnl->vtx[j];
        for ( k=0; k<3; k++ ) {
          //if ( x[k] > up[k] ) up[k] = x[k];
          up[k] = x[k] > up[k] ? x[k] : up[k];
          //if ( x[k] < lo[k] ) lo[k] = x[k];
          lo[k] = x[k] < lo[k] ? x[k] : lo[k];
        }
      }
    }
    cb->x[0] = 0.5*(lo[0] + up[0]);
    cb->x[1] = 0.5*(lo[1] + up[1]);
    cb->x[2] = 0.5*(lo[2] + up[2]);
    cb->eRad = 0.5*sqrt(SQR(up[0]-lo[0]) + SQR(up[1]-lo[1]) + SQR(up[2]-lo[2]));
  }

  for ( lev=maxLev-1; lev>=sys->height; lev-- ) {
    for ( cb=firstCbLst[lev]; cb!=NULL; cb=cb->next ) {
      up = cb->eBoxUp;
      lo = cb->eBoxLo;
      kid = cb->kids[0];
      for ( k=0; k<3; k++ ) {
        up[k] = kid->eBoxUp[k];
        lo[k] = kid->eBoxLo[k];
      }
      for ( i=1; i<cb->nKids; i++ ) {
        kid = cb->kids[i];
        for ( k=0; k<3; k++ ) {
          //if ( kid->eBoxUp[k] > up[k] ) up[k] = kid->eBoxUp[k];
          up[k] = kid->eBoxUp[k] > up[k] ? kid->eBoxUp[k] : up[k];
          //if ( kid->eBoxLo[k] < lo[k] ) lo[k] = kid->eBoxLo[k];
          lo[k] = kid->eBoxLo[k] < lo[k] ? kid->eBoxLo[k] : lo[k];
        }
      }
      cb->x[0] = 0.5*(lo[0] + up[0]);
      cb->x[1] = 0.5*(lo[1] + up[1]);
      cb->x[2] = 0.5*(lo[2] + up[2]);
      cb->eRad = 0.5*sqrt(SQR(up[0]-lo[0]) + SQR(up[1]-lo[1]) + SQR(up[2]-lo[2]));
    }
  }
} /* getEnclBoxs */




/*
 * returns 1 if c1 and c2 are first neighbors, but not equal
 *         2 if c1 == c2
 *         0 otherwise
 */
int near(cube *c1, cube *c2) {
  double dist, eta;

  if ( c1 == c2 ) return 2;
  dist = SQR(c1->x[0]-c2->x[0]) + SQR(c1->x[1]-c2->x[1]) + SQR(c1->x[2]-c2->x[2]);
  dist = sqrt(dist);

  /* this is for FMM */
  if ( dist<1e-15 ) return 1;

  eta = (c1->eRad + c2->eRad)/dist;
  if ( eta < eta0 ) {
    return 0;
  }
  else {
    return 1;
  }
} /* near */


/*
 * Get the first- and second nearest neighbors. All second-nearest neighbors
 * which are not first-nearest neighbors form the interaction list.
 * First-nearest neighbors are determined by the routine near().
 * In the coarsest level, the second nearest neighbors are all
 * cubes at that level. In the finer levels, the neighbors are determined
 * from the first-nearest neighbors at the parent's level.
 * In the neighborlist, the subject cube comes first, followed by the first
 * nearest and then the second nearest neighbors.
 *
 */
void getNbrs(ssystem *sys) {
  cube *cb, *nbr, *parent, *nbrParent;
  cube **nbrList;
  int i1, i2, nK, nP, lev, nNbrs, n2Nbrs, flag;
  int max1Nbrs=0;
  /*
   * coarsest level
   */
  for ( n2Nbrs=0, cb=sys->cubeList[sys->height]; cb!=NULL; cb=cb->next ) n2Nbrs++;

  for ( cb=sys->cubeList[sys->height]; cb!=NULL; cb=cb->next ) {
    CALLOC(nbrList, n2Nbrs, cube*);
    cb->nbrs = nbrList;
    for ( nNbrs=0, nbr=sys->cubeList[sys->height]; nbr!=NULL; nbr=nbr->next ) {
      if ( near(cb,nbr) ) nNbrs++;
    }
    for ( i1=1,i2=nNbrs,nbr=sys->cubeList[sys->height]; nbr!=NULL; nbr=nbr->next ) {
      flag = near(cb,nbr);
      if ( flag==1 ) {
        nbrList[i1++] = nbr;
      }
      else if ( flag==2 ) {
        nbrList[0] = nbr;
      }
      else {
        nbrList[i2++] = nbr;
      }
    }
    cb->nNbrs  = i1;
    cb->n2Nbrs = i2;
    max1Nbrs = MAX(max1Nbrs, i1);
  }

  /*
   * finer levels
   */
  for ( lev=sys->height+1; lev<=sys->depth; lev++ ) {
    for ( cb=sys->cubeList[lev]; cb!=NULL; cb=cb->next ) {
      parent = cb->parent;
      for ( n2Nbrs=nNbrs=nP=0; nP<parent->nNbrs; nP++ ) {
        nbrParent = parent->nbrs[nP];
        for ( nK=0; nK<nbrParent->nKids; nK++ ) {
          nbr = nbrParent->kids[nK];
          if ( near(nbr,cb) ) {
            nNbrs++;
            if ( near(nbr, cb)==2 ) cb->nthKid = nNbrs;
          }
          n2Nbrs++;
        }
      }
      CALLOC(nbrList, n2Nbrs, cube*);
      cb->nbrs = nbrList;
      for ( nP=0, i1=1, i2=nNbrs; nP<parent->nNbrs; nP++ ) {
        nbrParent = parent->nbrs[nP];
        for ( nK=0; nK<nbrParent->nKids; nK++ ) {
          nbr = nbrParent->kids[nK];
          flag = near(cb,nbr);
          if ( flag==1 ) {
            nbrList[i1++] = nbr;
          }
          else if ( flag==2 ) {
            nbrList[0] = nbr;
          }
          else {
            nbrList[i2++] = nbr;
          }
        }
      }
      cb->nNbrs  = i1;
      cb->n2Nbrs = i2;
      max1Nbrs = MAX(max1Nbrs, i1);
    }
  }
  sys->max1Nbrs = max1Nbrs;
} /* getNbrs */




/*  main driver in this file
 *  sets up the partitioning of space
 *  for the meaning of the parameters order and orderMom see getOrders()
 */
void gkInit(ssystem *sys, panel *pnlList, int order, int orderMom ) {
  int nPnls, i, j, k, l, lev, cnt;
  double area, dist2, minDist2;
  panel *pnl, *pnl1;
  cube *cb, *topCube, **nbrs;
  time_t time1, time2;

  maxLev = sys->depth;
  eta0 = sys->maxSepRatio;

  nPnls = loopPanels(sys, pnlList);
  eps = sideLen[0]*EPS;

  CALLOC(topCube, 1, cube);
  topCube->level = 0;
  topCube->i = topCube->j = topCube->k = 0;
  firstCube = NULL;

  /* allocate cubes in all levels
   * insert panels into finest-level cubes
   */
  for ( pnl = pnlList; pnl != NULL; pnl = pnl->next) {
    cb = goDown(pnl->x, topCube);

    /* append the panel to the end of cube */
    if ( (cb->pnls) == NULL ) {
      cb->pnls = pnl;
    } else {
      pnl1 = cb->pnls;
      while ( pnl1->nextC != NULL ) pnl1 = pnl1->nextC;
      pnl1->nextC = pnl;
    }
    cb->nPnls++;
    sys->nPnls++;
  }

  /* crate linked lists of cubes */
  CALLOC(lastCbLst, maxLev+1, cube*);
  CALLOC(firstCbLst, maxLev+1, cube*);
  sys->cubeList = firstCbLst;
  linkcubes(topCube);
  sys->pnlLst = sys->cubeList[maxLev]->pnls;
  sys->maxlevCudes = sys->maxlevnPnls = 0;
  for ( cb=sys->cubeList[maxLev]; cb!=NULL; cb=cb->next ) {
    sys->maxlevnPnls = cb->nPnls > sys->maxlevnPnls ? cb->nPnls : sys->maxlevnPnls;
    sys->maxlevCudes++;
  }
  //printf("%d\n", firstCbLst[3]->nPnls);
  //exit(0);
  //printf("The finest level contains %d cubes of %d imagine cubes\n",
  //       sys->maxlevCudes, (int)pow(8, maxLev));
  //countRecurse(topCube);
  for ( i=0, pnl=sys->pnlLst; pnl!=NULL; pnl=pnl->nextC ) pnl->idx = i++;

  getEnclBoxs(sys);
  getNbrs(sys);

  getOrders(sys, order, orderMom);
  initQuad(sys);
  initCalcMoments(sys);
  initExpan(sys);

} /* gkInit */
