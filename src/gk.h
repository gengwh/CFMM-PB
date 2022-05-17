/*
 * gk.h
 * copyright   Johannes Tausch
 */
struct panel {          /* panel */
  double vtx[3][3];     /* vertices */
  double nrm[3][3];     /* normal vector on vertices */
  double a[3][3];       /* sides */
  double x[3];          /* centroid */
  double normal[3];     /* normal, determined by right-hand rule */
  double area;          /* area  */
  //double area2;         /* 2*sqrt(area) */
  int shape;            /* 3=triangle, everything else generates an error */
  int nSurf;            /* number of face on surface (ith face) */
  int idx;              /* index in vector */
  struct cube *cube;    /* cube in which panel is located */
  struct panel *next;   /* next panel in input list. */
  struct panel *nextC;  /* next panel in linked list (Contiguous wn cubes) */
};
typedef struct panel panel;


struct edge {            /* edge */
  double v0[3];          /* startpoint */
  double v1[3];          /* endpoint */
};
typedef struct edge edge;


struct cube {               /* cube, actually a cluster of panels */
  int level;                /* 0 => root */
  int i, j, k;              /* cube is cubes[level][j][k][l] */
  double x[3];              /* Chebychev center */
  int nPnls;                /* number of panels in cube */
  panel *pnls;              /* linked list of panels with center in cube */
  int nNbrs;                /* Nbr of non-empty neighbors */
  int n2Nbrs;               /* Nbr of non-empty 2nd neighbors*/
  struct cube **nbrs;       /* Nbrs and 2nd nbrs with Panels */
  int nKids;                /* Number of kids */
  int nthKid;               /* nth kid of its parent */
  struct cube *kids[8];     /* Array of kids ptrs. */
  struct cube *parent;      /* parent cube */
  double *mom_pot;          /* moments for potential */
  double *mom_dpdn;         /* moments for potential derivative */
  double *lec_k1;           /* local expansion coefficients of kernel1 */
  double *lec_k2;           /* local expansion coefficients of kernel2 */
  double *lec_k3;           /* local expansion coefficients of kernel3 */
  double *lec_k4;           /* local expansion coefficients of kernel4 */
  double eBoxLo[3];         /* lower corner of the enclosing box */
  double eBoxUp[3];         /* upper corner of the enclosing box */
  double eRad;              /* half-diameter of enclosing box */
  struct cube *next;        /* Ptr to next nonempty cube with panels */
};
typedef struct cube cube;



struct ssystem {
  int depth;                /* # of levels of cubes. */
  int height;               /* highest level to throw out terms */
  int maxOrder;             /* maximal order */
  int nPnls;                /* nr of panels */
  int nChar;                /* nr of charges */
  int nVtxs;                /* nr of vertices */
//  int nSurf;                /* nr of surfaces */
  int nKerl;                /* nr of kernels */
  int layer;                /* single layer=0; double layer=1; adjoint=2 */
  int *ordMom;              /* order of moments to compute per level */
  int *ordM2L;              /* order used in M2L's */
  int *nMom;                /* number of moments (as a function of order) */
  int maxSngs;              /* max number of panels in a finest level cube */
  int max1Nbrs;             /* max number of frist neighbors that a cube can have */
  int maxQuadOrder;         /* maximal order for panel interactions */
  int maxlevCudes;          /* maximal number of finest lever cubes */
  int maxlevnPnls;          /* maximal number of nPnls the finest lever cubes have */
  int mesh_flag;            /* 1: MSMS 2: NanoShaper */
  double maxSepRatio;       /* maximal separation ratio for cubes to be neighbors */
  double *pos, *chr;        /* charge position and charges for rhs */
  panel *pnlLst;            /* linked list of panels (Contiguous wn cubes) */
  panel *pnlOLst;           /* linked list of original order panels */
  cube **cubeList;          /* heads of lists of cubes for each level */
};
typedef struct ssystem ssystem;
