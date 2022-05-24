/*
 * input.c
 *   various input/output routines for gk package
 *
 *   copyright by Johannes Tausch
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gkGlobal.h"
#include "gk.h"
#define MAXCOND 50
#define DENSITY 10    /* for mkMultiSpheres() */

int nSurf = 0;

/*
 * calculate area and normal
 */
double triangle_area(double v[3][3]){
  int i;
  double a[3], b[3], c[3], aa, bb, cc, ss, t_area;
  for (i=0;i<3;i++){
    a[i] = v[0][i]-v[1][i];
    b[i] = v[0][i]-v[2][i];
    c[i] = v[1][i]-v[2][i];
  }
  aa = sqrt(SQR(a[0])+SQR(a[1])+SQR(a[2]));
  bb = sqrt(SQR(b[0])+SQR(b[1])+SQR(b[2]));
  cc = sqrt(SQR(c[0])+SQR(c[1])+SQR(c[2]));
  ss = 0.5*(aa+bb+cc);
  t_area = sqrt(ss*(ss-aa)*(ss-bb)*(ss-cc));
  return(t_area);
}

/*
 * loadpanel returns a list of panel structs derived from passed data:
 * shape, vertices, and type.
 */
panel *loadPanel(char *panelfile, char *density, int *numSing, ssystem *sys) {
  int i, j, k, ii, shape, type, nSurf, mesh_flag=sys->mesh_flag;
  panel *pnlList, *pnl;
  char fpath[256], fname[256];
  FILE *fp, *wfp;

  char c,c1[10],c2[10],c3[10],c4[10],c5[10],c6[10];
  double a1,a2,a3,b1,b2,b3;//a_norm,r0_norm,v0_norm;
  int i1,i2,i3,j1,j2,j3,ierr,iface[3],jface[3],ialert;
  double den,prob_rds,xx[3],yy[3],face[3][3],tface[3][3],s_area;
  double **sptpos, **sptnrm;
  int **extr_v, **extr_f, *nvert;
  double dist_local, area_local, cpuf;
  int nspt, natm, nface;

  double *nrm, len;

  /* read in vertices */
  sys->nChar = 0;
  sprintf(fpath,"test_proteins/");
  sprintf(fname,"%s%s.pqr",fpath,panelfile);
  fp=fopen(fname,"r");
  sprintf(fname,"%s%s.xyzr",fpath,panelfile);
  wfp=fopen(fname,"w");
  /* new version of pqr file, 11 entries per line */
  while(fscanf(fp,"%s %s %s %s %s %lf %lf %lf %lf %lf %s",c1,c2,c3,
               c4,c5,&a1,&a2,&a3,&b1,&b2,c6) != EOF){
  /* old version of pqr file, 10 entries per line */
  /* while(fscanf(fp,"%s %s %s %s %s %lf %lf %lf %lf %lf",c1,c2,c3,
               c4,c5,&a1,&a2,&a3,&b1,&b2) != EOF){ */
    if (strcmp(c1,"ATOM")==0)
    {
      fprintf(wfp,"%f %f %f %f\n",a1,a2,a3,b2);
      sys->nChar++;
    }
  }
  printf("PDB ID = %s\n",panelfile);
  printf("#Atoms = %d\n",sys->nChar);
  fclose(fp);
  fclose(wfp);

  CALLOC(sys->pos, 3*sys->nChar, double);
  CALLOC(sys->chr, sys->nChar, double);
  sprintf(fname,"%s%s.pqr",fpath,panelfile);
  fp=fopen(fname,"r");
  for ( i=0; i<sys->nChar; i++ ){
    /* new version of pqr file, 11 entries per line */
    ierr=fscanf(fp,"%s %s %s %s %s %lf %lf %lf %lf %lf %s",c1,c2,c3,
                 c4,c5,&a1,&a2,&a3,&b1,&b2,c6);
    /* old version of pqr file, 10 entries per line */
    /*ierr=fscanf(fp,"%s %s %s %s %s %lf %lf %lf %lf %lf",c1,c2,c3,
                 c4,c5,&a1,&a2,&a3,&b1,&b2);*/
    if (strcmp(c1,"ATOM")==0)
    {
      sys->pos[3*i]=a1;
      sys->pos[3*i+1]=a2;
      sys->pos[3*i+2]=a3;
      sys->chr[i]=b1;
    }
  }
  fclose(fp);

  if ( mesh_flag == 1 ) {
  /* run msms */
    sprintf(fname,"msms -if %s%s.xyzr -prob 1.4 -de %s -of %s%s > msms.output",
                fpath, panelfile, density, fpath, panelfile);
    //printf("%s\n",fname);
    ierr=system(fname);
    sprintf(fname,"rm msms.output");
    ierr=system(fname);
  } else if ( mesh_flag == 2 ) {
    wfp = fopen("surfaceConfiguration.prm", "w");
    fprintf(wfp, "Grid_scale = %s\n", density);
    fprintf(wfp, "Grid_perfil = 90.0\n");
    fprintf(wfp, "XYZR_FileName = %s%s.xyzr\n", fpath, panelfile);
    fprintf(wfp, "Build_epsilon_maps = false\n");
    fprintf(wfp, "Build_status_map = false\n");
    fprintf(wfp, "Save_Mesh_MSMS_Format = true\n");
    fprintf(wfp, "Compute_Vertex_Normals = true\n");
    fprintf(wfp, "Surface = ses\n");

    fprintf(wfp, "Smooth_Mesh = true\n");
    fprintf(wfp, "Skin_Surface_Parameter = %f\n", 0.45);

    fprintf(wfp, "Cavity_Detection_Filling = false\n");
    fprintf(wfp, "Conditional_Volume_Filling_Value = %f\n", 11.4);
    fprintf(wfp, "Keep_Water_Shaped_Cavities = false\n");
    fprintf(wfp, "Probe_Radius = %f\n", 1.4);
    fprintf(wfp, "Accurate_Triangulation = true\n");
    fprintf(wfp, "Triangulation = true\n");
    fprintf(wfp, "Check_duplicated_vertices = true\n");
    fprintf(wfp, "Save_Status_map = false\n");
    fprintf(wfp, "Save_PovRay = false\n");
    fclose(wfp);

    ierr = system("NanoShaper >> nsout.txt");
    remove("nsout.txt");

    sprintf(fname, "%s%s.face", fpath, panelfile);
    rename("triangulatedSurf.face", fname);
    sprintf(fname, "%s%s.vert", fpath, panelfile);
    rename("triangulatedSurf.vert", fname);

    remove("stderror.txt");
    remove("surfaceConfiguration.prm");
    remove("triangleAreas.txt");
    remove("exposed.xyz");
    remove("exposedIndices.txt");
  }
  /*======================================================================*/

  /* read in vert */
  sprintf(fname, "%s%s.vert",fpath,panelfile);
  //printf("%s\n",fname);

  /* open the file and read through the first two rows */
  fp=fopen(fname,"r");
  for ( i=1; i<=2; i++ ) {
    while ( (c=getc(fp))!='\n' ){
   }
  }

  if ( mesh_flag == 1 ) {
    ierr=fscanf(fp,"%d %d %lf %lf ",&nspt,&natm,&den,&prob_rds);
    //printf("nspt=%d, natm=%d, den=%lf, prob=%lf\n", nspt,natm,den,prob_rds);
  } else if ( mesh_flag == 2 ) {
    ierr = fscanf(fp, "%d", &nspt);
  }

  /*allocate variables for vertices file*/
  CALLOC(sptpos, 3, double*);
  CALLOC(sptnrm, 3, double*);
  for ( i=0; i<3; i++ ) {
    CALLOC(sptpos[i], nspt, double);
    CALLOC(sptnrm[i], nspt, double);
  }

  for ( i=0; i<=nspt-1; i++ ) {
    ierr=fscanf(fp,"%lf %lf %lf %lf %lf %lf %d %d %d",&a1,&a2,&a3,&b1,&b2,&b3,&i1,&i2,&i3);

    sptpos[0][i]=a1;
    sptpos[1][i]=a2;
    sptpos[2][i]=a3;
    sptnrm[0][i]=b1;
    sptnrm[1][i]=b2;
    sptnrm[2][i]=b3;
  }
  fclose(fp);

  /* read in faces */
  ierr=sprintf(fname, "%s%s.face",fpath,panelfile);
  fp=fopen(fname,"r");
  for ( i=1; i<=2; i++ ) { while ((c=getc(fp))!='\n'){} }

  if ( mesh_flag == 1 ) {
    ierr=fscanf(fp,"%d %d %lf %lf ",&nface,&natm,&den,&prob_rds);
  //printf("nface=%d, natm=%d, den=%lf, prob=%lf\n", nface,natm,den,prob_rds);
  } else if ( mesh_flag == 2 ) {
    ierr = fscanf(fp, "%d", &nface);
  }

  /* allocate variables for vertices file */
  CALLOC(nvert, 3*nface, int);

  for ( i=0; i<=nface-1; i++ ) {
    ierr=fscanf(fp,"%d %d %d %d %d",&j1,&j2,&j3,&i1,&i2);
    nvert[3*i]=j1;
    nvert[3*i+1]=j2;
    nvert[3*i+2]=j3;
  }
  fclose(fp);

  /* we delete ill performence triangles */
  s_area = 0.0;
  printf("#ele=%d, ",nface);
  //printf("Number of vertices = %d\n",nspt);

  pnlList = NULL;
  *numSing = 0;
  for ( i=0; i<nface; i++ ) {

    for ( j=0; j<3; j++ ) {
      iface[j]=nvert[3*i+j];
      xx[j]=0.0;
    }
    for ( j=0; j<3; j++ ) {
      for ( k=0; k<3; k++ ) {
        face[k][j]=sptpos[j][iface[k]-1];
	      xx[j] += 1.0/3.0*(face[k][j]);
      }
    }
    /* compute the area of each triangule */
    area_local = triangle_area(face);
    /* if the point is too close with 10 points infront, delete it */
    ialert=0;
    for ( j=i-10; (j>=0&&j<i); j++ ) { /* like k=max(1,i-10), i-1 in fortran */
      for ( k=0; k<3; k++ ) {
        jface[k] = nvert[3*j+k];
        yy[k]=0.0;
      }
      for ( k=0; k<3; k++ ) {
        for ( ii=0; ii<3; ii++ ) {
          tface[ii][k] = sptpos[k][jface[ii]-1];
          yy[k] += 1.0/3.0*(tface[ii][k]);
        }
      }
      dist_local = 0.0;
      for( k=0; k<3; k++ ) dist_local+=(xx[k]-yy[k])*(xx[k]-yy[k]);/* dot_product */
      //dist_local=sqrt(dist_local);
      if ( dist_local<1.0e-5 ) {
        ialert=1;
        //printf("particles %d and %d are too close: %e\n", i,j,dist_local);
      }
    }

    /* allocate and fill the panel */
    if ( area_local>=1.0e-5 && ialert == 0 ) {
      (*numSing)++;
      if ( pnlList == NULL ) {
        CALLOC(pnlList, 1, panel);
        pnl = pnlList;
      }
      else{
        CALLOC(pnl->next, 1, panel);
        pnl = pnl->next;
      }
  /* Fill in corners. */
      for ( j=0; j<3; j++ ) {
        for ( k=0; k<3; k++ ) {
          pnl->vtx[j][k] = sptpos[k][iface[j]-1];
          pnl->nrm[j][k] = sptnrm[k][iface[j]-1];
        }
      }

      for ( j=0; j<3; j++ ) {
        pnl->a[0][j] = pnl->vtx[2][j] - pnl->vtx[1][j];
        pnl->a[1][j] = pnl->vtx[0][j] - pnl->vtx[2][j];
        pnl->a[2][j] = pnl->vtx[1][j] - pnl->vtx[0][j];
      }

      nrm = pnl->normal;
      nrm[0] = pnl->nrm[0][0]/2. + (pnl->nrm[1][0] + pnl->nrm[2][0])/4.;
      nrm[1] = pnl->nrm[0][1]/2. + (pnl->nrm[1][1] + pnl->nrm[2][1])/4.;
      nrm[2] = pnl->nrm[0][2]/2. + (pnl->nrm[1][2] + pnl->nrm[2][2])/4.;
      len = sqrt(SQR(nrm[0]) + SQR(nrm[1]) + SQR(nrm[2]));
      for ( j=0; j<3; j++ ) nrm[j] /= len;

      pnl->shape = 3;
      pnl->nSurf = *numSing;
      pnl->area = area_local;
      s_area += area_local;
    }
  }

  printf("Area=%f \n",s_area);
  //printf("%d ugly faces are deleted\n", nface-*numSing);

  sprintf(fname,"rm %s%s.xyzr",fpath,panelfile);
  ierr=system(fname);
  sprintf(fname,"rm %s%s.vert",fpath,panelfile);
  ierr=system(fname);
  sprintf(fname,"rm %s%s.face",fpath,panelfile);
  ierr=system(fname);

  for (i=0;i<3;i++){
    free(sptpos[i]);
    free(sptnrm[i]);
  }
  free(sptpos);
  free(sptnrm);
  free(nvert);

  return pnlList;
} /* loadPanel */
