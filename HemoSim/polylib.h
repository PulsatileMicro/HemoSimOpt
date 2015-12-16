/*
 *  LIBRARY ROUTINES FOR POLYNOMIAL CALCULUS AND INTERPOLATION
 */

#ifndef PLYLIB_H  
#define PLYLIB_H

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

/*-----------------------------------------------------------------------
                         M A I N     R O U T I N E S
  -----------------------------------------------------------------------*/

/* Points and weights */

void   zwgj    (double *, double *, int , double , double);
void   zwgrj   (double *, double *, int , double , double);
void   zwglj   (double *, double *, int , double , double);

/* Derivative operators */

void   dgj     (double **, double **, double *, int, double, double);
void   dgrj    (double **, double **, double *, int, double, double);
void   dglj    (double **, double **, double *, int, double, double);

/* Lagrangian interpolants */

double hgj     (int, double, double *, int, double, double);
double hgrj    (int, double, double *, int, double, double);
double hglj    (int, double, double *, int, double, double);

/* Interpolation operators */

void   igjm  (double**, double*, double*, int, int, double, double);
void   igrjm (double**, double*, double*, int, int, double, double);
void   igljm (double**, double*, double*, int, int, double, double);

/* Polynomial functions */

void jacobf (int, double *, double *, int, double, double);
void jacobd (int, double *, double *, int, double, double);

/*-----------------------------------------------------------------------
                         M A C R O S
  -----------------------------------------------------------------------*/

/* Points and weights */

#define  zwgl( z,w,np)   zwgj( z,w,np,0.0,0.0);
#define  zwgrl(z,w,np)   zwgrj(z,w,np,0.0,0.0);
#define  zwgll(z,w,np)   zwglj(z,w,np,0.0,0.0);

#define  zwgc( z,w,np)   zwgj( z,w,np,-0.5,-0.5);
#define  zwgrc(z,w,np)   zwgrj(z,w,np,-0.5,-0.5);
#define  zwglc(z,w,np)   zwglj(z,w,np,-0.5,-0.5);

/* Derivative operators */

#define dgl( d,dt,z,np)  dgj( d,dt,z,np,0.0,0.0);
#define dgrl(d,dt,z,np)  dgrj(d,dt,z,np,0.0,0.0);
#define dgll(d,dt,z,np)  dglj(d,dt,z,np,0.0,0.0);

#define dgc( d,dt,z,np)  dgj( d,dt,z,np,-0.5,-0.5);
#define dgrc(d,dt,z,np)  dgrj(d,dt,z,np,-0.5,-0.5);
#define dglc(d,dt,z,np)  dglj(d,dt,z,np,-0.5,-0.5);

/* Lagrangian interpolants */

#define hgl( i,z,zgj ,np)  hgj( i,z,zgj ,np,0.0,0.0);
#define hgrl(i,z,zgrj,np)  hgrj(i,z,zgrj,np,0.0,0.0);
#define hgll(i,z,zglj,np)  hglj(i,z,zglj,np,0.0,0.0);

#define hgc( i,z,zgj ,np)  hgj( i,z,zgj ,np,-0.5,-0.5);
#define hgrc(i,z,zgrj,np)  hgrj(i,z,zgrj,np,-0.5,-0.5);
#define hglc(i,z,zglj,np)  hglj(i,z,zglj,np,-0.5,-0.5);

/* Interpolation operators */

#define iglm( im12,zgl ,zm,nz,mz) igjm( im12,zgl ,zm,nz,mz,0.0,0.0)
#define igrlm(im12,zgrl,zm,nz,mz) igrjm(im12,zgrl,zm,nz,mz,0.0,0.0)
#define igllm(im12,zgll,zm,nz,mz) igljm(im12,zgll,zm,nz,mz,0.0,0.0)

#define igcm( im12,zgl ,zm,nz,mz) igjm( im12,zgl ,zm,nz,mz,-0.5,-0.5)
#define igrcm(im12,zgrl,zm,nz,mz) igrjm(im12,zgrl,zm,nz,mz,-0.5,-0.5)
#define iglcm(im12,zgll,zm,nz,mz) igljm(im12,zgll,zm,nz,mz,-0.5,-0.5)

#endif          /* END OF POLYLIB.H DECLARATIONS */
