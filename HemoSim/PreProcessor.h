#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Domain.h"

#define BUFSIZE 1024
#define MAXFIELDS 3
#define max(a,b) ( (b) < (a) ? (a) : (b) )
#define min(a,b) ( (b) > (a) ? (a) : (b) )
#define error_msg(a)  {fprintf(stderr,a); exit(-1);}

class Element;

class PreProcessor{
public:
  void    run(int argc, char *argv[]);

private:
  void    ReadParams(FILE *fp);
  char    *findSection(char *name, char *buf, FILE *fp);
  Element *Setup(FILE *fp, int *nel, Domain &omega_, int n);
  void    ReadBC(int Ndoms, int nvar, Domain *omega, FILE *fp);
  void    SortNode(Domain *omega);
  Element *CopyElmt(Element *U,int nel);
  void    ReadIC(int Ndoms, int nfields, Domain *omega, FILE *fp);
};

/* Function found in Basis.C */
void getzw (int q, double **z, double **w, char dir);

/* Matrix fcns */
void library_init(Domain *Omega, int Ndoms);
void setup_basis(void);
void get_basis (double ***g1, int L, int q);
void get_dbasis(double ***g1, int L, int q);
void setup_dinfo(void);
void getD(int q, double ***d1, double ***d1t);

/* Vector & Matrix manipulation */
int *ivector(int rmin, int rmax);  
float *svector(int rmin, int rmax);  
double *dvector(int rmin, int rmax);          
int **imatrix(int rmin, int rmax, int cmin, int cmax); 
float **smatrix(int rmin, int rmax, int cmin, int cmax); 
double **dmatrix(int rmin, int rmax, int cmin, int cmax); 
void free_dvector(double *v, int rmin);
void free_svector(float  *v, int rmin);
void free_ivector(int    *v, int rmin);
void free_dmatrix(double **mat, int rmin, int cmin);
void free_smatrix(float  **mat, int rmin, int cmin);
void free_imatrix(int    **mat, int rmin, int cmin);
void ifill(int n, int    a, int     *x, int incx);
void dfill(int n, double a, double  *x, int incx);
void dvmul(int n, double *x, int incx, double *y, int incy, double *z, int incz);
void dvvtvp(int n, double *w, int incw, double *x, int incx, double *y, int incy, double *z, int incz);
void dvvtvp_(int *np, double *w, int *iwp, double *x, int *ixp, double *y, int *iyp, double *z, int *izp);
void dvdiv(int n, double *x, int incx, double *y, int incy, double *z, int incz);
void dsvtvp(int n, double alpha, double *x, int incx, double *y, int incy, double *z, int incz);
void izero(int n, int *x, const int incx);
void dsmul(int n, double alpha, double *x, int incx, double *y, int incy);
void dvadd(int n, double *x, int incx, double *y, int incy, double *z, int incz);
void dsvvmt(int n, double alpha, double *x, int incx, double *y, int incy, double *z, int incz);
void dvsqrt(int n, double *x, int incx, double *y, int incy);

#endif
