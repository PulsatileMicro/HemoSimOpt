#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

double *dvector(int nl, int nh)
{
  double *v;

  v = (double*) malloc((nh-nl+1)*sizeof(double));
  for(int i=0; i<nh-nl+1; ++i)
    v[i] = 0.0;

  if(!v) fprintf(stderr, "out of memory for dvector");
  return v-nl;
}

float *svector(int nl, int nh)
{
  float *v;

  v = (float*) malloc((nh-nl+1)*sizeof(float));
  if(!v) fprintf(stderr, "out of memory for svector");
  return v-nl;
}

int *ivector(int nl, int nh)
{
  int *v;

  v = (int*) malloc((nh-nl+1)*sizeof(int));
  if(!v) fprintf(stderr, "out of memory for ivector");
  return v-nl;
}

void free_dvector(double *v, int nl)
{
  free(v+nl);
}

void free_svector(float *v, int nl)
{
  free(v+nl);
  return;
}

void free_ivector(int *v, int nl)
{
  free(v+nl);
  return;
}

void dvmul(int n, double *x, int incx, double *y, int incy,
           double *z, int incz)
{
  while( n-- ) {
    *z = (*x) * (*y);
    x += incx;
    y += incy;
    z += incz;
  }

  return;
}

void dvvtvp(int n, double *w, int incw, double *x, int incx, 
            double *y, int incy, double *z, int incz)
{
  while( n-- ) {
    *z = (*w) * (*x) + (*y);
    w += incw;
    x += incx;
    y += incy;
    z += incz;
  }
  return;
}

void dvvtvp_(int *np, double *w, int *iwp, double *x, int *ixp, 
             double *y, int *iyp, double *z, int *izp)
{
  register int n    = *np,
    incw = *iwp,
    incx = *ixp,
    incy = *iyp,
    incz = *izp;
  while( n-- ) {
    *z = (*w) * (*x) + (*y);
    w += incw;
    x += incx;
    y += incy;
    z += incz;
  }
  return;
}

void dvdiv(int n, double *x, int incx, double *y, int incy,
           double *z, int incz)
{
  while( n-- ) {
    *z = *x / *y;
    x += incx;
    y += incy;
    z += incz;
  }
  return;
}

void dsvtvp (int n, double alpha, 
             double *x, int incx, double *y, int incy, double *z, int incz)
{
  while (n--) {
    *z  = alpha * *x + *y;
    x  += incx;
    y  += incy;
    z  += incz;
  }
  return;
}

void izero(int n, int *x, const int incx)
{
  register int zero = 0;

  if(incx == 1)
    memset (x, '\0', n * sizeof(int));
  else
    while (n--) {
      *x = zero;
      x += incx;
    }

    return;
}

void dsmul (int n, double alpha, double *x, int incx, double *y, int incy)
{
  while (n--) {
    *y = alpha * (*x);
    x += incx;
    y += incy;
  }
  return;
}

void dvadd(int n, double *x, int incx, double *y, int incy,
           double *z, int incz)
{
  while( n-- ) {
    *z = *x + *y;
    x += incx;
    y += incy;
    z += incz;
  }
  return;
}

void dsvvmt (int n, double alpha, 
             double *x, int incx, double *y, int incy, double *z, int incz)
{
  while (n--) {
    *z = alpha * ( *x - *y );
    x += incx;
    y += incy;
    z += incz;
  }
  return;
}

void dvsqrt(int n, double *x, int incx, double *y, int incy)
{
  while (n--) {
    *y  = sqrt( *x );
    x  += incx;
    y  += incy;
  }
  return;
}