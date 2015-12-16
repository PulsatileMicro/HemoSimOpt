#include <stdlib.h>
#include <stdio.h>
#include "polylib.h"
#include <mkl_cblas.h>
#include "PreProcessor.h"

int LZero = 0;

/* zeros weights info */
typedef struct zwinfo {
	int     q;   /* number of quadrature points                       */
	char    dir; /* direction of zeros/weights either 'a', 'b' or 'c' */
	double *z;   /* zeros   */
	double *w;   /* weigths */
	struct zwinfo *next;
} ZWinfo;

/*-------------------------------------------------------------------*
*  link list for zeros and weights for jacobi polynomial in [0,1]   *
*  at the gauss points. Note: alpha > -1 , beta > -1                *
*-------------------------------------------------------------------*/

static ZWinfo *zwinf,*zwbase;
static ZWinfo *addzw(int, char);

void getzw(int q, double **z, double **w, char dir){
	/* check link list */

	for(zwinf = zwbase; zwinf; zwinf = zwinf->next)
		if(zwinf->q == q)
			if(zwinf->dir == dir){
				*z = zwinf->z;
				*w = zwinf->w;
				return;
			}

			/* else add new zeros and weights */

			zwinf  = zwbase;
			zwbase = addzw(q,dir);
			zwbase->next = zwinf;

			*z = zwbase->z;
			*w = zwbase->w;

			return;
}

/*--------------------------------------------------------------------------*
* This function gets the zeros and weights in [-1,1]. In the 'a' direction *
* alpha = 0, beta = 0. In the 'b' alpha = 1, beta = 0 and in the 'c'       *
* direction alpha = 2, beta = 0.                                           *
* Note: that the 'b','c' direction weights are scaled by (1/2)^alpha       *
* to be consistent with the jacobean d(r,s)/d(t,s)                         *
*--------------------------------------------------------------------------*/

static ZWinfo *addzw(int n, char dir){  
	int i;
	ZWinfo *zw = (ZWinfo *)calloc(1,sizeof(ZWinfo));

	zw->q     = n;
	zw->dir   = dir;

	zw->z = dvector(0,n-1);
	zw->w = dvector(0,n-1);
	switch(dir){
	case 'a':
		zwgll(zw->z,zw->w,n);
		break;
	case 'b':
		if(LZero){ // use Legendre weights 
			zwgrj(zw->z,zw->w,n,0.0,0.0);
			for(i = 0; i < n; ++i)
				zw->w[i] *= (1-zw->z[i])*0.5;
		}
		else{
			double a = 0.5;
			zwgrj(zw->z,zw->w,n,1.0,0.0);
			cblas_dscal(n,0.5,zw->w,1);
		}
		break;
	case 'c':
		if(LZero){
			zwgrj(zw->z,zw->w,n,0.0,0.0);
			for(i = 0; i < n; ++i){
				zw->w[i] *= (1-zw->z[i])*0.5;
				zw->w[i] *= (1-zw->z[i])*0.5;
			}
		}
		else{
			zwgrj(zw->z,zw->w,n,2.0,0.0);
			cblas_dscal(n,0.25,zw->w,1);
		}
		break;
	case 'g': /* just Gauss points */
		zwgl (zw->z,zw->w,n);
		break;
	case 'h': /* Gauss points with 1,0 weight */
		if(LZero){
			zwgj (zw->z,zw->w,n,0.0,0.0);
			for(i = 0; i < n; ++i)
				zw->w[i] *= (1-zw->z[i])*0.5;
		}
		else{
			zwgj (zw->z,zw->w,n,1.0,0.0);
			cblas_dscal(n,0.5,zw->w,1);
		}
		break;
	case 'i':
		if(LZero){
			zwgj (zw->z,zw->w,n,0.0,0.0);
			for(i = 0; i < n; ++i){
				zw->w[i] *= (1-zw->z[i])*0.5;
				zw->w[i] *= (1-zw->z[i])*0.5;
			}
		}
		else{
			zwgj (zw->z,zw->w,n,2.0,0.0);
			cblas_dscal(n,0.25,zw->w,1);
		}
		break;
	case 'j':
		zwgrj(zw->z,zw->w,n,0.0,0.0);
		break;
	default:
		fprintf(stderr,"addzw - incorrect direction specified");
		break;
	}

	return zw;
}
