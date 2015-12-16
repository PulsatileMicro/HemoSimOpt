#include "Domain.h"
#include "polylib.h"
#include <math.h>
#include <mkl.h>
#include "ModelParam.h"
#include "PreProcessor.h"

/**
\brief
It defines QGmax = max(QGmax, U->q) and LGmax = max(LGmax, U->L) taking into accound all the elements of
all the domains.
*/
void library_init(Domain *omega, int Ndoms){
	register int i;
	Element *U;

  ModelParam::oneD_MaxQ = ModelParam::oneD_MaxL = 0;

	/* set up info */
	for(i = 0; i < Ndoms; ++i)
		for(U=omega[i].U;U;U = U->next){
			ModelParam::oneD_MaxQ = max(ModelParam::oneD_MaxQ,U->q);
			ModelParam::oneD_MaxL = max(ModelParam::oneD_MaxL,U->L);
		}

		setup_basis();
		setup_dinfo();
}

typedef struct dinfo {
	double ***D1;     /* Differential matrix for a direction                  */
	double ***D1t;    /* Differential matrix for a direction (transposed)     */
} Dinfo;

static Dinfo Dbase;

void setup_dinfo(void){
	Dbase.D1  = (double ***)calloc(ModelParam::oneD_MaxQ,sizeof(double**))-1;
	Dbase.D1t = (double ***)calloc(ModelParam::oneD_MaxQ,sizeof(double**))-1;
}

#ifndef NODES
typedef struct basis {
	double ***g1;
	double ***dg1;
} Basis;

static Basis base;

void setup_basis(void){ /* modal setup */
	base.g1 = (double ***)calloc(ModelParam::oneD_MaxQ,sizeof(double **))-1;
	base.dg1 = (double ***)calloc(ModelParam::oneD_MaxQ,sizeof(double **))-1;
}

static void set_g1(Basis *be, int q);
static void set_g1_dis(Basis *be, int q);

void get_basis(double ***g1, int L, int q){
	if(!base.g1[q]){
		set_g1_dis(&base,q); 
	}

	*g1 = base.g1[q];
}

void get_dbasis(double ***dg1, int L, int q){

	if(!base.g1[q]){
		set_g1_dis(&base,q); 
	}
	*dg1 = base.dg1[q];
}

//static void set_g1(Basis *be, int q){
//	register int i;
//	int      L = q;
//	double  **d,**dt;
//	double  *z,*w,**s,**sd;
//
//	getzw(q,&z,&w,'a');
//
//	s  = be->g1 [q] = dmatrix(0,L-1,0,q-1); 
//	sd = be->dg1[q] = dmatrix(0,L-1,0,q-1); 
//
//	// dsadd(q,-1.0,z,1,s[0],1);            /*  (1-s)/2                     */
//	// dscal(q,-0.5,s[0],1);
//	// dsadd(q, 1.0,z,1,s[1],1);            /*  (1+s)/2                     */
//	// dscal(q, 0.5,s[1],1); 
//
//	if(L>2){
//		vdMul(q,s[0],s[1],s[2]);     /* (1-s)/2*(1+s)/2              */
//
//		for(i = 3; i < L; ++i){           
//			/* (1-s)/2*(1+s)/2 P^{1,1}_l(s) */
//#ifdef LEG
//			jacobf(q,z,s[i],i-2,0.0,0.0);
//#else
//#ifdef MDIAG
//			jacobf(q,z,s[i],i-2,2.0,2.0);
//#else 
//			jacobf(q,z,s[i],i-2,1.0,1.0);
//#endif
//#endif
//			vdMul (q,s[2],s[i],s[i]);
//		}
//	}
//
//	getD(q,&d,&dt);
//	for(i = 0; i < L; ++i)
//		cblas_dgemv(CblasColMajor,CblasTrans,q,q,1.0,*d,q,s[i],1,0.0,sd[i],1);
//
//}

static void set_g1_dis(Basis *be, int q){
	register int i;
	/* if q!=1, q is set to max(q,L+1) for ensuring stable quadrature integration (in setup module)
	if q==1, there's no need for quadrature integration, so L can be any value */
	int L;
	L = q;
	// L = g_L;

	double  **d,**dt;
	double  *z,*w,**s,**sd;

	getzw(q,&z,&w,'a');

	s  = be->g1 [q] = dmatrix(0,L-1,0,q-1);
	sd = be->dg1[q] = dmatrix(0,L-1,0,q-1);

	for(i = 0; i < L;  ++i){
		jacobf(q, z, s[i], i, 0.0, 0.0);
		/* normalise so that  (g,g)=1*/
		cblas_dscal(q,sqrt(0.5*(2.0*i+1.0)),s[i],1);
	}

	getD(q,&d,&dt);
	for(i = 0; i < L; ++i)
		cblas_dgemv(CblasColMajor,CblasTrans,q,q,1.0,*d,q,s[i],1,0.0,sd[i],1);
}

#else
// DEFINED NODES
typedef struct basis {
	double ****g1;
	double ****dg1;
} Basis;

static Basis base;

void setup_basis(void){ /* nodal setup */
	register int i;
	base.g1     = (double ****)malloc(LGmax*sizeof(double ***))-1;
	base.g1[1]  = (double *** )calloc(QGmax*LGmax,sizeof(double **))-1;

	base.dg1    = (double ****)malloc(LGmax*sizeof(double ***))-1;
	base.dg1[1] = (double *** )calloc(QGmax*LGmax,sizeof(double **))-1;

	for(i = 2; i <= LGmax; ++i){
		base.g1[i]  = base.g1 [i-1] + QGmax -1;
		base.dg1[i] = base.dg1[i-1] + QGmax -1;
	}    
}

static void set_g1(Basis *be, int L, int q);

void get_basis(double ***g1, int L, int q){
	if(!base.g1[L][q]) set_g1(&base,L,q); 
	*g1 = base.g1[L][q];
}

void get_dbasis(double ***dg1, int L, int q){
	if(!base.g1[L][q]) set_g1(&base,L,q); 
	*dg1 = base.dg1[L][q];
}

static void set_g1(Basis *be, int L, int q){
	register int i,j;
	double  **d,**dt;
	double  *z,*z1,*w,**s,**sd;

	getzw(L,&z ,&w,'a');
	getzw(q,&z1,&w,'a');

	s  = be-> g1[L][q] = dmatrix(0,L-1,0,q-1); 
	sd = be->dg1[L][q] = dmatrix(0,L-1,0,q-1); 

	for(i = 0; i < L; ++i)
		for(j = 0; j < q; ++j)
			s[i][j] = hgll(i,z1[j],z,L);

	getD(q,&d,&dt);
	for(i = 0; i < L; ++i)
		cblas_dgemv(CblasColMajor,CblasTrans,q,q,1.0,*d,q,s[i],1,0.0,sd[i],1);

}
#endif

static void set_D1(Dinfo *D, int q);

void getD(int q, double ***d1, double ***d1t){
	if(!Dbase.D1[q])
		set_D1(&Dbase,q);
	*d1  = Dbase.D1[q];
	*d1t = Dbase.D1t[q];
}

static void set_D1(Dinfo *D, int q){
	double *z,*w;

	D->D1 [q] = dmatrix(0,q-1,0,q-1);
	D->D1t[q] = dmatrix(0,q-1,0,q-1);

	getzw(q,&z,&w,'a');

	dgll(D->D1[q],D->D1t[q],z,q);
}