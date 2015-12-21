#include "Element.h"
#include "ModelParam.h"
#include <mkl.h>
#include <math.h>

/**
\brief
Pressure evaluation at the quadrature point "i" of the element.
*/
double Element::Get_P(int i){
	return eval_P(h[i],Ao[i],beta[i]);
}

/**
\brief
Pressure evaluation at the quadrature point "i" of the element.
*/
double Element::PtoA(int i, double p){
	return eval_A(p,Ao[i],beta[i]);
}

/**
\brief
Pressure evaluation at the quadrature point "i" of the element for a given input area "A".
*/
double Element::Get_P(int i, double A){
	return eval_P(A,Ao[i],beta[i]);
}

/**
\brief
It evaluates pressure at all the quadrature points of the element.
TODO(panqing): Get_P should be parallelized
*/
void Element::Get_P(double *p){
	/*int i;
	for(i = 0; i < q; ++i)
	p[i] = eval_P(h[i],Ao[i],beta[i]);*/

	double h_sqrt[128], ao_sqrt[128], tmp[128];
	// Laplace Law
	// Add nondinmensionalization scale coef
	// p = pext + beta*(sqrt(A)-sqrt(ao));		// beta--g_rho*scale_u0*scale_u0*scale_r0
	vdSqrt(q,h,h_sqrt);
	vdSqrt(q,Ao,ao_sqrt);
	vdSub(q,h_sqrt,ao_sqrt,tmp);
	vdMul(q,beta,tmp,p);
}

/**
\brief
Wave speed evaluation at the quadrature point "i" of the element.
*/
double Element::Get_c(int i){
	return eval_c(h[i],beta[i]);
}

/**
\brief
Wave speed evaluation at the quadrature point "i" of the element for a given input area "A".
*/
double Element::Get_c(int i, double A){
	return eval_c(A,beta[i]);
}

/**
\brief
It evaluates the wave speed at all the quadrature points of the element.
*/
void Element::Get_c(double *c){
	int i;
	for(i = 0; i < q; ++i)
		c[i] = eval_c(h[i],beta[i]);
}

/**
\brief
Linear wave speed evaluation at the quadrature point "i" of the element.
*/
double Element::Get_co(int i){
	return eval_c(Ao[i],beta[i]);
}

/**
\brief
It evaluates the linear wave speed at all the quadrature points of the element.
*/
void Element::Get_co(double *c){
	int i;
	for(i = 0; i < q; ++i)
		c[i] = eval_c(Ao[i],beta[i]);
}

/**
\brief
Area part of the Riemann invariant evaluation at the quadrature point "i" of the element.
*/
double Element::Get_W(int i){
	return eval_W(h[i],Ao[i],beta[i]);
}

/**
\brief
Area part of the Riemann invariant evaluation at the quadrature point "i" of the element
for a given input area "A".
*/
double Element::Get_W(int i, double A){
	return eval_W(A,Ao[i],beta[i]);
}

/**
\brief
It evaluates the area part of the Rieamann invariant at all the quadrature points of the element.
*/
void Element::Get_W(double *w){

	int i;
	for(i = 0; i < q; ++i)
		w[i] = eval_W(h[i],Ao[i],beta[i]);
}

/**
\brief
It evaluates pressure according to the tube law considered and as a function of A, Ao and the
material elastic properties.
*/
double Element::eval_P(double A, double ao, double beta){
	double p, pext;
	pext = 0.0; // Pa
	// Laplace Law
	// Add nondinmensionalization scale coef
	p = pext + beta*(sqrt(A)-sqrt(ao));		// beta--g_rho*scale_u0*scale_u0*scale_r0
	return p;
}

/**
\brief It evaluates area according to the tube law considered and as a
function of P, Ao and the material elastic properties.
*/
double Element::eval_A(double p, double ao, double beta){
	double a;
	a = (p/beta+sqrt(ao))*(p/beta+sqrt(ao));
	return a;
}

/**
\brief
It evaluates the wave speed according to the tube law considered and as a function of A and the
material elastic properties..
*/
double Element::eval_c(double A, double beta){
	double c;
	// Pow is much slower than sqrt
	if(ModelParam::nonDim)
		// c = sqrt(0.5*beta/ModelParam::rho)*pow(A,0.25)*sqrt(ModelParam::rho);
		c = sqrt(0.5*beta*sqrt(A));
	else
		// c = sqrt(0.5*beta/ModelParam::rho)*pow(A,0.25);
		c = sqrt(0.5*beta/ModelParam::rho*sqrt(A));

	return c;
}

/**
\brief
It evaluates the area part of the Riemann invariant according to the tube law considered and as a
function of A, Ao and the material elastic properties.
*/
double Element::eval_W(double A, double ao, double beta){
	double W;
	// Pow is much slower than sqrt
	if(ModelParam::nonDim)
		// W = 4.0*sqrt(0.5*beta/ModelParam::rho)*(pow(A,0.25) - pow(ao,0.25))*sqrt(ModelParam::rho);
		W = 4.0*(sqrt(0.5*beta*sqrt(A))-sqrt(0.5*beta*sqrt(ao)));
	else
		// W = 4.0*sqrt(0.5*beta/g_rho)*(pow(A,0.25) - pow(ao,0.25));
		W = 4.0*(sqrt(0.5*beta/ModelParam::rho*sqrt(A))-sqrt(0.5*beta/ModelParam::rho*sqrt(ao)));

	return W;
}