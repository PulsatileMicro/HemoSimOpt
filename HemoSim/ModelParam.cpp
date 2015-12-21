#include "ModelParam.h"
#include <stdio.h>
#include <math.h>

/*!
\brief Initialization of the static variables
*/
double ModelParam::dt=1e-3;
double ModelParam::rho=1050;
double ModelParam::oneD_alpha=1.26;
int ModelParam::nSteps=3200;
int ModelParam::hisSteps=1;
double ModelParam::relTol=1e-6;
double ModelParam::absTol=1e-6;
double ModelParam::riemannTol=1e-15;
double ModelParam::taperRate=0;
double ModelParam::scale_lamda=1;
double ModelParam::scale_u0=1;
double ModelParam::scale_r0=1;
int ModelParam::showCFL=0;
int ModelParam::showStepLapse=1;
int ModelParam::updVisc=0;
int ModelParam::riemann=1;
int ModelParam::nonDim=0;
int ModelParam::outMode=0;
int ModelParam::solverType=1;
int ModelParam::bifurloss=0;
int ModelParam::Ndoms=1;
int ModelParam::Nnodes=2;
int ModelParam::oneD_L=3;
int ModelParam::oneD_q=4;
int ModelParam::oneD_Je=2;
double *ModelParam::oneD_Uh=NULL;
double *ModelParam::oneD_Uhj=NULL;
double *ModelParam::oneD_Ah=NULL;
double *ModelParam::oneD_Ahj=NULL;
double *ModelParam::oneD_Hdh=NULL;
double *ModelParam::oneD_Hdhj=NULL;
double **ModelParam::oneD_Ufh=NULL;
double **ModelParam::oneD_Ufhj=NULL;
double **ModelParam::oneD_Afh=NULL;
double **ModelParam::oneD_Afhj=NULL;
int ModelParam::oneD_MaxQ=4;
int ModelParam::oneD_MaxL=3;
Domain *ModelParam::omega=NULL;
int ModelParam::argc=0;
char **ModelParam::argv=NULL;
int ModelParam::gammaI=0;
int ModelParam::gammaII=0;
double ModelParam::PI=4.0*atan(1.0); /* Number PI */
int ModelParam::sparseSolver=0;

ModelParam::ModelParam(){

}

ModelParam::~ModelParam(){

}
