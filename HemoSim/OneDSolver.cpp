#include "OneDSolver.h"
#include "ModelParam.h"
#include "PreProcessor.h"
#include "polylib.h"

double * OneDSolver::m_w=NULL;
double * OneDSolver::m_w_mat=NULL;
double **OneDSolver::m_g1=NULL;
double * OneDSolver::m_g1_mat=NULL;
double **OneDSolver::m_d=NULL;
double **OneDSolver::m_dt=NULL;
double * OneDSolver::ip_tmp_Umat=NULL;
double * OneDSolver::ip_tmp_Amat=NULL;
double * OneDSolver::unitDiag;
int *    OneDSolver::headerPrinted=NULL;

/*
CVODE Vars
*/
N_Vector OneDSolver::init_Y;
void *OneDSolver::cvode_mem=NULL;
realtype OneDSolver::tret;
int OneDSolver::flag;
int OneDSolver::reIntFlag;
double *OneDSolver::cvode_u=NULL;
double *OneDSolver::cvode_a=NULL;
double *OneDSolver::cvode_yu=NULL;
double *OneDSolver::cvode_ya=NULL;
double *OneDSolver::cvode_ydotu=NULL;
double *OneDSolver::cvode_ydota=NULL;
N_Vector OneDSolver::f1_Vec;
N_Vector OneDSolver::f2_Vec;

/*
File vars
*/
double OneDSolver::t_in;
double *OneDSolver::pff=NULL;
double *OneDSolver::pbb=NULL;
double *OneDSolver::uff=NULL;
double *OneDSolver::ubb=NULL;
double OneDSolver::Rt;
double *OneDSolver::q_in=NULL;
double *OneDSolver::A_star=NULL;
double *OneDSolver::u_star=NULL;
char OneDSolver::buf[BUFSIZE];
FILE **OneDSolver::fp_his;

/*
Timing vars
*/
clock_t OneDSolver::start_time;
clock_t OneDSolver::finish_time;
double OneDSolver::timelapse;
double OneDSolver::lastClock;
double OneDSolver::curTime;

/*
Energy loss at the bifurcations
*/
int 		nbif_div_right = 0; /**< Number of bifurcations with diverging flow to the right. */
int 		nbif_div_left = 0; /**< Number of bifurcations with diverging flow to the left. */
int 		nbif_conv_right = 0; /**< Number of bifurcations with converging flow to the right. */
int 		nbif_conv_left = 0; /**< Number of bifurcations with converging flow to the left. */
int 		nbif_ent = 0; /**< Number of bifurcations with flow entering the bifurcation through the side branch only. */
int 		nbif_leav = 0; /**< Number of bifurcations with flow leaving the bifurcation through the side branch only. */
int 		nbif_no_losses = 0; /**< Number of bifurcations with no energy losses because the combined flow is too small. */

void OneDSolver::initSolver(int solverType){
	// On the Macbook Air, VML functions do not work unless this is enabled
	// It might be a bug of MKL
	mkl_enable_instructions(MKL_SINGLE_PATH_ENABLE);

	// mkl_set_num_threads(4);
	/*
	Init the vars for integration and projection
	*/
	double *z;
	getD(ModelParam::oneD_q,&m_d,&m_dt);
	getzw(ModelParam::oneD_q,&z,&m_w,'a'); // Get the quadrature points, z, and the quadrature weights, w.
	m_w_mat = (double *)malloc(ModelParam::oneD_q*ModelParam::oneD_q*sizeof(double));
	for(int i=0; i<ModelParam::oneD_q; ++i)
	{
		for(int j=0; j<ModelParam::oneD_q; ++j)
		{
			if(i==j)
				m_w_mat[i*ModelParam::oneD_q+j] = m_w[j];
			else
				m_w_mat[i*ModelParam::oneD_q+j] = 0;
		}
	}

	get_basis(&m_g1,ModelParam::oneD_L,ModelParam::oneD_q); // Get the basis, g1, of the projected space (Legendre polynomials).
	m_g1_mat = (double *)malloc(ModelParam::oneD_L*ModelParam::oneD_q*sizeof(double));
	for(int i=0; i<ModelParam::oneD_q; ++i){
		for(int j=0; j<ModelParam::oneD_L; ++j){
			m_g1_mat[i*ModelParam::oneD_L+j] = m_g1[j][i];
			// printf("%d\n",i*ModelParam::oneD_L+j);
		}
	}

	// Init Unit diagnoal matrix
	unitDiag = (double *)malloc(ModelParam::oneD_q*ModelParam::oneD_q*sizeof(double));
	for(int i=0; i<ModelParam::oneD_q; ++i){
		for(int j=0; j<ModelParam::oneD_q; ++j){
			if(i==j)
				unitDiag[i*ModelParam::oneD_q+j] = 1;
			else
				unitDiag[i*ModelParam::oneD_q+j] = 0;
		}
	}

	// InnerProduct tmp var
	ip_tmp_Umat = (double *)malloc(ModelParam::Ndoms*ModelParam::oneD_q*sizeof(double));
	ip_tmp_Amat = (double *)malloc(ModelParam::Ndoms*ModelParam::oneD_q*sizeof(double));

	/*
	Init headerPrinted flag
	*/
	headerPrinted = (int *)malloc(ModelParam::Ndoms*sizeof(int));
	for (int i=0;i<ModelParam::Ndoms;++i){
		headerPrinted[i]=0;
	}

	// CVODE Initialization
	if(solverType==ModelParam::oneD_EXP)
		cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
	else if(solverType==ModelParam::oneD_IMP)
		cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	init_Y = N_VNew_Serial(ModelParam::oneD_q*2*ModelParam::Ndoms);
	cvode_u = NV_DATA_S(init_Y);
	cvode_a = cvode_u + ModelParam::oneD_q*ModelParam::Ndoms;

	f1_Vec=N_VNew_Serial(ModelParam::oneD_q*2*ModelParam::Ndoms);
	f2_Vec=N_VNew_Serial(ModelParam::oneD_q*2*ModelParam::Ndoms);
}

void OneDSolver::solve(){
	for (int i=0; i<ModelParam::nSteps; ++i){
		cblas_dcopy(ModelParam::oneD_q*ModelParam::Ndoms,ModelParam::oneD_Uh,1,cvode_u,1);
		cblas_dcopy(ModelParam::oneD_q*ModelParam::Ndoms,ModelParam::oneD_Ah,1,cvode_a,1);

		if(!i){
			CVodeInit(cvode_mem, CVODE_RHS, 0, init_Y);					// Init the ODE solver
			CVodeSStolerances(cvode_mem, ModelParam::relTol, ModelParam::absTol);		// Set the absolute & relative tolerance of the ODE solver
			// CVodeSetMaxOrd(cvode_mem, 2);										// Set the maximum order of the solver (1 to 5)
			CVodeSetMaxNumSteps(cvode_mem, 5000);								// Set the maximum iteration steps
		}
		else
			CVodeReInit(cvode_mem, 0, init_Y);

		CVLapackDense(cvode_mem, ModelParam::oneD_q*2*ModelParam::Ndoms);
		// CVDlsSetDenseJacFn(cvode_mem, Jac);
		while(1){
			if(CV_SUCCESS == CVode(cvode_mem,ModelParam::dt/ModelParam::scale_lamda*ModelParam::scale_u0,init_Y,&tret,CV_NORMAL)){
				reIntFlag=0;
				break;
			}
		}

		/* Copy the results back to omega */
		cblas_dcopy(ModelParam::oneD_q*ModelParam::Ndoms,cvode_u,1,ModelParam::oneD_Ufh[0],1);
		cblas_dcopy(ModelParam::oneD_q*ModelParam::Ndoms,cvode_a,1,ModelParam::oneD_Afh[0],1);

		// Update
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ModelParam::oneD_q,ModelParam::Ndoms,ModelParam::oneD_q,1,
			unitDiag,ModelParam::oneD_q,ModelParam::oneD_Ufh[0],ModelParam::oneD_q,0,ModelParam::oneD_Uh,ModelParam::oneD_q);
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ModelParam::oneD_q,ModelParam::Ndoms,ModelParam::oneD_q,1,
			unitDiag,ModelParam::oneD_q,ModelParam::oneD_Afh[0],ModelParam::oneD_q,0,ModelParam::oneD_Ah,ModelParam::oneD_q);

		if(ModelParam::updVisc){
			/* Update the viscosity and check the convergence (only when Hd is also a var) */
			UpdateVisc(ModelParam::omega);
		}

		curTime = (i+1)*ModelParam::dt;
		/* It is checked if the solution in the .his file can be written according to "hisstep". If it
		can, Write_history does the job. */
		if(ModelParam::hisSteps&&(!((i+1)%ModelParam::hisSteps))){
			Write_history(ModelParam::omega,ModelParam::argv[ModelParam::argc-1]);
			// BifurHd(omega);
			if(ModelParam::showCFL){
				/* The current time step, time and CFL number are printed on the screen. */
				fprintf(stdout,"Step %d: Time %lf: CFL = %lf \n",i+1,(i+1)*ModelParam::dt,CFL(ModelParam::omega));
			}
		}
	}
}

void OneDSolver::destroySolver(){
	CVodeFree(&cvode_mem);
}

int OneDSolver::CVODE_RHS(realtype t, N_Vector y, N_Vector ydot, void *user_data){
	int m;
	/* Get the pointer to the data of N_Vector struct */
	cvode_yu = NV_DATA_S(y);
	cvode_ya = cvode_yu + ModelParam::oneD_q*ModelParam::Ndoms;
	cvode_ydotu = NV_DATA_S(ydot);
	cvode_ydota = cvode_ydotu + ModelParam::oneD_q*ModelParam::Ndoms;

	cblas_dcopy(ModelParam::oneD_q*ModelParam::Ndoms,cvode_yu,1,ModelParam::oneD_Uh,1);
	cblas_dcopy(ModelParam::oneD_q*ModelParam::Ndoms,cvode_ya,1,ModelParam::oneD_Ah,1);

	Eval_RHS(ModelParam::argc, ModelParam::argv);
	reIntFlag=1;

	// cvode_ydotu,a store the right hand side of the ODE
	cblas_dcopy(ModelParam::oneD_q*ModelParam::Ndoms,ModelParam::oneD_Ufh[0],1,cvode_ydotu,1);
	cblas_dcopy(ModelParam::oneD_q*ModelParam::Ndoms,ModelParam::oneD_Afh[0],1,cvode_ydota,1);

	return 0;
}

/*!
\brief Evaluate the right hand side of the ODE
The form of the right hand side could be referred to Dr. Alastruey's doctoral thesis
*/
void OneDSolver::Eval_RHS(int argc, char *argv[]){
	BioAdv       (ModelParam::omega); /* Generate the space derivative of the flux terms and store them in the
									  the "h" component (physical space) of the elements Af and Uf. */
	BioSource    (ModelParam::omega); /* Generate the source term that takes into account the viscous dissipation
									  of the fluid due to the shear stress at the arterial wall and add the results to the
									  previous value of the "h" component (physical space) of the element Uf. */

	// InnerProduct
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ModelParam::oneD_q,ModelParam::Ndoms,ModelParam::oneD_q,1,m_w_mat,ModelParam::oneD_q,ModelParam::oneD_Ufh[0],ModelParam::oneD_q,0,ip_tmp_Umat,ModelParam::oneD_q);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ModelParam::oneD_q,ModelParam::Ndoms,ModelParam::oneD_q,1,m_w_mat,ModelParam::oneD_q,ModelParam::oneD_Afh[0],ModelParam::oneD_q,0,ip_tmp_Amat,ModelParam::oneD_q);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ModelParam::oneD_L,ModelParam::Ndoms,ModelParam::oneD_q,1,m_g1_mat,ModelParam::oneD_L,ip_tmp_Umat,ModelParam::oneD_q,0,ModelParam::oneD_Ufhj[0],ModelParam::oneD_L);
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ModelParam::oneD_L,ModelParam::Ndoms,ModelParam::oneD_q,1,m_g1_mat,ModelParam::oneD_L,ip_tmp_Amat,ModelParam::oneD_q,0,ModelParam::oneD_Afhj[0],ModelParam::oneD_L);

	Bifur        (ModelParam::omega); /* Set bc's at the bifurcations (both merging and splitting flow cases).
									  They are stored in bcval of each structure omega[d].*/
	BioFlux      (ModelParam::omega,argv[argc-1]); /* Determines the upwind variables at each interelement boundary by
												   solving a Riemann problem and performs the difference between the upwind
												   and the local fluxes, taking into account the bc's enforced at both sides
												   of each domain. The result is multiplied by the Legendre polynomial
												   basis and the inverse of the Jacobian. It is added to the previous value
												   of the "hj" component (projected space) of the elements Af and Uf. */
	// Transbwd
	cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,ModelParam::oneD_q,ModelParam::Ndoms,ModelParam::oneD_L,1,m_g1_mat,ModelParam::oneD_L,ModelParam::oneD_Ufhj[0],ModelParam::oneD_L,0,ModelParam::oneD_Ufh[0],ModelParam::oneD_q);
	cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,ModelParam::oneD_q,ModelParam::Ndoms,ModelParam::oneD_L,1,m_g1_mat,ModelParam::oneD_L,ModelParam::oneD_Afhj[0],ModelParam::oneD_L,0,ModelParam::oneD_Afh[0],ModelParam::oneD_q);
}

/** \brief
For each element \f$ i = 0, \cdots , nel \f$ of each domain \f$ n = 0, \cdots , Ndoms -1 \f$, it
generates the flux term
\f$
\mathbf{F} = \left[ \begin{array}{c}
F_{1} \\ F_{2}
\end{array} \right]
= \left[ \begin{array}{c}
uA \\ \frac{u^{2}}{2} + \frac{p}{\rho}
\end{array} \right] \; \; \; \;

\f$
and differenciates this flux with respect to \f$ x \f$.
\f$ \frac{\partial F_{1}}{\partial x} \f$  is stored in the "h" component of the element \f$ Af \f$ and
\f$ \frac{\partial F_{2}}{\partial x} \f$  is stored in the "h" component of the element \f$ Uf \f$ (both
in the physical space).
*/
void OneDSolver::BioAdv(Domain *omega){
	register int i, j, n, k;
	int nel;
	Element  *U, *A, **Uf, **Af;
	static double af[MAX_Q], uf[MAX_Q];
	double tmp_af[MAX_Q];

	for(n = 0; n < ModelParam::Ndoms; ++n){
		nel = omega[n].nel;
		U   = omega[n].U;
		A   = omega[n].A;
		Uf  = omega[n].Uf;
		Af  = omega[n].Af;

		// Generate flux terms and their space derivatives.
		for(i = 0; i < nel; ++i){
			// af = F1 = uA (Conservation of mass)
			vdmul(&ModelParam::oneD_q,U[i].h,A[i].h,af);
			// dvmul(ModelParam::oneD_q,U[i].h,1,A[i].h,1,af,1);
			cblas_dgemv(CblasColMajor,CblasTrans,ModelParam::oneD_q,ModelParam::oneD_q,-A[i].rx,*m_d,ModelParam::oneD_q,af,1,0.0,Af[0][i].h,1);

			// (Conservation of momentum)
			// uf = F2 = 0.5*u^2 + p/rho
			if(!ModelParam::gammaI && !ModelParam::gammaII){
				A[i].Get_P(uf); //uf = p
				/* In the nondimensionalized formula, no rho */
				if(ModelParam::nonDim)
					cblas_dscal(ModelParam::oneD_q,2.0,uf,1);  // uf = 2*p
				else
					cblas_dscal(ModelParam::oneD_q,2.0/ModelParam::rho,uf,1);  // uf = 2*p/rho

				dvvtvp(ModelParam::oneD_q,U[i].h,1,U[i].h,1,uf,1,uf,1); //uf = u^2 + 2*p/rho
				cblas_dscal(ModelParam::oneD_q,0.5,uf,1);  //  uf = uf/2
			}
			else{
				if(ModelParam::gammaI){  // New viscous wall model
					// uf = F2 = 0.5*u^2 + p/rho - gamma/rho/sqrt(A) * dQ/dx
					vdsqrt(&ModelParam::oneD_q,A[i].h,uf);	//uf=sqrt(A)
					cblas_dgemv(CblasColMajor,CblasTrans,ModelParam::oneD_q,ModelParam::oneD_q,-1.0,*m_d,ModelParam::oneD_q,af,1,0.0,tmp_af,1);	//af=-dQ/dx
					vdmul(&ModelParam::oneD_q,omega[n].gamma,tmp_af,tmp_af);	//af=af*gamma (-gamma * dQ/dx)
					// dvmul(ModelParam::oneD_q,omega[n].gamma,1,tmp_af,1,tmp_af,1);
					if(!ModelParam::nonDim)
						cblas_dscal(ModelParam::oneD_q,1.0/ModelParam::rho,tmp_af,1); //af=af/rho (-gamma/rho * dQ/dx)
					vddiv(&ModelParam::oneD_q,tmp_af,uf,af); //af = -gamma/rho * dQ/dx / sqrt(A)
					// dvdiv(ModelParam::oneD_q,tmp_af,1,uf,1,af,1);
					// for(k=0;k<ModelParam::oneD_q;k++) printf("%d, %lf, %lf\n", k, af[k]*1e6, uf[k]*1e6);

					A[i].Get_P(uf); //uf = p_elastic
					if(ModelParam::nonDim)
						cblas_dscal(ModelParam::oneD_q,2.0,uf,1);
					else
						cblas_dscal(ModelParam::oneD_q,2.0/ModelParam::rho,uf,1);
					dvvtvp(ModelParam::oneD_q,U[i].h,1,U[i].h,1,uf,1,uf,1); //uf = u^2 + 2*p/rho
					cblas_dscal(ModelParam::oneD_q,0.5,uf,1);
					// for(k=0;k<ModelParam::oneD_q;k++) printf("%d, %lf, %lf\n", k, af[k]*1e9, uf[k]*1e9);
					vdadd(&ModelParam::oneD_q,uf,af,uf); // uf = p_elastic+p_visc
					// dvadd(ModelParam::oneD_q,uf,1,af,1,uf,1);
				} 
				else if(ModelParam::gammaII){ // Old viscous wall model
					// uf = F2 = 0.5*u^2 + p/rho - gamma/rho * dQ/dx
					cblas_dgemv(CblasColMajor,CblasTrans,ModelParam::oneD_q,ModelParam::oneD_q,-1.0,*m_d,ModelParam::oneD_q,af,1,0.0,tmp_af,1);	//tmp_af=-dQ/dx
					vdmul(&ModelParam::oneD_q,omega[n].gamma,tmp_af,af);	//af=tmp_af*gamma (-gamma * dQ/dx)
					// dvmul(ModelParam::oneD_q,omega[n].gamma,1,tmp_af,1,af,1);
					if(!ModelParam::nonDim)
						cblas_dscal(ModelParam::oneD_q,1.0/ModelParam::rho,af,1); //af=af/rho (-gamma/rho * dQ/dx)

					A[i].Get_P(uf); //uf = p_elastic
					if(ModelParam::nonDim)
						cblas_dscal(ModelParam::oneD_q,2.0,uf,1);
					else
						cblas_dscal(ModelParam::oneD_q,2.0/ModelParam::rho,uf,1);
					dvvtvp(ModelParam::oneD_q,U[i].h,1,U[i].h,1,uf,1,uf,1); //uf = u^2 + 2*p/rho
					cblas_dscal(ModelParam::oneD_q,0.5,uf,1);
					// for(k=0;k<ModelParam::oneD_q;k++) printf("%d, %lf, %lf\n", k, af[k]*1e9, uf[k]*1e9);
					vdadd(&ModelParam::oneD_q,uf,af,uf); // uf = p_elastic+p_visc
					// dvadd(ModelParam::oneD_q,uf,1,af,1,uf,1);
				}
			}

			// Space derivative
			cblas_dgemv(CblasColMajor,CblasTrans,ModelParam::oneD_q,ModelParam::oneD_q,-U[i].rx,*m_d,ModelParam::oneD_q,uf,1,0.0,Uf[0][i].h,1);
		}
	}
}

/**
\brief
It adds the source term, which takes into account the dissipation due to the wall 
shear stress, to the "h" (physical space) component of the element \f$ Uf \f$.
*/
void OneDSolver::BioSource(Domain *omega){
	int i,j,n,nel;
	Element  *U, *A, **Uf;
	double   *temp_visc;
	static double tmp[MAX_Q], Kr[MAX_Q];

	if(ModelParam::oneD_alpha == 1.0) return;

	for(n = 0; n < ModelParam::Ndoms; ++n){
		nel = omega[n].nel;
		U   = omega[n].U;
		A   = omega[n].A;
		Uf  = omega[n].Uf;

		/* Different viscosity for each vessel & each quadrature point
		* Added by panqing. 2011-7-30 */
		temp_visc = omega[n].visc;
		if(ModelParam::nonDim){
			/* Nondimensionalization */
			double tmp_coef = -2.0*ModelParam::oneD_alpha*ModelParam::PI/ModelParam::rho/
				(ModelParam::oneD_alpha-1.0)*ModelParam::scale_lamda/ModelParam::scale_u0/(ModelParam::scale_r0*ModelParam::scale_r0);
			dsmul(ModelParam::oneD_q,tmp_coef,temp_visc,1,Kr,1);
		}
		else
			dsmul(ModelParam::oneD_q,-2.0*ModelParam::oneD_alpha*ModelParam::PI/ModelParam::rho/(ModelParam::oneD_alpha-1.0),temp_visc,1,Kr,1);

		for(i = 0; i < nel; ++i){
			// q = U[i].q;
			dvdiv(ModelParam::oneD_q,U[i].h,1,A[i].h,1,tmp,1); // tmp = u/A
			// vdDiv(ModelParam::oneD_q,U[i].h,A[i].h,tmp);
			dvvtvp(ModelParam::oneD_q,tmp,1,Kr,1,Uf[0][i].h,1,Uf[0][i].h,1);
		}
	}
}

void OneDSolver::Bifur(Domain *omega){
	int    d,nel_p,nel_l,d1,d_r,nel_d1,d2,nel_d2,temp_case=1,d_temp;
	double A_temp, u_temp, W_temp;
	double W[3],uu[3],Au[3],Hdu[3];
	double A_rate, Q_rate, K, K31, K32, K13, K23, energy_in, energy_out; // Energy losses variables
	double angle=45*ModelParam::PI/180;

	for(d = 0; d < ModelParam::Ndoms; ++d){
		// Splitting flow configuration.
		if(omega[d].bctype[3] == 'B'){ // Check end of domains to determine if it is a 'B' bc case (splitting flow bifurcation).
			// Parent vessel
			nel_p  = omega[d].nel;
			uu[0]  = omega[d].U[nel_p-1].h [omega[d].U[nel_p-1].q-1];
			Au[0]  = omega[d].A[nel_p-1].h [omega[d].A[nel_p-1].q-1];
			W[0]   = uu[0] + omega[d].A[nel_p-1].Get_W(omega[d].A[nel_p-1].q-1,Au[0]);

			// Daughter vessel 1
			d1     = omega[d].bifur[1][0];
			uu[1]  = omega[d1].U[0].h [0];
			Au[1]  = omega[d1].A[0].h [0];
			W[1]   = uu[1] - omega[d1].A[0].Get_W(0,Au[1]);

			// Daughter vessel 2
			d2     = omega[d].bifur[1][1];
			uu[2]  = omega[d2].U[0].h [0];
			Au[2]  = omega[d2].A[0].h [0];
			W[2]   = uu[2] - omega[d2].A[0].Get_W(0,Au[2]);

			if(ModelParam::bifurloss==0){
				// Riemann problem without considering energy losses (Bifurc_losses = 0).
				BifurRiem1to2(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,ModelParam::rho);
			}
			else{
				// Calculate the loss coefficients (Gardel, 1957 and Levin, 1958).
				/* First we solve the Riemann problem without considering energy losses to determine the main
				daughter vessel and to decide the flow directions. */
				BifurRiem1to2(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,ModelParam::rho);

				/* Decide which one is the main daughter vessel an assign it the domain location d1. */
				if (Au[2] > Au[1]){
					d_temp = d2;
					d2 = d1;
					d1 = d_temp;
					A_temp = Au[2];
					Au[2] = Au[1];
					Au[1] = A_temp;
					u_temp = uu[2];
					uu[2] = uu[1];
					uu[1] = u_temp;
					W_temp = W[2];
					W[2] = W[1];
					W[1] = W_temp;
					temp_case = 2;  // To undo the store location of the daughter vessels.
				}

				/* If the parent velocity is positive, 3 flow configurations are possible + any flow is leaving the
				configuration which is against the mass conservation statement. */
				if (uu[0] > 0.0) {
					if (uu[1] < 0.0 && uu[2] < 0.0){ // Any flow is leaving the bifurcation
						fprintf(stderr,"Error boundary condition 'B': no flow is leaving the bifurcation.\n");
						exit(1);
					}
					if (uu[1] > 0.0 && uu[2] > 0.0){ // Dividing flow to the right
						nbif_div_right = nbif_div_right + 1;
						if (fabs(uu[0]) < 1e-10){ /* If the velocity in the parent vessel is close to zero,
												  no energy losses are considered when solving the Riemann problem. */
							//fprintf(stdout,"No energy losses considered because u in the combined vessel is close to zero.\n");
							nbif_no_losses = nbif_no_losses + 1;
						}
						else{
							// Calculate area and flow rates.
							A_rate = Au[2]/Au[0];
							Q_rate = A_rate*uu[2]/uu[0];

							// Side leg
							K31 = 0.95*(1.0 - Q_rate)*(1.0 - Q_rate)
								+ Q_rate*Q_rate*(1.3/tan(0.5*(ModelParam::PI - angle)) - 0.3 + (0.4 - 0.1*A_rate)/A_rate/A_rate)
								+ 0.4*Q_rate*(1.0 - Q_rate)*(1.0 + 1.0/A_rate)/tan(0.5*(ModelParam::PI - angle));
							// Through leg
							K32 = 0.03*(1.0 - Q_rate)*(1.0 - Q_rate) + 0.35*Q_rate*Q_rate - 0.2*Q_rate*(1.0 - Q_rate);

							/* Solution of the Riemann problem (LU factorisation of the Jacobian matrix) considering
							energy losses K31 & K32. */
							BifurRiem1to2_losses_div_right(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,ModelParam::rho,K31,K32);
							//fprintf(stdout,"Dividing flow to the right.\n");
							//fprintf(stdout,"Loss coefficients K31: %lg and K32: %lg\n",K31,K32);
						}
					}
					else{
						if (uu[1] > 0.0 && uu[2] < 0.0){ // Combining flow to the right
							nbif_conv_right = nbif_conv_right + 1;
							if (fabs(uu[1]) < 1e-10){ /* If the velocity in the main daughter vessel is close to zero,
													  no energy losses are considered when solving the Riemann problem. */
								//BifurRiem1to2(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,rho);
								//fprintf(stdout,"No energy losses considered because u in the combined vessel is close to zero.\n");
								nbif_no_losses = nbif_no_losses + 1;
							}
							else{
								// Calculate area and flow rates.
								A_rate = Au[2]/Au[1];
								Q_rate = -A_rate*uu[2]/uu[1];

								K13 = - 0.92*(1.0 - Q_rate)*(1.0 - Q_rate)
									- Q_rate*Q_rate*(1.2*(cos(ModelParam::PI-angle)/A_rate - 1.0) + 0.8*(1.0 - 1.0/A_rate/A_rate)
									- (1.0 - A_rate)*cos(ModelParam::PI-angle)/A_rate) + (2.0 - A_rate)*Q_rate*(1.0 - Q_rate);
								// Parent leg
								K23 = 0.03*(1.0 - Q_rate)*(1.0 - Q_rate)
									- Q_rate*Q_rate*(1.0 + 1.62*(cos(ModelParam::PI-angle)/A_rate - 1.0) - 0.38*(1.0 - A_rate))
									+ (2.0 - A_rate)*Q_rate*(1.0 - Q_rate);

								/* Solution of the Riemann problem (LU factorisation of the Jacobian matrix) considering
								energy losses K13 & K23. */
								BifurRiem1to2_losses_comb_right(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,ModelParam::rho,K13,K23);
								//fprintf(stdout,"Combining flow to the right.\n");
								//fprintf(stdout,"Loss coefficients K13: %lg and k23: %lg\n",K13,K23);
							}
						}
						else{ // Flow leaving through the side branch only. The loss coefficients are only vaalid for angle=90, but we use them for any angle.
							nbif_leav = nbif_leav + 1;
							if (fabs(uu[2]) < 1e-10){ /* If the velocity in the smallest daughter vessel is close to zero,
													  no energy losses are considered when solving the Riemann problem. */
								//BifurRiem1to2(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,rho);
								//fprintf(stdout,"No energy losses considered because u in the combined vessel is close to zero.\n");
								nbif_no_losses = nbif_no_losses + 1;
							}
							else{
								// Calculate area and flow rates.
								A_rate = Au[0]/Au[2];
								Q_rate = A_rate*uu[0]/uu[2];

								// Parent leg
								K13 = 1.0 + 1.0/A_rate/A_rate + 3.0*Q_rate*(Q_rate - 1.0)/A_rate/A_rate;
								// Through leg
								K23 = K13;

								/* Solution of the Riemann problem (LU factorisation of the Jacobian matrix)
								considering energy losses K13 & K23. */
								BifurRiem1to2_losses_leav(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,ModelParam::rho,K13,K23);
								//fprintf(stdout,"Flow leaving through the side branch only.\n");
								//fprintf(stdout,"Loss coefficients K13: %lg and K23: %lg\n",K13,K23);
							}
						}
					}
				}
				/* If the parent velocity is negative, other 3 flow conditions are possible + all the flows leaving the
				configuration which is against the mass conservation statement. */
				else{
					if (uu[1] > 0.0 && uu[2] > 0.0){ // All the flows are leaving the bifurcation
						fprintf(stderr,"Error boundary condition 'B': all flows are leaving the bifurcation.\n");
						exit(1);
					}
					if (uu[1] < 0.0 && uu[2] > 0.0){ // Dividing flow to the left
						nbif_div_left = nbif_div_left + 1;
						if (fabs(uu[1]) < 1e-10){ /* If the velocity in the main daughter vessel is close to zero,
												  no energy losses are considered when solving the Riemann problem. */
							//BifurRiem1to2(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,rho);
							//fprintf(stdout,"No energy losses considered because u in the combined vessel is close to zero.\n");
							nbif_no_losses = nbif_no_losses + 1;
						}
						else{
							// Calculate area and flow rates.
							A_rate = Au[2]/Au[1];
							Q_rate = -A_rate*uu[2]/uu[1];

							// Side leg
							K31 = 0.95*(1.0 - Q_rate)*(1.0 - Q_rate)
								+ Q_rate*Q_rate*(1.3/tan(0.5*(angle)) - 0.3 + (0.4 - 0.1*A_rate)/A_rate/A_rate)
								+ 0.4*Q_rate*(1.0 - Q_rate)*(1.0 + 1.0/A_rate)/tan(0.5*(angle));
							// Parent leg
							K32 = 0.03*(1.0 - Q_rate)*(1.0 - Q_rate) + 0.35*Q_rate*Q_rate - 0.2*Q_rate*(1.0 - Q_rate);

							/* Solution of the Riemann problem (LU factorisation of the Jacobian matrix)
							considering enrgy losses K31 & K32. */
							BifurRiem1to2_losses_div_left(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,ModelParam::rho,K31,K32);
							//fprintf(stdout,"Dividing flow to the left.\n");
							//fprintf(stdout,"Loss coefficients K31: %lg and K32: %lg\n",K31,K32);
						}
					}
					else{
						if (uu[1] < 0.0 && uu[2] < 0.0){ // Combining flow to the left
							nbif_conv_left = nbif_conv_left + 1;
							if (fabs(uu[0]) < 1e-10){ /* If the velocity in the parent vessel is close to zero,
													  no energy losses are considered when solving the Riemann problem. */
								//BifurRiem1to2(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,rho);
								//fprintf(stdout,"No energy losses considered because u in the combined vessel is close to zero.\n");
								nbif_no_losses = nbif_no_losses + 1;
							}
							else{
								// Calculate area and flow rates.
								A_rate = Au[2]/Au[0];
								Q_rate = A_rate*uu[2]/uu[0];

								// Side leg
								K13 = - 0.92*(1.0 - Q_rate)*(1.0 - Q_rate)
									- Q_rate*Q_rate*(1.2*(cos(angle)/A_rate - 1.0) + 0.8*(1.0 - 1.0/A_rate/A_rate)
									- (1.0 - A_rate)*cos(angle)/A_rate) + (2.0 - A_rate)*Q_rate*(1.0 - Q_rate);
								// Through leg
								K23 = 0.03*(1.0 - Q_rate)*(1.0 - Q_rate)
									- Q_rate*Q_rate*(1.0 + 1.62*(cos(angle)/A_rate - 1.0) - 0.38*(1.0 - A_rate))
									+ (2.0 - A_rate)*Q_rate*(1.0 - Q_rate);

								/* Solution of the Riemann problem (LU factorisation of the Jacobian matrix)
								considering enrgy losses K13 & K23. */
								BifurRiem1to2_losses_comb_left(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,ModelParam::rho,K13,K23);
								//fprintf(stdout,"Combining flow to the left.\n");
								//fprintf(stdout,"Loss coefficients K13: %lg and K23: %lg\n",K13,K23);
							}
						}
						else{ // Flow entering through the side branch only. The loss coefficients are only vaalid for angle=90, but we use them for any angle.
							nbif_ent = nbif_ent + 1;
							if (fabs(uu[2]) < 1e-10){ /* If the velocity in the smallest daughter vessel is close to zero,
													  no energy losses are considered when solving the Riemann problem. */
								//BifurRiem1to2(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,rho);
								//fprintf(stdout,"No energy losses considered because u in the combined vessel is close to zero.\n");
								nbif_no_losses = nbif_no_losses + 1;
							}
							else{
								// Calculate area and flow rates.
								A_rate = Au[1]/Au[2];
								Q_rate = -A_rate*uu[1]/uu[2];

								// Through leg
								K31 = 1.0 + 0.3*Q_rate*Q_rate/A_rate/A_rate;
								// Parent leg
								K32 = 1.0 + 0.3*(1.0 - Q_rate)*(1.0 - Q_rate)/A_rate/A_rate;

								/* Solution of the Riemann problem (LU factorisation of the Jacobian matrix)
								considering enrgy losses K31 & K32. */
								BifurRiem1to2_losses_ent(nel_p,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,ModelParam::rho,K31,K32);
								//fprintf(stdout,"Flow entering through the side branch only.\n");
								//fprintf(stdout,"Loss coefficients K31: %lg and K32: %lg\n",K31,K32);
							}
						}
					}
				}

				/* If "temp_case == 2" we have to undo the store location of the daughter vessels. */
				if (temp_case == 2){
					d_temp = d2;
					d2 = d1;
					d1 = d_temp;
					A_temp = Au[2];
					Au[2] = Au[1];
					Au[1] = A_temp;
					u_temp = uu[2];
					uu[2] = uu[1];
					uu[1] = u_temp;
					W_temp = W[2];
					W[2] = W[1];
					W[1] = W_temp;
				}
			}

#ifdef DEBUG_1D
			// Check energy balance (energy per unit time IN > energy per unit time OUT)
			energy_in = 0.0;
			energy_out = 0.0;
			if(uu[0] > 0.0) energy_in = uu[0]*Au[0]*(0.5*ModelParam::rho*uu[0]*uu[0] + omega[d].A[nel_p-1].Get_P(omega[d].A[nel_p-1].q-1,Au[0]));
			else energy_out = fabs(uu[0])*Au[0]*(0.5*ModelParam::rho*uu[0]*uu[0] + omega[d].A[nel_p-1].Get_P(omega[d].A[nel_p-1].q-1,Au[0]));

			if(uu[1] > 0.0) energy_out = energy_out + uu[1]*Au[1]*(0.5*ModelParam::rho*uu[1]*uu[1] + omega[d1].A[0].Get_P(0,Au[1]));
			else energy_in = energy_in + fabs(uu[1])*Au[1]*(0.5*ModelParam::rho*uu[1]*uu[1] + omega[d1].A[0].Get_P(0,Au[1]));

			if(uu[2] > 0.0) energy_out = energy_out + uu[2]*Au[2]*(0.5*ModelParam::rho*uu[2]*uu[2] + omega[d2].A[0].Get_P(0,Au[2]));
			else energy_in = energy_in + fabs(uu[2])*Au[2]*(0.5*ModelParam::rho*uu[2]*uu[2] + omega[d2].A[0].Get_P(0,Au[2]));

			if(pow((energy_out - energy_in),2.0) > ModelParam::riemannTol){
				fprintf(stderr,"Error in boundary condition 'B': energy per unit time OUT > energy per unit time IN at bifurcacion with parent vessel %d. \n",d+1);
				exit(1);
			}
#endif

			/* The solution is assigned to the value of A & U enforced as bc's in the bifurcation sides of
			domains "d", "d1" and "d2". */
			omega[d].bcval[3] = Au[0];
			omega[d].bcval[4] = uu[0];

			omega[d1].bcval[0] = Au[1];
			omega[d1].bcval[1] = uu[1];

			omega[d2].bcval[0] = Au[2];
			omega[d2].bcval[1] = uu[2];

			//// Bifur of Hd, added by panqing
			//Hdu[0] = omega[d].Hd[nel_p-1].h[omega[d].Hd[nel_p-1].q-1];
			//Hdu[1] = omega[d1].Hd[0].h[0];
			//Hdu[2] = omega[d2].Hd[0].h[0];
			//
			//BifurHd1to2(uu,Au,Hdu);

			///* The solution is assigned to the value of Hd enforced as bc's in the bifurcation sides of
			// domains "d", "d1" and "d2". */
			//omega[d ].bcval[5] = Hdu[0];
			//omega[d1].bcval[2] = Hdu[1];
			//omega[d2].bcval[2] = Hdu[2];

			///* For fixed Hd mode, the Hd value could be set to all the quadrature points of the downstream vessels. */
			//dfill(omega[d1].Hd[0].q, Hdu[1], omega[d1].Hd[0].h, 1);
			//dfill(omega[d2].Hd[0].q, Hdu[2], omega[d2].Hd[0].h, 1);
		}

		// Merging flow configuration.
		if(omega[d].bctype[0] == 'C'){ // Check end of domains to determine if it is a 'C' bc case (merging flow bifurcation).
			// Parent vessel (Merging vessel)
			uu[0]   = omega[d].U[0].h [0];
			Au[0]   = omega[d].A[0].h [0];
			W[0]    = uu[0] - omega[d].A[0].Get_W(0,Au[0]);

			// Daughter vessel 1
			d1      = omega[d].bifur[0][0];
			nel_d1  = omega[d1].nel;
			uu[1]   = omega[d1].U[nel_d1-1].h [omega[d1].U[nel_d1-1].q-1];
			Au[1]   = omega[d1].A[nel_d1-1].h [omega[d1].A[nel_d1-1].q-1];
			W[1]    = uu[1] + omega[d1].A[nel_d1-1].Get_W(omega[d1].A[nel_d1-1].q-1,Au[1]);

			// Daughter vessel 2
			d2      = omega[d].bifur[0][1];
			nel_d2  = omega[d2].nel;
			uu[2]   = omega[d2].U[nel_d2-1].h [omega[d2].U[nel_d2-1].q-1];
			Au[2]   = omega[d2].A[nel_d2-1].h [omega[d2].A[nel_d2-1].q-1];
			W[2]    = uu[2] + omega[d2].A[nel_d2-1].Get_W(omega[d2].A[nel_d2-1].q-1,Au[2]);

			// Call the nonlinear bifurcation Riemann problem solver.
			BifurRiem2to1(nel_d1,nel_d2,omega[d].A,omega[d1].A,omega[d2].A,W,uu,Au,ModelParam::rho);

			/* The solution satisfies conservation of mass and continuity of the total pressure. It is assigned
			to the value of A & u enforced as bc's in the bifurcation sides of domains "d", "d1" and "d2". */
			omega[d].bcval[0] = Au[0];
			omega[d].bcval[1] = uu[0];

			omega[d1].bcval[3] = Au[1];
			omega[d1].bcval[4] = uu[1];

			omega[d2].bcval[3] = Au[2];
			omega[d2].bcval[4] = uu[2];

			//		// Converge of Hd, added by panqing
			//		Hdu[0]  = omega[d].Hd[0].h[0];
			//		Hdu[1]  = omega[d1].Hd[nel_d1-1].h[omega[d1].Hd[nel_d1-1].q-1];
			//		Hdu[2]  = omega[d2].Hd[nel_d2-1].h[omega[d2].Hd[nel_d2-1].q-1];

			//		// Call the nonlinear bifurcation Riemann problem solver.
			//		BifurHd2to1(uu,Au,Hdu);
			//		/* The solution satisfies conservation of mass and continuity of the total pressure. It is assigned
			//to the value of A & u enforced as bc's in the bifurcation sides of domains "d", "d1" and "d2". */
			//		omega[d ].bcval[5] = Hdu[0];
			//		omega[d1].bcval[2] = Hdu[1];
			//		omega[d2].bcval[2] = Hdu[2];

			//		/* For fixed Hd mode, the Hd value could be set to all the quadrature points of the downstream vessels. */
			//		dfill(omega[d].Hd[0].q, Hdu[0], omega[d].Hd[0].h, 1);
		}
		// Union configuration
		if(omega[d].bctype[3] == 'J'){ // Check end of domains to determine if it is a 'J' bc case (Union of two vessels).
			// Left vessel
			nel_l  = omega[d].nel;
			uu[0]  = omega[d].U[nel_l-1].h [omega[d].U[nel_l-1].q-1];
			Au[0]  = omega[d].A[nel_l-1].h [omega[d].A[nel_l-1].q-1];
			W[0]   = uu[0] + omega[d].A[nel_l-1].Get_W(omega[d].A[nel_l-1].q-1,Au[0]);

			// Right vessel
			d_r    = omega[d].bifur[1][0];
			uu[1]  = omega[d_r].U[0].h [0];
			Au[1]  = omega[d_r].A[0].h [0];
			W[1]   = uu[1] - omega[d_r].A[0].Get_W(0,Au[1]);

			if(ModelParam::bifurloss==0){
				/* Call the nonlinear junction Riemann problem solver */
				// Riemann problem without considering energy losses.
				JuncRiemann(nel_l,omega[d].A,omega[d_r].A,W,uu,Au,ModelParam::rho);
			}
			else{
				/* First we solve the Riemann problem without considering energy losses
				to determine the flow direction. */
				JuncRiemann(nel_l,omega[d].A,omega[d_r].A,W,uu,Au,ModelParam::rho);

				/* If the flow is from the to right, we can have a sudden expansion or a sudden
				contraction. */
				if (uu[0] > 0.0) {
					if (Au[0] < Au[1]) { // Sudden expansion to the right
						// Calculate area rate and the loss coefficient
						A_rate = Au[0]/Au[1];
						K = (1.0 - A_rate)*(1.0 - A_rate);

						// Solution of the Riemann problem (LU factorisation of the Jacobian matrix).
						JuncRiemann_losses_exp_right(nel_l,omega[d].A,omega[d_r].A,W,uu,Au,ModelParam::rho,K);
						//fprintf(stdout,"Sudden expansion to the right.\n");
						//fprintf(stdout,"Loss coefficient: %lg\n",K);
					}
					else { // Sudden contraction to the right
						// Calculate area rate
						A_rate = Au[1]/Au[0];
						// Loss coefficient depends upon area rate
						if (A_rate > 0.5776) K = (1.0 - A_rate)*(1.0 - A_rate);
						else K = 0.42*(1.0 - A_rate);

						// Solution of the Riemann problem (LU factorisation of the Jacobian matrix).
						JuncRiemann_losses_cont_right(nel_l,omega[d].A,omega[d_r].A,W,uu,Au,ModelParam::rho,K);
						//fprintf(stdout,"Sudden contraction to the right.\n");
						//fprintf(stdout,"Loss coefficient: %lg\n",K);
					}
				}
				else {   /* If the flow is from right to left, we can have a sudden expansion or a sudden contraction. */
					if (Au[0] > Au[1]) { // Sudden expansion to the left
						// Calculate area rate and the loss coefficient
						A_rate = Au[1]/Au[0];
						K = (1.0 - A_rate)*(1.0 - A_rate);

						// Solution of the Riemann problem (LU factorisation of the Jacobian matrix).
						JuncRiemann_losses_exp_left(nel_l,omega[d].A,omega[d_r].A,W,uu,Au,ModelParam::rho,K);
						//fprintf(stdout,"Sudden expansion to the left.\n");
						//fprintf(stdout,"Loss coefficient: %lg\n",K);
					}
					else { // Sudden contraction to the left
						// Calculate area rate
						A_rate = Au[0]/Au[1];

						// Loss coefficient depends upon area rate
						if (A_rate > 0.5776) K = (1.0 - A_rate)*(1.0 - A_rate);
						else K = 0.42*(1.0 - A_rate);

						// Solution of the Riemann problem (LU factorisation of the Jacobian matrix).
						JuncRiemann_losses_cont_left(nel_l,omega[d].A,omega[d_r].A,W,uu,Au,ModelParam::rho,K);
						//fprintf(stdout,"Sudden contraction to the left.\n");
						//fprintf(stdout,"Loss coefficient: %lg\n",K);
					}
				}
			}

#ifdef DEBUG_1D
			// Check flow rate balances
			if(pow((Au[0]*uu[0] - Au[1]*uu[1]),2.0) > ModelParam::riemannTol){
				fprintf(stderr,"Error in boundary condition 'J': Q does not balance at the junction with left element # %d.\n",d+1);
				exit(1);
			}
#endif
			/* The solution is assigned to the value of A & u enforced as bc's in the junction sides of
			domains "d" and "dr". */
			omega[d].bcval[3]  = Au[0];
			omega[d].bcval[4]  = uu[0];

			omega[d_r].bcval[0] = Au[1];
			omega[d_r].bcval[1] = uu[1];

			//// Junction for Hd, added by panqing
			//Hdu[0] = omega[d].Hd[nel_l-1].h [omega[d].Hd[nel_l-1].q-1];
			//Hdu[1] = omega[d_r].Hd[0].h[0];

			//JuncHd(uu,Au,Hdu);

			///* The solution is assigned to the value of A & u enforced as bc's in the junction sides of
			// domains "d" and "dr". */
			//omega[d  ].bcval[5] = Hdu[0];
			//omega[d_r].bcval[2] = Hdu[1];

			///* For fixed Hd mode, the Hd value could be set to all the quadrature points of the downstream vessels. */
			//dfill(omega[d_r].Hd[0].q, Hdu[1], omega[d_r].Hd[0].h, 1);
		}
	}
}

void OneDSolver::BioFlux(Domain *omega, char *name){
	int i, j, n;
	int nel;
	static double cl, cr;
	static double **uloc, **af, **aloc, **uf;
	Element *U, *A, *Uf, *Af; // Elements with the solution and the fluxes
	static double *bc[6]; // Value of A & u at each side of the domain
	static FILE  **fIN, **fOUT; // Points to the inlet and outlet of a domain at the boundary of the network
	static double co, Ao, W, beta;
	static double pinf = 0.0;
	double *gamma;

	if(!fIN){ // The first time the routine is called,
		fIN = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Allocate space for "fIN" to read from a .bcs file
	}
	if(!fOUT){ // The first time the routine is called,
		fOUT = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Allocate space for "fOUT" to read from a .bcs file
	}
	if(!uff){
		uff = (double *)malloc(ModelParam::Ndoms*sizeof(double));
	}
	if(!ubb){
		ubb = (double *)malloc(ModelParam::Ndoms*sizeof(double));
	}
	if(!pff){
		pff = (double *)malloc(ModelParam::Ndoms*sizeof(double));
	}
	if(!pbb){
		pbb = (double *)malloc(ModelParam::Ndoms*sizeof(double));
	}
	if(!q_in){
		q_in = (double *)malloc(ModelParam::Ndoms*sizeof(double));
	}
	if(!u_star){
		u_star = (double *)malloc(ModelParam::Ndoms*sizeof(double));
	}
	if(!A_star){
		A_star = (double *)malloc(ModelParam::Ndoms*sizeof(double));
	}

	for(n = 0; n < ModelParam::Ndoms; ++n){ // For each domain "n" (n=0,...,Ndoms),
		nel = omega[n].nel;   // Number of elements in the domain "n".
		U   = omega[n].U;
		A   = omega[n].A;
		Uf  = omega[n].Uf[0];
		Af  = omega[n].Af[0];
		gamma = omega[n].gamma;

		// Generate the local solution vectors at the ends of the elements of the domain "n".
		if(!uloc){  // The first time "BioFlux" is called (uloc not defined yet),
			aloc = dmatrix(-1,nel,0,1);  // Local area
			uloc = dmatrix(-1,nel,0,1);  // Local velocity
			af   = dmatrix(-1,nel,0,1);  // Local and upwind mass fluxs
			uf   = dmatrix(-1,nel,0,1);  // Local and upwind momentum fluxs
		}

		for(i = 0; i < nel; ++i){ /*  For each element "i" (i=0,...,nel) of the domain "n" (n=0,...,Ndoms),
								  the local area and the local velocity at the end of the element are defined. */
			aloc[i][0] = A[i].h[0];
			aloc[i][1] = A[i].h[A[i].q-1];
			uloc[i][0] = U[i].h[0];
			uloc[i][1] = U[i].h[U[i].q-1];
		}

		/* Set up boundary conditions at each side of the domain "n".
		The value of A and U on the l.h.s are stored in bc[0] and bc[1], respectively, and
		the value of A and U on the r.h.s are stored in bc[2] and bc[3], respectively. */

		/* Link the value of the A & U bc's at each side of the domain "n" with the local area and the
		local velocity at the end of the virtual element on the left of the first element and at the
		beggining of the virtual element on the right of the last element*/
		bc[0] = aloc[-1]+1;
		bc[1] = uloc[-1]+1;
		bc[3] = aloc[nel];
		bc[4] = uloc[nel];

		for(i = 0; i < 6; ++i){
			if(i != 2 && i != 5){
				switch(omega[n].bctype[i]){ // Find the A & U bc types at both sides of the domain "n".
				  case 'D': case 'B': case 'C': case 'J': case 'I':/** Dirichlet bc on upwind state. It takes into account the cases:
																   'D': Constant Dirichlet bc on upwind state.
																   'B': Bc from the solution of a Riemann problem at a
																   splitting flow bifurcation.
																   'C': Bc from the solution of a Riemann problem at a
																   merging flow bifurcation.
																   'J': Bc from the solution of a Riemann problem at a
																   union of two vessels.
																   'I': Bc from the solution of a Riemann problem at a
																   union of two vessels with a loss in energy to 
																   simulate an ischemic attack.
																   */
					  bc[i][0] = omega[n].bcval[i]; /** The bc is the value read in the input file for the case 'D' or
													the solution of the Riemann problem at the bifurcation for the
													cases 'B', 'C', 'J', and 'I'.
													*/
					  break;
				  case 'P': // Exact constant Dirichlet bc on p.
					  switch(i){
						  case 0:
							  double co; // Initial wave speed at the inlet
							  bc[i][0] = eval_A_from_P(omega[n].bcval[i],0.0,A[0].beta[0],A[0].Ao[0]); // Calculate A from p
							  bc[i][0] = eval_BC_A(bc[i][0],A[0].h[0]); // We enforce A_l that combined with A_r yields the desired area BC.
							  break;
						  case 1:
							  bc[i][0] = U[0].h[0];
							  break;
						  case 3:
							  bc[i][0] = eval_A_from_P(omega[n].bcval[i],0.0,A[nel-1].beta[A[nel-1].q-1],A[nel-1].Ao[A[nel-1].q-1]); // Calculate A from p
							  bc[i][0] = eval_BC_A(bc[i][0],A[nel-1].h[A[nel-1].q-1]); // We enforce A_l that combined with A_r yields the desired area BC.
							  break;
						  case 4:
							  bc[i][0] = U[nel-1].h[U[nel-1].q-1];
							  break;
					  }
					  break;
				  case 'p': // Exact Dirichlet bc on p.
					  switch(i){
						  case 0: // We let Ul = Ur and we impose bc through A evaluated from p(t) in the input file (ModelParam::rho, material properties = const).
							  if( (omega[n].bcval[i] == 2) || (omega[n].bcval[i] == 3)){ // The input data is read from a .bcs file				
								  if(!fIN[n]){ // Define the name of the .bcs file for domain n
									  if(ModelParam::Ndoms == 1) sprintf(buf,"%s_IN.bcs",strtok(name,"."));
									  else sprintf(buf,"%s_IN_%d.bcs",strtok(name,"."),n+1);
									  if(!(fIN[n]=fopen(buf,"rb"))){   /* If the input file .bcs is not found, print an error message and
																	   terminate executing the program. If found, open it to start reading it. */
										  fprintf(stderr,"Error in BCs: the file %s doesn't exist. \n",buf);
										  exit(-1);  /* Standard library function that terminates program execution.
													 It also calls "fclose" for each open output file in order to flush out any buffered output. */
									  }
								  }

								  if(!reIntFlag)
									  fread(&(pff[n]), sizeof(double), 1, fIN[n]);

								  bc[i][0] = eval_A_from_P(pff[n],0.0,A[0].beta[0],A[0].Ao[0]); // Calculate A from p
								  if(omega[n].bcval[i] == 3) bc[i][0] = eval_BC_A(bc[i][0],A[0].Ao[0]);  // We enforce an absorbing area BC.
								  else bc[i][0] = eval_BC_A(bc[i][0],A[0].h[0]);  // We enforce Al that combined with Ar yields the desired area BC.
							  }
							  else{ 
								  // The input pressure is read from the input file .in
							  }		
							  break;
						  case 1:
							  if((omega[n].bcval[i] == 1) || (omega[n].bcval[i] == 3)) bc[i][0] = 0.0;
							  else bc[i][0] = U[0].h[0];
							  break;
						  case 3: // We let Ul = Ur and we impose bc through A evaluated from p(t) in the input file (ModelParam::rho, material properties = const).
							  if( (omega[n].bcval[i] == 2) || (omega[n].bcval[i] == 3)){ // The input data is read from a .bcs file
								  if(!fOUT[n]){ // Define the name of the .bcs file for domain n
									  if(ModelParam::Ndoms == 1) sprintf(buf,"%s_OUT.bcs",strtok(name,"."));
									  else sprintf(buf,"%s_OUT_%d.bcs",strtok(name,"."),n+1);
									  if(!(fOUT[n]=fopen(buf,"rb"))){   /* If the input file .bcs is not found, print an error message and
																		terminate executing the program. If found, open it to start reading it. */
										  fprintf(stderr,"Error in BCs: the file %s doesn't exist. \n",buf);
										  exit(-1);  /* Standard library function that terminates program execution.
													 It also calls "fclose" for each open output file in order to flush out any buffered output. */
									  }					
								  }

								  if(!reIntFlag)
									  fread(&(pbb[n]), sizeof(double), 1, fIN[n]);

								  bc[i][0] = eval_A_from_P(pbb[n],0.0,A[nel-1].beta[A[nel-1].q-1],A[nel-1].Ao[A[nel-1].q-1]); // Calculate A from p
								  if(omega[n].bcval[i] == 3) bc[i][0] = eval_BC_A(bc[i][0],A[nel-1].Ao[A[nel-1].q-1]);  // We enforce an absorbing area BC.
								  else bc[i][0] = eval_BC_A(bc[i][0],A[nel-1].h[A[nel-1].q-1]);  // We enforce Ar that combined with Al yields the desired area BC.					
							  }		
							  else{
								  // The input is read from the input file .in	
							  }
							  break;
						  case 4:
							  if((omega[n].bcval[i] == 1) || (omega[n].bcval[i] == 3)) bc[i][0] = 0.0;
							  else bc[i][0] = U[nel-1].h[U[nel-1].q-1];
							  break;
				  }
				  break;
			  case 'u': // Exact Dirichlet bc on U.
				  switch(i){
					  case 0: // We let Al = Ar and we impose bc through U (ModelParam::rho, material properties = const).
						  if((omega[n].bcval[i] == 1) || (omega[n].bcval[i] == 3)) bc[i][0] = A[0].Ao[0];
						  else bc[i][0] = A[0].h[0];
						  break;
					  case 1:
						  if( (omega[n].bcval[i] == 2) || (omega[n].bcval[i] == 3)){ // The input data is read from a .bcs file				
							  if(!fIN[n]){ // Define the name of the .bcs file for domain n
								  if(ModelParam::Ndoms == 1) sprintf(buf,"%s_IN.bcs",strtok(name,"."));
								  else sprintf(buf,"%s_IN_%d.bcs",strtok(name,"."),n+1);
								  if(!(fIN[n]=fopen(buf,"rb"))){   /* If the input file .bcs is not found, print an error message and
																   terminate executing the program. If found, open it to start reading it. */
									  fprintf(stderr,"Error in BCs: the file %s doesn't exist. \n",buf);
									  exit(-1);  /* Standard library function that terminates program execution.
												 It also calls "fclose" for each open output file in order to flush out any buffered output. */
								  }
								  else{
									  // printf("%s opened\n", buf);
								  }	
							  }

							  if(!reIntFlag)
								  fread(&(uff[n]), sizeof(double), 1, fIN[n]);

							  if(omega[n].bcval[i] == 3) bc[i][0] = 2*uff[n];  // We enforce an absorbing velocity BC.
							  else bc[i][0] = 2*uff[n] - U[0].h[0];  // We enforce Ul that combined with Ur yields the desired velocity BC.
						  }		
						  else{
							  if(omega[n].bcval[i] == 1)	bc[i][0] = 2*bc[i][0]; // We enforce an absorbing velocity BC.
							  else bc[i][0] = 2*bc[i][0] - U[0].h[0]; // We enforce Ul that combined with Ur yields the desired velocity BC.
						  }
						  break;
					  case 3:
						  if((omega[n].bcval[i] == 1) || (omega[n].bcval[i] == 3)) bc[i][0] = A[nel-1].Ao[A[nel-1].q-1];
						  else bc[i][0] = A[nel-1].h[A[nel-1].q-1];
						  break;
					  case 4:
						  if( (omega[n].bcval[i] == 2) || (omega[n].bcval[i] == 3)){ // The input data is read from a .bcs file				
							  if(!fOUT[n]){ // Define the name of the .bcs file for domain n
								  if(ModelParam::Ndoms == 1) sprintf(buf,"%s_OUT.bcs",strtok(name,"."));
								  else sprintf(buf,"%s_OUT_%d.bcs",strtok(name,"."),n+1);
								  if(!(fOUT[n]=fopen(buf,"rb"))){   /* If the input file .bcs is not found, print an error message and
																	terminate executing the program. If found, open it to start reading it. */
									  fprintf(stderr,"Error in BCs: the file %s doesn't exist. \n",buf);
									  exit(-1);  /* Standard library function that terminates program execution.
												 It also calls "fclose" for each open output file in order to flush out any buffered output. */
								  }		
							  }

							  // fscanf(fOUT[n], "%lg %lg %lg %lg %lg %lg\n", &t_in, &pff, &pbb, &uff, &ubb, &Rt);
							  if(!reIntFlag){
								  fread(&(ubb[n]), sizeof(double), 1, fOUT[n]);
							  }
							  if(omega[n].bcval[i] == 3) bc[i][0] = 2*ubb[n];  // We enforce an absorbing velocity BC.
							  else bc[i][0] = 2*ubb[n] - U[nel-1].h[U[nel-1].q-1];  // We enforce Ur that combined with Ul yields the desired velocity BC.
						  }						
						  else{
							  // The input is read from the input file .in	
						  }				
						  break;
				  }
				  break;
			  case 'q': // Exact Dirichlet bc on Q.
				  switch(i){
					  case 0: // We let Al = Ar and we impose bc through U evaluated from Q(t) in the input file (ModelParam::rho, material properties = const).
						  bc[i][0] = A[0].h[0];
						  break;
					  case 1:
						  if( (omega[n].bcval[i] == 2) || (omega[n].bcval[i] == 3)){ // The input flow rate is read from a .bcs file
							  if(!fIN[n]){ // Define the name of the .bcs file for domain n
								  if(ModelParam::Ndoms == 1) sprintf(buf,"%s_IN.bcs",strtok(name,"."));
								  else sprintf(buf,"%s_IN_%d.bcs",strtok(name,"."),n+1);
								  if(!(fIN[n]=fopen(buf,"rb"))){   /* If the input file .bcs is not found, print an error message and
																   terminate executing the program. If found, open it to start reading it. */
									  fprintf(stderr,"Error in BCs: the file %s doesn't exist. \n",buf);
									  exit(-1);  /* Standard library function that terminates program execution.
												 It also calls "fclose" for each open output file in order to flush out any buffered output. */
								  }
							  }

							  if(!reIntFlag){
								  fread(&(q_in[n]), sizeof(double), 1, fIN[n]);
								  WKinflow(omega[n].A, U[0].h[0], A[0].h[0], q_in[n], &(u_star[n]), &(A_star[n]));
							  }

							  if(omega[n].bcval[i] == 3) bc[i][0] = 2.0*u_star[n]; // We enforce an absorbing velocity BC.
							  else bc[i][0] = 2.0*u_star[n] - U[0].h[0]; // We enforce Ul that combined with Ur yields the desired area BC.
						  }
						  else{ 
							  // The input is read from the input file .in	
						  }
						  break;
					  case 3:
						  bc[i][0] = A[nel-1].h[A[nel-1].q-1];
						  break;
					  case 4:
						  /*t = dparam("TIME");
						  vector_def("t",omega[n].bcstring[i]);
						  vector_set(1,&t,bc[i]);
						  bc[i][0]  = 2*bc[i][0]/A[nel-1].h[A[nel-1].q-1] -  U[nel-1].h[U[nel-1].q-1];*/
						  break;
				  }
				  break;
			  case 'T': // Terminal reflection coeficient bc
				  switch(i){
					  case 0: // We assume Al = Ar and we impose bc through U.
						  bc[i][0] = A[0].h[0];
						  //fprintf(stderr,"Terminal BC (T) not set up for l.h.s. \n");
						  break;
					  case 1:// We assume u0 = 0.
						  bc[i][0] = (1-omega[n].bcval[1])*(U[0].h[0]
						  - A[0].Get_W(0,A[0].h[0])); 
						  bc[i][0] = bc[i][0] - U[0].h[0];
						  break;
					  case 3: // We assume Ar = Al and we impose bc through U.
						  bc[i][0] = A[nel-1].h[A[nel-1].q-1];
						  break;
					  case 4:  // We assume u0 = 0.
						  bc[i][0] = (1-omega[n].bcval[4])*(U[nel-1].h[U[nel-1].q-1]
						  + A[nel-1].Get_W(A[nel-1].q-1,aloc[nel-1][1])); 
						  bc[i][0] = bc[i][0] - U[nel-1].h[U[nel-1].q-1];
						  //fprintf(fp[n],"%lg %lg %lg %lg\n", t, A[nel-1].Get_P(A[nel-1].q-1,bc[i-1][0]), bc[i][0], bc[i-1][0]);
						  break;
				  }
				  break;
			  case 't': // Terminal reflection coeficient bc enforced through R. Only for the r.h.s.
				  switch(i){
					  case 0:
						  bc[i][0] = 0;
						  fprintf(stderr,"BC type 't' not set up for l.h.s. \n");
						  break;
					  case 1:
						  bc[i][0] = 0;
						  fprintf(stderr,"BC type 't' not set up for l.h.s. \n");
						  break;
					  case 3: // We assume Ar = Al and we impose bc through U.
						  bc[i][0] = A[nel-1].h[A[nel-1].q-1];
						  break;
					  case 4:  // We assume u0 = 0.
						  double Rt, Wnl, R;

						  // We work out the terminal reflection coefficient from R and according to the linear theory
						  /* Linear wave speed and initial area*/
						  co = A[nel-1].Get_co(A[nel-1].q-1);
						  Ao = A[nel-1].Ao[A[nel-1].q-1];
						  beta = A[nel-1].beta[A[nel-1].q-1];

						  // Determine forward characteristic Wf from A,u on the left
						  W = uloc[nel-1][1] + co*aloc[nel-1][1]/Ao;

						  Wnl = U[nel-1].h[U[nel-1].q-1] + A[nel-1].Get_W(A[nel-1].q-1,aloc[nel-1][1])
							  - 4.0*sqrt(0.5*beta/ModelParam::rho*Ao)*sqrt(Ao);

						  // R from input file
						  R = omega[n].bcval[4];

						  // Terminal reflection coefficient
						  Rt = ( (R*Ao - ModelParam::rho*co)*0.5*W + pinf )/( (R*Ao + ModelParam::rho*co)*0.5*W);
						  //printf("Rt = %lf  %\lf \n", Rt, (R*Ao - rho*co)*0.5*W + pinf);
						  //printf("W_l = %lf    W_nl = %lf \n", W, co);


						  bc[i][0] = (1-Rt)*(U[nel-1].h[U[nel-1].q-1]
						  + A[nel-1].Get_W(A[nel-1].q-1,aloc[nel-1][1])); 

						  bc[i][0] = bc[i][0] - U[nel-1].h[U[nel-1].q-1];
						  break;
				  }
				  break;
			  case 'W': /** Windkessel 3-elements model (resistance + compliance + resistance) terminal bc.
						Only for the r.h.s. */
				  switch(i){
					  case 0:
						  bc[i][0]  = 0;
						  fprintf(stderr,"Terminal BC (W) not set up for l.h.s. \n");
						  break;
					  case 1:
						  bc[i][0]  = 0;
						  fprintf(stderr,"Terminal BC (W) not set up for l.h.s. \n");
						  break;
					  case 3: // We assume Al = Ar and we impose bc condition through U.
						  bc[i][0] =  A[nel-1].h[A[nel-1].q-1];
						  break;
					  case 4:
						  double p_star, p_left;
						  static double *p0D;
						  static double R2, C, R1;
						  double dt = ModelParam::dt;
						  double cl0;
						  // double time = dparam("TIME");
						  // static double CpT = dparam("Total_C_periph");
						  // static double RpT = dparam("Total_R_periph");
						  static int first, *init, *write;
						  char buf[BUFSIZ];
						  static FILE  **fp;
						  static int hisstep = ModelParam::hisSteps;

						  if(!fp){ // The first time the routine is called
							  fp = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Allocate space for "fp" to open an .lum file
						  }
						  if((!fp[n])){
							  /* Create a pointer to the .lum file. */
							  fp[n] = fopen(buf,"w");
						  }

						  // Create p0D and R1
						  if(!first){
							  init = new int [ModelParam::Ndoms]; izero(ModelParam::Ndoms,init,1); first = 1;
							  write = new int [ModelParam::Ndoms]; izero(ModelParam::Ndoms,write,1); p0D = new double [ModelParam::Ndoms];
						  }

						  // Header for the .lum file in each domain 'n'
						  if(init[n] == 0){
							  p0D[n] = 0.0;
							  init[n] = 1;
						  }

						  // We set R1 equal to the characteristic impedance
						  cl0 = A[nel-1].Get_c(A[nel-1].q-1,A[nel-1].Ao[A[nel-1].q-1]);
						  R1 = cl0*ModelParam::rho/A[nel-1].Ao[A[nel-1].q-1];

						  // Calculation of R2 = RT - R1
						  R2 = omega[n].bcval[4] - R1;

						  if (R2 < 0.0){
							  fprintf(stderr,"Error in in BC 'W': R2<0 in domain %d. \n", n+1); exit(-1);
						  }

						  // Set each peripheral compliance C
						  if(omega[n].bcval[3] == 0){
							  // Check if "Total_C_periph" is specified in the parameter list. If not, terminate with an error message.
						  }
						  else 	
							  C = omega[n].bcval[3];
						  //printf("Vessel %d CpT= %.3e  RpT= %.3e C= %.3e R= %.3e\n",n+1,CpT,RpT,C,omega[n].bcval[3]);
						  //printf("n = %d R1 = %lf C = %lg R2 = %lf \n", n+1, R1, C, R2);

						  p_left = A[nel-1].Get_P(A[nel-1].q-1,aloc[nel-1][1]);
						  write[n] += 1;
						  // Determine the upwind state (A_star,u_star).
						  WKoutflow(nel, omega[n].A, uloc[nel-1][1], aloc[nel-1][1], 0.0, p0D[n], R1, ModelParam::rho, &(u_star[n]), &(A_star[n]));

						  // Calculate the velocity enforced on the right to obtain u_star.
						  bc[i][0] = 2.0*u_star[n] - U[nel-1].h[U[nel-1].q-1];

						  // Update the value of p0D for the next time step.
						  p0D[n] += dt/C*(u_star[n]*A_star[n] - (p0D[n] - pinf)/R2);
						  break;
				  }
				  break;
			  case 'w': /** 2-element windkessel model (compliance + resistance) terminal bc.
						Only for the r.h.s. */
				  switch(i){
					  case 0:
						  bc[i][0]  = 0;
						  fprintf(stderr,"Terminal BC (w) not set up for l.h.s. \n");
						  break;
					  case 1:
						  bc[i][0]  = 0;
						  fprintf(stderr,"Terminal BC (w) not set up for l.h.s. \n");
						  break;
					  case 3: // We assume Al = Ar and we impose bc condition through U.
						  bc[i][0] =  A[nel-1].h[A[nel-1].q-1];
						  break;
					  case 4:
						  static double R, C, R_const;
						  static double *pC, *pout;
						  static int first, *init;

						  // C and R are taken from the input file .in
						  C = omega[n].bcval[3];
						  R = omega[n].bcval[4];
						  if(C < 0.0)
						  {fprintf(stderr,"Error in in BC 'w': C<0 in domain %d. \n", n+1); exit(-1);}
						  if(R < 0.0)
						  {fprintf(stderr,"Error in in BC 'w': R<0 in domain %d. \n", n+1); exit(-1);}

						  // Modified resistance
						  R_const = R*C/ModelParam::dt + 1.0;

						  // Create pC, qout and pout for each domain "n"
						  if(!first){
							  init = new int [ModelParam::Ndoms]; izero(ModelParam::Ndoms,init,1); first = 1; pC = new double [ModelParam::Ndoms]; pout = new double [ModelParam::Ndoms];
						  }

						  // Initialise pC and qout for the first time step.
						  if(init[n] == 0){
							  pC[n] = 0.0;
							  init[n] = 1;
							  if(!omega[n].bcval[3]) {
								  fprintf(stderr,"Error in BC 'Z': Resistance not set.\n");
								  exit(1);
							  }
						  }

						  pout[n] = pC[n]*R*C/ModelParam::dt/R_const + pinf/R_const;

						  // Determine the upwind state (A_star,u_star).
						  WKoutflow(nel, omega[n].A, uloc[nel-1][1], aloc[nel-1][1], 0.0, pout[n], R/R_const, ModelParam::rho, &(u_star[n]), &(A_star[n]));

						  // Calculate the value enforced on the right to obtain u_star.
						  bc[i][0] = 2.0*u_star[n] - U[nel-1].h[U[nel-1].q-1];

						  // Update the values of p0D and Q0D for the next call time step.
						  pC[n] = ( pC[n]*C/ModelParam::dt + A_star[n]*u_star[n])*R/R_const + pinf/R_const;
						  break;
				  }
				  break;
			  case 'R': // Resistance only terminal bc. Only for the r.h.s. It considers
				  // David & Moore's cerebral flow autoregulation model if total time > InitialTime.
				  switch(i){
					  case 0:
						  bc[i][0]  = 0;
						  fprintf(stderr,"BC type 'R' not set up for l.h.s. \n");
						  break;
					  case 1:
						  bc[i][0]  = 0;
						  fprintf(stderr,"BC type 'R' not set up for l.h.s. \n");
						  break;
					  case 3: // We assume Ar = Al and we impose bc condition through U.
						  bc[i][0] = A[nel-1].h[A[nel-1].q-1];
						  break;
					  case 4:
						  static double *R, *CBF;
						  /*double time = dparam("TIME");
						  static double dt = dparam("DT");
						  static int nsteps  = iparam("NSTEPS");
						  static int hisstep = iparam("HISSTEP");
						  static double InitialTime = dparam("InitialTime");
						  static double SetPointTime = dparam("SetPointTime");*/
						  static int first, *init, *write;
						  static FILE  **fp;

						  // Open the output files the first time the routine is called
						  if(!fp){ 
							  fp = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Allocate space for "fp" to open a .lum file
						  }
						  //if((!fp[n])&&(option("DLum"))){ /* Give the name of the .in file to the .lum file (plus the number of the domain if
						  //																more than one is defined). */
						  //	if(ModelParam::Ndoms == 1) sprintf(buf,"%s.lum",strtok(name,"."));
						  //	else sprintf(buf,"%s_%d.lum",strtok(name,"."),n+1);
						  //	fp[n] = fopen(buf,"w"); // Create a pointer to the .lum file
						  //}
						  // The first time this function is called:
						  if(!first){
							  init = new int [ModelParam::Ndoms]; izero(ModelParam::Ndoms,init,1); first = 1;
							  R = new double [ModelParam::Ndoms];
							  CBF = new double [ModelParam::Ndoms];
							  write = new int [ModelParam::Ndoms]; izero(ModelParam::Ndoms,write,1);
						  }

						  // Evaluate the outflow from the 1-D terminal vessel
						  CBF[n] = uloc[nel-1][1]*aloc[nel-1][1];
						  // The first time this function is called for the domain n
						  if(init[n] == 0){ 
							  R[n] = omega[n].bcval[4];
							  /*if(option("DLum")){
							  fprintf(fp[n],"# t  P1D   Q1D\n");
							  fprintf(fp[n],"#(s) (Pa) (m3/s)\n");
							  }*/
							  init[n] = 1;
						  }
						  // Dump the parameters of the terminal model in .lum taking into account the history step "HISSTEP" defined in the parameter list of the input file.
						  /*if((write[n]==g_hisstep)&&(option("DLum"))){
						  fprintf(fp[n],"%.4f %.4f %lg\n", time,CBF[n]*R[n]+pinf,CBF[n]);
						  write[n] = 0;
						  }*/

						  /*cl0 = A[nel-1].Get_co(A[nel-1].q-1);
						  double Z = cl0*ModelParam::rho/A[nel-1].Ao[A[nel-1].q-1];
						  printf("%d  %lf\n", n+1, (R[n]-Z)/(Z+R[n]));*/

						  write[n] += 1;

						  // Determine the upwind state (A_star,u_star) with the new R[n]
						  WKoutflow(nel, omega[n].A, uloc[nel-1][1], aloc[nel-1][1], 0.0, pinf, R[n], ModelParam::rho, &(u_star[n]), &(A_star[n]));

						  // Calculate the velocity enforced on the right to obtain u_star.
						  bc[i][0] = 2.0*u_star[n] - U[nel-1].h[U[nel-1].q-1];
						  break;
				    }
					break;
				}
			}
		}

		/* Generate upwind flux minus local flux for each interelement boundary "i" (i=0,...,nel+1) of the
		domain "n" (n=0,...,ModelParam::Ndoms). We consider the two virtual elements at the inlet and the outlet of the
		domain "n".
		* The upwind is based upon flux assuming subsonic. */
		double W[2], uu[2], Au[2]; // Riemann invariants, u & A at both sides of the interelement boundary "i".
		double Ao[2], co[2], Rf, dp[2], dP;
		double A_rate, K;

		// Calculate local wave speed and Riemann invariant at each side of each interelement boundary "i".
		for(i = 0; i < nel+1; ++i){
			if(ModelParam::riemann){
				// Left side
				if(!i) {  /* At the left boundary of the first element the information on the left comes from the
						  virtual element on the left of the domain. */
					cl   = A[0].Get_c(0,bc[0][0]);
					W[0] = bc[1][0] + A[0].Get_W(0,bc[0][0]);
				}
				else {
					cl = A[i-1].Get_c(A[i-1].q-1,aloc[i-1][1]);
					W[0] = uloc[i-1][1] + A[i-1].Get_W(A[i-1].q-1,aloc[i-1][1]);
				}

				// Right side
				if(i == nel) {  /* At the right boundary of the last element the information on the left comes
								from the virtual element on the right of the domain. */
					cr   = A[nel-1].Get_c(A[nel-1].q-1,bc[3][0]);
					W[1] = bc[4][0] - A[nel-1].Get_W(A[nel-1].q-1,bc[3][0]);
				}
				else {
					cr   = A[i].Get_c(0,aloc[i][0]);
					W[1] = uloc[i][0] - A[i].Get_W(0,aloc[i][0]);
				}

				// Check if flow is subsonic
				if((fabs(uloc[i-1][1]) > cl)||(fabs(uloc[i][0]) > cr)){
					fprintf(stderr,"Error in BioFlux: Flow is not subsonic.\n");
					exit(1);
				}
			}

			// A & U at each side of the interelement boundary.
			uu[0] = uloc[i-1][1];
			uu[1] = uloc[i][0];
			Au[0] = aloc[i-1][1];
			Au[1] = aloc[i][0];

			/* Call the nonlinear interelement Riemann problem solver
			at the element interface "i". The solution satisfies conservation
			of mass and continuity of the total pressure. Energy losses are
			included if "Junc_losses == 1". */

			// Couple the input
			if(ModelParam::riemann){
				// Riemann problem without considering energy losses.
				Riemann(i,nel,omega[n].A,W,uu,Au,ModelParam::rho);
			}
			else{
				Au[0] = (Au[0]+Au[1])/2;
				Au[1] = Au[0];
				uu[0] = (uu[0]+uu[1])/2;
				uu[1] = uu[0];
			}

			// Upwind fluxs at each interelement boundary "i" (i=0,...,nel+1)
			// Flow rate
			af[i-1][1] = Au[0]*uu[0]; // left of the interelement boundary "i"
			af[i][0]   = Au[1]*uu[1]; // right of the interelement boundary "i"

#ifdef DEBUG_1D
			if(pow((af[i-1][1] - af[i][0]),2.0) > ModelParam::riemannTol){
				fprintf(stderr,"Error in BioFlux: Q does not balance at the interelement boundary %d of domain %d.\n",i,n);
				exit(1);
			}
#endif
			// Energy
			// First element of domain
			if(!ModelParam::gammaI && !ModelParam::gammaII){
				if(!i){
					if(ModelParam::nonDim)
						uf[i-1][1] = 0.5*uu[0]*uu[0] + A[0].Get_P(0,Au[0]);
					else
						uf[i-1][1] = 0.5*uu[0]*uu[0] + A[0].Get_P(0,Au[0])/ModelParam::rho;
				}
				// Element whose boundaries are also elements
				else{
					if(ModelParam::nonDim)
						uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0]);
					else
						uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0])/ModelParam::rho;
				}
				// Beyond last element of domain
				if(i==nel){
					if(ModelParam::nonDim)
						uf[i][0]   = 0.5*uu[1]*uu[1] + A[nel-1].Get_P(A[nel-1].q-1,Au[1]);
					else
						uf[i][0]   = 0.5*uu[1]*uu[1] + A[nel-1].Get_P(A[nel-1].q-1,Au[1])/ModelParam::rho;
				}
				// Element whose boundaries are also elements
				else{
					if(ModelParam::nonDim)
						uf[i][0]   = 0.5*uu[1]*uu[1] + A[i].Get_P(0,Au[1]);
					else
						uf[i][0]   = 0.5*uu[1]*uu[1] + A[i].Get_P(0,Au[1])/ModelParam::rho;
				}
			}
			else{
				// Gamma
				if(ModelParam::gammaI){
					// First element of domain
					if(!i){
						if(ModelParam::nonDim){
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[0].Get_P(0,Au[0]) + gamma[0]*pow(A[0].h[0],-0.5)*Af[0].h[0];
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[0].Get_P(0,Au[1]) + gamma[0]*pow(A[0].h[0],-0.5)*Af[0].h[0];
						} 
						else{
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[0].Get_P(0,Au[0])/ModelParam::rho + gamma[0]*pow(A[0].h[0],-0.5)*Af[0].h[0]/ModelParam::rho;
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[0].Get_P(0,Au[1])/ModelParam::rho + gamma[0]*pow(A[0].h[0],-0.5)*Af[0].h[0]/ModelParam::rho;
						}

					}
					// Element whose boundaries are also elements
					if(i>0 && i<nel){
						if(ModelParam::nonDim){
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0]) + gamma[0]*0.5*(Af[i].h[0]*pow(A[i].h[0],-0.5) + Af[i-1].h[A[i-1].q-1]*pow(A[i-1].h[A[i-1].q-1],-0.5));
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[i].Get_P(0,Au[1]) + gamma[0]*0.5*(Af[i].h[0]*pow(A[i].h[0],-0.5) + Af[i-1].h[A[i-1].q-1]*pow(A[i-1].h[A[i-1].q-1],-0.5));
						}
						else{
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0])/ModelParam::rho + gamma[0]*0.5*(Af[i].h[0]*pow(A[i].h[0],-0.5) + Af[i-1].h[A[i-1].q-1]*pow(A[i-1].h[A[i-1].q-1],-0.5))/ModelParam::rho;
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[i].Get_P(0,Au[1])/ModelParam::rho + gamma[0]*0.5*(Af[i].h[0]*pow(A[i].h[0],-0.5) + Af[i-1].h[A[i-1].q-1]*pow(A[i-1].h[A[i-1].q-1],-0.5))/ModelParam::rho;
						}

					}
					// Last element of domain
					if(i==nel){
						if(ModelParam::nonDim){
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[nel-1].Get_P(A[nel-1].q-1,Au[1]) + gamma[0]*pow(A[nel-1].h[A[nel-1].q-1],-0.5)*Af[i-1].h[A[i-1].q-1];
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0]) + gamma[0]*pow(A[nel-1].h[A[nel-1].q-1],-0.5)*Af[i-1].h[A[i-1].q-1];
						}
						else{
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[nel-1].Get_P(A[nel-1].q-1,Au[1])/ModelParam::rho + gamma[0]*pow(A[nel-1].h[A[nel-1].q-1],-0.5)*Af[i-1].h[A[i-1].q-1]/ModelParam::rho;
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0])/ModelParam::rho + gamma[0]*pow(A[nel-1].h[A[nel-1].q-1],-0.5)*Af[i-1].h[A[i-1].q-1]/ModelParam::rho;
						}
					}		   
				}
				// GammaII
				else{
					// First element of domain
					if(!i){
						if(ModelParam::nonDim){
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[0].Get_P(0,Au[0]) + gamma[0]*Af[0].h[0];
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[0].Get_P(0,Au[1]) + gamma[0]*Af[0].h[0];
						}
						else{
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[0].Get_P(0,Au[0])/ModelParam::rho + gamma[0]*Af[0].h[0]/ModelParam::rho;
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[0].Get_P(0,Au[1])/ModelParam::rho + gamma[0]*Af[0].h[0]/ModelParam::rho;
						}
					}
					// Element whose boundaries are also elements
					if(i>0 && i<nel){
						if(ModelParam::nonDim){
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0]) + gamma[0]*0.5*(Af[i].h[0] + Af[i-1].h[A[i-1].q-1]);
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[i].Get_P(0,Au[1]) + gamma[0]*0.5*(Af[i].h[0] + Af[i-1].h[A[i-1].q-1]);
						}
						else{
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0])/ModelParam::rho + gamma[0]*0.5*(Af[i].h[0] + Af[i-1].h[A[i-1].q-1])/ModelParam::rho;
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[i].Get_P(0,Au[1])/ModelParam::rho + gamma[0]*0.5*(Af[i].h[0] + Af[i-1].h[A[i-1].q-1])/ModelParam::rho;
						} 
					}
					// Last element of domain
					if(i==nel){
						if(ModelParam::nonDim){
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[nel-1].Get_P(A[nel-1].q-1,Au[1]) + gamma[0]*Af[i-1].h[A[i-1].q-1];
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0]) + gamma[0]*Af[i-1].h[A[i-1].q-1];
						}
						else{
							uf[i][0]   = 0.5*uu[1]*uu[1] + A[nel-1].Get_P(A[nel-1].q-1,Au[1])/ModelParam::rho + gamma[0]*Af[i-1].h[A[i-1].q-1]/ModelParam::rho;
							uf[i-1][1] = 0.5*uu[0]*uu[0] + A[i-1].Get_P(A[i-1].q-1,Au[0])/ModelParam::rho + gamma[0]*Af[i-1].h[A[i-1].q-1]/ModelParam::rho;
						}
					}
				}
			}
#ifdef DEBUG_1D
			if(pow((uf[i-1][1] - uf[i][0]),2.0) > ModelParam::riemannTol){
				fprintf(stderr,"Error in BioFlux: energy does not balance at the interelement boundary %d of domain %d.\n",i,n);
				exit(1);
			}
#endif
		}

		// Subtract local fluxs from upwind fluxs at the inlet and the outlet of each element "i" (i=0,...,nel-1).
		if(!ModelParam::gammaI && !ModelParam::gammaII){
			for(i = 0; i < nel; ++i){
				// Inlet of the element 'i'
				af[i][0] -= aloc[i][0]*uloc[i][0];
				if(ModelParam::nonDim)
					uf[i][0] -= 0.5*uloc[i][0]*uloc[i][0] + A[i].Get_P(0,aloc[i][0]);
				else
					uf[i][0] -= 0.5*uloc[i][0]*uloc[i][0] + A[i].Get_P(0,aloc[i][0])/ModelParam::rho;
				// Outlet of the element 'i'
				af[i][1] -= aloc[i][1]*uloc[i][1];
				if(ModelParam::nonDim)
					uf[i][1] -= 0.5*uloc[i][1]*uloc[i][1] + A[i].Get_P(A[i].q-1,aloc[i][1]);
				else
					uf[i][1] -= 0.5*uloc[i][1]*uloc[i][1] + A[i].Get_P(A[i].q-1,aloc[i][1])/ModelParam::rho;
			}
		}
		else{
			for(i = 0; i < nel; ++i){
				if(ModelParam::gammaI){
					// Inlet of the element 'i'
					af[i][0] -= aloc[i][0]*uloc[i][0];
					if(ModelParam::nonDim)
						uf[i][0] -= 0.5*uloc[i][0]*uloc[i][0] + A[i].Get_P(0,aloc[i][0]) + gamma[0]/sqrt(A[i].h[0])*Af[i].h[0];
					else
						uf[i][0] -= 0.5*uloc[i][0]*uloc[i][0] + A[i].Get_P(0,aloc[i][0])/ModelParam::rho + gamma[0]/sqrt(A[i].h[0])*Af[i].h[0]/ModelParam::rho;
					// Outlet of the element 'i'
					af[i][1] -= aloc[i][1]*uloc[i][1];
					if(ModelParam::nonDim)
						uf[i][1] -= 0.5*uloc[i][1]*uloc[i][1] + A[i].Get_P(A[i].q-1,aloc[i][1]) + gamma[ModelParam::oneD_q-1]/sqrt(A[i].h[A[i].q-1])*Af[i].h[A[i].q-1];
					else
						uf[i][1] -= 0.5*uloc[i][1]*uloc[i][1] + A[i].Get_P(A[i].q-1,aloc[i][1])/ModelParam::rho + gamma[ModelParam::oneD_q-1]/sqrt(A[i].h[A[i].q-1])*Af[i].h[A[i].q-1]/ModelParam::rho;
				}
				else{
					// Inlet of the element 'i'
					af[i][0] -= aloc[i][0]*uloc[i][0];
					if(ModelParam::nonDim)
						uf[i][0] -= 0.5*uloc[i][0]*uloc[i][0] + A[i].Get_P(0,aloc[i][0]) + gamma[0]*Af[i].h[0];
					else
						uf[i][0] -= 0.5*uloc[i][0]*uloc[i][0] + A[i].Get_P(0,aloc[i][0])/ModelParam::rho + gamma[0]*Af[i].h[0]/ModelParam::rho;
					// Outlet of the element 'i'
					af[i][1] -= aloc[i][1]*uloc[i][1];
					if(ModelParam::nonDim)
						uf[i][1] -= 0.5*uloc[i][1]*uloc[i][1] + A[i].Get_P(A[i].q-1,aloc[i][1]) + gamma[ModelParam::oneD_q-1]*Af[i].h[A[i].q-1];
					else
						uf[i][1] -= 0.5*uloc[i][1]*uloc[i][1] + A[i].Get_P(A[i].q-1,aloc[i][1])/ModelParam::rho + gamma[ModelParam::oneD_q-1]*Af[i].h[A[i].q-1]/ModelParam::rho;
				}
			}
		}

		for(i = 0; i < nel; ++i){ // For every element of the domain "n",
			/* Add to the projected flux and source term, Af.hj and Uf.hj, the difference between the upwind and
			the local fluxes multiplied by the Legendre polynomial expansion basis and the inverse of the
			Jacobian. */
			for(j = 0; j < ModelParam::oneD_L; ++j){
				Af[i].hj[j] += m_g1[j][ 0 ]*af[i][0]/Af[i].jac;
				Af[i].hj[j] -= m_g1[j][ModelParam::oneD_q-1]*af[i][1]/Af[i].jac;

				Uf[i].hj[j] += m_g1[j][ 0 ]*uf[i][0]/Uf[i].jac;
				Uf[i].hj[j] -= m_g1[j][ModelParam::oneD_q-1]*uf[i][1]/Uf[i].jac;
			}
		}
	}
}

int OneDSolver::Riemann(int i, int nel, Element *A, double *W, double *uu, double *Au, double rho){
	double f[4], g[4], inv_J[4][4];
	double cl, cr, k, k1, k2;
	int proceed = 1, iter = 0;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector and the wave speed at each vessel.
		if(!i) {
			f[0] = uu[0] + A[0].Get_W(0,Au[0]) - W[0];
			f[1] = uu[1] - A[0].Get_W(0,Au[1]) - W[1];
			f[2] = uu[0]*Au[0] - uu[1]*Au[1];
			if(ModelParam::nonDim)
				f[3] = uu[0]*uu[0] + 2.0*A[0].Get_P(0,Au[0])
				- uu[1]*uu[1] - 2.0*A[0].Get_P(0,Au[1]);
			else
				f[3] = uu[0]*uu[0] + 2.0/rho*A[0].Get_P(0,Au[0])
				- uu[1]*uu[1] - 2.0/rho*A[0].Get_P(0,Au[1]);

			cl   = A[0].Get_c(0,Au[0]);
			cr   = A[0].Get_c(0,Au[1]);
		}
		else{
			if(i == nel) {
				f[0] = uu[0] + A[nel-1].Get_W(A[nel-1].q-1,Au[0]) - W[0];
				f[1] = uu[1] - A[nel-1].Get_W(A[nel-1].q-1,Au[1]) - W[1];
				f[2] = uu[0]*Au[0] - uu[1]*Au[1];
				if(ModelParam::nonDim)
					f[3] = uu[0]*uu[0] + 2.0*A[nel-1].Get_P(A[nel-1].q-1,Au[0])
					- uu[1]*uu[1] - 2.0*A[nel-1].Get_P(A[nel-1].q-1,Au[1]);
				else
					f[3] = uu[0]*uu[0] + 2.0/rho*A[nel-1].Get_P(A[nel-1].q-1,Au[0])
					- uu[1]*uu[1] - 2.0/rho*A[nel-1].Get_P(A[nel-1].q-1,Au[1]);
				cl   = A[nel-1].Get_c(A[nel-1].q-1,Au[0]);
				cr   = A[nel-1].Get_c(A[nel-1].q-1,Au[1]);
			}
			else {
				f[0] = uu[0] + A[i-1].Get_W(A[i-1].q-1,Au[0]) - W[0];
				f[1] = uu[1] - A[i].Get_W(0,Au[1]) - W[1];
				f[2] = uu[0]*Au[0] - uu[1]*Au[1];
				if(ModelParam::nonDim)
					f[3] = uu[0]*uu[0] + 2.0*A[i-1].Get_P(A[i-1].q-1,Au[0])
					- uu[1]*uu[1] - 2.0*A[i].Get_P(0,Au[1]);
				else
					f[3] = uu[0]*uu[0] + 2.0/rho*A[i-1].Get_P(A[i-1].q-1,Au[0])
					- uu[1]*uu[1] - 2.0/rho*A[i].Get_P(0,Au[1]);
				cl   = A[i-1].Get_c(A[i-1].q-1,Au[0]);
				cr   = A[i].Get_c(0,Au[1]);
			}
		}

		// Inverse Jacobian matrix dfdv
		k = (cl*Au[1]+Au[0]*cr);
		k1 = (cl-uu[0])*k;
		inv_J[0][0] = (Au[1]*cl*cl-cr*uu[0]*Au[0])/k1;
		inv_J[0][1] = Au[1]*(cr-uu[1])*cl/k1;
		inv_J[0][2] = cl*cr/k1;
		inv_J[0][3] = -0.5*cl*Au[1]/k1;

		k2 = (cr+uu[1])*k;
		inv_J[1][0] = Au[0]*(cl+uu[0])*cr/k2;
		inv_J[1][1] = (cl*uu[1]*Au[1]+cr*cr*Au[0])/k2;
		inv_J[1][2] = -cl*cr/k2;
		inv_J[1][3] = -0.5*Au[0]*cr/k2;

		inv_J[2][0] = Au[0]*(Au[0]*cr-uu[0]*Au[1])/k1;
		inv_J[2][1] = -Au[0]*Au[1]*(cr-uu[1])/k1;
		inv_J[2][2] = -Au[0]*cr/k1;
		inv_J[2][3] = 0.5*Au[1]*Au[0]/k1;

		inv_J[3][0] = Au[0]*Au[1]*(cl+uu[0])/k2;
		inv_J[3][1] = -Au[1]*(cl*Au[1]+uu[1]*Au[0])/k2;
		inv_J[3][2] = -cl*Au[1]/k2;
		inv_J[3][3] = -0.5*Au[1]*Au[0]/k2;

		/* Solve the linear system by inverting analyticaly the Jacobian: g = (dfdv)^(-1)*f */
		cblas_dgemv(CblasColMajor,CblasTrans, 4, 4, 1, *inv_J, 4, f, 1, 0, g, 1);

		// Update solution: x_new = x_old - dx
		uu[0] -= g[0];
		uu[1] -= g[1];
		Au[0] -= g[2];
		Au[1] -= g[3];

		// Check if the error of the solution is smaller than TOL.
		if((g[0]*g[0] + g[1]*g[1] + g[2]*g[2] + g[3]*g[3]) < ModelParam::riemannTol)      proceed = 0;
	}

	if(iter >= MAX_ITER){
		fprintf(stderr,"Error in Riemann: iteration failed to converge. \n");
		exit(-1);
	}

	return iter;
}

int OneDSolver::BifurRiem1to2(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho){
	double f[6], g[6], inv_J[6][6];
	double c1, c2, c3, k, k1, k2, k3;
	int proceed = 1, iter = 0, six, one, zero;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector
		f[0] = uu[0] + A_p[nel_p-1].Get_W(A_p[nel_p-1].q-1,Au[0]) - W[0]; // Continuity of the forward parent Riemann invariant.
		f[1] = uu[1] - A_d1[0].Get_W(0,Au[1]) - W[1]; // Continuity of the backward daughter 1 Riemann invariant.
		f[2] = uu[2] - A_d2[0].Get_W(0,Au[2]) - W[2]; // Continuity of the backward daughter 2 Riemann invariant.
		f[3] = uu[0]*Au[0] - uu[1]*Au[1] - uu[2]*Au[2]; // Conservation of mass.
		if(ModelParam::nonDim){
			f[4] = uu[0]*uu[0] + 2.0*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0])  /* Continuity of the total pressure
																				 between the parent and the daughter 1 vessels. */
																				 - uu[1]*uu[1] - 2.0*A_d1[0].Get_P(0,Au[1]);
			f[5] = uu[0]*uu[0] + 2.0*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0])  /* Continuity of the total pressure
																				 between the parent and the daughter 2 vessels. */
																				 - uu[2]*uu[2] - 2.0*A_d2[0].Get_P(0,Au[2]);
		}
		else{
			f[4] = uu[0]*uu[0] + 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0])  /* Continuity of the total pressure
																					 between the parent and the daughter 1 vessels. */
																					 - uu[1]*uu[1] - 2.0/rho*A_d1[0].Get_P(0,Au[1]);
			f[5] = uu[0]*uu[0] + 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0])  /* Continuity of the total pressure
																					 between the parent and the daughter 2 vessels. */
																					 - uu[2]*uu[2] - 2.0/rho*A_d2[0].Get_P(0,Au[2]);
		}

		// Wave speed at each vessel
		c1 = A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0]);
		c2 = A_d1[0].Get_c(0,Au[1]);
		c3 = A_d2[0].Get_c(0,Au[2]);

		// Inverse Jacobian matrix dfdv
		k = c1*Au[1]*c3+Au[0]*c3*c2+Au[2]*c1*c2;
		k1 = (c1-uu[0])*k;
		inv_J[0][0] = (-c2*uu[0]*c3*Au[0]+Au[2]*c2*c1*c1+Au[1]*c1*c1*c3)/k1;
		inv_J[0][1] = Au[1]*(c2-uu[1])*c1*c3/k1;
		inv_J[0][2] = Au[2]*(c3-uu[2])*c1*c2/k1;
		inv_J[0][3] = c1*c2*c3/k1;
		inv_J[0][4] = -0.5*c1*Au[1]*c3/k1;
		inv_J[0][5] = -0.5*Au[2]*c1*c2/k1;

		k2 = (c2+uu[1])*k;
		inv_J[1][0] = Au[0]*(c1+uu[0])*c2*c3/k2;
		inv_J[1][1] = (c1*uu[1]*c3*Au[1]+Au[2]*c1*c2*c2+c3*c2*c2*Au[0])/k2;
		inv_J[1][2] = -Au[2]*(c3-uu[2])*c1*c2/k2;
		inv_J[1][3] = -c1*c2*c3/k2;
		inv_J[1][4] = -0.5*(c1*Au[2]+Au[0]*c3)*c2/k2;
		inv_J[1][5] = 0.5*Au[2]*c1*c2/k2;

		k3 = (c3+uu[2])*k;
		inv_J[2][0] = Au[0]*(c1+uu[0])*c2*c3/k3;
		inv_J[2][1] = -Au[1]*(c2-uu[1])*c1*c3/k3;
		inv_J[2][2] = (c1*c2*uu[2]*Au[2]+c1*Au[1]*c3*c3+c2*c3*c3*Au[0])/k3;
		inv_J[2][3] = -c1*c2*c3/k3;
		inv_J[2][4] = 0.5*c1*Au[1]*c3/k3;
		inv_J[2][5] = -0.5*(Au[1]*c1+c2*Au[0])*c3/k3;

		inv_J[3][0] = Au[0]*(Au[0]*c3*c2-uu[0]*c3*Au[1]-uu[0]*c2*Au[2])/k1;
		inv_J[3][1] = -Au[0]*Au[1]*(c2-uu[1])*c3/k1;
		inv_J[3][2] = -Au[0]*Au[2]*(c3-uu[2])*c2/k1;
		inv_J[3][3] = -Au[0]*c3*c2/k1;
		inv_J[3][4] = 0.5*Au[0]*Au[1]*c3/k1;
		inv_J[3][5] = 0.5*Au[0]*c2*Au[2]/k1;

		inv_J[4][0] = Au[0]*Au[1]*(c1+uu[0])*c3/k2;
		inv_J[4][1] = -Au[1]*(c1*Au[1]*c3+c1*uu[1]*Au[2]+c3*uu[1]*Au[0])/k2;
		inv_J[4][2] = -Au[2]*Au[1]*(c3-uu[2])*c1/k2;
		inv_J[4][3] = -c1*Au[1]*c3/k2;
		inv_J[4][4] = -0.5*Au[1]*(c1*Au[2]+Au[0]*c3)/k2;
		inv_J[4][5] = 0.5*Au[2]*Au[1]*c1/k2;

		inv_J[5][0] = Au[0]*Au[2]*(c1+uu[0])*c2/k3;
		inv_J[5][1] = -Au[2]*Au[1]*(c2-uu[1])*c1/k3;
		inv_J[5][2] = -Au[2]*(Au[2]*c1*c2+c1*uu[2]*Au[1]+c2*uu[2]*Au[0])/k3;
		inv_J[5][3] = -Au[2]*c1*c2/k3;
		inv_J[5][4] = 0.5*Au[2]*Au[1]*c1/k3;
		inv_J[5][5] = -0.5*Au[2]*(Au[1]*c1+c2*Au[0])/k3;

		one = 1;
		six = 6;
		zero = 0;

		/* Solve the linear system by inverting analyticaly the Jacobian: g = (dfdv)^(-1)*f */
		cblas_dgemv(CblasColMajor,CblasTrans, six, six, one, *inv_J, six, f, one, zero, g, one);

		// Update solution: x_new = x_old - dx
		uu[0] -= g[0];
		uu[1] -= g[1];
		uu[2] -= g[2];
		Au[0] -= g[3];
		Au[1] -= g[4];
		Au[2] -= g[5];

		// Check if the error of the solution is smaller than TOL.
		if((g[0]*g[0] + g[1]*g[1] + g[2]*g[2] + g[3]*g[3] +	g[4]*g[4] + g[5]*g[5]) < ModelParam::riemannTol)      proceed = 0;
	}

	/*if(iter > 3){
	printf("Bifur1to2 iter = %d, parent id=%d\n", iter, d);
	}*/

	if(iter >= MAX_ITER){
		fprintf(stderr,"Error in BifurRiem1to2: iteration failed to converge. \n");
		exit(-1);
	}

	return iter;
}

int OneDSolver::BifurRiem1to2_losses_div_right(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K31, double K32){
	double f[6],dfdv[6][6];
	int proceed = 1,iter = 0,ipiv[6],six,one,info;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector.
		f[0] = uu[0] + A_p[nel_p-1].Get_W(A_p[nel_p-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_d1[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[2] - A_d2[0].Get_W(0,Au[2]) - W[2];
		f[3] = uu[0]*Au[0] - uu[1]*Au[1] - uu[2]*Au[2];
		f[4] = uu[0]*uu[0]*(1.0 - K32) + 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0])
			- uu[1]*uu[1] - 2.0/rho*A_d1[0].Get_P(0,Au[1]);
		f[5] = uu[0]*uu[0]*(1.0 - K31) + 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0])
			- uu[2]*uu[2] - 2.0/rho*A_d2[0].Get_P(0,Au[2]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  0;
		dfdv[0][3] =  A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0])/Au[0];
		dfdv[0][4] =  0;
		dfdv[0][5] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  0;
		dfdv[1][4] =  -A_d1[0].Get_c(0,Au[1])/Au[1];
		dfdv[1][5] =  0;

		dfdv[2][0] =  0;
		dfdv[2][1] =  0;
		dfdv[2][2] =  1;
		dfdv[2][3] =  0;
		dfdv[2][4] =  0;
		dfdv[2][5] =  -A_d2[0].Get_c(0,Au[2])/Au[2];

		dfdv[3][0] =  Au[0];
		dfdv[3][1] = -Au[1];
		dfdv[3][2] = -Au[2];
		dfdv[3][3] =  uu[0];
		dfdv[3][4] = -uu[1];
		dfdv[3][5] = -uu[2];

		dfdv[4][0] =  2*uu[0]*(1.0 - K32);
		dfdv[4][1] = -2*uu[1];
		dfdv[4][2] =  0;
		dfdv[4][3] =  2.0*pow(A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0]),2)/Au[0];
		dfdv[4][4] = -2.0*pow(A_d1[0].Get_c(0,Au[1]),2)/Au[1];
		dfdv[4][5] =  0;

		dfdv[5][0] =  2*uu[0]*(1.0 - K31);
		dfdv[5][1] =  0;
		dfdv[5][2] = -2*uu[2];
		dfdv[5][3] =  2.0*pow(A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0]),2)/Au[0];
		dfdv[5][4] =  0;
		dfdv[5][5] = -2.0*pow(A_d2[0].Get_c(0,Au[2]),2)/Au[2];

		one = 1;
		six = 6;

		/*  Solve the linear system by factorising dfdv into LU by dgetrf and solving the systems Ly=f and
		U(dx)=y	by dgetrs. dx is stored in f. */
		dgetrf(&six, &six, *dfdv, &six, ipiv, &info); // LU factorisation
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_div_right: dgetrf failed at iter %d. \n",iter); exit(1);}

		dgetrs("T", &six, &one, *dfdv, &six, ipiv, f, &six, &info); // Ly=f and U(dx)=y solution
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_div_right: dgetrs failed at iter %d. \n",iter); exit(1);}

		// Update solution: x_new = x_old - dx
		uu[0] -= f[0];
		uu[1] -= f[1];
		uu[2] -= f[2];
		Au[0] -= f[3];
		Au[1] -= f[4];
		Au[2] -= f[5];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3] + f[4]*f[4] + f[5]*f[5]) < ModelParam::riemannTol)
			proceed = 0;

	}
	if(iter >= MAX_ITER)
		fprintf(stderr,"Error in BifurRiem1to2_losses_div_right: iteration failed to converge. \n");

	return iter;
}

int OneDSolver::BifurRiem1to2_losses_div_left(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K31, double K32){
	double f[6],dfdv[6][6];
	int proceed = 1,iter = 0,ipiv[6],six,one,info;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector.
		f[0] = uu[0] + A_p[nel_p-1].Get_W(A_p[nel_p-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_d1[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[2] - A_d2[0].Get_W(0,Au[2]) - W[2];
		f[3] = uu[0]*Au[0] - uu[1]*Au[1] - uu[2]*Au[2];
		f[4] = uu[1]*uu[1]*(1.0 - K32) + 2.0/rho*A_d1[0].Get_P(0,Au[1])
			- uu[0]*uu[0] - 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0]);
		f[5] = uu[1]*uu[1]*(1.0 - K31) + 2.0/rho*A_d1[0].Get_P(0,Au[1])
			- uu[2]*uu[2] - 2.0/rho*A_d2[0].Get_P(0,Au[2]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  0;
		dfdv[0][3] =  A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0])/Au[0];
		dfdv[0][4] =  0;
		dfdv[0][5] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  0;
		dfdv[1][4] =  -A_d1[0].Get_c(0,Au[1])/Au[1];
		dfdv[1][5] =  0;

		dfdv[2][0] =  0;
		dfdv[2][1] =  0;
		dfdv[2][2] =  1;
		dfdv[2][3] =  0;
		dfdv[2][4] =  0;
		dfdv[2][5] =  -A_d2[0].Get_c(0,Au[2])/Au[2];

		dfdv[3][0] =  Au[0];
		dfdv[3][1] = -Au[1];
		dfdv[3][2] = -Au[2];
		dfdv[3][3] =  uu[0];
		dfdv[3][4] = -uu[1];
		dfdv[3][5] = -uu[2];

		dfdv[4][0] = -2*uu[0];
		dfdv[4][1] =  2*uu[1]*(1.0 - K32);
		dfdv[4][2] =  0;
		dfdv[4][3] = -2.0*pow(A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0]),2)/Au[0];
		dfdv[4][4] =  2.0*pow(A_d1[0].Get_c(0,Au[1]),2)/Au[1];
		dfdv[4][5] =  0;

		dfdv[5][0] =  0;
		dfdv[5][1] =  2*uu[1]*(1.0 - K31);
		dfdv[5][2] = -2*uu[2];
		dfdv[5][3] =  0;
		dfdv[5][4] =  2.0*pow(A_d1[0].Get_c(0,Au[1]),2)/Au[1];
		dfdv[5][5] = -2.0*pow(A_d2[0].Get_c(0,Au[2]),2)/Au[2];

		one = 1;
		six = 6;

		/*  Solve the linear system by factorising dfdv into LU by dgetrf and solving the systems Ly=f and
		U(dx)=y	by dgetrs. dx is stored in f. */
		dgetrf(&six, &six, *dfdv, &six, ipiv, &info); // LU factorisation
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_div_left: dgetrf failed at iter %d. \n",iter); exit(1);}

		dgetrs("T", &six, &one, *dfdv, &six, ipiv, f, &six, &info); // Ly=f and U(dx)=y solution
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_div_left: dgetrs failed at iter %d. \n",iter); exit(1);}

		// Update solution: x_new = x_old - dx
		uu[0] -= f[0];
		uu[1] -= f[1];
		uu[2] -= f[2];
		Au[0] -= f[3];
		Au[1] -= f[4];
		Au[2] -= f[5];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3] + f[4]*f[4] + f[5]*f[5]) < ModelParam::riemannTol)      
			proceed = 0;

	}
	if(iter >= MAX_ITER)
		fprintf(stderr,"Error in BifurRiem1to2_losses_div_left: iteration failed to converge. \n");

	return iter;
}

int OneDSolver::BifurRiem1to2_losses_comb_left(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K13, double K23){
	double f[6],dfdv[6][6];
	int proceed = 1,iter = 0,ipiv[6],six,one,info;

	while((proceed)&&(iter++ < MAX_ITER)){

		// Calculate constraint vector.
		f[0] = uu[0] + A_p[nel_p-1].Get_W(A_p[nel_p-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_d1[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[2] - A_d2[0].Get_W(0,Au[2]) - W[2];
		f[3] = uu[0]*Au[0] - uu[1]*Au[1] - uu[2]*Au[2];
		f[4] = uu[1]*uu[1]*(1.0 + K23) + 2.0/rho*A_d1[0].Get_P(0,Au[1])
			- uu[0]*uu[0] - 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0]);
		f[5] = uu[1]*uu[1]*(1.0 + K13) + 2.0/rho*A_d1[0].Get_P(0,Au[1])
			- uu[2]*uu[2] - 2.0/rho*A_d2[0].Get_P(0,Au[2]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  0;
		dfdv[0][3] =  A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0])/Au[0];
		dfdv[0][4] =  0;
		dfdv[0][5] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  0;
		dfdv[1][4] =  -A_d1[0].Get_c(0,Au[1])/Au[1];
		dfdv[1][5] =  0;

		dfdv[2][0] =  0;
		dfdv[2][1] =  0;
		dfdv[2][2] =  1;
		dfdv[2][3] =  0;
		dfdv[2][4] =  0;
		dfdv[2][5] =  -A_d2[0].Get_c(0,Au[2])/Au[2];

		dfdv[3][0] =  Au[0];
		dfdv[3][1] = -Au[1];
		dfdv[3][2] = -Au[2];
		dfdv[3][3] =  uu[0];
		dfdv[3][4] = -uu[1];
		dfdv[3][5] = -uu[2];

		dfdv[4][0] = -2*uu[0];
		dfdv[4][1] =  2*uu[1]*(1.0 + K23);
		dfdv[4][2] =  0;
		dfdv[4][3] = -2.0*pow(A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0]),2)/Au[0];
		dfdv[4][4] =  2.0*pow(A_d1[0].Get_c(0,Au[1]),2)/Au[1];
		dfdv[4][5] =  0;

		dfdv[5][0] =  0;
		dfdv[5][1] =  2*uu[1]*(1.0 + K13);
		dfdv[5][2] = -2*uu[2];
		dfdv[5][3] =  0;
		dfdv[5][4] =  2.0*pow(A_d1[0].Get_c(0,Au[1]),2)/Au[1];
		dfdv[5][5] = -2.0*pow(A_d2[0].Get_c(0,Au[2]),2)/Au[2];

		one = 1;
		six = 6;

		/*  Solve the linear system by factorising dfdv into LU by dgetrf and solving the systems Ly=f and
		U(dx)=y	by dgetrs. dx is stored in f. */
		dgetrf(&six, &six, *dfdv, &six, ipiv, &info); // LU factorisation
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_comb_right: dgetrf failed at iter %d. \n",iter); exit(1);}

		dgetrs("T", &six, &one, *dfdv, &six, ipiv, f, &six, &info); // Ly=f and U(dx)=y solution
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_comb_right: dgetrs failed at iter %d. \n",iter); exit(1);}

		// Update solution: x_new = x_old - dx
		uu[0] -= f[0];
		uu[1] -= f[1];
		uu[2] -= f[2];
		Au[0] -= f[3];
		Au[1] -= f[4];
		Au[2] -= f[5];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3] + f[4]*f[4] + f[5]*f[5]) < ModelParam::riemannTol)      
			proceed = 0;
	}
	if(iter >= MAX_ITER)
		fprintf(stderr,"Error in BifurRiem1to2_losses_comb_right: iteration failed to converge. \n");

	return iter;
}

int OneDSolver::BifurRiem1to2_losses_comb_right(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K13, double K23){
	double f[6],dfdv[6][6];
	int proceed = 1,iter = 0,ipiv[6],six,one,info;

	while((proceed)&&(iter++ < MAX_ITER)){

		// Calculate constraint vector.
		f[0] = uu[0] + A_p[nel_p-1].Get_W(A_p[nel_p-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_d1[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[2] - A_d2[0].Get_W(0,Au[2]) - W[2];
		f[3] = uu[0]*Au[0] - uu[1]*Au[1] - uu[2]*Au[2];
		f[4] = uu[0]*uu[0]*(1.0 + K23) + 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0])
			- uu[1]*uu[1] - 2.0/rho*A_d1[0].Get_P(0,Au[1]);
		f[5] = uu[0]*uu[0]*(1.0 + K13) + 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0])
			- uu[2]*uu[2] - 2.0/rho*A_d2[0].Get_P(0,Au[2]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  0;
		dfdv[0][3] =  A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0])/Au[0];
		dfdv[0][4] =  0;
		dfdv[0][5] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  0;
		dfdv[1][4] =  -A_d1[0].Get_c(0,Au[1])/Au[1];
		dfdv[1][5] =  0;

		dfdv[2][0] =  0;
		dfdv[2][1] =  0;
		dfdv[2][2] =  1;
		dfdv[2][3] =  0;
		dfdv[2][4] =  0;
		dfdv[2][5] =  -A_d2[0].Get_c(0,Au[2])/Au[2];

		dfdv[3][0] =  Au[0];
		dfdv[3][1] = -Au[1];
		dfdv[3][2] = -Au[2];
		dfdv[3][3] =  uu[0];
		dfdv[3][4] = -uu[1];
		dfdv[3][5] = -uu[2];

		dfdv[4][0] =  2*uu[0]*(1.0 + K23);
		dfdv[4][1] = -2*uu[1];
		dfdv[4][2] =  0;
		dfdv[4][3] =  2.0*pow(A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0]),2)/Au[0];
		dfdv[4][4] = -2.0*pow(A_d1[0].Get_c(0,Au[1]),2)/Au[1];
		dfdv[4][5] =  0;

		dfdv[5][0] =  2*uu[0]*(1.0 + K13);
		dfdv[5][1] =  0;
		dfdv[5][2] = -2*uu[2];
		dfdv[5][3] =  2.0*pow(A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0]),2)/Au[0];
		dfdv[5][4] =  0;
		dfdv[5][5] = -2.0*pow(A_d2[0].Get_c(0,Au[2]),2)/Au[2];

		one = 1;
		six = 6;

		/*  Solve the linear system by factorising dfdv into LU by dgetrf and solving the systems Ly=f and
		U(dx)=y	by dgetrs. dx is stored in f. */
		dgetrf(&six, &six, *dfdv, &six, ipiv, &info); // LU factorisation
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_comb_left: dgetrf failed at iter %d. \n",iter); exit(1);}

		dgetrs("T", &six, &one, *dfdv, &six, ipiv, f, &six, &info); // Ly=f and U(dx)=y solution
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_comb_left: dgetrs failed at iter %d. \n",iter); exit(1);}

		// Update solution: x_new = x_old - dx
		uu[0] -= f[0];
		uu[1] -= f[1];
		uu[2] -= f[2];
		Au[0] -= f[3];
		Au[1] -= f[4];
		Au[2] -= f[5];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3] + f[4]*f[4] + f[5]*f[5]) < ModelParam::riemannTol)
			proceed = 0;

	}
	if(iter >= MAX_ITER)
		fprintf(stderr,"Error in BifurRiem1to2_losses_comb_left: iteration failed to converge. \n");

	return iter;
}

int OneDSolver::BifurRiem1to2_losses_ent(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K31, double K32){
	double f[6],dfdv[6][6];
	int proceed = 1,iter = 0,ipiv[6],six,one,info;

	while((proceed)&&(iter++ < MAX_ITER)){

		// Calculate constraint vector.
		f[0] = uu[0] + A_p[nel_p-1].Get_W(A_p[nel_p-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_d1[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[2] - A_d2[0].Get_W(0,Au[2]) - W[2];
		f[3] = uu[0]*Au[0] - uu[1]*Au[1] - uu[2]*Au[2];
		f[4] = uu[2]*uu[2]*(1.0 - K32) + 2.0/rho*A_d2[0].Get_P(0,Au[2])
			- uu[0]*uu[0] - 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0]);
		f[5] = uu[2]*uu[2]*(1.0 - K31) + 2.0/rho*A_d2[0].Get_P(0,Au[2])
			- uu[1]*uu[1] - 2.0/rho*A_d1[0].Get_P(0,Au[1]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  0;
		dfdv[0][3] =  A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0])/Au[0];
		dfdv[0][4] =  0;
		dfdv[0][5] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  0;
		dfdv[1][4] =  -A_d1[0].Get_c(0,Au[1])/Au[1];
		dfdv[1][5] =  0;

		dfdv[2][0] =  0;
		dfdv[2][1] =  0;
		dfdv[2][2] =  1;
		dfdv[2][3] =  0;
		dfdv[2][4] =  0;
		dfdv[2][5] =  -A_d2[0].Get_c(0,Au[2])/Au[2];

		dfdv[3][0] =  Au[0];
		dfdv[3][1] = -Au[1];
		dfdv[3][2] = -Au[2];
		dfdv[3][3] =  uu[0];
		dfdv[3][4] = -uu[1];
		dfdv[3][5] = -uu[2];

		dfdv[4][0] = -2*uu[0];
		dfdv[4][1] =  0;
		dfdv[4][2] =  2*uu[2]*(1.0 - K32);
		dfdv[4][3] = -2.0*pow(A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0]),2)/Au[0];
		dfdv[4][4] =  0;
		dfdv[4][5] =  2.0*pow(A_d2[0].Get_c(0,Au[2]),2)/Au[2];

		dfdv[5][0] =  0;
		dfdv[5][1] = -2*uu[1];
		dfdv[5][2] =  2*uu[2]*(1.0 - K31);
		dfdv[5][3] =  0;
		dfdv[5][4] = -2.0*pow(A_d1[0].Get_c(0,Au[1]),2)/Au[1];
		dfdv[5][5] =  2.0*pow(A_d2[0].Get_c(0,Au[2]),2)/Au[2];

		one = 1;
		six = 6;

		/*  Solve the linear system by factorising dfdv into LU by dgetrf and solving the systems Ly=f and
		U(dx)=y	by dgetrs. dx is stored in f. */
		dgetrf(&six, &six, *dfdv, &six, ipiv, &info); // LU factorisation
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_ent: dgetrf failed at iter %d. \n",iter); exit(1);}

		dgetrs("T", &six, &one, *dfdv, &six, ipiv, f, &six, &info); // Ly=f and U(dx)=y solution
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_ent: dgetrs failed at iter %d. \n",iter); exit(1);}

		// Update solution: x_new = x_old - dx
		uu[0] -= f[0];
		uu[1] -= f[1];
		uu[2] -= f[2];
		Au[0] -= f[3];
		Au[1] -= f[4];
		Au[2] -= f[5];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3] + f[4]*f[4] + f[5]*f[5]) < ModelParam::riemannTol)
			proceed = 0;

	}
	if(iter >= MAX_ITER)
		fprintf(stderr,"Error in BifurRiem1to2_losses_ent: iteration failed to converge. \n");

	return iter;
}

int OneDSolver::BifurRiem1to2_losses_leav(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K13, double K23){
	double f[6],dfdv[6][6];
	int proceed = 1,iter = 0,ipiv[6],six,one,info;

	while((proceed)&&(iter++ < MAX_ITER)){

		// Calculate constraint vector.
		f[0] = uu[0] + A_p[nel_p-1].Get_W(A_p[nel_p-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_d1[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[2] - A_d2[0].Get_W(0,Au[2]) - W[2];
		f[3] = uu[0]*Au[0] - uu[1]*Au[1] - uu[2]*Au[2];
		f[4] = uu[2]*uu[2]*(1.0 + K13) + 2.0/rho*A_d2[0].Get_P(0,Au[2])
			- uu[0]*uu[0] - 2.0/rho*A_p[nel_p-1].Get_P(A_p[nel_p-1].q-1,Au[0]);
		f[5] = uu[2]*uu[2]*(1.0 + K23) + 2.0/rho*A_d2[0].Get_P(0,Au[2])
			- uu[1]*uu[1] - 2.0/rho*A_d1[0].Get_P(0,Au[1]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  0;
		dfdv[0][3] =  A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0])/Au[0];
		dfdv[0][4] =  0;
		dfdv[0][5] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  0;
		dfdv[1][4] =  -A_d1[0].Get_c(0,Au[1])/Au[1];
		dfdv[1][5] =  0;

		dfdv[2][0] =  0;
		dfdv[2][1] =  0;
		dfdv[2][2] =  1;
		dfdv[2][3] =  0;
		dfdv[2][4] =  0;
		dfdv[2][5] =  -A_d2[0].Get_c(0,Au[2])/Au[2];

		dfdv[3][0] =  Au[0];
		dfdv[3][1] = -Au[1];
		dfdv[3][2] = -Au[2];
		dfdv[3][3] =  uu[0];
		dfdv[3][4] = -uu[1];
		dfdv[3][5] = -uu[2];

		dfdv[4][0] = -2*uu[0];
		dfdv[4][1] =  0;
		dfdv[4][2] =  2*uu[2]*(1.0 + K13);
		dfdv[4][3] = -2.0*pow(A_p[nel_p-1].Get_c(A_p[nel_p-1].q-1,Au[0]),2)/Au[0];
		dfdv[4][4] =  0;
		dfdv[4][5] =  2.0*pow(A_d2[0].Get_c(0,Au[2]),2)/Au[2];

		dfdv[5][0] =  0;
		dfdv[5][1] = -2*uu[1];
		dfdv[5][2] =  2*uu[2]*(1.0 + K23);
		dfdv[5][3] =  0;
		dfdv[5][4] = -2.0*pow(A_d1[0].Get_c(0,Au[1]),2)/Au[1];
		dfdv[5][5] =  2.0*pow(A_d2[0].Get_c(0,Au[2]),2)/Au[2];

		one = 1;
		six = 6;

		/*  Solve the linear system by factorising dfdv into LU by dgetrf and solving the systems Ly=f and
		U(dx)=y	by dgetrs. dx is stored in f. */
		dgetrf(&six, &six, *dfdv, &six, ipiv, &info); // LU factorisation
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_leav: dgetrf failed at iter %d. \n",iter); exit(1);}

		dgetrs("T", &six, &one, *dfdv, &six, ipiv, f, &six, &info); // Ly=f and U(dx)=y solution
		if(info){fprintf(stderr,"Error in BifurRiem1to2_losses_leav: dgetrs failed at iter %d. \n",iter); exit(1);}

		// Update solution: x_new = x_old - dx
		uu[0] -= f[0];
		uu[1] -= f[1];
		uu[2] -= f[2];
		Au[0] -= f[3];
		Au[1] -= f[4];
		Au[2] -= f[5];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3] +	f[4]*f[4] + f[5]*f[5]) < ModelParam::riemannTol)
			proceed = 0;

	}
	if(iter >= MAX_ITER)
		fprintf(stderr,"Error in BifurRiem1to2_losses_leav: iteration failed to converge. \n");

	return iter;
}

int OneDSolver::BifurRiem2to1(int nel_d1, int nel_d2, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho){
	double f[6], g[6], inv_J[6][6];
	double c1, c2, c3, k, k1, k2, k3;
	int proceed = 1, iter = 0, six, one, zero;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector.
		f[0] = uu[0] - A_p[0].Get_W(0,Au[0]) - W[0];
		f[1] = uu[1] + A_d1[nel_d1-1].Get_W(A_d1[nel_d1-1].q-1,Au[1]) - W[1];
		f[2] = uu[2] + A_d2[nel_d2-1].Get_W(A_d2[nel_d2-1].q-1,Au[2]) - W[2];
		f[3] = uu[0]*Au[0] - uu[1]*Au[1] - uu[2]*Au[2];
		if(ModelParam::nonDim){
			f[4] = uu[0]*uu[0] + 2.0*A_p[0].Get_P(0,Au[0])
				- uu[1]*uu[1] - 2.0*A_d1[nel_d1-1].Get_P(A_d1[nel_d1-1].q-1,Au[1]);
			f[5] = uu[0]*uu[0] + 2.0*A_p[0].Get_P(0,Au[0])
				- uu[2]*uu[2] - 2.0*A_d2[nel_d2-1].Get_P(A_d2[nel_d2-1].q-1,Au[2]);
		}
		else{
			f[4] = uu[0]*uu[0] + 2.0/rho*A_p[0].Get_P(0,Au[0])
				- uu[1]*uu[1] - 2.0/rho*A_d1[nel_d1-1].Get_P(A_d1[nel_d1-1].q-1,Au[1]);
			f[5] = uu[0]*uu[0] + 2.0/rho*A_p[0].Get_P(0,Au[0])
				- uu[2]*uu[2] - 2.0/rho*A_d2[nel_d2-1].Get_P(A_d2[nel_d2-1].q-1,Au[2]);
		}

		// Wave speed at each vessel
		c1 = A_p[0].Get_c(0,Au[0]);
		c2 = A_d1[nel_d1-1].Get_c(A_d1[nel_d1-1].q-1,Au[1]);
		c3 = A_d2[nel_d2-1].Get_c(A_d2[nel_d2-1].q-1,Au[2]);

		// Inverse Jacobian matrix dfdv
		k = c1*Au[1]*c3+Au[0]*c3*c2+Au[2]*c1*c2;
		k1 = (c1+uu[0])*k;
		inv_J[0][0] = (c2*uu[0]*c3*Au[0]+Au[2]*c2*c1*c1+Au[1]*c1*c1*c3)/k1;
		inv_J[0][1] = Au[1]*(c2+uu[1])*c1*c3/k1;
		inv_J[0][2] = Au[2]*(c3+uu[2])*c1*c2/k1;
		inv_J[0][3] = c1*c2*c3/k1;
		inv_J[0][4] = 0.5*Au[1]*c1*c3/k1;
		inv_J[0][5] = 0.5*Au[2]*c1*c2/k1;

		k2 = (c2-uu[1])*k;
		inv_J[1][0] = Au[0]*(c1-uu[0])*c2*c3/k2;
		inv_J[1][1] = (-c1*uu[1]*c3*Au[1]+Au[2]*c1*c2*c2+c3*c2*c2*Au[0])/k2;
		inv_J[1][2] = -Au[2]*(c3+uu[2])*c1*c2/k2;
		inv_J[1][3] = -c1*c2*c3/k2;
		inv_J[1][4] = 0.5*(c1*Au[2]+Au[0]*c3)*c2/k2;
		inv_J[1][5] = -0.5*Au[2]*c1*c2/k2;

		k3 = (c3-uu[2])*k;
		inv_J[2][0] = Au[0]*(c1-uu[0])*c2*c3/k3;
		inv_J[2][1] = -Au[1]*(c2+uu[1])*c1*c3/k3;
		inv_J[2][2] = -(c1*uu[2]*c2*Au[2]-Au[1]*c1*c3*c3-c2*c3*c3*Au[0])/k3;
		inv_J[2][3] = -c1*c2*c3/k3;
		inv_J[2][4] = -0.5*Au[1]*c1*c3/k3;
		inv_J[2][5] = 0.5*(Au[1]*c1+Au[0]*c2)*c3/k3;

		inv_J[3][0] = -Au[0]*(Au[0]*c3*c2+uu[0]*c3*Au[1]+uu[0]*c2*Au[2])/k1;
		inv_J[3][1] = Au[0]*Au[1]*(c2+uu[1])*c3/k1;
		inv_J[3][2] = Au[0]*Au[2]*(c3+uu[2])*c2/k1;
		inv_J[3][3] = Au[0]*c3*c2/k1;
		inv_J[3][4] = 0.5*Au[0]*Au[1]*c3/k1;
		inv_J[3][5] = 0.5*Au[0]*c2*Au[2]/k1;

		inv_J[4][0] = -Au[0]*Au[1]*(c1-uu[0])*c3/k2;
		inv_J[4][1] = Au[1]*(Au[1]*c1*c3-c1*uu[1]*Au[2]-c3*uu[1]*Au[0])/k2;
		inv_J[4][2] = Au[2]*Au[1]*(c3+uu[2])*c1/k2;
		inv_J[4][3] = Au[1]*c1*c3/k2;
		inv_J[4][4] = -0.5*Au[1]*(c1*Au[2]+Au[0]*c3)/k2;
		inv_J[4][5] = 0.5*Au[2]*Au[1]*c1/k2;

		inv_J[5][0] = -Au[0]*Au[2]*(c1-uu[0])*c2/k3;
		inv_J[5][1] = Au[2]*Au[1]*(c2+uu[1])*c1/k3;
		inv_J[5][2] = Au[2]*(Au[2]*c1*c2-c1*uu[2]*Au[1]-c2*uu[2]*Au[0])/k3;
		inv_J[5][3] = Au[2]*c1*c2/k3;
		inv_J[5][4] = 0.5*Au[2]*Au[1]*c1/k3;
		inv_J[5][5] = -0.5*Au[2]*(Au[1]*c1+Au[0]*c2)/k3;

		one = 1;
		six = 6;
		zero = 0;

		/* Solve the linear system by inverting analyticaly the Jacobian: g = (dfdv)^(-1)*f */
		cblas_dgemv(CblasColMajor,CblasTrans, six, six, one, *inv_J, six, f, one, zero, g, one);

		// Update solution: x_new = x_old - dx
		uu[0] -= g[0];
		uu[1] -= g[1];
		uu[2] -= g[2];
		Au[0] -= g[3];
		Au[1] -= g[4];
		Au[2] -= g[5];

		// Check if the error of the solution is smaller than TOL.
		if((g[0]*g[0] + g[1]*g[1] + g[2]*g[2] + g[3]*g[3] +
			g[4]*g[4] + g[5]*g[5]) < ModelParam::riemannTol)      proceed = 0;
	}

	if(iter >= MAX_ITER)
		fprintf(stderr,"Error in BifurRiem2to1: iteration failed to converge. \n");

	return iter;
}

int OneDSolver::JuncRiemann(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho){
	double f[4], g[4], inv_J[4][4];
	double cl, cr, k, k1, k2;
	int proceed = 1, iter = 0;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector and the wave speed at each vessel.
		f[0] = uu[0] + A_l[nel_l-1].Get_W(A_l[nel_l-1].q-1,Au[0]) - W[0]; // Continuity of the forward left domain Riemann invariant.
		f[1] = uu[1] - A_r[0].Get_W(0,Au[1]) - W[1]; // Continuity of the backward right domain Riemann invariant.
		f[2] = uu[0]*Au[0] - uu[1]*Au[1]; // Conservation of mass.
		if(ModelParam::nonDim){
			f[3] = uu[0]*uu[0] + 2.0*A_l[nel_l-1].Get_P(A_l[nel_l-1].q-1,Au[0]) /* Continuity of the total pressure. */
				- uu[1]*uu[1] - 2.0*A_r[0].Get_P(0,Au[1]);
		}
		else{
			f[3] = uu[0]*uu[0] + 2.0/rho*A_l[nel_l-1].Get_P(A_l[nel_l-1].q-1,Au[0]) /* Continuity of the total pressure. */
				- uu[1]*uu[1] - 2.0/rho*A_r[0].Get_P(0,Au[1]);
		}

		// Wave speed at each vessel
		cl   = A_l[nel_l-1].Get_c(A_l[nel_l-1].q-1,Au[0]);
		cr   = A_r[0].Get_c(0,Au[1]);

		// Inverse Jacobian matrix dfdv
		k = (cl*Au[1]+Au[0]*cr);
		k1 = (cl-uu[0])*k;
		inv_J[0][0] = (Au[1]*cl*cl-cr*uu[0]*Au[0])/k1;
		inv_J[0][1] = Au[1]*(cr-uu[1])*cl/k1;
		inv_J[0][2] = cl*cr/k1;
		inv_J[0][3] = -0.5*cl*Au[1]/k1;

		k2 = (cr+uu[1])*k;
		inv_J[1][0] = Au[0]*(cl+uu[0])*cr/k2;
		inv_J[1][1] = (cl*uu[1]*Au[1]+cr*cr*Au[0])/k2;
		inv_J[1][2] = -cl*cr/k2;
		inv_J[1][3] = -0.5*Au[0]*cr/k2;

		inv_J[2][0] = Au[0]*(Au[0]*cr-uu[0]*Au[1])/k1;
		inv_J[2][1] = -Au[0]*Au[1]*(cr-uu[1])/k1;
		inv_J[2][2] = -Au[0]*cr/k1;
		inv_J[2][3] = 0.5*Au[1]*Au[0]/k1;

		inv_J[3][0] = Au[0]*Au[1]*(cl+uu[0])/k2;
		inv_J[3][1] = -Au[1]*(cl*Au[1]+uu[1]*Au[0])/k2;
		inv_J[3][2] = -cl*Au[1]/k2;
		inv_J[3][3] = -0.5*Au[1]*Au[0]/k2;

		/* Solve the linear system by inverting analyticaly the Jacobian: g = (dfdv)^(-1)*f */
		cblas_dgemv(CblasColMajor,CblasTrans, 4, 4, 1, *inv_J, 4, f, 1, 0, g, 1);

		// Update solution: x_new = x_old - dx
		uu[0] -= g[0];
		uu[1] -= g[1];
		Au[0] -= g[2];
		Au[1] -= g[3];

		// Check if the error of the solution is smaller than TOL.
		if((g[0]*g[0] + g[1]*g[1] + g[2]*g[2] + g[3]*g[3]) < ModelParam::riemannTol)      proceed = 0;
	}

	if(iter >= MAX_ITER){
		fprintf(stderr,"Error in JuncRiemann: iteration failed to converge. \n");
		exit(-1);
	}

	return iter;
}

int OneDSolver::JuncRiemann_losses_cont_left(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho, double K){
	double f[4], dfdv[4][4];
	int proceed = 1, iter = 0, ipiv[4], info;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector
		f[0] = uu[0] + A_l[nel_l-1].Get_W(A_l[nel_l-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_r[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[0]*Au[0] - uu[1]*Au[1];
		f[3] = uu[0]*uu[0]*(1.0 + K) + 2.0/rho*A_l[nel_l-1].Get_P(A_l[nel_l-1].q-1,Au[0])
			- uu[1]*uu[1] - 2.0/rho*A_r[0].Get_P(0,Au[1]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  A_l[nel_l-1].Get_c(A_l[nel_l-1].q-1,Au[0])/Au[0];
		dfdv[0][3] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  -A_r[0].Get_c(0,Au[1])/Au[1];

		dfdv[2][0] =  Au[0];
		dfdv[2][1] = -Au[1];
		dfdv[2][2] =  uu[0];
		dfdv[2][3] = -uu[1];

		dfdv[3][0] =  2*uu[0]*(1.0 + K);
		dfdv[3][1] = -2*uu[1];
		dfdv[3][2] =  2.0*pow(A_l[nel_l-1].Get_c(A_l[nel_l-1].q-1,Au[0]),2)/Au[0];
		dfdv[3][3] = -2.0*pow(A_r[0].Get_c(0,Au[1]),2)/Au[1];

		int four = 4;
		int one  = 1;
		// invert
		dgetrf(&four, &four, *dfdv, &four, ipiv, &info);
		if(info){
			fprintf(stderr,"Error in JuncRiemann_losses_cont_left: dgetrf failed at iter %d \n",iter);
			exit(1);
		}
		dgetrs("T", &four, &one, *dfdv, &four, ipiv, f, &four, &info);
		if(info){
			fprintf(stderr,"Error in JuncRiemann_losses_cont_left: dgetrs failed at iter %d \n",iter);
			exit(1);
		}

		// Update solution
		uu[0] -= f[0];
		uu[1] -= f[1];
		Au[0] -= f[2];
		Au[1] -= f[3];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3]) < ModelParam::riemannTol)
			proceed = 0;
	}

	if(iter >= MAX_ITER)
		fprintf(stderr,"Warning JuncRiemann_losses_cont_left: iteration failed to converge\n");

	return iter;
}

int OneDSolver::JuncRiemann_losses_cont_right(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho, double K){
	double f[4], dfdv[4][4];
	int proceed = 1, iter = 0, ipiv[4], info;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector
		f[0] = uu[0] + A_l[nel_l-1].Get_W(A_l[nel_l-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_r[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[0]*Au[0] - uu[1]*Au[1];
		f[3] = uu[0]*uu[0] + 2.0/rho*A_l[nel_l-1].Get_P(A_l[nel_l-1].q-1,Au[0])
			- uu[1]*uu[1]*(1.0 + K) - 2.0/rho*A_r[0].Get_P(0,Au[1]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  A_l[nel_l-1].Get_c(A_l[nel_l-1].q-1,Au[0])/Au[0];
		dfdv[0][3] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  -A_r[0].Get_c(0,Au[1])/Au[1];

		dfdv[2][0] =  Au[0];
		dfdv[2][1] = -Au[1];
		dfdv[2][2] =  uu[0];
		dfdv[2][3] = -uu[1];

		dfdv[3][0] =  2*uu[0];
		dfdv[3][1] = -2*uu[1]*(1.0 + K);
		dfdv[3][2] =  2.0*pow(A_l[nel_l-1].Get_c(A_l[nel_l-1].q-1,Au[0]),2)/Au[0];
		dfdv[3][3] = -2.0*pow(A_r[0].Get_c(0,Au[1]),2)/Au[1];

		int four = 4;
		int one  = 1;

		// invert
		dgetrf(&four, &four, *dfdv, &four, ipiv, &info);
		if(info){
			fprintf(stderr,"Error in JuncRiemann_losses_cont_right: dgetrf failed at iter %d \n",iter);
			exit(1);
		}
		dgetrs("T", &four, &one, *dfdv, &four, ipiv, f, &four, &info);
		if(info){
			fprintf(stderr,"Error in JuncRiemann_losses_cont_right: dgetrs failed at iter %d \n",iter);
			exit(1);
		}

		// Update solution
		uu[0] -= f[0];
		uu[1] -= f[1];
		Au[0] -= f[2];
		Au[1] -= f[3];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3]) < ModelParam::riemannTol)
			proceed = 0;
	}

	if(iter >= MAX_ITER)
		fprintf(stderr,"Warning JuncRiemann_losses_cont_right: iteration failed to converge\n");

	return iter;
}

int OneDSolver::JuncRiemann_losses_exp_left(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho, double K){
	double f[4], dfdv[4][4];
	int proceed = 1, iter = 0, ipiv[4], info;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector
		f[0] = uu[0] + A_l[nel_l-1].Get_W(A_l[nel_l-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_r[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[0]*Au[0] - uu[1]*Au[1];
		f[3] = uu[0]*uu[0] + 2.0/rho*A_l[nel_l-1].Get_P(A_l[nel_l-1].q-1,Au[0])
			- uu[1]*uu[1]*(1.0 - K) - 2.0/rho*A_r[0].Get_P(0,Au[1]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  A_l[nel_l-1].Get_c(A_l[nel_l-1].q-1,Au[0])/Au[0];
		dfdv[0][3] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  -A_r[0].Get_c(0,Au[1])/Au[1];

		dfdv[2][0] =  Au[0];
		dfdv[2][1] = -Au[1];
		dfdv[2][2] =  uu[0];
		dfdv[2][3] = -uu[1];

		dfdv[3][0] =  2*uu[0];
		dfdv[3][1] = -2*uu[1]*(1.0 - K);
		dfdv[3][2] =  2.0*pow(A_l[nel_l-1].Get_c(A_l[nel_l-1].q-1,Au[0]),2)/Au[0];
		dfdv[3][3] = -2.0*pow(A_r[0].Get_c(0,Au[1]),2)/Au[1];

		int four = 4;
		int one  = 1;

		// invert
		dgetrf(&four, &four, *dfdv, &four, ipiv, &info);
		if(info){
			fprintf(stderr,"Error in JuncRiemann_losses_exp_left: dgetrf failed at iter %d \n",iter);
			exit(1);
		}
		dgetrs("T", &four, &one, *dfdv, &four, ipiv, f, &four, &info);
		if(info){
			fprintf(stderr,"Error in JuncRiemann_losses_exp_left: dgetrs failed at iter %d \n",iter);
			exit(1);
		}

		// Update solution
		uu[0] -= f[0];
		uu[1] -= f[1];
		Au[0] -= f[2];
		Au[1] -= f[3];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3]) < ModelParam::riemannTol)
			proceed = 0;
	}

	if(iter >= MAX_ITER)
		fprintf(stderr,"Warning JuncRiemann_losses_exp_left: iteration failed to converge\n");

	return iter;
}

int OneDSolver::JuncRiemann_losses_exp_right(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho, double K){
	double f[4], dfdv[4][4];
	int proceed = 1, iter = 0, ipiv[4], info;

	while((proceed)&&(iter++ < MAX_ITER)){
		// Calculate constraint vector
		f[0] = uu[0] + A_l[nel_l-1].Get_W(A_l[nel_l-1].q-1,Au[0]) - W[0];
		f[1] = uu[1] - A_r[0].Get_W(0,Au[1]) - W[1];
		f[2] = uu[0]*Au[0] - uu[1]*Au[1];
		f[3] = uu[0]*uu[0]*(1.0 - K) + 2.0/rho*A_l[nel_l-1].Get_P(A_l[nel_l-1].q-1,Au[0])
			- uu[1]*uu[1] - 2.0/rho*A_r[0].Get_P(0,Au[1]);

		dfdv[0][0] =  1;
		dfdv[0][1] =  0;
		dfdv[0][2] =  A_l[nel_l-1].Get_c(A_l[nel_l-1].q-1,Au[0])/Au[0];
		dfdv[0][3] =  0;

		dfdv[1][0] =  0;
		dfdv[1][1] =  1;
		dfdv[1][2] =  0;
		dfdv[1][3] =  -A_r[0].Get_c(0,Au[1])/Au[1];

		dfdv[2][0] =  Au[0];
		dfdv[2][1] = -Au[1];
		dfdv[2][2] =  uu[0];
		dfdv[2][3] = -uu[1];

		dfdv[3][0] =  2*uu[0]*(1.0 - K);
		dfdv[3][1] = -2*uu[1];
		dfdv[3][2] =  2.0*pow(A_l[nel_l-1].Get_c(A_l[nel_l-1].q-1,Au[0]),2)/Au[0];
		dfdv[3][3] = -2.0*pow(A_r[0].Get_c(0,Au[1]),2)/Au[1];

		int four = 4;
		int one  = 1;

		// invert
		dgetrf(&four, &four, *dfdv, &four, ipiv, &info);
		if(info){
			fprintf(stderr,"Error in JuncRiemann_losses_exp_right: dgetrf failed at iter %d \n",iter);
			exit(1);
		}
		dgetrs("T", &four, &one, *dfdv, &four, ipiv, f, &four, &info);
		if(info){
			fprintf(stderr,"Error in JuncRiemann_losses_exp_right: dgetrs failed at iter %d \n",iter);
			exit(1);
		}

		// Update solution
		uu[0] -= f[0];
		uu[1] -= f[1];
		Au[0] -= f[2];
		Au[1] -= f[3];

		if((f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3]) < ModelParam::riemannTol)
			proceed = 0;
	}

	if(iter >= MAX_ITER)
		fprintf(stderr,"Warning JuncRiemann_losses_exp_right: iteration failed to converge\n");

	return iter;
}

/**
\brief
According to the tube law considered, it evaluates the area bc in order to obtain the desire "A_star"
taking into account "A" at the other side of the bc.
*/
double OneDSolver::eval_BC_A(double A_star, double A){
	double BC_A,tmp;
	// BC_A = pow(2.0*pow(A_star,0.25) - pow(A,0.25),4);
	tmp = 2*sqrt(sqrt(A_star))-sqrt(sqrt(A));
	BC_A = tmp*tmp*tmp*tmp;
	return BC_A;
}

/**
\brief
According to the tube law considered, it calculates the area corresponding to the input pressure.
*/
double OneDSolver::eval_A_from_P(double p, double po, double beta, double Ao){
	double A, tmp;
	// A = pow((p-po)/beta + pow(Ao,0.5),2);
	tmp = (p-po)/beta+sqrt(Ao);
	A = tmp*tmp;
	return A;
}

int OneDSolver::WKinflow(Element *A, double u_right, double A_right, double Q, double *uu, double *Au){
	int proceed = 1, iter = 0;
	double W, fa, dfa, A_calc, delta_A_calc;

	/* Tolerances for the algorithm */
	double      tol = 1.0e-10;
	double    SMALL = 1.0e-200;

	/* Riemann invariant W2(Ar,ur). */
	W = u_right - A[0].Get_W(0,A_right);

	/* Ensure that flow is subsonic. */
	if(u_right > A[0].Get_c(0,A_right)){
		fprintf(stderr,"Error in WKinflow: flow is not subsonic.\n");
		return iter = -1;
	}

	/* Newton's iteration (Area only). */
	A_calc = A_right;

	while( (proceed) && (iter++ < MAX_ITER) ){
		fa  =  Q - W*A_calc - A_calc*A[0].Get_W(0,A_calc);

		dfa =  -W - A[0].Get_W(0,A_calc) - A[0].Get_c(0,A_calc);

		if( fabs(dfa) < SMALL ){
			fprintf(stderr,"Error in WKinflow: zero derivative in the Newton's iteration.\n");
			return iter = -2;
		}

		delta_A_calc = fa/dfa;
		if( fabs(delta_A_calc) < tol ) proceed = 0;
		A_calc -= delta_A_calc;
		// fprintf(stdout,"iter: %d fa: %lf\n",iter, fa);
		// fprintf(stdout,"A: %lf  dA: %lf\n",A_calc,delta_A_calc);
	}

	/* Obtain u from W2 and A_calc. */
	*uu = W + A[0].Get_W(0,A_calc);
	*Au = A_calc;

	//printf("uu = %lg A_calc = %lg dfa = %lg dfa = %lg\n",u,A_calc,fa,dfa);

	if(iter >= MAX_ITER)
		fprintf(stderr,"Error in WKinflow: iteration failed to converge.\n");

	return iter;
}

int OneDSolver::WKoutflow(int nel, Element *A, double u_left, double A_left, double p_ext, double p_inf, double R, double rho, double *uu, double *Au){
	int proceed = 1, iter = 0;
	double W, fa, dfa, A_calc, delta_A_calc;

	/* Tolerances for the algorithm */
	double      tol = 1.0e-10;
	double    SMALL = 1.0e-200;

	/* Riemann invariant W1(Al,ul). */
	W = u_left + A[nel-1].Get_W(A[nel-1].q-1,A_left);

	/* Ensure that it is outflow. */
	if(u_left > A[nel-1].Get_c(A[nel-1].q-1,A_left)){
		fprintf(stderr,"Error in WKoutflow: flow is not subsonic.\n");
		return iter = -1;
	}

	/* Newton's iteration (Area only). */
	A_calc = A_left;

	while( (proceed) && (iter++ < MAX_ITER) ){
		fa  =  R*(W*A_calc - A_calc*A[nel-1].Get_W(A[nel-1].q-1,A_calc))
			- A[nel-1].Get_P(A[nel-1].q-1,A_calc) - p_ext + p_inf;

		dfa =  R*(W - A[nel-1].Get_W(A[nel-1].q-1,A_calc) - A[nel-1].Get_c(A[nel-1].q-1,A_calc))
			- rho*pow(A[nel-1].Get_c(A[nel-1].q-1,A_calc),2)/A_calc;

		if( fabs(dfa) < SMALL ){
			fprintf(stderr,"Error in WKoutflow: zero derivative in the Newton's iteration.\n");
			return iter = -2;
		}

		delta_A_calc = fa/dfa;
		// if( fabs(delta_A_calc) < tol ) proceed = 0;
		if( fabs(delta_A_calc) < ModelParam::riemannTol ) proceed = 0;
		A_calc -= delta_A_calc;
		// fprintf(stdout,"iter: %d fa: %lf\n",iter, fa);
		// fprintf(stdout,"A: %lf  dA: %lf\n",A_calc,delta_A_calc);
	}

	/* Obtain u from W1 and A_calc. */
	*uu = W - A[nel-1].Get_W(A[nel-1].q-1,A_calc);
	*Au = A_calc;

	//printf("uu = %lg A_calc = %lg dfa = %lg dfa = %lg\n",u,A_calc,fa,dfa);

	if(iter >= MAX_ITER)
		fprintf(stderr,"Error in WKoutflow: iteration failed to converge.\n");

	return iter;
}

/**
\brief 
Update the viscosity according to diameter of the current step and the fixed Hd
*/
void OneDSolver::UpdateVisc(Domain *omega){
	int i,j,k,n;
	Element *A, *U, *Hd;
	double *visc, *oldAhis, *oldUhis, *oldHdhis;

	for(n=0; n<ModelParam::Ndoms; ++n){
		A = omega[n].A;
		U = omega[n].U;
		Hd = omega[n].Hd;
		visc = omega[n].visc;

		if(ModelParam::nonDim){
			for(i=0; i<A[0].q; ++i){
				if(ModelParam::updVisc==1)
					visc[i] = Eta_vitro(Hd[0].h[i], 2*sqrt(A[0].h[i]*ModelParam::scale_r0*ModelParam::scale_r0/ModelParam::PI)*1e6)/1000;
				else if(ModelParam::updVisc==2)
					visc[i] = Eta_vivo_ESL(Hd[0].h[i], 2*sqrt(A[0].h[i]*ModelParam::scale_r0*ModelParam::scale_r0/ModelParam::PI)*1e6)/1000;
				else if(ModelParam::updVisc==3)
					visc[i] = Eta_vivo(Hd[0].h[i], 2*sqrt(A[0].h[i]*ModelParam::scale_r0*ModelParam::scale_r0/ModelParam::PI)*1e6, 0.85)/1000;
			}
		}
		else{
			for(i=0; i<A[0].q; ++i){
				if(ModelParam::updVisc==1)
					visc[i] = Eta_vitro(Hd[0].h[i], 2*sqrt(A[0].h[i]/ModelParam::PI)*1e6)/1000;
				else if(ModelParam::updVisc==2)
					visc[i] = Eta_vivo_ESL(Hd[0].h[i], 2*sqrt(A[0].h[i]/ModelParam::PI)*1e6)/1000;
				else if(ModelParam::updVisc==3)
					visc[i] = Eta_vivo(Hd[0].h[i], 2*sqrt(A[0].h[i]/ModelParam::PI)*1e6, 0.85)/1000;
			}
		}
	}
}

/**
\brief
Ref. 1994 Resistance to blood flow in microvessels in vivo
*/
double OneDSolver::Eta_vitro(double Hd, double diam){
	double C = Eta_C(diam);
	double visc_045 = Eta_045(diam);
	double Eta_vitro = 1+(visc_045-1)*(pow(1-Hd,C)-1)/(pow(1-0.45,C)-1);
	return Eta_vitro;
}

double OneDSolver::Eta_045(double diam){
	return 220*exp(-1.3*diam)+3.2-2.44*exp(-0.06*pow(diam,0.645));
}

double OneDSolver::Eta_C(double diam){
	return (0.8+exp(-0.075*diam))*(-1+1/(1+pow(double(10),-11)+pow(diam,12)))+1/(1+pow(double(10),-11)+pow(diam,12));
}

double OneDSolver::Eta_Was(double diam){
	double Doff = 2.4;	// um
	double D50 = 100;		// um
	double Wmax = 2.6;	// um, depends on the data

	if(diam > Doff)
		return (diam-Doff)/(diam+D50-2*Doff)*Wmax;
	else
		return 0;
}

double OneDSolver::Eta_Wpeak(double diam){
	double Doff = 2.4;	// um
	double Dcrit = 10.5	;	// um
	double Eamp = 1.1;
	double Ewidth = 0.03;
	double Ehd = 1.18;

	if(diam>Dcrit)
		return Eamp*exp(-Ewidth*(diam-Dcrit));
	else if(diam>Doff && diam<=Dcrit)
		return Eamp*(diam-Doff)/(Dcrit-Doff);
	else
		return 0;
}

/*
Ref. 1994 Resistance to blood flow in microvessels in vivo
*/
double OneDSolver::Eta_vivo(double Hd, double diam, double W){
	double C=(0.8+exp(-0.075*diam))*(-1+1/(1+10e-11*pow(diam,12)))+1/(1+10e-11*pow(diam,12));
	double eta_45=6*exp(-0.085*diam)+3.2-2.44*exp(-0.06*pow(diam,0.645));
	return (1+(eta_45-1)*(pow(1-Hd,C)-1)/(pow(1-0.45,C)-1)*pow(diam/(diam-W),2))*pow(diam/(diam-W),2);
}

/*
Ref. 2005 Microvascular blood viscosity in vivo and the endothelial surface layer
*/
double OneDSolver::Eta_vivo_ESL(double Hd, double diam){
	double Deff = Eta_Deff(Hd, diam);
	double Dph = Eta_Dph(diam);
	double eta_vitro = Eta_vitro(Hd, Dph);
	return eta_vitro*pow(diam/Deff, 4);
}

double OneDSolver::Eta_Deff(double Hd, double diam){
	double Ehd = 1.18;
	double Weff = Eta_Was(diam)+Eta_Wpeak(diam)*(1+Hd*Ehd);
	return diam-2*Weff;
}

double OneDSolver::Eta_Dph(double diam){
	double Epeak = 0.6;
	double Wph = Eta_Was(diam)+Eta_Wpeak(diam)*Epeak;
	return diam-2*Wph;
}

void OneDSolver::Write_history(Domain *omega, char *name){
	register int i,j,k,n;
	int      nel, eid;
	char     buf[BUFSIZ];
	double   *z,*w;
	Element  *U,*A,**Af;
	static   double ***His_interp;
	static   int   **His_elmt;
	double   zint;
	double   W1,W2,W,c,Ahis,Uhis,phis,dQdxhis,Hdhis,Qhis,Vischis;
	double   uavg,aavg,pavg,len; // Average values
	double   tmp_1[MAX_Q], tmp_2[MAX_Q], tmp_3[MAX_Q];

#ifdef _WIN32
	double now = clock();
	if(ModelParam::showStepLapse)
		printf("step %f, time=%f\n",curTime, (now-lastClock)/1000);
	lastClock = now;
#endif

	if(!fp_his) { // The first time the routine is called, allocate space for "fp", "His_elmt" and "His_interp".
		fp_his     = (FILE **)   calloc(ModelParam::Ndoms,sizeof(FILE *)); // Open a .his file.
		His_elmt   = (int **)    calloc(ModelParam::Ndoms,sizeof(int *));
		His_interp = (double ***)calloc(ModelParam::Ndoms,sizeof(double **));
	}

	for(n=0; n < ModelParam::Ndoms; ++n){ // For each domain,
		if(omega[n].hispts){		/* Check if in the domain "n" a history point solution is required
									in the input list. */
			nel = omega[n].nel;
			U   = omega[n].U;
			A   = omega[n].A;

			if(!fp_his[n]){ /* Give the name of the .in file to the .his file (plus the number of the domain if
							there are more than one). */
				if(ModelParam::Ndoms == 1) sprintf(buf,"%s.his",strtok(name,"."));
				else			 sprintf(buf,"%s_%d.his",strtok(name,"."),n+1);

				//// Create a pointer to the .out file.
				//if(ModelParam::outMode)
				//  fp_his[n] = fopen(buf,"wb");
				//else
				//  fp_his[n] = fopen(buf,"w");

				// Set up the interpolation matrices.
				His_interp[n] = (double **) malloc(omega[n].hispts*sizeof(double *));
				His_elmt[n]   = ivector(0,omega[n].hispts-1);
				ifill(omega[n].hispts,-1,His_elmt[n],1);

				for(i = 0; i < omega[n].hispts; ++i){ /* For all the history points of the domain "n", */
					for(j = 0; j < nel; ++j) // For all the elements of the domain "n",
						if((omega[n].hisptsx[i] >= U[j].x[0])&&(omega[n].hisptsx[i] <= U[j].x[1])){ // Check if the history point location is inside the element "j"
							// q = U[j].q;
							His_elmt  [n][i] = j;
							His_interp[n][i] = dvector(0,ModelParam::oneD_q-1);
							zint = 2*(omega[n].hisptsx[i]-U[j].x[0])/(U[j].x[1]-U[j].x[0])-1;
							getzw(ModelParam::oneD_q,&z,&w,'a');
							for(k = 0; k < ModelParam::oneD_q; ++k)  His_interp[n][i][k] = hgll(k,zint,z,ModelParam::oneD_q);
							break;
						}

						// If the history point is not inside the element "j" write an error message.
						if(His_elmt[n][i] == -1){
							fprintf(stderr,"Error in Write_history: History point %d (x=%lf) is not in domain %d\n", i+1, omega[n].hisptsx[i], n+1);
							--i;
							--omega[n].hispts;
							exit(-1);
						}
				}

				if(!ModelParam::outMode){
					if(!headerPrinted[n]){
						fp_his[n] = fopen(buf,"w");
						// Print the header of the .his file.
						fprintf(fp_his[n],"# 1D nonlinear hp code outfile \n");
						fprintf(fp_his[n],"# History points: %d \n", omega[n].hispts);
						for(i = 0; i < omega[n].hispts; ++i) // For all the history points of the domain "n",
							fprintf(fp_his[n],"# Point %d: x = %lf in elmt %d \n", i+1, omega[n].hisptsx[i], His_elmt[n][i]+1);

						fprintf(fp_his[n],"# t, P(x,t), U(x,t), Q(x,t), A(x,t), Visc(x,t), Wf, Wb, P forw, P backw, U forw, U backw # point\n");
						headerPrinted[n] = 1;
					}
					else
						fp_his[n] = fopen(buf,"a");
				}
			}

			// Calculate average values
			len = uavg = aavg = pavg = 0.0;
			for(j = 0; j < nel; ++j){
				// q = U[j].q;
				getzw(ModelParam::oneD_q,&z,&w,'a');
				uavg += cblas_ddot(ModelParam::oneD_q,w,1,U[j].h,1)*U[j].jac;
				aavg += cblas_ddot(ModelParam::oneD_q,w,1,A[j].h,1)*A[j].jac;
				for(i = 0; i < ModelParam::oneD_q; ++i) pavg += w[i]*A[j].jac*A[j].Get_P(i);
				len += 2*U[j].jac;
			}
			uavg /= len;
			aavg /= len;
			pavg /= len;

			double *tmp_1 = dvector(0,ModelParam::oneD_MaxQ-1);
			double *tmp_2 = dvector(0,ModelParam::oneD_MaxQ-1);

			for(i = 0; i < omega[n].hispts; ++i){ // For all the history points of the domain "n",
				eid = His_elmt[n][i];

				// Evaluate A & u at the history point by means of the interpolation
				if(ModelParam::nonDim){
					Ahis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,omega[n].A[eid].h,1)*(ModelParam::scale_r0*ModelParam::scale_r0);
					Uhis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,omega[n].U[eid].h,1)*ModelParam::scale_u0;
				}
				else{
					Ahis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,omega[n].A[eid].h,1);
					Uhis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,omega[n].U[eid].h,1);
				}
				Hdhis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,omega[n].Hd[eid].h,1);
				Vischis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,omega[n].visc,1);
				Qhis = Ahis*Uhis;

				if(ModelParam::gammaI){
					// Evaluate the viscous part of pressure
					vdmul(&ModelParam::oneD_q,omega[n].U[eid].h,omega[n].A[eid].h,tmp_1); //tmp_1 = A*u
					// dvmul(ModelParam::oneD_q,omega[n].U[eid].h,1,omega[n].A[eid].h,1,tmp_1,1);
					cblas_dgemv(CblasColMajor,CblasTrans,ModelParam::oneD_q,ModelParam::oneD_q,1,*m_d,ModelParam::oneD_q,tmp_1,1,0.0,tmp_3,1); //tmp_3 = d(tmp_1)/dx
					vdmul(&ModelParam::oneD_q,tmp_3,omega[n].gamma,tmp_3);	//tmp_3 = tmp_3*gamma
					// dvmul(ModelParam::oneD_q,tmp_3,1,omega[n].gamma,1,tmp_3,1);
					vdsqrt(&ModelParam::oneD_q,A[eid].h,tmp_2); //tmp_2 = sqrt(A)
					// dvsqrt(ModelParam::oneD_q,A[eid].h,1,tmp_2,1);
					vddiv(&ModelParam::oneD_q,tmp_3,tmp_2,tmp_1); //tmp_1 = tmp_3 / tmp_2
					// dvdiv(ModelParam::oneD_q,tmp_3,1,tmp_2,1,tmp_1,1);
					dQdxhis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,tmp_1,1);
				}
				else if(ModelParam::gammaI){
					// Evaluate the viscous part of pressure
					vdmul(&ModelParam::oneD_q,omega[n].U[eid].h,omega[n].A[eid].h,tmp_1); //tmp_1 = A*u
					// dvmul(ModelParam::oneD_q,omega[n].U[eid].h,1,omega[n].A[eid].h,1,tmp_1,1);
					cblas_dgemv(CblasColMajor,CblasTrans,ModelParam::oneD_q,ModelParam::oneD_q,omega[n].gamma[0],*m_d,ModelParam::oneD_q,tmp_1,1,0.0,tmp_3,1); //tmp_3 = d(tmp_1)/dx
					// vdmul(&ModelParam::oneD_q,tmp_3,omega[n].gamma,tmp_3);	//tmp_3 = tmp_3*gamma
					dQdxhis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,tmp_3,1);
				}

				// Evaluate W at the history point by means of the interpolation
				A[eid].Get_W(tmp_1);
				W = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,tmp_1,1);
				// Work out the forward and backward Riemann invariants
				W1 = Uhis + W;
				W2 = Uhis - W;
				// Evaluate c at the history point by means of the interpolation
				A[eid].Get_c(tmp_1);
				c = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,tmp_1,1);
				// Evaluate p at the history point by means of the interpolation
				A[eid].Get_P(tmp_1);
				if(ModelParam::nonDim)
					phis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,tmp_1,1)*ModelParam::rho*ModelParam::scale_u0*ModelParam::scale_u0;
				else
					phis = cblas_ddot(ModelParam::oneD_q,His_interp[n][i],1,tmp_1,1);
				// Dump the solution in the .his file
				if(!ModelParam::gammaI && !ModelParam::gammaII){
					if(ModelParam::outMode){
						// Write binary data
						int curPos = i+1;
						fwrite(&curTime, sizeof(double), 1, fp_his[n]);
						fwrite(&phis, sizeof(double), 1, fp_his[n]);
						fwrite(&Uhis, sizeof(double), 1, fp_his[n]);
						fwrite(&Qhis, sizeof(double), 1, fp_his[n]);
						fwrite(&Ahis, sizeof(double), 1, fp_his[n]);
						fwrite(&curPos, sizeof(int), 1, fp_his[n]);
					}
					else{
						// Write txt data
						fprintf(fp_his[n],"%lg %lg %lg %lg %lg %lg %d\n", 
							curTime, phis, Uhis, Ahis*Uhis, Ahis, Vischis, i+1);
						/*fprintf(fp_his[n],"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d\n", 
						curTime, phis, Uhis, Ahis*Uhis, Ahis, Vischis, W1, W2, g_rho*c*W1/2.0, -g_rho*c*W2/2.0, W1/2.0, W2/2.0, i+1);*/
					}
				}
				else{
					fprintf(fp_his[n],"%lg %lg %lg %lg %lg %lg %d\n", 
						curTime, phis-dQdxhis, Uhis, Ahis*Uhis, Ahis, Vischis, i+1);
					/*fprintf(fp_his[n],"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d\n", 
					curTime, phis-dQdxhis, Uhis, Ahis*Uhis, Ahis, Vischis, W1, W2, g_rho*c*W1/2.0, -g_rho*c*W2/2.0, W1/2.0, W2/2.0, i+1);*/
					// printf("phis=%lf, dQdxhis=%lf\n", phis, dQdxhis);
				}
			}
		}
		// fflush(fp_his[n]);
		fclose(fp_his[n]);
		fp_his[n]=0;
	}
}

double OneDSolver::CFL(Domain *omega){
	register int i,j,n;
	int    nel;
	double cfl = 0.0, c, dx;
	double *z,*w;
	Element *U, *A;

	for(n = 0; n < ModelParam::Ndoms; ++n){   // For all domains
		nel = omega[n].nel;
		U   = omega[n].U;
		A   = omega[n].A;

		for(i = 0; i < nel; ++i){
			getzw(ModelParam::oneD_q,&z,&w,'a'); /* Determine the location of the quadrature points of each
												 element, z, to later worked out dx. */
			for(j = 0; j < U[i].q; ++j){

				c = A[i].Get_c(j);
				if(j) dx = (U[i].x[1]-U[i].x[0])*(z[j+1]-z[j])/2;
				else  dx = (U[i].x[1]-U[i].x[0])*(z[j]-z[j-1])/2;
				cfl = max(cfl,ModelParam::dt*(fabs(U[i].h[j]) + c)/dx);
			}
		}
	}

	return cfl;
}

/*
Compute the Jacobian matrix of the right-hand side
TODO: For 1D model, it takes very long time to converge
*/
int OneDSolver::Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, 
					N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	int n=(int)N;
	/* precisions for jacobi_matrix calculation */
	double eps = 1e-15;
	/* jacobi-matrix solver handle */
	_JACOBIMATRIX_HANDLE_t handle;
	/* controls of rci cycle */
	MKL_INT successful, rci_request, i;

	double *f1, *f2;

	f1=NV_DATA_S(f1_Vec);
	f2=NV_DATA_S(f2_Vec);

	if (djacobi_init (&handle, &n, &n, NV_DATA_S(y), J->data, &eps) != TR_SUCCESS){
		/* if function does not complete successfully then print error message */
		printf ("\n#fail: error in djacobi_init\n"); fflush (0);
		MKL_FreeBuffers();
		return 1;
	}

	/* set initial rci cycle variables */
	rci_request = 0;
	successful  = 0;
	/* rci cycle */
	while (successful == 0) {
		/* call solver */
		if (djacobi_solve (&handle, f1, f2, &rci_request) != TR_SUCCESS){
			/* if function does not complete successfully then print error message */
			printf ("\n#fail: error in djacobi_solve\n"); fflush (0);
			MKL_FreeBuffers();
			return 1;
		}
		if (rci_request == 1) {
			/* calculate function value f1 = f(x+eps) */
			// f (&m, &n, y->content, f1);
			CVODE_RHS(t, y, f1_Vec, NULL);
		} else if (rci_request == 2) {
			/* calculate function value f1 = f(x-eps) */
			// f (&m, &n, y->content, f2);
			CVODE_RHS(t, y, f2_Vec, NULL);
		} else if (rci_request == 0)
			/* exit rci cycle */
			successful = 1;
	} /* rci cycle */
	/* free handle memory */
	if (djacobi_delete (&handle) != TR_SUCCESS) {
		/* if function does not complete successfully then print error message */
		printf ("\n#fail: error in djacobi_delete\n"); fflush (0);
		MKL_FreeBuffers();
		return 1;
	}
	MKL_FreeBuffers();
	return 0;
}