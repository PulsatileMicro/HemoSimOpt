#include "ZeroDSolver.h"
#include "ModelParam.h"
#include "polylib.h"
#include "PreProcessor.h"
#include <string.h>
#include <fstream>
#include <iostream>
using namespace std;

// CVODE vars
N_Vector ZeroDSolver::init_Y;
void *ZeroDSolver::cvode_mem=NULL;
realtype ZeroDSolver::tret;
int ZeroDSolver::flag;
int ZeroDSolver::reIntFlag;
double *ZeroDSolver::cvode_p=NULL;
double *ZeroDSolver::cvode_q=NULL;
double *ZeroDSolver::cvode_yp=NULL;
double *ZeroDSolver::cvode_yq=NULL;
double *ZeroDSolver::cvode_ydotp=NULL;
double *ZeroDSolver::cvode_ydotq=NULL;

char ZeroDSolver::buf[BUFSIZE];

FILE **ZeroDSolver::fp_his;
clock_t ZeroDSolver::start_time;
clock_t ZeroDSolver::finish;
double ZeroDSolver::duration;
double ZeroDSolver::lastClock;
double ZeroDSolver::curTime;

int *ZeroDSolver::headerPrinted=NULL;

N_Vector ZeroDSolver::f1_Vec;
N_Vector ZeroDSolver::f2_Vec;

/*!
\brief Initialize the ODE Solver (CVODE module)
Solver Type:
1. Explicit: Adams-Bashforth
2. Implicit: BDF
*/
void ZeroDSolver::initSolver(int solverType){
	if(solverType==ModelParam::RLC_EXP)
		cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
	else if(solverType==ModelParam::RLC_IMP)
		cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	init_Y = N_VNew_Serial(ModelParam::Ndoms+ModelParam::Nnodes);
	cvode_q = NV_DATA_S(init_Y);
	cvode_p = cvode_q + ModelParam::Ndoms;

	headerPrinted = (int *)malloc(ModelParam::Ndoms*sizeof(int));
	for (int i=0;i<ModelParam::Ndoms;++i){
		headerPrinted[i]=0;
	}

	f1_Vec = N_VNew_Serial(ModelParam::Nnodes+ModelParam::Ndoms);
	f2_Vec = N_VNew_Serial(ModelParam::Nnodes+ModelParam::Ndoms);
}

/*!
\brief Solve the ODE
*/
void ZeroDSolver::solve(){
	for (int i=0; i<ModelParam::nSteps; ++i){
		if(!i){
			CVodeInit(cvode_mem, CVODE_RHS, 0, init_Y);					// Init the ODE solver
			CVodeSStolerances(cvode_mem, ModelParam::relTol, ModelParam::absTol);		// Set the absolute & relative tolerance of the ODE solver
			// CVodeSetMaxNumSteps(cvode_mem, 50000);								// Set the maximum iteration steps
			// CVodeSetMaxOrd(cvode_mem, 5);
		}
		else
			CVodeReInit(cvode_mem, 0, init_Y);

		CVLapackDense(cvode_mem, ModelParam::Ndoms+ModelParam::Nnodes);
		// CVDlsSetDenseJacFn(cvode_mem, Jac);

		while(1){
			if(CV_SUCCESS == CVode(cvode_mem,ModelParam::dt/ModelParam::scale_lamda*ModelParam::scale_u0,init_Y,&tret,CV_NORMAL)){
				reIntFlag=0;
				break;
			}
		}

		curTime = (i+1)*ModelParam::dt;
		/* It is checked if the solution in the .his file can be written according to "hisstep". If it
		can, Write_history does the job. */
		if(ModelParam::hisSteps&&(!((i+1)%ModelParam::hisSteps))){
			Write_history(ModelParam::omega,ModelParam::argv[ModelParam::argc-1]);
		}
	}
}

/*!
\brief Destroy the ODE Solver
Free CVODE memory
*/
void ZeroDSolver::destroySolver(){
	CVodeFree(&cvode_mem);
}

int ZeroDSolver::CVODE_RHS(realtype t, N_Vector y, N_Vector ydot, void *user_data){
	/* Get the pointer to the data of N_Vector struct */
	cvode_yq = NV_DATA_S(y);
	cvode_yp = cvode_yq + ModelParam::Ndoms;
	cvode_ydotq = NV_DATA_S(ydot);
	cvode_ydotp = cvode_ydotq + ModelParam::Ndoms;

	Eval_RHS(ModelParam::argc, ModelParam::argv); // use y to compute ydot
	reIntFlag=1;

	return 0;
}

/*!
\brief Evaluate the right hand side of the ODE of the Bond Graph Model
*/
void ZeroDSolver::Eval_RHS(int argc, char *argv[]){
	int i, j, n;
	int nel;
	int start, end, from, to;
	int d1, d2, d3, d4;
	double Capacitor, Resistance;
	static FILE  **fIN, **fOUT; // Points to the inlet and outlet of a domain at the boundary of the network
	static double *q_in;
	char *name = argv[argc-1];

	if(!fIN){ // The first time the routine is called,
		fIN = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Allocate space for "fIN" to read from a .bcs file
	}
	if(!fOUT){ // The first time the routine is called,
		fOUT = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Allocate space for "fOUT" to read from a .bcs file
	}
	if(!q_in){
		q_in = new double[ModelParam::Ndoms];// (double *)calloc(ModelParam::Ndoms,sizeof(double));
	}

	for(n = 0; n < ModelParam::Ndoms; ++n){ // For each domain "n" (n=0,...,Ndoms),
		start = ModelParam::omega[n].nodes[0][0];     // start node number
		end   = ModelParam::omega[n].nodes[1][0];     // end node number
		if(ModelParam::Ndoms>1){
			d1    = ModelParam::omega[n].bifur[0][0];
			d2    = ModelParam::omega[n].bifur[0][1];
			d3    = ModelParam::omega[n].bifur[1][0];
			d4    = ModelParam::omega[n].bifur[1][1];
		}

		// Process the flow
		// If the domain is not a boundary segment
		// If the domain is an outlet
		if(ModelParam::omega[n].bctype[3] == 'R'){
			Resistance = ModelParam::omega[n].bcval[4];
			cvode_ydotq[n] = (cvode_yp[start]-cvode_yq[n]*ModelParam::omega[n].ZeroD_R-Resistance*cvode_yq[n])/ModelParam::omega[n].ZeroD_I;
		}
		else if(ModelParam::omega[n].bctype[3] == 'T'){
			Resistance = 0;
			cvode_ydotq[n] = (cvode_yp[start]-cvode_yq[n]*ModelParam::omega[n].ZeroD_R-Resistance*cvode_yq[n])/ModelParam::omega[n].ZeroD_I;
		}
		else{
			// dQ/dt=(Pin-Pout-RQ)/I
			cvode_ydotq[n] = (cvode_yp[start]-cvode_yp[end]-cvode_yq[n]*ModelParam::omega[n].ZeroD_R)/ModelParam::omega[n].ZeroD_I;
		}

		// Process the pressure
		// dP/dt=(Qin-Qout)/C
		// Start node
		switch(ModelParam::omega[n].bctype[0]){
		  case 'B': // If the start node is a 'B' node, the current segment is a daughter segment in the bifurcation
			  if(ModelParam::omega[d1].nodes[0][0] == start){
				  // If the start node of domain d1 == the start node of n, d1 is the daughter segment, d2 is the parent segment
				  Capacitor = ModelParam::omega[n].ZeroD_C + ModelParam::omega[d1].ZeroD_C;
				  cvode_ydotp[start] = (cvode_yq[d2]-cvode_yq[n]-cvode_yq[d1])/Capacitor;
			  }
			  else if(ModelParam::omega[d2].nodes[0][0] == start){
				  Capacitor = ModelParam::omega[n].ZeroD_C + ModelParam::omega[d2].ZeroD_C;
				  cvode_ydotp[start] = (cvode_yq[d1]-cvode_yq[n]-cvode_yq[d2])/Capacitor;
			  }
			  else{
				  printf("error\n");
			  }
			  break;
		  case 'C':
			  Capacitor = ModelParam::omega[n].ZeroD_C;
			  cvode_ydotp[start] = (cvode_yq[d1]+cvode_yq[d2]-cvode_yq[n])/Capacitor;
			  break;
		  case 'J':
			  Capacitor = ModelParam::omega[n].ZeroD_C;
			  cvode_ydotp[start] = (cvode_yq[d1]-cvode_yq[n])/Capacitor;
			  break;
		  case 'q':     
			  if(!fIN[n]){ // Define the name of the .bcs file for domain n
				  if(ModelParam::Ndoms == 1) 
					  sprintf(buf,"%s_IN.bcs",strtok(name,"."));
				  else 
					  sprintf(buf,"%s_IN_%d.bcs",strtok(name,"."),n+1);

				  if(!(fIN[n]=fopen(buf,"rb"))){   /* If the input file .bcs is not found, print an error message and
												   terminate executing the program. If found, open it to start reading it. */
					  fprintf(stderr,"Error in BCs: the file %s doesn't exist. \n",buf);
					  exit(-1);  /* Standard library function that terminates program execution.
								 It also calls "fclose" for each open output file in order to flush out any buffered output. */
				  }
			  }

			  if(!reIntFlag){
				  fread(&(q_in[n]), sizeof(double), 1, fIN[n]);
			  }

			  Capacitor = ModelParam::omega[n].ZeroD_C;
			  cvode_ydotp[start] = (q_in[n]-cvode_yq[n])/Capacitor;
			  break;
		  default:
			  printf("start default\n");
			  break;
		}

		// End node
		switch (ModelParam::omega[n].bctype[3]){
		  case 'B':
			  Capacitor = ModelParam::omega[d3].ZeroD_C + ModelParam::omega[d4].ZeroD_C;
			  cvode_ydotp[end] = (cvode_yq[n]-cvode_yq[d3]-cvode_yq[d4])/Capacitor;
			  break;
		  case 'C':
			  if(ModelParam::omega[d3].nodes[1][0] == end){
				  // If the end node of domain d3 == the end node of n, d3 is the daughter segment, d4 is the parent segment
				  Capacitor = ModelParam::omega[d4].ZeroD_C;
				  cvode_ydotp[end] = (cvode_yq[n]+cvode_yq[d3]-cvode_yq[d4])/Capacitor;
			  }
			  else if(ModelParam::omega[d4].nodes[1][0] == end){
				  Capacitor = ModelParam::omega[d3].ZeroD_C;
				  cvode_ydotp[end] = (cvode_yq[n]+cvode_yq[d4]-cvode_yq[d3])/Capacitor;
			  }
			  else{
				  printf("error\n");
			  }
			  break;
		  case 'J':
			  Capacitor = ModelParam::omega[d3].ZeroD_C;
			  cvode_ydotp[end] = (cvode_yq[n]-cvode_yq[d3])/Capacitor;
			  break;
		  default:
			  // printf("end default\n");
			  break;
		}
	}
}

void ZeroDSolver::Write_history(Domain *omega, char *name){
	register int i,j,k,n;
	int      nel, eid;
	char     buf[BUFSIZ];

#ifdef _WIN32
	double now = clock();
	if(ModelParam::showStepLapse)
		printf("step %f, time=%f\n",curTime, (now-lastClock)/1000);
	lastClock = now;
#endif

	if(!fp_his) { // The first time the routine is called, allocate space for "fp", "His_elmt" and "His_interp".
		fp_his     = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Open a .his file.
	}

	for(n=0; n<ModelParam::Ndoms; ++n){
		if(!fp_his[n]){ /* Give the name of the .in file to the .his file (plus the number of the domain if
						there are more than one). */
			if(ModelParam::Ndoms==1)
				sprintf(buf,"%s.his",strtok(name,"."));
			else
				sprintf(buf,"%s_%d.his",strtok(name,"."),n+1);

			// Create a pointer to the .out file.
			if(!headerPrinted[n]){
				fp_his[n] = fopen(buf,"w");
				fprintf(fp_his[n],"# 0D blood flow outfile \n");
				fprintf(fp_his[n],"# t, Q(t), P_start(t), P_end(t)\n");
				headerPrinted[n] = 1;
			}
			else
				fp_his[n] = fopen(buf,"a");
		}

		if (omega[n].nodes[1][0]<0){ // boundary node
			//fprintf(fp_his[n],"%lg %lg %lg %lg\n", curTime, cvode_q[n], cvode_p[omega[n].nodes[0][0]], 0);
			fprintf(fp_his[n],"%lg %lg %lg %lg\n", curTime, cvode_q[n], cvode_p[omega[n].nodes[0][0]], cvode_p[omega[n].nodes[1][0]]);
		}
		else
			fprintf(fp_his[n],"%lg %lg %lg %lg\n", curTime, cvode_q[n], cvode_p[omega[n].nodes[0][0]], cvode_p[omega[n].nodes[1][0]]);

		fclose(fp_his[n]);
		fp_his[n]=0;
	}
}

int ZeroDSolver::Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, 
					 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	int n=(int)N;
	/* precisions for jacobi_matrix calculation */
	double eps = 0.00000001;
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