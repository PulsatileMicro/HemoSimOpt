#include "RCSolver.h"
#include "ModelParam.h"
#include "polylib.h"
#include "PreProcessor.h"
#include <string.h>

// CVODE vars
N_Vector RCSolver::init_Y;
void *RCSolver::cvode_mem=NULL;
realtype RCSolver::tret;
int RCSolver::flag;
int RCSolver::reIntFlag;
double *RCSolver::cvode_p=NULL;
double *RCSolver::cvode_yp=NULL;
double *RCSolver::cvode_ydotp=NULL;
double *RCSolver::q=NULL;

char RCSolver::buf[BUFSIZE];
double *RCSolver::q_in=NULL;

FILE **RCSolver::fp_his;
clock_t RCSolver::start_time;
clock_t RCSolver::finish;
double RCSolver::duration;
double RCSolver::lastClock;
double RCSolver::curTime;

int *RCSolver::headerPrinted=NULL;

N_Vector RCSolver::f1_Vec;
N_Vector RCSolver::f2_Vec;

/*!
  \brief Initialize the ODE Solver (CVODE module)
  Solver Type:
    1. Explicit: Adams-Bashforth
    2. Implicit: BDF
*/
void RCSolver::initSolver(int solverType){
	if(solverType==ModelParam::RC_EXP)
		cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
	else if(solverType==ModelParam::RC_IMP)
		cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	init_Y = N_VNew_Serial(ModelParam::Nnodes);
	cvode_p = NV_DATA_S(init_Y);

	headerPrinted = (int *)malloc(ModelParam::Ndoms*sizeof(int));
	for (int i=0;i<ModelParam::Ndoms;++i){
		headerPrinted[i]=0;
	}

  q = (double *)malloc(ModelParam::Ndoms*sizeof(double));

  f1_Vec = N_VNew_Serial(ModelParam::Nnodes);
  f2_Vec = N_VNew_Serial(ModelParam::Nnodes);
}

void RCSolver::solve(){
  for (int i=0; i<ModelParam::nSteps; ++i){
    if(!i){
      CVodeInit(cvode_mem, CVODE_RHS, 0, init_Y);					// Init the ODE solver
      CVodeSStolerances(cvode_mem, ModelParam::relTol, ModelParam::absTol);		// Set the absolute & relative tolerance of the ODE solver
      // CVodeSetMaxNumSteps(cvode_mem, 50000);								// Set the maximum iteration steps
      // CVodeSetMaxOrd(cvode_mem, 5);
    }
    else
      CVodeReInit(cvode_mem, 0, init_Y);

    CVLapackDense(cvode_mem, ModelParam::Nnodes);
    CVDlsSetDenseJacFn(cvode_mem, Jac);
    // CVDense(cvode_mem, ModelParam::Nnodes);
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

void RCSolver::destroySolver(){
  CVodeFree(&cvode_mem);
}

int RCSolver::CVODE_RHS(realtype t, N_Vector y, N_Vector ydot, void *user_data){
  /* Get the pointer to the data of N_Vector struct */
  cvode_yp = NV_DATA_S(y);
  cvode_ydotp = NV_DATA_S(ydot);

  Eval_RHS(ModelParam::argc, ModelParam::argv); // use y to compute ydot
  reIntFlag=1;

  return 0;
}

/*!
  \brief Evaluate the right hand side of the ODE of the Bond Graph Model
*/
void RCSolver::Eval_RHS(int argc, char *argv[]){
  int i, j, n;
  int nel;
  int start, end, from, to;
  int d1, d2, d3, d4;
  double Capacitor, Resistance;
  static FILE  **fIN, **fOUT; // Points to the inlet and outlet of a domain at the boundary of the network
  char *name = argv[argc-1];

  if(!fIN){ // The first time the routine is called,
    fIN = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Allocate space for "fIN" to read from a .bcs file
  }
  if(!fOUT){ // The first time the routine is called,
    fOUT = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Allocate space for "fOUT" to read from a .bcs file
  }
  if(!q_in){
    q_in = (double *)malloc(ModelParam::Ndoms*sizeof(double));
  }

  for(n = 0; n < ModelParam::Ndoms; ++n){ // For each domain "n" (n=0,...,Ndoms),
    start = ModelParam::omega[n].nodes[0][0];     // start node number
    end   = ModelParam::omega[n].nodes[1][0];     // end node number
    if(ModelParam::Ndoms>1){ // Not a single vessel
      d1    = ModelParam::omega[n].bifur[0][0];
      d2    = ModelParam::omega[n].bifur[0][1];
      d3    = ModelParam::omega[n].bifur[1][0];
      d4    = ModelParam::omega[n].bifur[1][1];
    }

    // Process the flow
    // If the domain is an outlet
    if(ModelParam::omega[n].bctype[3] == 'R'){
      Resistance = ModelParam::omega[n].bcval[4]; // Boundary Resistance
      // Q=Pin/(Rp+R)
      q[n] = cvode_yp[start]/(Resistance+ModelParam::omega[n].ZeroD_R);
    }
    else if(ModelParam::omega[n].bctype[3] == 'T'){
      q[n] = cvode_yp[start]/(ModelParam::omega[n].ZeroD_R);
    }
    else{
      // Q=(Pin-Pout)/R
      q[n] = (cvode_yp[start]-cvode_yp[end])/ModelParam::omega[n].ZeroD_R;
    }

    // Process the pressure
    // dP/dt=(Qin-Qout)/C
    // Start node
    switch(ModelParam::omega[n].bctype[0]){
      case 'B': // If the start node is a 'B' node, the current segment is a daughter segment in the bifurcation
        if(ModelParam::omega[d1].nodes[0][0] == start){
          // If the start node of domain d1 == the start node of n, d1 is the daughter segment, d2 is the parent segment
          Capacitor = ModelParam::omega[n].ZeroD_C + ModelParam::omega[d1].ZeroD_C;
          cvode_ydotp[start] = (q[d2]-q[n]-q[d1])/Capacitor;
        }
        else if(ModelParam::omega[d2].nodes[0][0] == start){
          Capacitor = ModelParam::omega[n].ZeroD_C + ModelParam::omega[d2].ZeroD_C;
          cvode_ydotp[start] = (q[d1]-q[n]-q[d2])/Capacitor;
        }
        else{
          printf("error\n");
        }
        break;
      case 'C':
        Capacitor = ModelParam::omega[n].ZeroD_C;
        cvode_ydotp[start] = (q[d1]+q[d2]-q[n])/Capacitor;
        break;
      case 'J':
        Capacitor = ModelParam::omega[n].ZeroD_C;
        cvode_ydotp[start] = (q[d1]-q[n])/Capacitor;
        break;
      case 'q':
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
          fread(&(q_in[n]), sizeof(double), 1, fIN[n]);

        Capacitor = ModelParam::omega[n].ZeroD_C;
        cvode_ydotp[start] = (q_in[n]-q[n])/Capacitor;
        break;
      default:
        printf("start default\n");
        break;
    }

    // End node
    switch (ModelParam::omega[n].bctype[3]){
      case 'B':
        Capacitor = ModelParam::omega[d3].ZeroD_C + ModelParam::omega[d4].ZeroD_C;
        cvode_ydotp[end] = (q[n]-q[d3]-q[d4])/Capacitor;
        break;
      case 'C':
        if(ModelParam::omega[d3].nodes[1][0] == end){
          // If the end node of domain d3 == the end node of n, d3 is the daughter segment, d4 is the parent segment
          Capacitor = ModelParam::omega[d4].ZeroD_C;
          cvode_ydotp[end] = (q[n]+q[d3]-q[d4])/Capacitor;
        }
        else if(ModelParam::omega[d4].nodes[1][0] == end){
          Capacitor = ModelParam::omega[d3].ZeroD_C;
          cvode_ydotp[end] = (q[n]+q[d4]-q[d3])/Capacitor;
        }
        else{
          printf("error\n");
        }
        break;
      case 'J':
        Capacitor = ModelParam::omega[d3].ZeroD_C;
        cvode_ydotp[end] = (q[n]-q[d3])/Capacitor;
        break;
      default:
        // printf("end default\n");
        break;
    }
  }
}

void RCSolver::Write_history(Domain *omega, char *name){
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
      fprintf(fp_his[n],"%lg %lg %lg %lg\n", curTime, q[n], cvode_p[omega[n].nodes[0][0]], 0);
    }
    else
      fprintf(fp_his[n],"%lg %lg %lg %lg\n", curTime, q[n], cvode_p[omega[n].nodes[0][0]], cvode_p[omega[n].nodes[1][0]]);

    fclose(fp_his[n]);
    fp_his[n]=0;
  }
}

int RCSolver::Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, 
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