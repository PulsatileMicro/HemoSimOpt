#ifndef ZERODSOLVER_H
#define ZERODSOLVER_H

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <cvode/cvode_lapack.h>      /* prototype for CVLapackDense */
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., and macros */
#include "Domain.h"
#include <time.h>

#define MAX_Q 128
#define BUFSIZE 1024

class ZeroDSolver{
public:
  static void initSolver(int solverType);
  static void solve();
  static void destroySolver();

private:
  static int  CVODE_RHS(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  static void Eval_RHS(int argc, char *argv[]);
  static void Write_history(Domain *omega, char *name);
  static int Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, 
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  // CVODE vars
  static N_Vector init_Y;
  static void *cvode_mem;
  static realtype tret;
  static int flag, reIntFlag;
  static double *cvode_p, *cvode_q, *cvode_yp, *cvode_yq, *cvode_ydotp, *cvode_ydotq;

  static char buf[BUFSIZE];

  static FILE **fp_his;
  static clock_t start_time, finish;   
  static double duration, lastClock, curTime;

  static int *headerPrinted;

  static N_Vector f1_Vec;
  static N_Vector f2_Vec;
};

#endif
