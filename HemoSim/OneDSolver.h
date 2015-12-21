#ifndef ONEDSOLVER_H
#define ONEDSOLVER_H

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <cvode/cvode_lapack.h>      /* prototype for CVLapackDense */
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., and macros */
#include "Domain.h"
#include "ModelParam.h"
#include <time.h>

#define MAX_Q 128
#define MAX_ITER 500
#define BUFSIZE 1024

class OneDSolver{
public:
	static void initSolver(int solverType);
	static void solve();
	static void destroySolver();

private:
	static void BioAdv(Domain *omega);
	static void BioFlux(Domain *omega, char *name);
	static void Bifur(Domain *omega);
	static void BioSource(Domain *omega);
	static void UpdateVisc(Domain *omega);
	static int  Riemann(int i, int nel, Element *A, double *W, double *uu, double *Au, double rho);
	static int  BifurRiem2to1(int nel_d1, int nel_d2, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho);
	static int  BifurRiem1to2(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho);
	static int  BifurRiem1to2_losses_div_right(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K31, double K32);
	static int  BifurRiem1to2_losses_div_left(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K31, double K32);
	static int  BifurRiem1to2_losses_comb_right(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K13, double K23);
	static int  BifurRiem1to2_losses_comb_left(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K13, double K23);
	static int  BifurRiem1to2_losses_leav(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K13, double K23);
	static int  BifurRiem1to2_losses_ent(int nel_p, Element *A_p, Element *A_d1, Element *A_d2, double *W, double *uu, double *Au, double rho, double K31, double K32);
	static int  JuncRiemann(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho);
	static int  JuncRiemann_losses_exp_right(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho, double K);
	static int  JuncRiemann_losses_exp_left(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho, double K);
	static int  JuncRiemann_losses_cont_right(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho, double K);
	static int  JuncRiemann_losses_cont_left(int nel_l, Element *A_l, Element *A_r, double *W, double *uu, double *Au, double rho, double K);
	static double eval_A_from_P(double p, double po, double beta, double Ao);
	static double eval_BC_A(double A_star, double A);
	static int  WKinflow(Element *A, double u_right, double A_right, double Q, double *uu, double *Au);
	static int  WKoutflow(int nel, Element *A, double u_left, double A_left, double p_ext, double p_inf, double R, double rho, double *uu, double *Au);
	static int  CVODE_RHS(realtype t, N_Vector y, N_Vector ydot, void *user_data);
	static void Eval_RHS(int argc, char *argv[]);
	static void Write_history(Domain *omega, char *name);
	static double CFL(Domain *omega);
	static int Jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, 
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

	// Fahraeus-Lindqvist effects functions
	static double Eta_vitro(double Hd, double diam);
	static double Eta_vivo(double Hd, double diam, double W);
	static double Eta_vivo_ESL(double Hd, double diam);
	static double Eta_045(double diam);
	static double Eta_C(double diam);
	static double Eta_Was(double diam);
	static double Eta_Wpeak(double diam);
	static double Eta_Deff(double Hd, double diam);
	static double Eta_Dph(double diam);

	static char buf[BUFSIZE];

	static double *m_w, *m_w_mat;	/**< used for getzw function */
	static double **m_g1, *m_g1_mat;		/**< used for get_basis function */
	static double **m_d, **m_dt;		/**< used for getD function */
	static double *ip_tmp_Umat, *ip_tmp_Amat;
	static double *unitDiag;       // diagonal matrix
	static int *headerPrinted;

	// CVODE vars
	static N_Vector init_Y;
	static void *cvode_mem;
	static realtype tret;
	static int flag, reIntFlag;
	static double *cvode_u, *cvode_a, *cvode_yu, *cvode_ya, *cvode_ydotu, *cvode_ydota;
	static N_Vector f1_Vec;
	static N_Vector f2_Vec;

	// Used for reading from input file
	static double t_in, *pff, *pbb, *uff, *ubb, Rt, *q_in, *A_star, *u_star;

	static FILE **fp_his;
	static clock_t start_time, finish_time;   
	static double timelapse, lastClock, curTime;
};

#endif
