#ifndef SPARSSSOLVER_H
#define SPARSSSOLVER_H

#include "Domain.h"
#include <stdio.h>
#include <mkl.h>

class SparSSSolver{
public:
	static void initSolver();
	static void solve();
	static void destroySolver();

private:
	static void Write_history(Domain *omega, char *name);
	static void initPardiso();

	static MKL_INT *ia;
	static MKL_INT *ja;
	static double *sparse_a;
	static double *results;
	static double **JMat;
	static double *NodeFlow;
	static double *SS_Press;
	static double *SS_DeltaP;
	static double *MeanP;
	static double *MeanQ;

	static double *q_in;
	static FILE **fp_his;

	// For PARDISO
	static void *pt[64];
	static MKL_INT iparm[64];
	static MKL_INT maxfct, mnum, error, msglvl;
	static MKL_INT intel_mtype;		/* Real asymmetric matrix */
	static MKL_INT phase;
	static double ddum;			/* Double dummy */
	static MKL_INT idum;			/* Integer dummy. */
	static MKL_INT nrhs;
	static int nonzero_cnt;
};

#endif
