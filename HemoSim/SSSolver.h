#ifndef SSSOLVER_H
#define SSSOLVER_H

#include "Domain.h"
#include <stdio.h>

class SSSolver{
public:
  static void initSolver();
  static void solve();
  static void destroySolver();

private:
  static void Write_history(Domain *omega, char *name);

  static double **JMat;
  static double *NodeFlow;
  static double *SS_Press;
  static double *SS_DeltaP;
  static double *MeanP;
  static double *MeanQ;

  static double *q_in;
  static FILE **fp_his;
};

#endif
