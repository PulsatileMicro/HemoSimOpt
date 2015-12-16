#include "PreProcessor.h"
#include "ModelParam.h"
#include "OneDSolver.h"
#include "ZeroDSolver.h"
#include "SSSolver.h"
#include "RCSolver.h"
#include "SparSSSolver.h"
#include "AdapSSWall.h"
#include <iostream>
using namespace std;

// #include <vld.h>

/*!
  Print a message on the screen. 
*/
void Display_screen(Domain *omega){
  fprintf(stdout,"Hemodynamics Simulation Program. \n");
  switch(ModelParam::solverType){
    case ModelParam::oneD_EXP:
      fprintf(stdout,"Solver Type: 1D Explicit Solver\n");
      break;
    case ModelParam::oneD_IMP:
      fprintf(stdout,"Solver Type: 1D Implicit Solver\n");
      break;
    case ModelParam::RLC_EXP:
      fprintf(stdout,"Solver Type: RLC Circuit Explicit Solver\n");
      break;
    case ModelParam::RLC_IMP:
      fprintf(stdout,"Solver Type: RLC Circuit Implicit Solver\n");
      break;
    case ModelParam::SS:
      fprintf(stdout,"Steady State Solver\n");
      break;
    case ModelParam::SS_Sparse:
      fprintf(stdout,"Sparse Steady State Solver\n");
      break;
    case ModelParam::Womer_Tree:
      fprintf(stdout,"Womersley Solver for arterial tree\n");
      break;
    case ModelParam::Womer_Net:
      fprintf(stdout,"Womersley Solver for network\n");
      break;
    case ModelParam::RC_EXP:
      fprintf(stdout,"RC Circuit Explicit Solver\n");
      break;
    case ModelParam::RC_IMP:
      fprintf(stdout,"RC Circuit Implicit Solver\n");
      break;
    case ModelParam::Adap_SS_Wall:
      cout << "Steady state simulation for structural adaptation" << endl;
      break;
    default:
      break;
  }
  fprintf(stdout,"dt:                     %lg\n",ModelParam::dt);
  fprintf(stdout,"Number of time steps:   %d\n", ModelParam::nSteps);
  fprintf(stdout,"Number of domains:      %d\n", ModelParam::Ndoms);
  fprintf(stdout,"Quadrature order:       %d\n", ModelParam::oneD_q);
  fprintf(stdout,"Polynomial order:       %d\n", ModelParam::oneD_L);
  fprintf(stdout,"Nondimensionalized:     %d\n", ModelParam::nonDim);
  fprintf(stdout,"UpdVisc strategy:       %d\n", ModelParam::updVisc);
  fprintf(stdout,"Riemann solver:         %d\n", ModelParam::riemann);
  fprintf(stdout,"Binary output:          %d\n", ModelParam::outMode);
}

int main(int argc, char *argv[]){
  ModelParam::argc=argc;
  ModelParam::argv=argv;
  /* 预处理，从.in文件读取数据 */
  PreProcessor preProcessor;
  preProcessor.run(argc, argv);

  Display_screen(ModelParam::omega);

  switch (ModelParam::solverType){
    /*
      1D Model
      Explicit & Implicit Solvers
    */
    case ModelParam::oneD_EXP:
    case ModelParam::oneD_IMP:
      OneDSolver::initSolver(ModelParam::solverType);
      OneDSolver::solve();
      OneDSolver::destroySolver();
      break;
    /*
      0D Model (Bond Graph/Circuit method)
      Explicit & Implicit Solvers
    */
    case ModelParam::RLC_EXP:
    case ModelParam::RLC_IMP:
      ZeroDSolver::initSolver(ModelParam::solverType);
      ZeroDSolver::solve();
      ZeroDSolver::destroySolver();
      break;
    /*
      Steady state Model
    */
    case ModelParam::SS:
      SSSolver::initSolver();
      SSSolver::solve();
      SSSolver::destroySolver();
      break;
    case ModelParam::SS_Sparse:
      SparSSSolver::initSolver();
      SparSSSolver::solve();
      SparSSSolver::destroySolver();
      break;
      /*
      Womersley Solution for Tree-like vasculature
      TODO
    */
    case ModelParam::Womer_Tree:
      break;
    /*
      Womersley Solution for general network
      TODO
    */
    case ModelParam::Womer_Net:
      break;
    /*
      0D RC circuit Model
      Explicit & Implicit Solvers
    */
    case ModelParam::RC_EXP:
    case ModelParam::RC_IMP:
      RCSolver::initSolver(ModelParam::solverType);
      RCSolver::solve();
      RCSolver::destroySolver();
      break;
    /*
      Steady state model for simulating structural adaptation
    */
    case ModelParam::Adap_SS_NoWall:
    case ModelParam::Adap_SS_Wall:
      Adap_SS_Solver::initSolver();
      Adap_SS_Solver::solve();
      Adap_SS_Solver::destroySolver();
      break;
    default:
      break;
  }
}
