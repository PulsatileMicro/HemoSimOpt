#ifndef ADAP_SS_SOLVER_H
#define ADAP_SS_SOLVER_H

#include "Domain.h"
#include <stdio.h>
#include <vector>
#include <algorithm>
using namespace std;

class Adap_SS_Solver{
public:
	static void initSolver();
	static void solve();
	static void destroySolver();

	// For optimization loop
	static int AdapObjFunc();
	// For debug and output
	static vector<double> Debug_MeanP, Debug_MeanQ;

private:
	// General
	static void hemoRun();
	static void Write_history(Domain *omega, char *name);
	// For viscosity loop
	static void adjustFlowDir();
	static void adjustHdCalcPOrder();
	static void calcHd();
	static void UpdateVisc(Domain *omega);
	// For metabolic adaptation loop
	static void adjustScCalcOrder();
	static void calcPO2();
	static void calcSm();
	static void calcSc();
	static double Hill(double SO2);
	static double mean_abs(vector<double> vec);
	// For multiple call the AdapObjFunc()
	static void resetInitState();

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
	static double PS_Effect(double Df,double Da,double Db,double Hd,double FQb);
	static double logit(double x);

	static double **JMat;
	static double *NodeFlow;
	static double *SS_Press;
	static double *SS_DeltaP;
	static double *MeanP;
	static int *Pmarker, *Nmarker;
	static vector<double> MeanQ;
	static vector<double> lastVisc;
	static vector<double> Dm,dDm,nDm;
	static vector<double> Aw,dAw,nAw;
	static vector<double> MPO;
	static vector<double> tau;      // 剪切力
	static vector<double> o;        // 周向应力
	static vector<double> Stm,Som,Stot;
	static int sortCnt;
	static double dDmTOL,dAwTOL,StotTOL;
	static int logerr;

	// For debug
	static vector<double> Debug_Hd;
	static vector<double> Debug_Sm, Debug_Sc, Debug_Sp, Debug_Stau;
	static vector<double> Debug_SO2, Debug_PO2;
	static vector<double> Debug_Diam;

	// Store the original structural parameters
	static vector<double> Org_Diam, Org_WallTh, Org_Visc;
	static int ***org_bifur;
	static char **org_bctype;
	static int ***org_nodes;

	static double *q_in;
	static FILE **fp_his;
};

#endif
