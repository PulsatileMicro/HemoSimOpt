#ifndef MODELPARAM_H
#define MODELPARAM_H

#include "Domain.h"

class ModelParam{
public:
  ModelParam();
  ~ModelParam();

  // ����
  static double dt;
  // ѪҺ�ܶ�
  static double rho;
  // 1Dģ���е�alpha
  static double oneD_alpha;
  // ���沽��
  static int nSteps;
  // 1Dģ���У�ÿ����hisSteps������¼һ��
  static int hisSteps;
  // ODE Solver ��Ծ���
  static double relTol;
  // ODE Solver ���Ծ���
  static double absTol;
  // Riemann Solver ����
  static double riemannTol;
  // Ѫ��׶��
  static double taperRate;
  // ȥ���ٻ� - lamda
  static double ModelParam::scale_lamda;
  // ȥ���ٻ� - u0
  static double ModelParam::scale_u0;
  // ȥ���ٻ� - r0
  static double ModelParam::scale_r0;
  // �Ƿ���ʾCFLֵ
  static int showCFL;
  // �Ƿ���ʾÿ����ʱ
  static int showStepLapse;
  // �Ƿ�����F-LЧӦ����Viscosity
  static int updVisc;
  // DG�����Ƿ�ʹ��Riemann Solver
  static int riemann;
  // �Ƿ�ȥ���ٻ�
  static int nonDim;
  // ����ļ���ʽ
  static int outMode;
  // ���������
  static int solverType;
  // �Ƿ��Ƿֲ洦������ʧ
  static int bifurloss;
  // Ѫ�ܶ���
  static int Ndoms;
  // �ڵ���
  static int Nnodes;
  // ճ����I,II
  static int gammaI;  // New Viscoelastic Model
  static int gammaII; // Old Viscoelastic Model

  // 1Dģ��L
  static int oneD_L;
  // 1Dģ��q
  static int oneD_q;
  // 1Dģ��Adams-Bashforth��������
  static int oneD_Je;
  // 1Dģ������ָ��
  static double *oneD_Uh, *oneD_Uhj, *oneD_Ah, *oneD_Ahj, *oneD_Hdh, *oneD_Hdhj, **oneD_Ufh, **oneD_Ufhj, **oneD_Afh, **oneD_Afhj;
  // 1Dģ�������q,L
  static int oneD_MaxQ, oneD_MaxL;

  // ��̬ģ��
  static int sparseSolver;

  // Domain����
  static Domain *omega;

  // �������������
  static int argc;
  static char **argv;

  enum SolverType{
    oneD_EXP=0,
    oneD_IMP=1,
    RLC_EXP=2,
    RLC_IMP=3,
    SS=4,
    Womer_Tree=5,  // Only support tree-like vasculature
    Womer_Net=6,
	  RC_EXP=7,
	  RC_IMP=8,
    SS_Sparse=9,
    Adap_SS_Wall=10,
    Adap_Pulse_Wall=11,
    Adap_SS_NoWall=12,
  };

  // Բ����
  static double PI;
};

#endif
