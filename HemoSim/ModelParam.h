#ifndef MODELPARAM_H
#define MODELPARAM_H

#include "Domain.h"

class ModelParam{
public:
  ModelParam();
  ~ModelParam();

  // 步长
  static double dt;
  // 血液密度
  static double rho;
  // 1D模型中的alpha
  static double oneD_alpha;
  // 仿真步数
  static int nSteps;
  // 1D模型中，每仿真hisSteps步，记录一次
  static int hisSteps;
  // ODE Solver 相对精度
  static double relTol;
  // ODE Solver 绝对精度
  static double absTol;
  // Riemann Solver 精度
  static double riemannTol;
  // 血管锥度
  static double taperRate;
  // 去量纲化 - lamda
  static double ModelParam::scale_lamda;
  // 去量纲化 - u0
  static double ModelParam::scale_u0;
  // 去量纲化 - r0
  static double ModelParam::scale_r0;
  // 是否显示CFL值
  static int showCFL;
  // 是否显示每步耗时
  static int showStepLapse;
  // 是否依据F-L效应更新Viscosity
  static int updVisc;
  // DG法中是否使用Riemann Solver
  static int riemann;
  // 是否去量纲化
  static int nonDim;
  // 输出文件格式
  static int outMode;
  // 求解器类型
  static int solverType;
  // 是否考虑分叉处能量损失
  static int bifurloss;
  // 血管段数
  static int Ndoms;
  // 节点数
  static int Nnodes;
  // 粘弹性I,II
  static int gammaI;  // New Viscoelastic Model
  static int gammaII; // Old Viscoelastic Model

  // 1D模型L
  static int oneD_L;
  // 1D模型q
  static int oneD_q;
  // 1D模型Adams-Bashforth方法阶数
  static int oneD_Je;
  // 1D模型数据指针
  static double *oneD_Uh, *oneD_Uhj, *oneD_Ah, *oneD_Ahj, *oneD_Hdh, *oneD_Hdhj, **oneD_Ufh, **oneD_Ufhj, **oneD_Afh, **oneD_Afhj;
  // 1D模型中最大q,L
  static int oneD_MaxQ, oneD_MaxL;

  // 稳态模型
  static int sparseSolver;

  // Domain数组
  static Domain *omega;

  // 主程序输入参数
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

  // 圆周率
  static double PI;
};

#endif
