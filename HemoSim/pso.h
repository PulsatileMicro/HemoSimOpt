#ifndef PSO_H
#define PSO_H

#define PSO_Dim   8     // 每个粒子的维度 13 or 8
#define PSO_PNum  20    // 粒子数量
#define PSO_N     20    // 迭代次数

typedef struct PARTICLE{
  double X[PSO_Dim];
  double P[PSO_Dim];
  double V[PSO_Dim];
  double LastX[PSO_Dim];    // 在更新前，保存上次的X和V，万一更新的粒子不能算出目标函数值，则需要回到上一次的值，重新计算新的X和V
  double LastV[PSO_Dim];
  double Fitness;
  int newFlag;              // 初始为0，如果该粒子是淘汰后新引入的，则计数加1
}particle;

typedef struct SWARM{  
  particle Particle[PSO_PNum];
  int GBestIndex;  
  double GBestFitness;
  double GBest[PSO_Dim];
  double W;
  double C1;
  double C2;
  double Xup[PSO_Dim];
  double Xdown[PSO_Dim];
  double Vmax[PSO_Dim];
}swarm;  

void    PSO_RandInitofSwarm(void);
void    PSO_PriesInitofSwarm(void); 
void    PSO_ComputFitofSwarm(int, int);
void    PSO_FirstComputPandGbest(void);  
void    PSO_UpdateofVandX(void);  
void    PSO_ModifyVandX(int i);
void    PSO_UpdatePandGbest(void);  

double  PSO_InertWeight(void);  
int     PSO_CriteriaofStop(void);
double  PSO_ComputAFitness(double X[]);

extern swarm *s;
extern particle *p;

#endif