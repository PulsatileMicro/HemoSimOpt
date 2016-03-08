#ifndef PSO_H
#define PSO_H

#define PSO_Dim   8     // 每个粒子的维度 13 or 8
#define PSO_PNum  20    // 粒子数量
#define PSO_N     5000    // 迭代次数, 经验值1500 enough！！
//#define NOADAPERR TRUE	// 每个粒子必须计算出结果

typedef struct PARTICLE{
	double X[PSO_Dim];
	double PBest[PSO_Dim];
	double V[PSO_Dim];
	double LastX[PSO_Dim];    // 在更新前，保存上次的X和V，万一更新的粒子不能算出目标函数值，则需要回到上一次的值，重新计算新的X和V
	double LastV[PSO_Dim];
	double Fitness;
	double PBestFitness;
	int newFlag;              // 初始为0，如果该粒子是淘汰后新引入的，则计数加1
}particle;

typedef struct SWARM{  
	particle Particle[PSO_PNum];
	int	consist_convergence_time;
	int total_convergence_time;
	double GBestFitness;
	double GBest[PSO_Dim];
	double W;
	double C1;
	double C2;
	double Xup[PSO_Dim];
	double Xdown[PSO_Dim];
	double Vmax[PSO_Dim];
	double Alpha;
	double MeanX[PSO_Dim];
}swarm;  

void    PSO_RandInitofSwarm(void);
void    PSO_PriesInitofSwarm(void); 
void    PSO_ComputFitofSwarm(int);
void    PSO_FirstComputPandGbest(void);  
void    PSO_UpdateofVandX(void);
void	PSO_UpdateofVandX_CompressMutation(void);
void	PSO_UpdateofVandX_QuantumBehavior(int);
void	PSO_UpdateofVandX_SecondBehavior(void);
void    PSO_ModifyVandX(int i);
void    PSO_UpdatePandGbest(int);
void    PSO_UpdatePandGbest_SelectBehavior(void);  

double  PSO_InertWeight(void);  
int     PSO_CriteriaofStop(void);
double  PSO_ComputAFitness(double X[]);
void	PSO_QuickSort(particle[], int, int);
int		PSO_QuickSortPartion(particle[], int, int);
extern swarm *s;
extern particle *p;

#endif