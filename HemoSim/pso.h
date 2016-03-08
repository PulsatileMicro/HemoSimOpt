#ifndef PSO_H
#define PSO_H

#define PSO_Dim   8     // ÿ�����ӵ�ά�� 13 or 8
#define PSO_PNum  20    // ��������
#define PSO_N     5000    // ��������, ����ֵ1500 enough����
//#define NOADAPERR TRUE	// ÿ�����ӱ����������

typedef struct PARTICLE{
	double X[PSO_Dim];
	double PBest[PSO_Dim];
	double V[PSO_Dim];
	double LastX[PSO_Dim];    // �ڸ���ǰ�������ϴε�X��V����һ���µ����Ӳ������Ŀ�꺯��ֵ������Ҫ�ص���һ�ε�ֵ�����¼����µ�X��V
	double LastV[PSO_Dim];
	double Fitness;
	double PBestFitness;
	int newFlag;              // ��ʼΪ0���������������̭��������ģ��������1
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