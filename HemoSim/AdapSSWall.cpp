#include "AdapSSWall.h"
#include "ModelParam.h"
#include "AdapParam.h"
#include <stdlib.h>
#include <mkl.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <string>
#include "pso.h"
#include "PreProcessor.h"
using namespace std;

//#ifndef _DEBUG
//#	define _DEBUG
//#endif

double  **Adap_SS_Solver::JMat=NULL;
double  *Adap_SS_Solver::NodeFlow=NULL;
double  *Adap_SS_Solver::SS_Press=NULL;
double  *Adap_SS_Solver::SS_DeltaP=NULL;
double  *Adap_SS_Solver::MeanP=NULL;
double  *Adap_SS_Solver::q_in=NULL;
int     *Adap_SS_Solver::Pmarker=NULL;
int     *Adap_SS_Solver::Nmarker=NULL;
int   ***Adap_SS_Solver::org_nodes=NULL;
char   **Adap_SS_Solver::org_bctype=NULL;
int   ***Adap_SS_Solver::org_bifur=NULL;
FILE    **Adap_SS_Solver::fp_his;
int     Adap_SS_Solver::sortCnt=1;
double  Adap_SS_Solver::dDmTOL=1e-3; /*��������Ҫ��*/
double  Adap_SS_Solver::dAwTOL=1e-3;
double  Adap_SS_Solver::StotTOL=1e-5;
int     Adap_SS_Solver::logerr=0;

vector<double> Adap_SS_Solver::Dm;
vector<double> Adap_SS_Solver::dDm;
vector<double> Adap_SS_Solver::nDm;
vector<double> Adap_SS_Solver::Aw;
vector<double> Adap_SS_Solver::dAw;
vector<double> Adap_SS_Solver::nAw;
vector<double> Adap_SS_Solver::tau;      // ������
vector<double> Adap_SS_Solver::o;        // ����Ӧ��
vector<double> Adap_SS_Solver::Stm;
vector<double> Adap_SS_Solver::Som;
vector<double> Adap_SS_Solver::Stot;
vector<double> Adap_SS_Solver::Debug_Hd;
vector<double> Adap_SS_Solver::Debug_Sm;
vector<double> Adap_SS_Solver::Debug_Sc;
vector<double> Adap_SS_Solver::Debug_Sp;
vector<double> Adap_SS_Solver::Debug_Stau;
vector<double> Adap_SS_Solver::Debug_SO2;
vector<double> Adap_SS_Solver::Debug_PO2;
vector<double> Adap_SS_Solver::Debug_MeanP;
vector<double> Adap_SS_Solver::Debug_Diam;
vector<double> Adap_SS_Solver::Debug_MeanQ;
vector<double> Adap_SS_Solver::MPO;
vector<double> Adap_SS_Solver::lastVisc;
vector<double> Adap_SS_Solver::MeanQ;
vector<double> Adap_SS_Solver::Org_Diam;
vector<double> Adap_SS_Solver::Org_Visc;
vector<double> Adap_SS_Solver::Org_WallTh;

void Adap_SS_Solver::initSolver(){
	int i,j,k;
	/*
	Create the JMat matrix and initialize it with zeros
	*/
	JMat=(double **)malloc(ModelParam::Nnodes*sizeof(double*));
	for(i=0; i<ModelParam::Nnodes; ++i){
		JMat[i]=(double *)malloc(ModelParam::Nnodes*sizeof(double));
	}
	for(i=0; i<ModelParam::Nnodes; ++i)
		for (j=0; j<ModelParam::Nnodes; ++j)
			JMat[i][j]=0;

	/*!
	Create the original bctype, bifur, and nodes structures
	*/
	org_bctype=(char **)malloc(ModelParam::Ndoms*sizeof(char *));
	for(i=0;i<ModelParam::Ndoms;i++){
		org_bctype[i]=(char *)malloc(6*sizeof(char)); // �����������ڵ�bctype
	}
	org_bifur=(int ***)malloc(ModelParam::Ndoms*sizeof(int **));
	org_nodes=(int ***)malloc(ModelParam::Ndoms*sizeof(int **));
	for(i=0;i<ModelParam::Ndoms;i++){
		org_bifur[i]=imatrix(0,1,0,1);
		org_nodes[i]=imatrix(0,1,0,3);
	}

	/*
	Create the boundary vector and initialize it with zeros
	*/
	NodeFlow=(double *)malloc(ModelParam::Nnodes*sizeof(double));
	for(i=0; i<ModelParam::Nnodes; ++i)
		NodeFlow[i]=0;

	SS_Press  =(double *)malloc(ModelParam::Nnodes*sizeof(double));
	SS_DeltaP =(double *)malloc(ModelParam::Ndoms*sizeof(double));
	MeanP     =(double *)malloc(ModelParam::Ndoms*sizeof(double));
	Pmarker   =(int *)malloc(ModelParam::Ndoms*sizeof(int));
	Nmarker   =(int *)malloc(ModelParam::Ndoms*sizeof(int));
	for(i=0;i<ModelParam::Ndoms;i++){
		Pmarker[i]=0;
		Nmarker[i]=0;
	}

	for(i=0;i<ModelParam::Ndoms;i++){
		lastVisc.push_back(0.0);
		MeanQ.push_back(0.0);
		Dm.push_back(0.0);
		Aw.push_back(0.0);
		tau.push_back(0.0);
		o.push_back(0.0);
		MPO.push_back(0.0);
		Debug_Hd.push_back(0.0);
		Debug_Sm.push_back(0.0);
		Debug_Sc.push_back(0.0);
		Debug_Sp.push_back(0.0);
		Debug_Stau.push_back(0.0);
		Debug_SO2.push_back(0.0);
		Debug_PO2.push_back(0.0);
		Debug_Diam.push_back(0.0);
		Debug_MeanP.push_back(0.0);
		Debug_MeanQ.push_back(0.0);
		Stm.push_back(0.0);
		Som.push_back(0.0);
		Stot.push_back(0.0);
		dDm.push_back(0.0);
		dAw.push_back(0.0);
		nDm.push_back(0.0);
		nAw.push_back(0.0);
		Pmarker[i]=0;
		Nmarker[i]=0;

		// �����ʼ�Ľṹ�����˲�������resetInitState()����ʹ��
		double tempDiam=sqrt(4*ModelParam::omega[i].A[0].h[0]/ModelParam::PI);
		Org_Diam.push_back(tempDiam);
		double tempVisc=ModelParam::omega[i].visc[0];
		Org_Visc.push_back(tempVisc);
		double tempWallTh=ModelParam::omega[i].WallTh;
		Org_WallTh.push_back(tempWallTh);

		// ����bifur�ĳ�ʼ״̬
		for(j=0;j<2;j++){
			for(k=0;k<2;k++){
				org_bifur[i][j][k]=ModelParam::omega[i].bifur[j][k];
			}
		}
		// ����bctype�ĳ�ʼ״̬
		for(j=0;j<6;j++){
			org_bctype[i][j]=ModelParam::omega[i].bctype[j];    // 0,1,2 ���, 3,4,5 ����
		}
		// ����nodes�ĳ�ʼ״̬
		for(j=0;j<2;j++){
			for(k=0;k<4;k++){
				org_nodes[i][j][k]=ModelParam::omega[i].nodes[j][k];
			}
		}  
	}

	AdapParam::optCate=AdapParam::PSO;
	if (!strncmp(ModelParam::argv[ModelParam::argc-2], "STDPSO",6))
	{
		AdapParam::optMethod = AdapParam::STDPSO;
	} 
	else if (!strncmp(ModelParam::argv[ModelParam::argc-2], "CMPPSO",6))
	{
		AdapParam::optMethod = AdapParam::CMPPSO;
	}
	else if (!strncmp(ModelParam::argv[ModelParam::argc-2], "QUAPSO",6))
	{
		AdapParam::optMethod = AdapParam::QUAPSO;
	}
	else{
		cout << "connot figure out adapParam::optType!" << endl;
		exit(EXIT_FAILURE);
	}
	AdapParam::adapLogFile.open("adap_results.txt");
	if(ModelParam::solverType==ModelParam::Adap_SS_Wall){
		AdapParam::adapLogFile << setw(15) << "Particle_ID" << setw(15) << "ErrorV" << setw(15) << "ErrorD" 
			<< setw(15) << "ErrorQ" << setw(15) << "PO2Ref" << setw(15) << "kc" << setw(15) << "kmd" 
			<< setw(15) << "kmg" << setw(15) << "ksd" << setw(15) << "ksg" << setw(15) << "kwtau" 
			<< setw(15) << "kwsigma" << setw(15) << "tauRef" << setw(15) << "sigmaRef" << setw(15) << "wRef" 
			<< setw(15) << "J0" << setw(15) << "LRef" << endl;
	}
	else if(ModelParam::solverType==ModelParam::Adap_SS_NoWall){
		AdapParam::adapLogFile << setw(15) << "Particle_ID" << setw(15) << "ErrorV" << setw(15) << "ErrorD" 
			<< setw(15) << "ErrorQ" << setw(15) << "kc" << setw(15) << "kp" << setw(15) << "km" 
			<< setw(15) << "ks" << setw(15) << "J0" << setw(15) << "LRef" << setw(15) << "tauRef" 
			<< setw(15) << "QRef" << endl;
	}
}

/*!
\brief Optimize the parameters of the structural adaptation model
*/
void Adap_SS_Solver::solve(){
	int n=0;
	AdapParam::AdapErrCnt=0;
	switch(AdapParam::optCate){
	case AdapParam::PSO:
		AdapParam::set_model_coeff_bounds(5, 0.2);   // ��ʼ������Ӧ����
		PSO_RandInitofSwarm();              // ��ʼ��PSO����
		// PSO_PriesInitofSwarm();

		// Start Iteration
		while(n++!=PSO_N)
		{  
			AdapParam::adapLogFile << "Iteration " << n << endl;
			printf("The %dth time to calculate .\n", n);
			//printf("Updated of the swarm's Fitness:\n");  
			PSO_ComputFitofSwarm(n);
			//printf("Replaced of P and Gbest:\n\n");  
			PSO_UpdatePandGbest(n);
			//printf("Updated of the swarm:\n\n"); 
			switch(AdapParam::optMethod){
				case AdapParam::STDPSO:
					PSO_UpdateofVandX();
					break;
				case AdapParam::CMPPSO:
					PSO_UpdateofVandX_CompressMutation();
					break;
				case AdapParam::QUAPSO:
					PSO_UpdateofVandX_QuantumBehavior(n);
					break;
				case AdapParam::SECPSO:
					PSO_UpdateofVandX_SecondBehavior();
					break;
			}
			printf("AdapErrCnt=%d\n", AdapParam::AdapErrCnt);
		}  
		break;
	case AdapParam::DOWNHILL:
		break;
	/*case AdapParam::NO_OPT:
		AdapParam::initPriesAdapParam();
		AdapObjFunc();
		break;*/
	default:
		break;
	}
}

/*!
\brief Solve the structural adaptation process
*/
int Adap_SS_Solver::AdapObjFunc()
{
	int n,i,j,k;
	int loop1_cnt=0,loop1_visc_num=30;      // visc loop������
	int loop2_cnt=0,loop2_adap_num=1000;    // adap loop������
	double sumViscErr=0.0;
	double maxViscErr=0.0;
	double viscErr;
	double diam, tempP, tauP;
	double mean_Stot, last_mean_Stot=0.0;

	// Parameters of the structural adaptation model
	// Pries et al. (2005) Hypertension
	double ktd=1,kog=1,ee=1e-9;
	double Rtd,Rod,Rtg,Rog,Rw;

	// Reset the structural parameters and connection matrix to the original state
	resetInitState();

	AdapParam::AdapErrCnt=0;
	AdapParam::errFlag=AdapParam::NO_ADAP_ERR;
	while(++loop2_cnt<loop2_adap_num){
		loop1_cnt=0;
		while(++loop1_cnt<loop1_visc_num){
			// Loop 1: compute the hemodynamic parameters and update the blood apparent viscosity
			maxViscErr=0.0;
			hemoRun();
			if(AdapParam::errFlag){
				break;
			}
			for(n=0;n<ModelParam::Ndoms;n++){
				Debug_MeanQ[n]=MeanQ[n];
				Debug_MeanP[n]=MeanP[n];
			}
			// �������¼������Ѫ������ѧ����������Ѫ������
			adjustFlowDir();
			if(AdapParam::errFlag){
				break;
			}
			adjustHdCalcPOrder();
			if(AdapParam::errFlag){
				// AdapParam::AdapErrCnt++;
				break;
			}

			// ����Hd
			calcHd();
			// ���ݸ��µ�Hd����viscosity
			UpdateVisc(ModelParam::omega);
			// ����SOR��������viscosity, ���ж��Ƿ�����
			if(lastVisc[0]!=0){ // �ж��Ƿ��ǵ�һ�μ���viscosity
				// ����viscosity
				for(n=0;n<ModelParam::Ndoms;++n){
					for(i=0;i<ModelParam::oneD_q-1;++i)
						ModelParam::omega[n].visc[i]=(ModelParam::omega[n].visc[i]+lastVisc[n])*0.5;
				}
				// �ж��Ƿ�����
				for(n=0;n<ModelParam::Ndoms;++n){
					viscErr=(ModelParam::omega[n].visc[0]-lastVisc[n])/lastVisc[n];
					if(fabs(viscErr)>maxViscErr)
						maxViscErr=fabs(viscErr);
				}
				// cout << "ViscLoop=" << loop1_cnt << " MaxViscErr=" << maxViscErr << endl;
				// ��������
				if(maxViscErr<1e-3)
					break;  // ����loop1
			}

			// �����µ�visc������lastVisc����
			// TODO: ʹ��MATLAB�����е�modifyViscosity����
			for(n=0;n<ModelParam::Ndoms;++n)
				lastVisc[n]=ModelParam::omega[n].visc[0];
#ifdef _DEBUG
			cout << "ViscLoop=" << loop1_cnt << " MaxViscErr=" << maxViscErr;
#endif // _DEBUG
		}

		if(AdapParam::errFlag)
			break;

		// ��ʼ��������Ӧ����
		for(n=0;n<ModelParam::Ndoms;++n){
			// For wall adaptation
			diam=sqrt(4*ModelParam::omega[n].A[0].h[0]/ModelParam::PI)*1e3;   // unit:mm
			Dm[n]=diam+ModelParam::omega[n].WallTh*1e3/2;                     // unit:mm
			Aw[n]=ModelParam::omega[n].WallTh*1e3*2*ModelParam::PI*Dm[n];     // unit:mm2
			tau[n]=32*ModelParam::omega[n].visc[0]*MeanQ[n]/60*1e-12/(ModelParam::PI*diam*diam*diam)*1e9*10;
			o[n]=133*MeanP[n]*diam/(2*ModelParam::omega[n].WallTh*1e3)*10;
			// For no wall adaptation
			ModelParam::omega[n].Stau=log10(fabs(tau[n])+AdapParam::tauRef2[0]);
			Debug_Stau[n]=ModelParam::omega[n].Stau;
			tempP=MeanP[n];
			if(tempP<10)
				tempP=10;
			tauP=100-86*exp(-5000*pow(log10(log10(tempP)),5.4));
			ModelParam::omega[n].Sp=-log10(tauP);
			Debug_Sp[n]=ModelParam::omega[n].Sp;
		}

		// LrefVec
		// TODO

		// ����������Ӧ�źż���
		calcPO2();
		calcSm();
		// ����Sc����˳��
		adjustScCalcOrder();
		if(AdapParam::errFlag)
			break;

		calcSc();

		if(ModelParam::solverType==ModelParam::Adap_SS_Wall){
			// �����ܴ̼��ź�
			// ����Pries et al.(2005), R��ֵ�����������
			// ���￼��: �ܾ��仯��Ҫ��tau���𣬱ں�仯��Ҫ��sigma����
			Rtd=0.8;
			Rod=0.05;
			Rtg=0.05;
			Rog=0.8;
			Rw=0.5;

			for(n=0;n<ModelParam::Ndoms;++n){
				double WallTh=ModelParam::omega[n].WallTh*1e6;    // um
				// ���tau, sigma�Ȳ���С��0�������Ѫ������ѧ�����������
				if(tau[n]<=0 || WallTh<=0 || o[n]<=0){
					AdapParam::errFlag=AdapParam::ERR_HEMO;
					break;
				}
				Stm[n]=ktd*log10(tau[n]/AdapParam::tauRef[0]+ee)/(1+AdapParam::kwtau[0]*log10(WallTh/AdapParam::wRef[0]+ee))
					+AdapParam::kmd[0]*(ModelParam::omega[n].Sm+AdapParam::kc[0]*ModelParam::omega[n].Sc)-AdapParam::ksd[0];
				Som[n]=kog*log10(o[n]/AdapParam::sigmaRef[0]+ee)/(1+AdapParam::kwsigma[0]*log10(WallTh/AdapParam::wRef[0]+ee))
					+AdapParam::kmg[0]*(ModelParam::omega[n].Sm+AdapParam::kc[0]*ModelParam::omega[n].Sc)-AdapParam::ksg[0];
				dDm[n]=(Rtd*Stm[n]+Rod*Som[n])*Dm[n]*AdapParam::dt;
				dAw[n]=Rw*(Rtg*Stm[n]+Rog*Som[n])*Aw[n]*AdapParam::dt;
				// ���deltaD��deltaA����1������ֲ���������
				if(fabs(dDm[n])>1 || fabs(dAw[n])>1){
					AdapParam::errFlag=AdapParam::NO_CONV;
					break;
				}
				nDm[n]=Dm[n]+dDm[n];
				nAw[n]=Aw[n]+dAw[n];
				WallTh=nAw[n]/(2*ModelParam::PI*nDm[n]);
				ModelParam::omega[n].WallTh=WallTh/1e3;

				for(i=0;i<ModelParam::oneD_q-1;++i){
					double tempDiam=nDm[n]-WallTh/2;
					// �ܾ����Ϊ������˴�����Ӧʧ��
					if(tempDiam<0){
						AdapParam::errFlag=AdapParam::NEGDIAM;
						break;
					}
					else{
						ModelParam::omega[n].A[0].h[i]=0.25*ModelParam::PI*tempDiam*tempDiam/1e6;
						Debug_Diam[n]=tempDiam*1e3;
					}
				}
				if(AdapParam::errFlag)
					break;
			}
			/* 
			����Ӧѭ������ģ��
			�ж��Ƿ����
			*/
			if(!AdapParam::errFlag){
				// û�г���������������ж��Ƿ����
				double mean_dDm=mean_abs(dDm);
				double mean_dAw=mean_abs(dAw);
				// cout << "dDm=" << setprecision(6) << setw(12) << mean_dDm  << " dAw=" << setprecision(6) << setw(12) << mean_dAw << " visc iter=" << loop1_cnt << endl;
				if(mean_dDm<dDmTOL && mean_dAw<dAwTOL){
					// �ﵽ���ȣ���������Ӧѭ��
					break;
				}
			}
			else{
				AdapParam::AdapErrCnt++;
				// ��������Ӧ��������
				break;
			}
		}
		else if(ModelParam::solverType==ModelParam::Adap_SS_NoWall){
			// ��Ѫ�ܱ�����Ӧ
			for(n=0;n<ModelParam::Ndoms;++n){
				diam=sqrt(4*ModelParam::omega[n].A[0].h[0]/ModelParam::PI)*1e3;   // unit:mm
				// �����ܴ̼�
				Stot[n]=AdapParam::kp2[0]*ModelParam::omega[n].Sp+ModelParam::omega[n].Stau
					+AdapParam::km2[0]*(ModelParam::omega[n].Sm+AdapParam::kc2[0]*ModelParam::omega[n].Sc)-AdapParam::ks2[0];
				diam=diam+Stot[n]*diam*AdapParam::dt; // Update the diameter

				if(diam<=0){
					AdapParam::errFlag=AdapParam::NEGDIAM;
					break;
				}
				else if(diam>200){
					AdapParam::errFlag=AdapParam::NO_CONV;
					break;
				}
				else{
					ModelParam::omega[n].A[0].h[0]=0.25*ModelParam::PI*diam*diam/1e6;
					Debug_Diam[n]=diam*1e3;
				}
			}

			/* 
			����Ӧѭ������ģ��
			�ж��Ƿ����
			*/
			if(!AdapParam::errFlag){
				// û�г���������������ж��Ƿ����
				double sum_Stot=0.0;
				for(i=0;i<ModelParam::Ndoms;i++){
					sum_Stot+=Stot[i];
				}
				mean_Stot=sum_Stot/ModelParam::Ndoms;
				// cout << "Stot=" << setprecision(6) << setw(12) << mean_Stot << " visc iter=" << loop1_cnt << " adap iter=" << loop2_cnt << endl;
				if(fabs(mean_Stot)<StotTOL){
					// if(last_mean_Stot*mean_Stot<0 && fabs(last_mean_Stot+mean_Stot)<StotTOL){   // �����෴�Һ�С����ֵ
					// �ﵽ���ȣ���������Ӧѭ��
					break;
				}
				else
					last_mean_Stot=mean_Stot;
			}
			else{
				AdapParam::AdapErrCnt++;
				// ��������Ӧ��������
				break;
			}
		}
	}
	// ����Ӧѭ������
#ifdef _DEBUG
	cout << " AdapLoop=" << loop2_cnt << " last_mean_Stot=" << last_mean_Stot << endl;
#endif // _DEBUG

	if(loop2_cnt==loop2_adap_num){
		AdapParam::errFlag=AdapParam::NO_CONV;
	}

	/* 
	����Ŀ�꺯��
	����������� ErrorV
	����ܾ���� ErrorD
	*/
	double adapDiam, diam_err_sum=0.0;
	double adapVel,  vel_err_sum=0.0;
	double adapFlow, orgFlow, flow_err_sum=0.0;

	for(n=0;n<ModelParam::Ndoms;n++){
		orgFlow=fabs(0.25*Org_Diam[n]*Org_Diam[n]*ModelParam::PI*ModelParam::omega[n].MesVel);
		adapDiam=sqrt(4*ModelParam::omega[n].A[0].h[0]/ModelParam::PI);
		adapVel=(MeanQ[n]/60/1e12/ModelParam::omega[n].A[0].h[0])*1e3;    // mm/s
		adapFlow=0.25*adapDiam*adapDiam*ModelParam::PI*adapVel;
		diam_err_sum+=(Org_Diam[n]-adapDiam)*(Org_Diam[n]-adapDiam)/(Org_Diam[n]*Org_Diam[n]);
		vel_err_sum+=4*(fabs(ModelParam::omega[n].MesVel)-adapVel)*(fabs(ModelParam::omega[n].MesVel)-adapVel)/((fabs(ModelParam::omega[n].MesVel)+adapVel)*(fabs(ModelParam::omega[n].MesVel)+adapVel));
		flow_err_sum+=(adapFlow-orgFlow)*(adapFlow-orgFlow)/(orgFlow*orgFlow);
	}
	AdapParam::ErrorD=sqrt(diam_err_sum/ModelParam::Ndoms);
	AdapParam::ErrorV=sqrt(vel_err_sum/ModelParam::Ndoms);
	// AdapParam::ErrorQ=AdapParam::ErrorD+AdapParam::ErrorV;
	AdapParam::ErrorQ=AdapParam::ErrorV;

	if(AdapParam::errFlag){
		string errStr;
		switch(AdapParam::errFlag){
		  case AdapParam::FLOWDIR:
			  errStr="Flow Direction Error";
			  break;
		  case AdapParam::PORDER:
			  errStr="Positive Calculation Order Error";
			  break;
		  case AdapParam::NORDER:
			  errStr="Negative Calculation Order Error";
			  break;
		  case AdapParam::NEGDIAM:
			  errStr="Negative Diameter Error";
			  break;
		  case AdapParam::LINEQU:
			  errStr="Linear Equation Solver Error";
			  break;
		  case AdapParam::NO_CONV:
			  errStr="Not Converged";
			  break;
		  case AdapParam::ERR_HEMO:
			  errStr="Hemodynamic Parameter Error";
			  break;
		  default:
			  break;
		}
		cout << "error code=" << errStr << endl;
		/*AdapParam::ErrorQ=AdapParam::ErrorQ*10;
		AdapParam::ErrorV=AdapParam::ErrorV*10;*/
		AdapParam::ErrorV = AdapParam::ErrorQ = 100.00;
	}

	// Write_history(ModelParam::omega,ModelParam::argv[ModelParam::argc-1]);

	return AdapParam::errFlag;
}

/*!
\brief Evaluate the X in the equ. AX=B in the steady state model
*/
void Adap_SS_Solver::hemoRun(){
	int i, j, n;
	int start, end, from, to;
	int d1, d2, d3, d4;
	double radius, len;
	int From[2], To[2];

	if(!q_in){
		q_in = (double *)malloc(ModelParam::Ndoms*sizeof(double));
	}

	for(n = 0; n < ModelParam::Ndoms; ++n){         // For each domain "n" (n=0,...,Ndoms),
		start = ModelParam::omega[n].nodes[0][0];     // start node number
		end   = ModelParam::omega[n].nodes[1][0];     // end node number
		d1    = ModelParam::omega[n].bifur[0][0];     // start node�����ӵ�����Ѫ�ܵ�id
		d2    = ModelParam::omega[n].bifur[0][1];
		d3    = ModelParam::omega[n].bifur[1][0];     // end node�����ӵ�����Ѫ�ܵ�id
		d4    = ModelParam::omega[n].bifur[1][1];

		radius = sqrt(ModelParam::omega[n].A[0].h[0]/ModelParam::PI);
		len = ModelParam::omega[n].A[0].x[1]-ModelParam::omega[n].A[0].x[0];
		ModelParam::omega[n].J=ModelParam::PI*radius*radius*radius*radius/(8*ModelParam::omega[n].visc[0]*len);

		// start node
		switch(ModelParam::omega[n].bctype[0]){
			/*
			���start node������ΪB����ö�Ѫ��ֻ������Ѫ��
			d1,d2�����ʹ������������d1ΪĸѪ�ܣ�d2Ϊ��Ѫ�ܣ������෴��
			*/
		  case 'B':
			  if(d1<0 || d2<0 || start<0){
				  AdapParam::errFlag=AdapParam::FLOWDIR;
				  return;
			  }
			  else if(ModelParam::omega[d1].nodes[0][0] == start){
				  /*
				  ���d1������뵱ǰ����Ѫ�ܶε������ͬ��˵��d2��ĸѪ�ܣ�d1��n����Ѫ��
				  From[0]ΪĸѪ��d2�����
				  To[0]Ϊ��ǰѪ�ܶε��յ�
				  To[1]Ϊ��һ��Ѫ��d1���յ�
				  */
				  From[0]=ModelParam::omega[d2].nodes[0][0];
				  To[0]=end;
				  To[1]=ModelParam::omega[d1].nodes[1][0];
				  JMat[start][start]=-(ModelParam::omega[n].J+ModelParam::omega[d1].J+ModelParam::omega[d2].J);
				  JMat[start][From[0]]=ModelParam::omega[d2].J;

				  if(To[0]==To[1]){
					  JMat[start][To[0]]=ModelParam::omega[d1].J+ModelParam::omega[n].J;
				  }
				  else if(To[0]==From[0]){
					  JMat[start][From[0]]=ModelParam::omega[d2].J+ModelParam::omega[n].J;
					  JMat[start][To[1]]=ModelParam::omega[d1].J;
				  }
				  else if(To[1]==From[0]){
					  JMat[start][From[0]]=ModelParam::omega[d2].J+ModelParam::omega[d1].J;
					  JMat[start][To[0]]=ModelParam::omega[n].J;
				  }
				  else if(0){
					  // ѹ���߽�Ϊ�ֲ�����(CAM����)����δʵ��

				  }
				  else{
					  JMat[start][To[0]]=ModelParam::omega[n].J;
					  JMat[start][To[1]]=ModelParam::omega[d1].J;
				  }
			  }
			  else if(ModelParam::omega[d2].nodes[0][0] == start){
				  /*
				  ���d2������뵱ǰ����Ѫ�ܶε������ͬ��˵��d1��ĸѪ�ܣ�d2��n����Ѫ��
				  From[0]ΪĸѪ��d1�����
				  To[0]Ϊ��ǰѪ�ܶε��յ�
				  To[1]Ϊ��һ��Ѫ��d2���յ�
				  */
				  From[0] = ModelParam::omega[d1].nodes[0][0];
				  To[0]=end;
				  To[1]=ModelParam::omega[d2].nodes[1][0];
				  JMat[start][start]=-(ModelParam::omega[n].J+ModelParam::omega[d1].J+ModelParam::omega[d2].J);
				  JMat[start][From[0]]=ModelParam::omega[d1].J;
				  if(To[0]==To[1]){
					  JMat[start][To[0]]=ModelParam::omega[n].J+ModelParam::omega[d2].J;
				  }
				  else if(To[0]==From[0]){
					  JMat[start][From[0]]=ModelParam::omega[d1].J+ModelParam::omega[n].J;
					  JMat[start][To[1]]=ModelParam::omega[d2].J;
				  }
				  else if(To[1]==From[0]){
					  JMat[start][From[0]]=ModelParam::omega[d1].J+ModelParam::omega[d2].J;
					  JMat[start][To[0]]=ModelParam::omega[n].J;
				  }
				  else if(0){
					  // ѹ���߽�Ϊ�ֲ�����(CAM����)

				  }
				  else{
					  JMat[start][To[0]]=ModelParam::omega[n].J;
					  JMat[start][To[1]]=ModelParam::omega[d2].J;
				  }
			  }
			  break;
			  /*
			  ��ڵ�ΪC����ʱ��ֻ��һ�����������ǰѪ��nΪĸѪ�ܣ�d1,d2Ϊ��Ѫ��
			  */   
		  case 'C':
			  /*
			  From[0]Ϊ��Ѫ��d1�����
			  From[1]Ϊ��Ѫ��d2�����
			  To[0]ΪĸѪ��n���յ�
			  */
			  if(d1<0 || d2<0 || start<0){
				  AdapParam::errFlag=AdapParam::FLOWDIR;
				  return;
			  }
			  From[0]=ModelParam::omega[d1].nodes[0][0];
			  From[1]=ModelParam::omega[d2].nodes[0][0];
			  To[0]=end;
			  // AdapParam::adapLogFile << "error occur!\t" << "start=" << start << "\tn=" << n << "\td1=" << d1 << "\td2=" << d2 << endl;
			  JMat[start][start]=-(ModelParam::omega[n].J+ModelParam::omega[d1].J+ModelParam::omega[d2].J);
			  JMat[start][To[0]]=ModelParam::omega[n].J;
			  if (From[0]==From[1]){
				  JMat[start][From[0]]=ModelParam::omega[d1].J+ModelParam::omega[d2].J;
			  }
			  else if(From[0]==To[0]){
				  JMat[start][To[0]]=ModelParam::omega[n].J+ModelParam::omega[d1].J;
				  JMat[start][From[1]]=ModelParam::omega[d2].J;
			  }
			  else if(From[1]==To[0]){
				  JMat[start][To[0]]=ModelParam::omega[n].J+ModelParam::omega[d2].J;
				  JMat[start][From[1]]=ModelParam::omega[d1].J;
			  }
			  else{
				  JMat[start][From[0]]=ModelParam::omega[d1].J;
				  JMat[start][From[1]]=ModelParam::omega[d2].J;
			  }
			  break;
		  case 'J':
			  From[0]=ModelParam::omega[d1].nodes[0][0];
			  To[0]=end;
			  JMat[start][start]=-(ModelParam::omega[d1].J+ModelParam::omega[n].J);
			  JMat[start][From[0]]=ModelParam::omega[d1].J;
			  JMat[start][end]=ModelParam::omega[n].J;
			  break;
		  case 'q':
			  q_in[n]=ModelParam::omega[n].bcval[1];
			  JMat[start][start]=ModelParam::omega[n].J;
			  JMat[start][end]=-ModelParam::omega[n].J;
			  NodeFlow[start]=q_in[n];
			  break;
		  case 'p':
			  // TODO
			  break;
		  default:
			  break;
		}

		// End node
		switch (ModelParam::omega[n].bctype[3]){
			/* ĩ��ΪB����ʱ��ֻ��һ�������nΪĸѪ�ܣ�d3,d4Ϊ��Ѫ�� */
		  case 'B':
			  if(d3<0 || d4<0 || end<0){
				  AdapParam::errFlag=AdapParam::FLOWDIR;
				  return;
			  }
			  From[0]=start;
			  To[0]=ModelParam::omega[d3].nodes[1][0];
			  To[1]=ModelParam::omega[d4].nodes[1][0];
			  JMat[end][end]=-(ModelParam::omega[n].J+ModelParam::omega[d3].J+ModelParam::omega[d4].J);
			  JMat[end][From[0]]=ModelParam::omega[n].J;
			  if (To[0]==To[1]){
				  JMat[end][To[0]]=ModelParam::omega[d3].J+ModelParam::omega[d4].J;
			  }
			  else if(To[0]==From[0]){
				  JMat[end][From[0]]=ModelParam::omega[n].J+ModelParam::omega[d3].J;
				  JMat[end][To[0]]=ModelParam::omega[d4].J;
			  }
			  else if(To[1]==From[0]){
				  JMat[end][From[0]]=ModelParam::omega[n].J+ModelParam::omega[d4].J;
				  JMat[end][To[1]]=ModelParam::omega[d3].J;
			  }
			  else if(0){
				  // ѹ���߽�Ϊ�ֲ�����(CAM����)
			  }
			  else{
				  JMat[end][To[0]]=ModelParam::omega[d3].J;
				  JMat[end][To[1]]=ModelParam::omega[d4].J;
			  }
			  break;
			  /*
			  ĩ��ΪC����ʱ�������ΪB�������ƣ���ֱ�����
			  */
		  case 'C':
			  if(d3<0 || d4<0 || end<0){
				  AdapParam::errFlag=AdapParam::FLOWDIR;
				  return;
			  }
			  else if(ModelParam::omega[d3].nodes[1][0] == end){
				  From[0]=start;
				  From[1]=ModelParam::omega[d3].nodes[0][0];
				  To[0]=ModelParam::omega[d4].nodes[1][0];
				  JMat[end][end]=-(ModelParam::omega[n].J+ModelParam::omega[d3].J+ModelParam::omega[d4].J);
				  JMat[end][To[0]]=ModelParam::omega[d4].J;
				  if(From[0]==From[1]){
					  JMat[end][From[0]]=ModelParam::omega[d3].J+ModelParam::omega[n].J;
				  }
				  else if(From[0]==To[0]){
					  JMat[end][To[0]]=ModelParam::omega[d4].J+ModelParam::omega[n].J;
					  JMat[end][From[1]]=ModelParam::omega[d3].J;
				  }
				  else if(From[1]==To[0]){
					  JMat[end][To[0]]=ModelParam::omega[d4].J+ModelParam::omega[d3].J;
					  JMat[end][From[0]]=ModelParam::omega[n].J;
				  }
				  else{
					  JMat[end][From[0]]=ModelParam::omega[n].J;
					  JMat[end][From[1]]=ModelParam::omega[d3].J;
				  }
			  }
			  else if(ModelParam::omega[d4].nodes[1][0] == end){
				  From[0]=start;
				  From[1]=ModelParam::omega[d4].nodes[0][0];
				  To[0]=ModelParam::omega[d3].nodes[1][0];
				  JMat[end][end]=-(ModelParam::omega[n].J+ModelParam::omega[d3].J+ModelParam::omega[d4].J);
				  JMat[end][To[0]]=ModelParam::omega[d3].J;
				  if (From[0]==From[1]){
					  JMat[end][From[0]]=ModelParam::omega[d4].J+ModelParam::omega[n].J;
				  }
				  else if(From[0]==To[0]){
					  JMat[end][To[0]]=ModelParam::omega[d3].J+ModelParam::omega[n].J;
					  JMat[end][From[1]]=ModelParam::omega[d4].J;
				  }
				  else if(From[1]==To[0]){
					  JMat[end][To[0]]=ModelParam::omega[d3].J+ModelParam::omega[d4].J;
					  JMat[end][From[1]]=ModelParam::omega[n].J;
				  }
				  else{
					  JMat[end][From[0]]=ModelParam::omega[n].J;
					  JMat[end][From[1]]=ModelParam::omega[d4].J;
				  }
			  }
			  break;
		  case 'J':
			  From[0]=start;
			  To[0]=ModelParam::omega[d3].nodes[1][0];
			  JMat[end][end]=-(ModelParam::omega[n].J+ModelParam::omega[d3].J);
			  JMat[end][From[0]]=ModelParam::omega[n].J;
			  JMat[end][To[0]]=ModelParam::omega[d3].J;
			  break;
		  case 'q':
			  q_in[n]=ModelParam::omega[n].bcval[4];
			  JMat[end][end]=ModelParam::omega[n].J;
			  JMat[end][start]=-ModelParam::omega[n].J;
			  NodeFlow[end]=q_in[n];
			  break;
		  case 'R':
		  case 'p':
			  // Main output
			  NodeFlow[start]=-ModelParam::omega[n].bcval[4]*133*ModelParam::omega[n].J;
			  break;
		  default:
			  // printf("end default\n");
			  break;
		}
	}

	int lapack_n=ModelParam::Nnodes;
	int one=1;
	double *JMatArray=(double *)malloc(ModelParam::Nnodes*ModelParam::Nnodes*sizeof(double));
	for(int i=0; i<ModelParam::Nnodes; ++i){
		for(int j=0; j<ModelParam::Nnodes; ++j){
			JMatArray[i*ModelParam::Nnodes+j]=JMat[j][i];   // JMatArray������˳��ڷ�Ԫ��
		}
	}

	int ok;   // dgesv״̬
	int *pivot=(int *)malloc(ModelParam::Nnodes*sizeof(int));
	dgesv(&lapack_n,&one,JMatArray,&lapack_n,pivot,NodeFlow,&lapack_n,&ok);

	if(ok!=0){
		AdapParam::errFlag=AdapParam::LINEQU;
		// Clear the state
		free(JMatArray);
		free(pivot);

		for(int i=0; i<ModelParam::Nnodes; ++i)
			for (int j=0; j<ModelParam::Nnodes; ++j)
				JMat[i][j]=0;

		for(int i=0; i<ModelParam::Nnodes; ++i){
			NodeFlow[i]=0;
		}
		return;
	}

	// NodeFlow stores the X(Pressure) in AX=B
	for(n=0; n<ModelParam::Ndoms; ++n){
		start = ModelParam::omega[n].nodes[0][0];     // start node number
		end   = ModelParam::omega[n].nodes[1][0];     // end node number

		if (ModelParam::omega[n].bctype[3]=='R'){
			SS_DeltaP[n]=NodeFlow[start]-ModelParam::omega[n].bcval[4]*133;
			MeanP[n]=NodeFlow[start]-SS_DeltaP[n]/2;
		}
		else{
			SS_DeltaP[n]=NodeFlow[start]-NodeFlow[end];
			MeanP[n]=NodeFlow[start]-SS_DeltaP[n]/2;
		}
	}

	for(n=0; n<ModelParam::Ndoms; ++n){
		MeanQ[n]=60*1e12*ModelParam::omega[n].J*SS_DeltaP[n]; // nL/min
		SS_DeltaP[n]=SS_DeltaP[n]/133;
		MeanP[n]=MeanP[n]/133;                                // mmHg
	}

	// Clear the state
	free(JMatArray);
	free(pivot);

	for(int i=0; i<ModelParam::Nnodes; ++i)
		for (int j=0; j<ModelParam::Nnodes; ++j)
			JMat[i][j]=0;

	for(int i=0; i<ModelParam::Nnodes; ++i){
		NodeFlow[i]=0;
	}
}

void Adap_SS_Solver::destroySolver(){
	// TODO: Everything initialized in the initSolver() should be destroyed.
}

void Adap_SS_Solver::Write_history(Domain *omega, char *name){
	register int n;
	char     buf[BUFSIZ];
	int      start,end;

	if(!fp_his) { // The first time the routine is called, allocate space for "fp", "His_elmt" and "His_interp".
		fp_his     = (FILE **)calloc(ModelParam::Ndoms,sizeof(FILE *)); // Open a .his file.
	}

	for(n=0;n<ModelParam::Ndoms;++n){
		if(!fp_his[n]){ /* Give the name of the .in file to the .his file (plus the number of the domain if
						there are more than one). */
			if(ModelParam::Ndoms==1)
				sprintf(buf,"%s.his",strtok(name,"."));
			else
				sprintf(buf,"%s_%d.his",strtok(name,"."),n+1);

			// Create a pointer to the .out file.
			fp_his[n] = fopen(buf,"w");
			fprintf(fp_his[n],"# Steady state blood flow outfile \n");
			fprintf(fp_his[n],"# Flow, Pressure\n");

			start = omega[n].nodes[0][0];
			end   = omega[n].nodes[1][0];
			fprintf(fp_his[n],"%lg %lg\n", MeanQ[n], MeanP[n]);
			fclose(fp_his[n]);
		}
	}
}

/*!
\brief Adjust the flow directions according to the hemodynamic simulation results
*/
void Adap_SS_Solver::adjustFlowDir()
{
	int i,j,n;
	int start,end,d1,d2,d3,d4;
	int tempNode;
	for (n=0;n<ModelParam::Ndoms;++n){
		if(MeanQ[n]<0){
			// ��n��Ѫ�ܵ�start��end����(nodes[0/1][0])���ڵ������ӵ�Ѫ�ܶα�Ż���(nodes[0/1][1,2,3])
			for(j=0;j<4;j++){
				tempNode = ModelParam::omega[n].nodes[0][j];
				ModelParam::omega[n].nodes[0][j] = ModelParam::omega[n].nodes[1][j];
				ModelParam::omega[n].nodes[1][j] = tempNode;
			}

			// start��end�����ӵ�Ѫ�ܶλ���(bifur[0/1][0,1])
			for(j=0;j<2;j++){
				tempNode = ModelParam::omega[n].bifur[0][j];
				ModelParam::omega[n].bifur[0][j] = ModelParam::omega[n].bifur[1][j];
				ModelParam::omega[n].bifur[1][j] = tempNode;
			}

			start = ModelParam::omega[n].nodes[0][0];     // start node number
			end   = ModelParam::omega[n].nodes[1][0];     // end node number
			d1    = ModelParam::omega[n].bifur[0][0];     // start node�����ӵ�����Ѫ�ܵ�id
			d2    = ModelParam::omega[n].bifur[0][1];
			d3    = ModelParam::omega[n].bifur[1][0];     // end node�����ӵ�����Ѫ�ܵ�id
			d4    = ModelParam::omega[n].bifur[1][1];

			// ���ΪB��, TODO:���б�������?
			// cout << "segid=" << n << "\tstart type: " << ModelParam::omega[n].bctype[0] << "\tend type: " << ModelParam::omega[n].bctype[3] << endl;
			if(d1<0 || d2<0){
				AdapParam::errFlag=AdapParam::FLOWDIR;
				return;
			}
			// ���d1�������start��ͬ��˵��d1���ΪB�ͣ�d2�յ�ΪB��
			else if(ModelParam::omega[d1].nodes[0][0] == start){
				for(i=0;i<3;i++){
					ModelParam::omega[d1].bctype[i]='B';
					ModelParam::omega[d2].bctype[i+3]='B';
				}
			}
			else if(ModelParam::omega[d2].nodes[0][0] == start){
				for(i=0;i<3;i++){
					ModelParam::omega[d2].bctype[i]='B';
					ModelParam::omega[d1].bctype[i+3]='B';
				}
			}

			// �յ�ΪC��, TODO
			// ���d3�յ���end��ͬ��˵��d3�յ�ΪC�ͺţ�d4���ΪC��
			if(d3<0 || d4<0){
				AdapParam::errFlag=AdapParam::FLOWDIR;
				return;
			}
			else if(ModelParam::omega[d3].nodes[1][0] == end){
				for(i=0;i<3;i++){
					ModelParam::omega[d4].bctype[i]='C';
					ModelParam::omega[d3].bctype[i+3]='C';
				}
			}
			else if(ModelParam::omega[d4].nodes[1][0] == end){
				for(i=0;i<3;i++){
					ModelParam::omega[d3].bctype[i]='C';
					ModelParam::omega[d4].bctype[i+3]='C';
				}
			}

			// ������ֵ����Ϊ��ֵ
			MeanQ[n]      = -MeanQ[n];
			SS_DeltaP[n]  = -SS_DeltaP[n];
		}
	}

	// cout << "adjust finished" << endl;
}

/*!
\brief ����Hd����˳��
ɨ������Ѫ�ܣ�������ϼ�Ѫ���Ѿ���������¼�Ѫ������
���򣬵ȴ���һ��ɨ��
ֱ������Ѫ�ܶ��������
*/
void Adap_SS_Solver::adjustHdCalcPOrder()
{
	int n, cnt=0;
	int start, end, d1, d2, d3, d4;

	// ÿ�ζ����PosCalcOrder
	if(!AdapParam::PosCalcOrder.empty())
		AdapParam::PosCalcOrder.clear();

	for(n=0;n<ModelParam::Ndoms;n++)
		Pmarker[n]=0;

	while(AdapParam::PosCalcOrder.size()<ModelParam::Ndoms){
		if(++cnt>ModelParam::Ndoms/5){
			AdapParam::errFlag=AdapParam::PORDER;
			return;
		}
		for(n=0;n<ModelParam::Ndoms;n++){
			start = ModelParam::omega[n].nodes[0][0];     // start node number
			end   = ModelParam::omega[n].nodes[1][0];     // end node number
			d1    = ModelParam::omega[n].bifur[0][0];     // start node�����ӵ�����Ѫ�ܵ�id
			d2    = ModelParam::omega[n].bifur[0][1];
			d3    = ModelParam::omega[n].bifur[1][0];     // end node�����ӵ�����Ѫ�ܵ�id
			d4    = ModelParam::omega[n].bifur[1][1];

			// start node
			switch(ModelParam::omega[n].bctype[0]){
				/*
				��ڵ�ΪC����ʱ, �ö�Ѫ��һ��Ϊ��۽ڵ��ĸѪ��
				*/   
			  case 'C':
				  if(d1<0 || d1>=ModelParam::Ndoms || d2<0 || d2>=ModelParam::Ndoms){
					  AdapParam::errFlag=AdapParam::FLOWDIR;
					  return;
				  }
				  else if(Pmarker[d1]==1 && Pmarker[d2]==1){
					  // ���������Ѫ������һ�λ�δ������������
					  if(Pmarker[n]==0){
						  Pmarker[n]=1;
						  AdapParam::PosCalcOrder.push_back(n);
					  }
				  }
				  break;
			  case 'q':
			  case 'p':
				  /*
				  ���start node����Ϊq or p, ˵�������Ѫ��, ���Ϊ���ȼ���
				  */
				  if(Pmarker[n]==0){
					  Pmarker[n]=1;
					  AdapParam::PosCalcOrder.push_back(n);
				  }
				  break;
			  default:
				  break;
			}

			// End node
			switch (ModelParam::omega[n].bctype[3]){
				/*
				ĩ��ΪB����ʱ��ֻ��һ�������nΪĸѪ�ܣ�d3,d4Ϊ��Ѫ��
				��ʱ��1. ���n�ļ���˳�򣨲�Ҫ�ظ���ǣ� 2. �����������Ѫ�ܵļ���˳�򣨲�Ҫ�ظ���ǣ�
				*/
			  case 'B':
				  if(d3<0 || d3>=ModelParam::Ndoms || d4<0 || d4>=ModelParam::Ndoms){
					  AdapParam::errFlag=AdapParam::FLOWDIR;
					  return;
				  }
				  if(Pmarker[n]==1){
					  if(Pmarker[d3]==0){
						  Pmarker[d3]=1;
						  AdapParam::PosCalcOrder.push_back(d3);
					  }
					  if(Pmarker[d4]==0){
						  Pmarker[d4]=1;
						  AdapParam::PosCalcOrder.push_back(d4);
					  }
				  }
				  break;
			  case 'J':
				  if(d3<0 || d3>=ModelParam::Ndoms){
					  AdapParam::errFlag=AdapParam::FLOWDIR;
					  return;
				  }
				  /* ���n�ѱ�Ǽ���˳������d3 */
				  if(Pmarker[n]==1){
					  if(Pmarker[d3]==0){
						  Pmarker[d3]=1;
						  AdapParam::PosCalcOrder.push_back(d3);
					  }
				  }
				  break;
				  /*case 'q':
				  case 'p':
				  case 'R':
				  if(Pmarker[n]==0){
				  Pmarker[n]=1;
				  AdapParam::PosCalcOrder.push_back(n);
				  }
				  break;*/
			  default:
				  // printf("end default\n");
				  break;
			}
		}
	}
}

/*!
\brief ����Sc����˳��
Sc�ļ���˳��Ϊ����˳��
ɨ������Ѫ�ܣ�������¼�Ѫ���Ѿ���������ϼ�Ѫ������
���򣬵ȴ���һ��ɨ�裬ֱ������Ѫ�ܶ��������
*/
void Adap_SS_Solver::adjustScCalcOrder()
{
	int n, cnt=0;
	int start, end, d1, d2, d3, d4;

	if(!AdapParam::NegCalcOrder.empty())
		AdapParam::NegCalcOrder.clear();

	for(n=0;n<ModelParam::Ndoms;n++)
		Nmarker[n]=0;

	while(AdapParam::NegCalcOrder.size()<ModelParam::Ndoms){
		if(++cnt>ModelParam::Ndoms/5){
			AdapParam::errFlag=AdapParam::NORDER;
			break;
		}
		for(n=ModelParam::Ndoms-1;n>=0;n--){
			start = ModelParam::omega[n].nodes[0][0];     // start node number
			end   = ModelParam::omega[n].nodes[1][0];     // end node number
			d1    = ModelParam::omega[n].bifur[0][0];     // start node�����ӵ�����Ѫ�ܵ�id
			d2    = ModelParam::omega[n].bifur[0][1];
			d3    = ModelParam::omega[n].bifur[1][0];     // end node�����ӵ�����Ѫ�ܵ�id
			d4    = ModelParam::omega[n].bifur[1][1];

			// End node
			switch (ModelParam::omega[n].bctype[3]){
				/* ĩ��ΪB����ʱ������ΪC�� */
			  case 'B':
				  if(Nmarker[d3]==1 && Nmarker[d4]==1){
					  if(Nmarker[n]==0){
						  Nmarker[n]=1;
						  AdapParam::NegCalcOrder.push_back(n);
					  }
				  }
				  break;
				  /* ĩ��Ϊq,p,R��Ϊ�߽� */
			  case 'q':
			  case 'p':
			  case 'R':
				  if(Nmarker[n]==0){
					  Nmarker[n]=1;
					  AdapParam::NegCalcOrder.push_back(n);
				  }
				  break;
			  default:
				  break;
			}

			// start node
			switch(ModelParam::omega[n].bctype[0]){
				/* ���start node������ΪC, �������������ԣ�Ϊ�ֲ� */
			  case 'C':
				  if(Nmarker[n]==1){
					  if(Nmarker[d1]==0){
						  Nmarker[d1]=1;
						  AdapParam::NegCalcOrder.push_back(d1);
					  }
					  if(Nmarker[d2]==0){
						  Nmarker[d2]=1;
						  AdapParam::NegCalcOrder.push_back(d2);
					  }
				  }
				  break;
			  case 'J':
				  // �������Ѫ���Ѿ�����˼���˳������ǰ�εļ���˳��
				  if(Nmarker[n]==1){
					  if(Nmarker[d1]==0){
						  Nmarker[d1]=1;
						  AdapParam::NegCalcOrder.push_back(d1);
					  }
				  }
				  break;
				  /*case 'q':
				  case 'p':
				  if(Nmarker[n]==0){
				  Nmarker[n]=1;
				  AdapParam::NegCalcOrder.push_back(n);
				  }
				  break;*/
			  default:
				  break;
			}
		}
	}
}

/*
����Hd�ļ���˳�򣬼���Hd
*/
void Adap_SS_Solver::calcHd()
{
	int i,j,n;
	int start,end,d1,d2,d3,d4;
	double FlowRatio1, FlowRatio2, FQe1, FQe2;
	double Df,Da,Db;

	for(n=0;n<ModelParam::Ndoms;n++){
		for(i=0;i<ModelParam::oneD_q-1;i++)
			ModelParam::omega[n].Hd[0].h[i]=ModelParam::omega[n].BHd;
	}

	for (i=0;i<ModelParam::Ndoms;i++){
		// ����˳������n��Ѫ��
		n=AdapParam::PosCalcOrder[i];
		start = ModelParam::omega[n].nodes[0][0];     // start node number
		end   = ModelParam::omega[n].nodes[1][0];     // end node number
		d1    = ModelParam::omega[n].bifur[0][0];     // start node�����ӵ�����Ѫ�ܵ�id
		d2    = ModelParam::omega[n].bifur[0][1];
		d3    = ModelParam::omega[n].bifur[1][0];     // end node�����ӵ�����Ѫ�ܵ�id
		d4    = ModelParam::omega[n].bifur[1][1];

		// start node
		switch(ModelParam::omega[n].bctype[0]){
		  case 'q':
		  case 'p':
			  /* ��������Ѫ�� */
			  for(j=0;j<ModelParam::oneD_q-1;j++)
				  ModelParam::omega[n].Hd[0].h[j]=ModelParam::omega[n].BHd;
			  break; 
		  case 'C':
			  for(j=0;j<ModelParam::oneD_q-1;j++)
				  ModelParam::omega[n].Hd[0].h[j]=(ModelParam::omega[d1].Hd[0].h[j]*MeanQ[d1]+ModelParam::omega[d2].Hd[0].h[j]*MeanQ[d2])/MeanQ[n];
			  break;
		  default:
			  break;
		}

		// end node
		switch(ModelParam::omega[n].bctype[3]){
			/*
			ĩ��ΪB����ʱ��ֻ��һ�������nΪĸѪ�ܣ�d3,d4Ϊ��Ѫ��
			��ʹ��Phase separation effect����Hd
			*/
		  case 'B':
			  FlowRatio1=MeanQ[d3]/MeanQ[n];
			  FlowRatio2=MeanQ[d4]/MeanQ[n];
			  Df=sqrt(4*ModelParam::omega[n].A[0].h[0]/ModelParam::PI)*1e6;
			  Da=sqrt(4*ModelParam::omega[d3].A[0].h[0]/ModelParam::PI)*1e6;
			  Db=sqrt(4*ModelParam::omega[d4].A[0].h[0]/ModelParam::PI)*1e6;
			  FQe1=PS_Effect(Df,Da,Db,ModelParam::omega[n].Hd[0].h[0],FlowRatio1);
			  FQe2=PS_Effect(Df,Db,Da,ModelParam::omega[n].Hd[0].h[0],FlowRatio2);
			  // �ж�ʵ��
			  if(logerr){
				  if(FlowRatio1>FlowRatio2){
					  FQe1=1;
					  FQe2=0;
				  }
				  else{
					  FQe1=0;
					  FQe2=1;
				  }
				  logerr=0;
			  }
			  for(j=0;j<ModelParam::oneD_q-1;j++){
				  ModelParam::omega[d3].Hd[0].h[j]=ModelParam::omega[n].Hd[0].h[0]*FQe1/FlowRatio1;
				  ModelParam::omega[d4].Hd[0].h[j]=ModelParam::omega[n].Hd[0].h[0]*FQe2/FlowRatio2;
			  }
			  break;
		  case 'J':
			  for(j=0;j<ModelParam::oneD_q-1;j++){
				  ModelParam::omega[d3].Hd[0].h[j]=ModelParam::omega[n].Hd[0].h[j];
			  }
			  break;
		  default:
			  break;
		}
	}

	for(n=0;n<ModelParam::Ndoms;n++)
		Debug_Hd[n]=ModelParam::omega[n].Hd[0].h[0];
}

/**
\brief
Ref. 1994 Resistance to blood flow in microvessels in vivo
*/
double Adap_SS_Solver::Eta_vitro(double Hd, double diam){
	double C = Eta_C(diam);
	double visc_045 = Eta_045(diam);
	double Eta_vitro = 1+(visc_045-1)*(pow(1-Hd,C)-1)/(pow(1-0.45,C)-1);
	return Eta_vitro;
}

double Adap_SS_Solver::Eta_045(double diam){
	return 220*exp(-1.3*diam)+3.2-2.44*exp(-0.06*pow(diam,0.645));
}

double Adap_SS_Solver::Eta_C(double diam){
	return (0.8+exp(-0.075*diam))*(-1+1/(1+pow(double(10),-11)+pow(diam,12)))+1/(1+pow(double(10),-11)+pow(diam,12));
}

double Adap_SS_Solver::Eta_Was(double diam){
	double Doff = 2.4;	// um
	double D50 = 100;		// um
	double Wmax = 2.6;	// um, depends on the data

	if(diam > Doff)
		return (diam-Doff)/(diam+D50-2*Doff)*Wmax;
	else
		return 0;
}

double Adap_SS_Solver::Eta_Wpeak(double diam){
	double Doff = 2.4;	// um
	double Dcrit = 10.5	;	// um
	double Eamp = 1.1;
	double Ewidth = 0.03;
	double Ehd = 1.18;

	if(diam>Dcrit)
		return Eamp*exp(-Ewidth*(diam-Dcrit));
	else if(diam>Doff && diam<=Dcrit)
		return Eamp*(diam-Doff)/(Dcrit-Doff);
	else
		return 0;
}

/*
Ref. 1994 Resistance to blood flow in microvessels in vivo
*/
double Adap_SS_Solver::Eta_vivo(double Hd, double diam, double W){
	double C=(0.8+exp(-0.075*diam))*(-1+1/(1+10e-11*pow(diam,12)))+1/(1+10e-11*pow(diam,12));
	double eta_45=6*exp(-0.085*diam)+3.2-2.44*exp(-0.06*pow(diam,0.645));
	return (1+(eta_45-1)*(pow(1-Hd,C)-1)/(pow(1-0.45,C)-1)*pow(diam/(diam-W),2))*pow(diam/(diam-W),2);
}

/*
Ref. 2005 Microvascular blood viscosity in vivo and the endothelial surface layer
*/
double Adap_SS_Solver::Eta_vivo_ESL(double Hd, double diam){
	double Deff = Eta_Deff(Hd, diam);
	double Dph = Eta_Dph(diam);
	double eta_vitro = Eta_vitro(Hd, Dph);
	return eta_vitro*pow(diam/Deff, 4);
}

double Adap_SS_Solver::Eta_Deff(double Hd, double diam){
	double Ehd = 1.18;
	double Weff = Eta_Was(diam)+Eta_Wpeak(diam)*(1+Hd*Ehd);
	return diam-2*Weff;
}

double Adap_SS_Solver::Eta_Dph(double diam){
	double Epeak = 0.6;
	double Wph = Eta_Was(diam)+Eta_Wpeak(diam)*Epeak;
	return diam-2*Wph;
}

/*
\brief Phase separation effect
*/
double Adap_SS_Solver::PS_Effect(double Df,double Da,double Db,double Hd,double FQb){
	double A,B,X0,FQe;
	A=-13.29*((Da*Da-Db*Db)/(Da*Da+Db*Db))*(1-Hd)/Df;
	B=1+6.98*(1-Hd)/Df;
	X0=0.964*(1-Hd)/Df;
	FQe=1/(1+exp(-(A+B*logit((FQb-X0)/(1-2*X0)))));
	return FQe;
}

double Adap_SS_Solver::logit(double x){
	double y;
	if(x/(1-x)>0)
		y=log(x/(1-x));
	else{
		y=0;
		logerr=1;   // log�����־��Ϊ1
	}
	return y;
}

/**
\brief 
Update the viscosity according to diameter of the current step and the fixed Hd
*/
void Adap_SS_Solver::UpdateVisc(Domain *omega){
	int i,j,k,n;
	Element *A, *U, *Hd;
	double *visc, *oldAhis, *oldUhis, *oldHdhis;

	for(n=0; n<ModelParam::Ndoms; ++n){
		A = omega[n].A;
		U = omega[n].U;
		Hd = omega[n].Hd;
		visc = omega[n].visc;

		if(ModelParam::nonDim){
			for(i=0; i<A[0].q; ++i){
				if(ModelParam::updVisc==1)
					visc[i] = Eta_vitro(Hd[0].h[i], 2*sqrt(A[0].h[i]*ModelParam::scale_r0*ModelParam::scale_r0/ModelParam::PI)*1e6)/1000;
				else if(ModelParam::updVisc==2)
					visc[i] = Eta_vivo_ESL(Hd[0].h[i], 2*sqrt(A[0].h[i]*ModelParam::scale_r0*ModelParam::scale_r0/ModelParam::PI)*1e6)/1000;
				else if(ModelParam::updVisc==3)
					visc[i] = Eta_vivo(Hd[0].h[i], 2*sqrt(A[0].h[i]*ModelParam::scale_r0*ModelParam::scale_r0/ModelParam::PI)*1e6, 0.85)/1000;
			}
		}
		else{
			for(i=0; i<A[0].q; ++i){
				if(ModelParam::updVisc==1)
					visc[i] = Eta_vitro(Hd[0].h[i], 2*sqrt(A[0].h[i]/ModelParam::PI)*1e6)/1000;
				else if(ModelParam::updVisc==2)
					visc[i] = Eta_vivo_ESL(Hd[0].h[i], 2*sqrt(A[0].h[i]/ModelParam::PI)*1e6)/1000;
				else if(ModelParam::updVisc==3)
					visc[i] = Eta_vivo(Hd[0].h[i], 2*sqrt(A[0].h[i]/ModelParam::PI)*1e6, 0.85)/1000;
			}
		}
	}
}

/*
\brief Calculate PO2
*/
void Adap_SS_Solver::calcPO2(){
	int n,i,j,k;
	double Len;
	vector<double> MeanFlow(ModelParam::Ndoms);
	int start,end,d1,d2,d3,d4;
	double MPO2=4e-11;  // ƽ������

	// ���ݼ���˳�����SO2
	for (n=0;n<ModelParam::Ndoms;n++){
		// ���߽�ֵBSO2��ֵ��SO2
		ModelParam::omega[n].SO2[0]=ModelParam::omega[n].BSO2;
		MeanFlow[n]=MeanQ[n]/60/1e6;    // ��λת��Ϊcm3/s
		Len=(ModelParam::omega[n].U->x[1]-ModelParam::omega[n].U->x[0])*1e3;     // mm
		MPO[n]=MPO2*1e4/60*Len/10;  // ���ĵ�λת��O2cm3.cm3/s

		// ���ڱ߽�Ѫ�ܣ�����JconMin, JconMmid, JconMout
		if(ModelParam::omega[n].BSO2>0){
			ModelParam::omega[n].JconM[0]=0.5*MeanFlow[n]*ModelParam::omega[n].Hd[0].h[0]*ModelParam::omega[n].SO2[0];
			ModelParam::omega[n].JconM[1]=ModelParam::omega[n].JconM[0]-MPO[n]/2;
			ModelParam::omega[n].JconM[2]=ModelParam::omega[n].JconM[0]-MPO[n];
			ModelParam::omega[n].SO2[1]=ModelParam::omega[n].JconM[1]/(0.5*MeanFlow[n]*ModelParam::omega[n].Hd[0].h[0]);
			ModelParam::omega[n].SO2[2]=ModelParam::omega[n].JconM[2]/(0.5*MeanFlow[n]*ModelParam::omega[n].Hd[0].h[0]);
			if(ModelParam::omega[n].SO2[1]<0){
				ModelParam::omega[n].SO2[1]=0;
				ModelParam::omega[n].SO2[2]=0;
				ModelParam::omega[n].PO2=0;
			}
			else{
				ModelParam::omega[n].PO2=Hill(ModelParam::omega[n].SO2[1]);
			}
		}
	}

	// ���ݼ���˳�����SO2
	for (i=0;i<ModelParam::Ndoms;i++){
		// ����˳������n��Ѫ��
		n=AdapParam::PosCalcOrder[i];

		start = ModelParam::omega[n].nodes[0][0];     // start node number
		end   = ModelParam::omega[n].nodes[1][0];     // end node number
		d1    = ModelParam::omega[n].bifur[0][0];     // start node�����ӵ�����Ѫ�ܵ�id
		d2    = ModelParam::omega[n].bifur[0][1];
		d3    = ModelParam::omega[n].bifur[1][0];     // end node�����ӵ�����Ѫ�ܵ�id
		d4    = ModelParam::omega[n].bifur[1][1];

		// start node
		switch(ModelParam::omega[n].bctype[0]){
		  case 'C':
			  ModelParam::omega[n].JconM[0]=MeanFlow[d1]*0.5*ModelParam::omega[d1].Hd[0].h[0]*ModelParam::omega[d1].SO2[2]+
				  MeanFlow[d2]*0.5*ModelParam::omega[d2].Hd[0].h[0]*ModelParam::omega[d2].SO2[2];
			  ModelParam::omega[n].JconM[1]=ModelParam::omega[n].JconM[0]-MPO[n]/2;
			  ModelParam::omega[n].JconM[2]=ModelParam::omega[n].JconM[0]-MPO[n];
			  ModelParam::omega[n].SO2[0]=ModelParam::omega[n].JconM[0]/(MeanFlow[n]*0.5*ModelParam::omega[n].Hd[0].h[0]);
			  ModelParam::omega[n].SO2[1]=ModelParam::omega[n].JconM[1]/(MeanFlow[n]*0.5*ModelParam::omega[n].Hd[0].h[0]);
			  ModelParam::omega[n].SO2[2]=ModelParam::omega[n].JconM[2]/(MeanFlow[n]*0.5*ModelParam::omega[n].Hd[0].h[0]);
			  if(ModelParam::omega[n].Hd[0].h[0]==0){
				  ModelParam::omega[n].SO2[0]=0;
				  ModelParam::omega[n].SO2[2]=0;
				  ModelParam::omega[n].PO2=0;
			  }
			  else if(ModelParam::omega[n].SO2[1]<0){
				  ModelParam::omega[n].SO2[1]=0;
				  ModelParam::omega[n].SO2[2]=0;
				  ModelParam::omega[n].PO2=0;
			  }
			  else if(ModelParam::omega[n].SO2[1]>1){
				  ModelParam::omega[n].PO2=95;      // mmHg
			  }
			  else{
				  ModelParam::omega[n].PO2=Hill(ModelParam::omega[n].SO2[1]);
			  }
			  break;
		  default:
			  break;
		}

		// end node
		switch(ModelParam::omega[n].bctype[3]){
			/* ĩ��ΪB����ʱ��ֻ��һ�������nΪĸѪ�ܣ�d3,d4Ϊ��Ѫ�� */
		  case 'B':
			  ModelParam::omega[d3].SO2[0]=ModelParam::omega[n].SO2[2];
			  ModelParam::omega[d4].SO2[0]=ModelParam::omega[n].SO2[2];
			  ModelParam::omega[d3].JconM[0]=MeanFlow[d3]*0.5*ModelParam::omega[d3].Hd[0].h[0]*ModelParam::omega[d3].SO2[0];
			  ModelParam::omega[d4].JconM[0]=MeanFlow[d4]*0.5*ModelParam::omega[d4].Hd[0].h[0]*ModelParam::omega[d4].SO2[0];
			  ModelParam::omega[d3].JconM[1]=ModelParam::omega[d3].JconM[0]-MPO[d3]/2;
			  ModelParam::omega[d4].JconM[1]=ModelParam::omega[d4].JconM[0]-MPO[d4]/2;
			  ModelParam::omega[d3].JconM[2]=ModelParam::omega[d3].JconM[0]-MPO[d3];
			  ModelParam::omega[d4].JconM[2]=ModelParam::omega[d4].JconM[0]-MPO[d4];
			  ModelParam::omega[d3].SO2[1]=ModelParam::omega[d3].JconM[1]/(MeanFlow[d3]*0.5*ModelParam::omega[d3].Hd[0].h[0]);
			  ModelParam::omega[d4].SO2[1]=ModelParam::omega[d4].JconM[1]/(MeanFlow[d4]*0.5*ModelParam::omega[d4].Hd[0].h[0]);
			  ModelParam::omega[d3].SO2[2]=ModelParam::omega[d3].JconM[2]/(MeanFlow[d3]*0.5*ModelParam::omega[d3].Hd[0].h[0]);
			  ModelParam::omega[d4].SO2[2]=ModelParam::omega[d4].JconM[2]/(MeanFlow[d4]*0.5*ModelParam::omega[d4].Hd[0].h[0]);
			  if(ModelParam::omega[d3].Hd[0].h[0]==0){
				  ModelParam::omega[d3].SO2[1]=0;
				  ModelParam::omega[d3].SO2[2]=0;
			  }
			  if(ModelParam::omega[d4].Hd[0].h[0]==0){
				  ModelParam::omega[d4].SO2[1]=0;
				  ModelParam::omega[d4].SO2[2]=0;
			  }
			  ModelParam::omega[d3].PO2=Hill(ModelParam::omega[d3].SO2[1]);
			  ModelParam::omega[d4].PO2=Hill(ModelParam::omega[d4].SO2[1]);

			  if(ModelParam::omega[d3].SO2[1]<0){
				  ModelParam::omega[d3].SO2[2]=0;
				  ModelParam::omega[d3].PO2=0;
			  }
			  if(ModelParam::omega[d4].SO2[1]<0){
				  ModelParam::omega[d4].SO2[2]=0;
				  ModelParam::omega[d4].PO2=0;
			  }
			  break;
		  case 'J':
			  ModelParam::omega[d3].SO2[0]=ModelParam::omega[n].SO2[2];
			  ModelParam::omega[d3].JconM[0]=MeanFlow[d3]*0.5*ModelParam::omega[d3].Hd[0].h[0]*ModelParam::omega[d3].SO2[0];
			  ModelParam::omega[d3].JconM[1]=ModelParam::omega[d3].JconM[0]-MPO[d3]/2;
			  ModelParam::omega[d3].JconM[2]=ModelParam::omega[d3].JconM[0]-MPO[d3];
			  ModelParam::omega[d3].SO2[1]=ModelParam::omega[d3].JconM[1]/(MeanFlow[d3]*0.5*ModelParam::omega[d3].Hd[0].h[0]);
			  ModelParam::omega[d3].SO2[2]=ModelParam::omega[d3].JconM[2]/(MeanFlow[d3]*0.5*ModelParam::omega[d3].Hd[0].h[0]);
			  ModelParam::omega[d3].PO2=Hill(ModelParam::omega[d3].SO2[1]);
			  if(ModelParam::omega[n].SO2[1]<0){
				  ModelParam::omega[n].SO2[1]=0;
				  ModelParam::omega[n].SO2[2]=0;
				  ModelParam::omega[n].PO2=0;
			  }
			  else if(ModelParam::omega[n].SO2[1]>1)
				  ModelParam::omega[n].PO2=95;
			  else
				  ModelParam::omega[n].PO2=Hill(ModelParam::omega[n].SO2[1]);
			  break;
		  default:
			  break;
		}

		// printf("n=%5d, SO2=%f, PO2=%f\n", n+1, ModelParam::omega[n].SO2[1], ModelParam::omega[n].PO2);
	}

	for(n=0;n<ModelParam::Ndoms;n++){
		Debug_SO2[n]=ModelParam::omega[n].SO2[1];
		Debug_PO2[n]=ModelParam::omega[n].PO2;
	}
}

/*
\brief Calculate Sm
*/
void Adap_SS_Solver::calcSm(){
	int n,i,j,k;
	int start,end,d1,d2,d3,d4;
	double PO2;
	vector<double> Jm(ModelParam::Ndoms);
	vector<double> Len(ModelParam::Ndoms);
	double M0=1000;   // TODO: ����ʲô?

	// ���ݼ���˳�����SO2
	for (n=0;n<ModelParam::Ndoms;n++){
		PO2 = ModelParam::omega[n].PO2;
		Jm[n]=1e-10;  // ���ü�С�ķ����ֵ
		Len[n]=(ModelParam::omega[n].U->x[1]-ModelParam::omega[n].U->x[0])*1e6;

		// ���ڱ߽�Ѫ�ܣ�����Jm
		if(ModelParam::omega[n].BJm>0){
			if(PO2<AdapParam::PO2Ref[0]){
				Jm[n]=(ModelParam::omega[n].BJm+Len[n]*(1-PO2/AdapParam::PO2Ref[0]))*exp(-AdapParam::dt/M0);
				ModelParam::omega[n].Jmmid=ModelParam::omega[n].BJm+0.5*Len[n]*(1-PO2/AdapParam::PO2Ref[0]);
			}
			else{
				Jm[n]=ModelParam::omega[n].BJm*exp(-AdapParam::dt/M0);
				ModelParam::omega[n].Jmmid=ModelParam::omega[n].BJm;
			}
		}
	}

	// ���ݼ���˳�����SO2
	for (i=0;i<ModelParam::Ndoms;i++){
		// ����˳������n��Ѫ��
		n=AdapParam::PosCalcOrder[i];

		start = ModelParam::omega[n].nodes[0][0];     // start node number
		end   = ModelParam::omega[n].nodes[1][0];     // end node number
		d1    = ModelParam::omega[n].bifur[0][0];     // start node�����ӵ�����Ѫ�ܵ�id
		d2    = ModelParam::omega[n].bifur[0][1];
		d3    = ModelParam::omega[n].bifur[1][0];     // end node�����ӵ�����Ѫ�ܵ�id
		d4    = ModelParam::omega[n].bifur[1][1];

		PO2   = ModelParam::omega[n].PO2;

		// start node
		switch(ModelParam::omega[n].bctype[0]){
		  case 'C':
			  if(PO2<AdapParam::PO2Ref[0]){
				  Jm[n]=Jm[d1]+Jm[d2]+Len[n]*(1-PO2/AdapParam::PO2Ref[0])*exp(-AdapParam::dt/M0);
				  ModelParam::omega[n].Jmmid=Jm[d1]+Jm[d2]+0.5*Len[n]*(1-PO2/AdapParam::PO2Ref[0]);
			  }
			  else{
				  Jm[n]=(Jm[d1]+Jm[d2])*exp(-AdapParam::dt/M0);
				  ModelParam::omega[n].Jmmid=Jm[d1]+Jm[d2];
			  }
			  break;
		  default:
			  break;
		}

		// end node
		switch(ModelParam::omega[n].bctype[3]){
			/* ĩ��ΪB����ʱ��ֻ��һ�������nΪĸѪ�ܣ�d3,d4Ϊ��Ѫ�� */
		  case 'B':
			  if(ModelParam::omega[d3].PO2>AdapParam::PO2Ref[0]){
				  Jm[d3]=Jm[n]*MeanQ[d3]/MeanQ[n]*exp(-AdapParam::dt/M0);
				  ModelParam::omega[d3].Jmmid=Jm[n]*MeanQ[d3]/MeanQ[n];
			  }
			  else{
				  Jm[d3]=(Jm[n]*MeanQ[d3]/MeanQ[n]+Len[d3]*(1-ModelParam::omega[d3].PO2/AdapParam::PO2Ref[0]))*exp(-AdapParam::dt/M0);
				  ModelParam::omega[d3].Jmmid=Jm[n]*MeanQ[d3]/MeanQ[n]+0.5*Len[d3]*(1-ModelParam::omega[d3].PO2/AdapParam::PO2Ref[0]);
			  }

			  if(ModelParam::omega[d4].PO2>AdapParam::PO2Ref[0]){
				  Jm[d4]=Jm[n]*MeanQ[d4]/MeanQ[n]*exp(-AdapParam::dt/M0);
				  ModelParam::omega[d4].Jmmid=Jm[n]*MeanQ[d4]/MeanQ[n];
			  }
			  else{
				  Jm[d4]=(Jm[n]*MeanQ[d4]/MeanQ[n]+Len[d4]*(1-ModelParam::omega[d4].PO2/AdapParam::PO2Ref[0]))*exp(-AdapParam::dt/M0);
				  ModelParam::omega[d4].Jmmid=Jm[n]*MeanQ[d4]/MeanQ[n]+0.5*Len[d4]*(1-ModelParam::omega[d4].PO2/AdapParam::PO2Ref[0]);
			  }
			  break;
		  case 'J':
			  if(ModelParam::omega[d3].PO2<AdapParam::PO2Ref[0]){
				  Jm[d3]=(Jm[n]+Len[d3]*(1-ModelParam::omega[d3].PO2/AdapParam::PO2Ref[0]))*exp(-AdapParam::dt/M0);
				  ModelParam::omega[d3].Jmmid=Jm[n]+0.5*Len[d3]*(1-ModelParam::omega[d3].PO2/AdapParam::PO2Ref[0]);
			  }
			  else{
				  Jm[d3]=Jm[n]*exp(-AdapParam::dt/M0);
				  ModelParam::omega[d3].Jmmid=Jm[n];
			  }
			  break;
		  default:
			  break;
		}
	}

	// Calculate Sm from Jm
	for(n=0;n<ModelParam::Ndoms;++n){
		ModelParam::omega[n].Sm=log10(1+ModelParam::omega[n].Jmmid/(MeanQ[n]+AdapParam::QRef));
		Debug_Sm[n]=ModelParam::omega[n].Sm;
	}
}

/*
\brief Calculate Sc
*/
void Adap_SS_Solver::calcSc(){
	int n,i,j,k;
	int start,end,d1,d2,d3,d4;
	vector<double> Jc(ModelParam::Ndoms);
	vector<double> Len(ModelParam::Ndoms);

	// ����Jc�߽�
	for (n=0;n<ModelParam::Ndoms;n++){
		Jc[n]=1e-10;  // ���ü�С�ķ����ֵ
		Len[n]=(ModelParam::omega[n].U->x[1]-ModelParam::omega[n].U->x[0])*1e6;   // um

		// ���ڱ߽�Ѫ�ܣ�����Jc
		if(ModelParam::omega[n].BJc>0){
			Jc[n]=(ModelParam::omega[n].BJc+Len[n]*ModelParam::omega[n].Sm)*exp(-Len[n]/AdapParam::LRef[0]);
			ModelParam::omega[n].Jcmid=ModelParam::omega[n].BJc*exp(-Len[n]/AdapParam::LRef[0])+0.5*Len[n]*ModelParam::omega[n].Sm*exp(-Len[n]/AdapParam::LRef[0]);
		}
	}

	// ���ݼ���˳�����Sc
	for (i=0;i<ModelParam::Ndoms;i++){
		// ����˳������n��Ѫ��
		n=AdapParam::NegCalcOrder[i];

		start = ModelParam::omega[n].nodes[0][0];     // start node number
		end   = ModelParam::omega[n].nodes[1][0];     // end node number
		d1    = ModelParam::omega[n].bifur[0][0];     // start node�����ӵ�����Ѫ�ܵ�id
		d2    = ModelParam::omega[n].bifur[0][1];
		d3    = ModelParam::omega[n].bifur[1][0];     // end node�����ӵ�����Ѫ�ܵ�id
		d4    = ModelParam::omega[n].bifur[1][1];

		// NOTICE: �����������ĩ�㣬�������
		// end node
		switch(ModelParam::omega[n].bctype[3]){
			/* ĩ��ΪB����ʱ������Sc����Ϊ��۽ڵ� */
		  case 'B':
			  Jc[n]=(Jc[d3]+Jc[d4]+Len[n]*ModelParam::omega[n].Sm)*exp(-Len[n]/AdapParam::LRef[0]);
			  ModelParam::omega[n].Jcmid=(Jc[d3]+Jc[d4])*exp(-0.5*Len[n]/AdapParam::LRef[0])+0.5*Len[n]*ModelParam::omega[n].Sm*exp(-0.5*Len[n]/AdapParam::LRef[0]);
			  break;
		  default:
			  break;
		}

		// start node
		switch(ModelParam::omega[n].bctype[0]){
			// ���ΪC���ͣ�����Sc�������Ϊ�ֲ�ڵ�
		  case 'C':
			  Jc[d1]=(Jc[n]/2+Len[d1]*ModelParam::omega[d1].Sm)*exp(-Len[d1]/AdapParam::LRef[0]);
			  ModelParam::omega[d1].Jcmid=Jc[n]/2*exp(-0.5*Len[d1]/AdapParam::LRef[0])+0.5*Len[d1]*ModelParam::omega[d1].Sm*exp(-0.5*Len[d1]/AdapParam::LRef[0]);
			  Jc[d2]=(Jc[n]/2+Len[d2]*ModelParam::omega[d2].Sm)*exp(-Len[d2]/AdapParam::LRef[0]);
			  ModelParam::omega[d2].Jcmid=Jc[n]/2*exp(-0.5*Len[d2]/AdapParam::LRef[0])+0.5*Len[d2]*ModelParam::omega[d2].Sm*exp(-0.5*Len[d2]/AdapParam::LRef[0]);
			  break;
		  case 'J':
			  Jc[d1]=(Jc[n]+Len[d1]*ModelParam::omega[d1].Sm)*exp(-Len[d1]/AdapParam::LRef[0]);
			  ModelParam::omega[d1].Jcmid=Jc[n]*exp(-0.5*Len[d1]/AdapParam::LRef[0])+0.5*Len[d1]*ModelParam::omega[d1].Sm*exp(-0.5*Len[d1]/AdapParam::LRef[0]);
			  break;
		  default:
			  break;
		}
	}

	// Calculate Sc from Jc
	for(n=0;n<ModelParam::Ndoms;++n){
		ModelParam::omega[n].Sc=ModelParam::omega[n].Jcmid/(ModelParam::omega[n].Jcmid+AdapParam::J0[0]);
		Debug_Sc[n]=ModelParam::omega[n].Sc;
	}
}

double Adap_SS_Solver::Hill(double SO2){
	double P50=38;  // P50=38 mmHg
	int N=3;
	return P50*pow(SO2/(1-SO2),1.0/N);
}

double Adap_SS_Solver::mean_abs(vector<double> vec){
	double total=0.0;
	int n;
	for(n=0;n<vec.size();++n){
		total+=abs(vec[n]);
	}
	return total/vec.size();
}

void Adap_SS_Solver::resetInitState()
{
	int n,i,j,k;
	// Reset the diameter, wall thickness and viscosity to the input values
	for(n=0;n<ModelParam::Ndoms;n++){
		for(i=0;i<ModelParam::oneD_q-1;i++){
			ModelParam::omega[n].A[0].h[i]=0.25*ModelParam::PI*Org_Diam[n]*Org_Diam[n];
			ModelParam::omega[n].visc[i]=Org_Visc[n];
		}
		ModelParam::omega[n].WallTh=Org_WallTh[n];

		// Reset the connection matrix
		// ����bifur�ĳ�ʼ״̬
		for(j=0;j<2;j++){
			for(k=0;k<2;k++){
				ModelParam::omega[n].bifur[j][k]=org_bifur[n][j][k];
			}
		}
		// ����bctype�ĳ�ʼ״̬
		for(j=0;j<6;j++){
			ModelParam::omega[n].bctype[j]=org_bctype[n][j];    // 0,1,2 ���, 3,4,5 ����
		}
		// ����nodes�ĳ�ʼ״̬
		for(j=0;j<2;j++){
			for(k=0;k<4;k++){
				ModelParam::omega[n].nodes[j][k]=org_nodes[n][j][k];
			}
		} 

		lastVisc[n]=0.0;
	}
}