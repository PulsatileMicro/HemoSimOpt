#include "pso.h"
#include "AdapParam.h"
#include "AdapSSWall.h"
#include "ModelParam.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <iomanip>
#include <windows.h>
#include <iostream>
using namespace std;

swarm *s=(swarm *)malloc(sizeof(swarm));
particle *p=(particle *)malloc(sizeof(particle));

void PSO_RandInitofSwarm(void){
	int i,j;  

	s->W=1.4;  
	s->C1=2.05;  
	s->C2=2.05;
	s->Alpha=0.5;
	s->consist_convergence_time = 0;
	s->total_convergence_time = 0;
	for(j=0;j<PSO_Dim;j++)
	{
		if(ModelParam::solverType==ModelParam::Adap_SS_Wall){
			switch(j){
	  case 0: // PO2Ref
		  s->Xup[j]    = AdapParam::PO2Ref[1];
		  s->Xdown[j]  = AdapParam::PO2Ref[2];
		  break;
	  case 1: // kc
		  s->Xup[j]    = AdapParam::kc[1];
		  s->Xdown[j]  = AdapParam::kc[2];
		  break;
	  case 2: // kmd
		  s->Xup[j]    = AdapParam::kmd[1];
		  s->Xdown[j]  = AdapParam::kmd[2];
		  break;
	  case 3: // kmg
		  s->Xup[j]    = AdapParam::kmg[1];
		  s->Xdown[j]  = AdapParam::kmg[2];
		  break;
	  case 4: // ksd
		  s->Xup[j]    = AdapParam::ksd[1];
		  s->Xdown[j]  = AdapParam::ksd[2];
		  break;
	  case 5: // ksg
		  s->Xup[j]    = AdapParam::ksg[1];
		  s->Xdown[j]  = AdapParam::ksg[2];
		  break;
	  case 6: // kwtau
		  s->Xup[j]    = AdapParam::kwtau[1];
		  s->Xdown[j]  = AdapParam::kwtau[2];
		  break;
	  case 7: // kwsigma
		  s->Xup[j]    = AdapParam::kwsigma[1];
		  s->Xdown[j]  = AdapParam::kwsigma[2];
		  break;
	  case 8:
		  s->Xup[j]    = AdapParam::tauRef[1];
		  s->Xdown[j]  = AdapParam::tauRef[2];
		  break;
	  case 9:
		  s->Xup[j]    = AdapParam::sigmaRef[1];
		  s->Xdown[j]  = AdapParam::sigmaRef[2];
		  break;
	  case 10:
		  s->Xup[j]    = AdapParam::wRef[1];
		  s->Xdown[j]  = AdapParam::wRef[2];
		  break;
	  case 11:
		  s->Xup[j]    = AdapParam::J0[1];
		  s->Xdown[j]  = AdapParam::J0[2];
		  break;
	  case 12:
		  s->Xup[j]    = AdapParam::LRef[1];
		  s->Xdown[j]  = AdapParam::LRef[2];
		  break;
	  default:
		  break;
			}
		}
		else if(ModelParam::solverType==ModelParam::Adap_SS_NoWall){
			switch(j){
		case 0: // kc
			s->Xup[j]    = AdapParam::kc2[1];
			s->Xdown[j]  = AdapParam::kc2[2];
			break;
		case 1: // kp
			s->Xup[j]    = AdapParam::kp2[1];
			s->Xdown[j]  = AdapParam::kp2[2];
			break;
		case 2: // km
			s->Xup[j]    = AdapParam::km2[1];
			s->Xdown[j]  = AdapParam::km2[2];
			break;
		case 3: // ks
			s->Xup[j]    = AdapParam::ks2[1];
			s->Xdown[j]  = AdapParam::ks2[2];
			break;
		case 4: // J0
			s->Xup[j]    = AdapParam::J02[1];
			s->Xdown[j]  = AdapParam::J02[2];
			break;
		case 5: // LRef
			s->Xup[j]    = AdapParam::LRef2[1];
			s->Xdown[j]  = AdapParam::LRef2[2];
			break;
		case 6: // tauRef
			s->Xup[j]    = AdapParam::tauRef2[1];
			s->Xdown[j]  = AdapParam::tauRef2[2];
			break;
		case 7: // QRef
			s->Xup[j]    = AdapParam::QRef2[1];
			s->Xdown[j]  = AdapParam::QRef2[2];
			break;
		default:
			break;
			}
		}
		s->Vmax[j] = (s->Xup[j]+s->Xdown[j])/10;
	}  

	for(i=0; i<PSO_PNum; i++){  
		//printf(" The %dth of X is: ",i);  
		for(j=0; j<PSO_Dim; j++){
			s->Particle[i].X[j] = i*(s->Xup[j]-s->Xdown[j])/PSO_PNum+s->Xdown[j];
			s->Particle[i].V[j] = 0.3*(s->Xup[j]-s->Xdown[j]);
			cout << setprecision(6) << setw(15) << s->Particle[i].X[j];
			if (AdapParam::optMethod == AdapParam::SECPSO)
			{
				s->Particle[i].LastX[j] = i*(s->Xup[j]-s->Xdown[j])/PSO_PNum+s->Xdown[j];
				s->Particle[i].LastV[j] = 0.3*(s->Xup[j]-s->Xdown[j]);
				cout << setprecision(6) << setw(15) << s->Particle[i].LastV[j];
			}
		}
	}  
}

void PSO_PriesInitofSwarm(void){
	int i,j;  

	s->W=1.4;  
	s->C1=2.0;  
	s->C2=2.0;
	s->consist_convergence_time = 0;
	s->total_convergence_time = 0;
	for(j=0;j<PSO_Dim;j++)
	{
		if(ModelParam::solverType==ModelParam::Adap_SS_Wall){
			switch(j){
		case 0: // PO2Ref
			s->Xup[j]    = AdapParam::PO2Ref[0];
			s->Xdown[j]  = AdapParam::PO2Ref[0];
			break;
		case 1: // kc
			s->Xup[j]    = AdapParam::kc[0];
			s->Xdown[j]  = AdapParam::kc[0];
			break;
		case 2: // kmd
			s->Xup[j]    = AdapParam::kmd[0];
			s->Xdown[j]  = AdapParam::kmd[0];
			break;
		case 3: // kmg
			s->Xup[j]    = AdapParam::kmg[0];
			s->Xdown[j]  = AdapParam::kmg[0];
			break;
		case 4: // ksd
			s->Xup[j]    = AdapParam::ksd[0];
			s->Xdown[j]  = AdapParam::ksd[0];
			break;
		case 5: // ksg
			s->Xup[j]    = AdapParam::ksg[0];
			s->Xdown[j]  = AdapParam::ksg[0];
			break;
		case 6: // kwtau
			s->Xup[j]    = AdapParam::kwtau[0];
			s->Xdown[j]  = AdapParam::kwtau[0];
			break;
		case 7: // kwsigma
			s->Xup[j]    = AdapParam::kwsigma[0];
			s->Xdown[j]  = AdapParam::kwsigma[0];
			break;
		case 8:
			s->Xup[j]    = AdapParam::tauRef[0];
			s->Xdown[j]  = AdapParam::tauRef[0];
			break;
		case 9:
			s->Xup[j]    = AdapParam::sigmaRef[0];
			s->Xdown[j]  = AdapParam::sigmaRef[0];
			break;
		case 10:
			s->Xup[j]    = AdapParam::wRef[0];
			s->Xdown[j]  = AdapParam::wRef[0];
			break;
		case 11:
			s->Xup[j]    = AdapParam::J0[0];
			s->Xdown[j]  = AdapParam::J0[0];
			break;
		case 12:
			s->Xup[j]    = AdapParam::LRef[0];
			s->Xdown[j]  = AdapParam::LRef[0];
			break;
		default:
			break;
			}
		}
		else if(ModelParam::solverType==ModelParam::Adap_SS_NoWall){
			switch(j){
		case 0: // kc
			s->Xup[j]    = AdapParam::kc2[0];
			s->Xdown[j]  = AdapParam::kc2[0];
			break;
		case 1: // kp
			s->Xup[j]    = AdapParam::kp2[0];
			s->Xdown[j]  = AdapParam::kp2[0];
			break;
		case 2: // km
			s->Xup[j]    = AdapParam::km2[0];
			s->Xdown[j]  = AdapParam::km2[0];
			break;
		case 3: // ks
			s->Xup[j]    = AdapParam::ks2[0];
			s->Xdown[j]  = AdapParam::ks2[0];
			break;
		case 4: // J0
			s->Xup[j]    = AdapParam::J02[0];
			s->Xdown[j]  = AdapParam::J02[0];
			break;
		case 5: // LRef
			s->Xup[j]    = AdapParam::LRef2[0];
			s->Xdown[j]  = AdapParam::LRef2[0];
			break;
		case 6: // tauRef
			s->Xup[j]    = AdapParam::tauRef2[0];
			s->Xdown[j]  = AdapParam::tauRef2[0];
			break;
		case 7: // QRef
			s->Xup[j]    = AdapParam::QRef2[0];
			s->Xdown[j]  = AdapParam::QRef2[0];
			break;
		default:
			break;
			}
		}
		s->Vmax[j] = (s->Xup[j]+s->Xdown[j])/10;
	}  

	srand((unsigned)time(NULL));
	for(i=0; i<PSO_PNum; i++){  
		//printf(" The %dth of X is: ",i);  
		for(j=0; j<PSO_Dim; j++){
			s->Particle[i].X[j] = rand()/(double)RAND_MAX*(s->Xup[j]-s->Xdown[j])+s->Xdown[j];
			s->Particle[i].V[j] = rand()/(double)RAND_MAX*s->Vmax[j]*2-s->Vmax[j];
			//printf(" %.2f \n ",s->Particle[i].X[j]);  
		}  
	}
}

void PSO_ComputFitofSwarm(int nIter){
	int i,j;
	LARGE_INTEGER tick;
	LARGE_INTEGER timestamp;
	double ratio;
	char str[50];
	int loop_cnt=0;

	for(i=0; i<PSO_PNum; i++)
	{ 
		/*QueryPerformanceFrequency(&tick);
		QueryPerformanceCounter(&timestamp);
		unsigned int us=(unsigned int)(timestamp.QuadPart % tick.QuadPart)*1E6/tick.QuadPart;
		srand(us);*/
		
		AdapParam::ErrorV=PSO_ComputAFitness(s->Particle[i].X);
		
#ifdef NOADAPERR
		// 如果更新后的粒子计算不出目标值，则淘汰该粒子，引入新粒子，直至能够计算出目标值
		while(AdapParam::errFlag != AdapParam::NO_ADAP_ERR){
		// 随机化粒子及速度
			for(j=0; j<PSO_Dim; j++){
				s->Particle[i].X[j] = rand()/(double)RAND_MAX*(s->Xup[j]-s->Xdown[j])+s->Xdown[j];
				s->Particle[i].V[j] = rand()/(double)RAND_MAX*s->Vmax[j]*2-s->Vmax[j];
			}
			AdapParam::ErrorV = PSO_ComputAFitness(s->Particle[i].X);
		}
		if (AdapParam::errFlag != AdapParam::NO_ADAP_ERR && AdapParam::ErrorV < 0)
		{
			AdapParam::ErrorV = -AdapParam::ErrorV;
		}
#endif

		//printf("Particle, The Fitness of %dth Particle: ",i+1);
		// s->Particle[i].Fitness = AdapParam::ErrorD;
		s->Particle[i].Fitness = AdapParam::ErrorV;
		if (nIter == 1 && i == 0)
		{
			s->GBestFitness = s->Particle[i].Fitness;
		}
		//cout << setprecision(6) << setw(15) << s->Particle[i].Fitness << endl;

		if (AdapParam::errFlag == AdapParam::NO_ADAP_ERR)
		{
			// 将Good Particle及其对应的Fitness value写入文件记录
			AdapParam::adapLogFile << setw(15) << i << setw(15) << AdapParam::errFlag
				<< setprecision(6) << setw(15) << AdapParam::ErrorV
				<< setprecision(6) << setw(15) << AdapParam::ErrorD 
				<< setprecision(6) << setw(15) << AdapParam::ErrorQ;
			for(j=0;j<PSO_Dim;j++)
				AdapParam::adapLogFile << setprecision(6) << setw(15) << s->Particle[i].X[j];
			AdapParam::adapLogFile << endl;

			// 对第nIter次迭代的第i个粒子，记录仿真结果
			sprintf(str,"AdapResult_Iter%02d_Par%02d.dat",nIter,i+1);
			AdapParam::adapHemoFile.open(str);
			for(j=0;j<ModelParam::Ndoms;j++){
				AdapParam::adapHemoFile << setw(12) << j << setw(12) << sqrt(4*ModelParam::omega[j].A[0].h[0]/ModelParam::PI)*1e6 
					<< setw(12) << ModelParam::omega[j].WallTh*1e6 << setw(12) << Adap_SS_Solver::Debug_MeanP[j] 
				<< setw(12) << Adap_SS_Solver::Debug_MeanQ[j] << endl;
			}
			AdapParam::adapHemoFile.close();
		}
	}
}  

void PSO_UpdateofVandX(){
	int i,j;  
	srand((unsigned)time(NULL));  
	for(i=0; i<PSO_PNum; i++){  
		//printf(" The %dth of X is: ",i);  
		for(j=0; j<PSO_Dim; j++)
			s->Particle[i].LastV[j] = s->Particle[i].V[j];    // 先保存上一次的值
		s->Particle[i].V[j] = s->W*s->Particle[i].V[j]+
			(double)rand()/(RAND_MAX + 1)*s->C1*(s->Particle[i].PBest[j] - s->Particle[i].X[j])+  
			(double)rand()/(RAND_MAX + 1)*s->C2*(s->GBest[j] - s->Particle[i].X[j]);  
		for(j=0; j<PSO_Dim; j++){  
			if(s->Particle[i].V[j]>s->Vmax[j])
				s->Particle[i].V[j] = s->Vmax[j];  
			if(s->Particle[i].V[j]<-s->Vmax[j])   
				s->Particle[i].V[j] = -s->Vmax[j];
		}  

		for(j=0; j<PSO_Dim; j++){  
			s->Particle[i].LastX[j] = s->Particle[i].X[j];
			s->Particle[i].X[j] += s->Particle[i].V[j];  
			if(s->Particle[i].X[j]>s->Xup[j])  
				s->Particle[i].X[j]=s->Xup[j];  
			if(s->Particle[i].X[j]<s->Xdown[j])  
				s->Particle[i].X[j]=s->Xdown[j];  
			cout << setw(15) << setprecision(6) << s->Particle[i].X[j];
		}   
	}  
}


void PSO_UpdateofVandX_CompressMutation(){
	int i,j; 
	srand((unsigned)time(NULL));  
	for(i=0; i<PSO_PNum; i++){  
		//printf(" The %dth of X is: ",i);  
		for(j=0; j<PSO_Dim; j++)
			s->Particle[i].LastV[j] = s->Particle[i].V[j];    // 先保存上一次的值
		double phi = s->C1 + s->C2;
		double ksi = 2/abs(2-phi-sqrt(phi*phi - 4*phi));
		s->Particle[i].V[j] = ksi*(s->Particle[i].V[j]+
			(double)rand()/(RAND_MAX + 1)*s->C1*(s->Particle[i].PBest[j] - s->Particle[i].X[j])+  
			(double)rand()/(RAND_MAX + 1)*s->C2*(s->GBest[j] - s->Particle[i].X[j]));  
		for(j=0; j<PSO_Dim; j++){  
			double vb = (s->Xup-s->Xdown)*0.3;
			if(s->Particle[i].V[j]>vb)
				s->Particle[i].V[j] = vb - (double)rand()/(RAND_MAX + 1)*2*vb;  
			if(s->Particle[i].V[j]<-vb)   
				s->Particle[i].V[j] = -vb + (double)rand()/(RAND_MAX + 1)*2*vb;
		}  

		for(j=0; j<PSO_Dim; j++){  
			s->Particle[i].LastX[j] = s->Particle[i].X[j];
			s->Particle[i].X[j] += s->Particle[i].V[j];  
			if(s->Particle[i].X[j]>s->Xup[j])  
				s->Particle[i].X[j]=s->Xup[j] - (double)rand()/(RAND_MAX + 1)*0.8*(s->Xup[j] - s->Xdown[j]);  
			if(s->Particle[i].X[j]<s->Xdown[j])  
				s->Particle[i].X[j]=s->Xdown[j] + (double)rand()/(RAND_MAX + 1)*0.8*(s->Xup[j] - s->Xdown[j]);  
			cout << setw(15) << setprecision(6) << s->Particle[i].X[j];
		}   
	}  
}

void PSO_UpdateofVandX_QuantumBehavior(int nIter){
	int i,j;  
	// double Beita =  (1 - s->Alpha)*(PSO_N - nIter + 1)/PSO_N + s->Alpha;
	double Beita =  0.7 + (double)rand()/(RAND_MAX+1);
	srand((unsigned)time(NULL));
	for(i=0; i<PSO_PNum; i++){  
		//printf(" The %dth of X is: ",i);  
		for(j=0; j<PSO_Dim; j++){
			s->Particle[i].LastX[j] = s->Particle[i].X[j];    // 先保存上一次的值
			double pn = (double)rand()/(RAND_MAX + 1);
			s->Particle[i].PBest[j] = pn*s->Particle[i].PBest[j] + (1-pn)*s->GBest[j];
			pn = (double)rand()/(RAND_MAX + 1);

			/*if rand < 0.5
			x(i,d) = p(i,d)+beita*abs(mbest(d) - x(i,d))*log(1/fy);
			else
			x(i,d) = p(i,d)-beita*abs(mbest(d) - x(i,d))*log(1/fy);
			end*/

			if (pn < 0.5)
			{
				s->Particle[i].X[j] = s->Particle[i].PBest[j] + 
					Beita*abs(s->MeanX[j] - s->Particle[i].X[j])*log(1/pn);
			} 
			else
			{
				s->Particle[i].X[j] = s->Particle[i].PBest[j] -
					Beita*abs(s->MeanX[j] - s->Particle[i].X[j])*log(1/pn);
			}

			if(s->Particle[i].X[j]>s->Xup[j])  
				s->Particle[i].X[j]=s->Xup[j] - (double)rand()/(RAND_MAX + 1)*0.8*(s->Xup[j] - s->Xdown[j]);  
			if(s->Particle[i].X[j]<s->Xdown[j])  
				s->Particle[i].X[j]=s->Xdown[j] + (double)rand()/(RAND_MAX + 1)*0.8*(s->Xup[j] - s->Xdown[j]);  
			//cout << setw(15) << setprecision(6) << s->Particle[i].X[j];
		}   
	}
}

void PSO_UpdateofVandX_SecondBehavior(){
	int i,j;  
	srand((unsigned)time(NULL));  
	for(i=0; i<PSO_PNum; i++){  
		//printf(" The %dth of X is: ",i);  
		for(j=0; j<PSO_Dim; j++)
			s->Particle[i].LastV[j] = s->Particle[i].V[j];    // 先保存上一次的值
		s->Particle[i].V[j] = s->W*s->Particle[i].V[j]+
			(double)rand()/(RAND_MAX + 1)*s->C1*(s->Particle[i].PBest[j] - 2*s->Particle[i].X[j] + s->Particle[i].LastX[j])+  
			(double)rand()/(RAND_MAX + 1)*s->C2*(s->GBest[j] - 2*s->Particle[i].X[j] + s->Particle[i].LastX[j]);  
		for(j=0; j<PSO_Dim; j++){  
			double vb = (s->Xup-s->Xdown)*0.3;
			if(s->Particle[i].V[j]>vb)
				s->Particle[i].V[j] = vb - (double)rand()/(RAND_MAX + 1)*2*vb;  
			if(s->Particle[i].V[j]<-vb)   
				s->Particle[i].V[j] = -vb + (double)rand()/(RAND_MAX + 1)*2*vb;
		}  

		for(j=0; j<PSO_Dim; j++){  
			s->Particle[i].LastX[j] = s->Particle[i].X[j];
			s->Particle[i].X[j] += s->Particle[i].V[j];  
			if(s->Particle[i].X[j]>s->Xup[j])  
				s->Particle[i].X[j]=s->Xup[j] - (double)rand()/(RAND_MAX + 1)*0.8*(s->Xup[j] - s->Xdown[j]);  
			if(s->Particle[i].X[j]<s->Xdown[j])  
				s->Particle[i].X[j]=s->Xdown[j] + (double)rand()/(RAND_MAX + 1)*0.8*(s->Xup[j] - s->Xdown[j]);  
			cout << setw(15) << setprecision(6) << s->Particle[i].X[j];
		}   
	}
}

void PSO_ModifyVandX(int i){
	int j;
	// 重新更新一次Velocity
	for(j=0; j<PSO_Dim; j++)
		s->Particle[i].V[j] = s->W*s->Particle[i].LastV[j]+
		rand()/(double)RAND_MAX*s->C1*(s->Particle[i].PBest[j] - s->Particle[i].LastX[j])+  
		rand()/(double)RAND_MAX*s->C2*(s->GBest[j] - s->Particle[i].LastX[j]);  
	for(j=0; j<PSO_Dim; j++){  
		if(s->Particle[i].V[j]>s->Vmax[j])
			s->Particle[i].V[j] = s->Vmax[j];  
		if(s->Particle[i].V[j]<-s->Vmax[j])   
			s->Particle[i].V[j] = -s->Vmax[j];  
	}  

	// 重新用Velocity计算一次X
	for(j=0; j<PSO_Dim; j++){
		s->Particle[i].X[j] = s->Particle[i].LastX[j]+s->Particle[i].V[j];  
		if(s->Particle[i].X[j]>s->Xup[j])  
			s->Particle[i].X[j]=s->Xup[j];  
		if(s->Particle[i].X[j]<s->Xdown[j])  
			s->Particle[i].X[j]=s->Xdown[j];
	}  
}

void PSO_UpdatePandGbest( int nIter){
	int i,j;
	if ( nIter == 1 ) {
		//P=X;  
		for(i=0; i<PSO_PNum; i++){  
			for(j=0;j<PSO_Dim;j++){  
				s->Particle[i].PBest[j]=s->Particle[i].X[j];  
			}  
			s->Particle[i].PBestFitness = s->Particle[i].Fitness;
		}
	}
	else
	{
		//update of P if the X is bigger than current P  
		for (i = 0; i < PSO_PNum; i++){  
			//printf(" The %dth of P is: ",i);  
			if (!(s->Particle[i].Fitness > s->Particle[i].PBestFitness)){  
				for(j=0;j<PSO_Dim;j++){  
					s->Particle[i].PBest[j] = s->Particle[i].X[j];  
				}  
				s->Particle[i].PBestFitness = s->Particle[i].Fitness;
			}  
			//printf(" %.2f %.2f \n",s->Particle[i].P[0],s->Particle[i].P[1]);  
		}  
	}
	
	//update of GBest 
	double gbest_temp = s->GBestFitness;
	for(i=0; i<PSO_PNum; i++)
		// if(PSO_ComputAFitness(s->Particle[i].P) <= s->Particle[s->GBestIndex].Fitness)
		if(!(s->Particle[i].PBestFitness > s->GBestFitness)){
			s->GBestFitness = s->Particle[i].PBestFitness;
			for(j=0;j<PSO_Dim;j++){
				s->GBest[j]=s->Particle[i].PBest[j];
			}
		}
	if (s->GBestFitness < 1)
	{
		s->total_convergence_time++;
		if(fabs(gbest_temp - s->GBestFitness) < 1e-7)
		{
			s->consist_convergence_time++;
		}
		else
		{
			s->consist_convergence_time = 1;
		}
	}
	//cout<< "Fitness of GBest: " << setprecision(6) << setw(12) << s->GBestFitness <<endl;
	if ( nIter == 1)
	{
		AdapParam::adapGBestFile.open("Global_Best_Record.dat",ios::app);
		AdapParam::adapGBestFile << setw(12) << "GBest=" << s->GBestFitness << endl;
	} 
	else
	{
		AdapParam::adapGBestFile << setw(12) << "GBest=" << s->GBestFitness << endl;
	}
	//calc MeanX
	if (AdapParam::optMethod == AdapParam::QUAPSO)
	{
		for (j=0; j<PSO_Dim; j++)
		{
			double temp = 0;
			for (i=0; i<PSO_PNum; i++)
			{
				temp += s->Particle[i].PBest[j];
			}
			s->MeanX[j] = temp/PSO_PNum;
		}

		if(nIter != 1)
		{
			constexpr std::size_t pso_size_t = sizeof s->Particle / sizeof *s->Particle;

			std::qsort(s->Particle, pso_size_t, sizeof *s->Particle, [](const void* x, const void* y)
			{
				particle arg1 = *static_cast<const struct PARTICLE*>(x);
				particle arg2 = *static_cast<const struct PARTICLE*>(y);

				if (arg1.Fitness - arg2.Fitness < 1e-7)
					return -1;
				else if (arg1.Fitness - arg2.Fitness > 1e-7)
					return 1;
				else
					return 0;
			});

			for (i = 0; i < PSO_PNum/2; i++)
			{
				for (j = 0; j < PSO_Dim; j++)
				{
					s->Particle[PSO_PNum - 1 - i].X[j] = s->Particle[i].X[j];
				}
			}
		}
	}
	else if (AdapParam::optMethod == AdapParam::SELPSO)
	{
		PSO_QuickSort(s->Particle, 0, PSO_PNum);
		for (i=0; i<PSO_PNum/2; i++)
		{
			for (j=0; j<PSO_Dim; j++)
			{
				s->Particle[PSO_PNum - 1 - i].X[j] = s->Particle[i].X[j];
			}
		}
	}
}

static double PSO_ComputAFitness(double X[])  
{  
	AdapParam::setPara2PSO(X);
	int errType=Adap_SS_Solver::AdapObjFunc();
	// return AdapParam::ErrorD;
	return AdapParam::ErrorV;
	// return AdapParam::ErrorQ;
}  


void PSO_QuickSort(particle A[], int p,int q)
{
	int r;
	if(p<q)
	{
		r=PSO_QuickSortPartion(A, p,q);
		PSO_QuickSort(A,p,r);  
		PSO_QuickSort(A,r+1,q);
	}
}


int PSO_QuickSortPartion(particle A[], int p,int q)
{
	double x= A[p].Fitness;
	int i=p;
	int j;

	for(j=p+1; j<q; j++)
	{
		if(A[j].Fitness <= x)
		{
			i=i+1;
			swap(A[i],A[j]);
		}

	}

	swap(A[i],A[p]);
	return i;
}
