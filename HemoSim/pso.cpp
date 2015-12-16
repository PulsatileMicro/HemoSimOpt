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
  s->C1=2.0;  
  s->C2=2.0;
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

void PSO_PriesInitofSwarm(void){
  int i,j;  

  s->W=1.4;  
  s->C1=2.0;  
  s->C2=2.0;
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

void PSO_ComputFitofSwarm(int isFirst, int nIter){
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
    if(isFirst){
      // 尝试随机化的粒子，直到能够计算出适应度函数
      do{
        // 随机化粒子及速度
        for(j=0; j<PSO_Dim; j++){
          s->Particle[i].X[j] = rand()/(double)RAND_MAX*(s->Xup[j]-s->Xdown[j])+s->Xdown[j];
          s->Particle[i].V[j] = rand()/(double)RAND_MAX*s->Vmax[j]*2-s->Vmax[j];
        }
        // printf("Bad particle\n");
        // Debug
        /*for(j=0;j<PSO_Dim;j++)
        AdapParam::adapLogFile << setprecision(6) << setw(15) << s->Particle[i].X[j];
        AdapParam::adapLogFile << endl;*/

        // 尝试计算适应度函数
		AdapParam::ErrorV=PSO_ComputAFitness(s->Particle[i].X);
      }while(AdapParam::errFlag!=AdapParam::NO_ADAP_ERR);
	  if (i == 0)
	  {
		  s->GBestFitness = s->Particle[i].Fitness;
	  }
    }
    else{
      // 不是第一次，则尝试计算适应度函数
      AdapParam::ErrorV=PSO_ComputAFitness(s->Particle[i].X);
      // 如果更新后的粒子计算不出目标值，则淘汰该粒子，引入新粒子，直至能够计算出目标值
      while(AdapParam::errFlag!=AdapParam::NO_ADAP_ERR){
        // 随机化粒子及速度
        for(j=0; j<PSO_Dim; j++){
          s->Particle[i].X[j] = rand()/(double)RAND_MAX*(s->Xup[j]-s->Xdown[j])+s->Xdown[j];
          s->Particle[i].V[j] = rand()/(double)RAND_MAX*s->Vmax[j]*2-s->Vmax[j];
        }
        AdapParam::ErrorV=PSO_ComputAFitness(s->Particle[i].X);
      }
    }

    printf("Good Particle, The Fitness of %dth Particle: ",i);
    // s->Particle[i].Fitness = AdapParam::ErrorD;
    s->Particle[i].Fitness = AdapParam::ErrorV;
    printf(" %.2f\n",s->Particle[i].Fitness);

    // 将Good Particle及其对应的Fitness value写入文件记录
    AdapParam::adapLogFile << setw(15) << i << setprecision(6) << setw(15) << AdapParam::ErrorV
      << setprecision(6) << setw(15) << AdapParam::ErrorD 
      << setprecision(6) << setw(15) << AdapParam::ErrorQ;
    for(j=0;j<PSO_Dim;j++)
      AdapParam::adapLogFile << setprecision(6) << setw(15) << s->Particle[i].X[j];
    AdapParam::adapLogFile << endl;

    // 对第nIter次迭代的第i个粒子，记录仿真结果
    sprintf(str,"AdapResult_Iter%02d_Par%02d.dat",nIter,i);
    AdapParam::adapHemoFile.open(str);
    for(j=0;j<ModelParam::Ndoms;j++){
      AdapParam::adapHemoFile << setw(12) << j << setw(12) << sqrt(4*ModelParam::omega[j].A[0].h[0]/ModelParam::PI)*1e6 
        << setw(12) << ModelParam::omega[j].WallTh*1e6 << setw(12) << Adap_SS_Solver::Debug_MeanP[j] 
      << setw(12) << Adap_SS_Solver::Debug_MeanQ[j] << endl;
    }
     AdapParam::adapHemoFile.close();
  }
}

void PSO_FirstComputPandGbest(void)  
{  
  int i,j;  
  //P=X;  
  for(i=0; i<PSO_PNum; i++){  
    for(j=0;j<PSO_Dim;j++){  
      s->Particle[i].PBest[j]=s->Particle[i].X[j];  
    }  
	s->Particle[i].PBestFitness = s->Particle[i].Fitness;
  }  
  //Computation of GBest  
  for(i=0; i<PSO_PNum; i++)  
    if(s->Particle[i].Fitness <= s->GBestFitness){  
      s->GBestFitness = s->Particle[i].Fitness;
      for(j=0;j<PSO_Dim;j++){  
        s->GBest[j]=s->Particle[i].X[j];
      }
    }
  
  printf("Fitness of GBest:%d ,%.2f \n",s->GBestFitness);

  AdapParam::adapGBestFile.open("Global_Best_Record.dat",ios::app);
  AdapParam::adapGBestFile << setw(12) << "GBest=" << s->GBestFitness << endl;
}  

void PSO_UpdateofVandX(){
  int i,j;  
  // srand((unsigned)time(NULL));  
  for(i=0; i<PSO_PNum; i++){  
    //printf(" The %dth of X is: ",i);  
    for(j=0; j<PSO_Dim; j++)
      s->Particle[i].LastV[j] = s->Particle[i].V[j];    // 先保存上一次的值
    s->Particle[i].V[j] = s->W*s->Particle[i].V[j]+
      rand()/(double)RAND_MAX*s->C1*(s->Particle[i].PBest[j] - s->Particle[i].X[j])+  
      rand()/(double)RAND_MAX*s->C2*(s->GBest[j] - s->Particle[i].X[j]);  
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
    }  
    //printf(" %.2f %.2f \n",s->Particle[i].X[0],s->Particle[i].X[1]);  
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

void PSO_UpdatePandGbest(){
  int i,j;
  //update of P if the X is bigger than current P  
  for (i = 0; i < PSO_PNum; i++){  
    //printf(" The %dth of P is: ",i);  
    if (s->Particle[i].Fitness <= s->Particle[i].PBestFitness){  
      for(j=0;j<PSO_Dim;j++){  
        s->Particle[i].PBest[j] = s->Particle[i].X[j];  
      }  
	  s->Particle[i].PBestFitness = s->Particle[i].Fitness;
    }  
    //printf(" %.2f %.2f \n",s->Particle[i].P[0],s->Particle[i].P[1]);  
  }  
  for (i = 0; i < PSO_PNum; i++){  
    //printf("The %dth of P's Fitness is : %.2f  \n",i,ComputAFitness(s->Particle[i].P));  
  }  
  //update of GBest  
  for(i=0; i<PSO_PNum; i++)
    // if(PSO_ComputAFitness(s->Particle[i].P) <= s->Particle[s->GBestIndex].Fitness)
    if(s->Particle[i].PBestFitness <= s->GBestFitness){
      s->GBestFitness = s->Particle[i].PBestFitness;
      for(j=0;j<PSO_Dim;j++){
        s->GBest[j]=s->Particle[i].PBest[j];
      }
    }
  printf("Fitness of GBest:%.2f \n", s->GBestFitness); 
  AdapParam::adapGBestFile << setw(12) << "GBest=" << s->GBestFitness << endl;
}

static double PSO_ComputAFitness(double X[])  
{  
  AdapParam::setPara2PSO(X);
  int errType=Adap_SS_Solver::AdapObjFunc();
  // return AdapParam::ErrorD;
  return AdapParam::ErrorV;
  // return AdapParam::ErrorQ;
}  