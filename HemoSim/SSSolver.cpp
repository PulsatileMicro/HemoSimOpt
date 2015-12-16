#include "SSSolver.h"
#include "ModelParam.h"
#include <stdlib.h>
#include <mkl.h>
#include <math.h>
#include <string.h>

double **SSSolver::JMat=NULL;
double *SSSolver::NodeFlow=NULL;
double *SSSolver::SS_Press=NULL;
double *SSSolver::SS_DeltaP=NULL;
double *SSSolver::MeanP=NULL;
double *SSSolver::MeanQ=NULL;
double *SSSolver::q_in=NULL;
FILE **SSSolver::fp_his;

void SSSolver::initSolver(){
  /*
    Create the JMat matrix and initialize it with zeros
  */
  JMat=(double **)malloc(ModelParam::Nnodes*sizeof(double*));
  for(int i=0; i<ModelParam::Nnodes; ++i){
    JMat[i]=(double *)malloc(ModelParam::Nnodes*sizeof(double));
  }
  for(int i=0; i<ModelParam::Nnodes; ++i)
    for (int j=0; j<ModelParam::Nnodes; ++j)
      JMat[i][j]=0;
    
  /*
    Create the boundary vector and initialize it with zeros
  */
  NodeFlow=(double *)malloc(ModelParam::Nnodes*sizeof(double));
  for(int i=0; i<ModelParam::Nnodes; ++i)
    NodeFlow[i]=0;

  SS_Press  =(double *)malloc(ModelParam::Nnodes*sizeof(double));
  SS_DeltaP =(double *)malloc(ModelParam::Ndoms*sizeof(double));
  MeanP     =(double *)malloc(ModelParam::Ndoms*sizeof(double));
  MeanQ     =(double *)malloc(ModelParam::Ndoms*sizeof(double));
}

/*!
  \brief Evaluate the X in the equ. AX=B in the steady state model
*/
void SSSolver::solve(){
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
    d1    = ModelParam::omega[n].bifur[0][0];     // start node所连接的两条血管的id
    d2    = ModelParam::omega[n].bifur[0][1];
    d3    = ModelParam::omega[n].bifur[1][0];     // end node所连接的两条血管的id
    d4    = ModelParam::omega[n].bifur[1][1];

    radius = sqrt(ModelParam::omega[n].A[0].h[0]/ModelParam::PI);
    len = ModelParam::omega[n].A[0].x[1]-ModelParam::omega[n].A[0].x[0];
    ModelParam::omega[n].J=ModelParam::PI*radius*radius*radius*radius/(8*ModelParam::omega[n].visc[0]*len);

    switch(ModelParam::omega[n].bctype[0]){
      /*
        如果start node的类型为B，则该段血管只能是子血管
        d1,d2的类型存在两种情况：d1为母血管，d2为子血管；或者相反。
      */
      case 'B':
        if(ModelParam::omega[d1].nodes[0][0] == start){
          /*
            如果d1的起点与当前分析血管段的起点相同，说明d2是母血管，d1和n是子血管
              From[0]为母血管d2的起点
              To[0]为当前血管段的终点
              To[1]为另一子血管d1的终点
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
            // 压力边界为分叉的情况(CAM网络)，暂未实现

          }
          else{
            JMat[start][To[0]]=ModelParam::omega[n].J;
            JMat[start][To[1]]=ModelParam::omega[d1].J;
          }
        }
        else if(ModelParam::omega[d2].nodes[0][0] == start){
          /*
            如果d2的起点与当前分析血管段的起点相同，说明d1是母血管，d2和n是子血管
              From[0]为母血管d1的起点
              To[0]为当前血管段的终点
              To[1]为另一子血管d2的终点
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
            // 压力边界为分叉的情况(CAM网络)

          }
          else{
            JMat[start][To[0]]=ModelParam::omega[n].J;
            JMat[start][To[1]]=ModelParam::omega[d2].J;
          }
        }
        break;
        /*
          入口点为C类型时，只有一种情况，即当前血管n为母血管，d1,d2为子血管
        */   
      case 'C':
        /*
          From[0]为子血管d1的起点
          From[1]为子血管d2的起点
          To[0]为母血管n的终点
        */
        From[0]=ModelParam::omega[d1].nodes[0][0];
        From[1]=ModelParam::omega[d2].nodes[0][0];
        To[0]=end;
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
        /*
          末点为B类型时，只有一种情况：n为母血管，d3,d4为子血管
        */
      case 'B':
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
          // 压力边界为分叉的情况(CAM网络)
        }
        else{
          JMat[end][To[0]]=ModelParam::omega[d3].J;
          JMat[end][To[1]]=ModelParam::omega[d4].J;
        }
        break;
        /*
          末点为C类型时，与起点为B类型相似，需分别讨论
        */
      case 'C':
        if(ModelParam::omega[d3].nodes[1][0] == end){
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
      JMatArray[i*ModelParam::Nnodes+j]=JMat[j][i];   // JMatArray按照列顺序摆放元素
    }
  }
  
  int ok;   // dgesv状态
  int *pivot=(int *)malloc(ModelParam::Nnodes*sizeof(int));
  dgesv(&lapack_n,&one,JMatArray,&lapack_n,pivot,NodeFlow,&lapack_n,&ok);

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
    MeanQ[n]=60*1e12*ModelParam::omega[n].J*SS_DeltaP[n];
    SS_DeltaP[n]=SS_DeltaP[n]/133;
    MeanP[n]=MeanP[n]/133;
  }

  Write_history(ModelParam::omega,ModelParam::argv[ModelParam::argc-1]);
}

void SSSolver::destroySolver(){
  // TODO: Everything initialized in the initSolver() should be destroyed.
}

void SSSolver::Write_history(Domain *omega, char *name){
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
