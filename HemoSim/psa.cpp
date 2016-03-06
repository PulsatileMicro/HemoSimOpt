#include "psa.h"
#include "AdapParam.h"
#include "AdapSSWall.h"
#include "ModelParam.h"
#include <math.h>
#include <iomanip>
#include <iostream>

PSA::PSA(void)
{
}

PSA::~PSA(void)
{
}

void PSA::set_step_len(int step_len)
{
	this->step_len = step_len;
}

int PSA::get_step_len()
{
	return this->step_len;
}


void PSA::init_particle_temp()
{
	for (int j=0; j<PSO_Dim; j++)
	{
		if(ModelParam::solverType==ModelParam::Adap_SS_NoWall){
			switch(j){
			case 0: // kc
				this->particle_temp[0].X[j] = AdapParam::kc2[0];
				this->particle_temp[1].X[j] = AdapParam::kc2[1];
				this->particle_temp[2].X[j] = AdapParam::kc2[2];
				//cout << "KC2:" << setw(15) << AdapParam::kc2[0] << setw(15) << AdapParam::kc2[1] << setw(15) << AdapParam::kc2[2] << endl;
				break;
			case 1: // kp
				this->particle_temp[0].X[j] = AdapParam::kp2[0];
				this->particle_temp[1].X[j] = AdapParam::kp2[1];
				this->particle_temp[2].X[j] = AdapParam::kp2[2];
				//cout << "kp2:" << setw(15) << AdapParam::kp2[0] << setw(15) << AdapParam::kp2[1] << setw(15) << AdapParam::kp2[2] << endl;
				break;
			case 2: // km
				this->particle_temp[0].X[j] = AdapParam::km2[0];
				this->particle_temp[1].X[j] = AdapParam::km2[1];
				this->particle_temp[2].X[j] = AdapParam::km2[2];
				//cout << "km2:" << setw(15) << AdapParam::km2[0] << setw(15) << AdapParam::km2[1] << setw(15) << AdapParam::km2[2] << endl;
				break;
			case 3: // ks
				this->particle_temp[0].X[j] = AdapParam::ks2[0];
				this->particle_temp[1].X[j] = AdapParam::ks2[1];
				this->particle_temp[2].X[j] = AdapParam::ks2[2];
				//cout << "KC2:" << setw(15) << AdapParam::kc2[0] << setw(15) << AdapParam::kc2[1] << setw(15) << AdapParam::kc2[2] << endl;
				break;
			case 4: // J0
				this->particle_temp[0].X[j] = AdapParam::J02[0];
				this->particle_temp[1].X[j] = AdapParam::J02[1];
				this->particle_temp[2].X[j] = AdapParam::J02[2];
				//cout << "J02:" << setw(15) << AdapParam::J02[0] << setw(15) << AdapParam::J02[1] << setw(15) << AdapParam::J02[2] << endl;
				break;
			case 5: // LRef
				this->particle_temp[0].X[j] = AdapParam::LRef2[0];
				this->particle_temp[1].X[j] = AdapParam::LRef2[1];
				this->particle_temp[2].X[j] = AdapParam::LRef2[2];
				//cout << "LRef2:" << setw(15) << AdapParam::LRef2[0] << setw(15) << AdapParam::LRef2[1] << setw(15) << AdapParam::LRef2[2] << endl;
				break;
			case 6: // tauRef
				this->particle_temp[0].X[j] = AdapParam::tauRef2[0];
				this->particle_temp[1].X[j] = AdapParam::tauRef2[1];
				this->particle_temp[2].X[j] = AdapParam::tauRef2[2];
				//cout << "tauRef2:" << setw(15) << AdapParam::tauRef2[0] << setw(15) << AdapParam::tauRef2[1] << setw(15) << AdapParam::tauRef2[2] << endl;
				break;
			case 7: // QRef
				this->particle_temp[0].X[j] = AdapParam::QRef2[0];
				this->particle_temp[1].X[j] = AdapParam::QRef2[1];
				this->particle_temp[2].X[j] = AdapParam::QRef2[2];
				//cout << "QRef2:" << setw(15) << AdapParam::QRef2[0] << setw(15) << AdapParam::QRef2[1] << setw(15) << AdapParam::QRef2[2] << endl;
				break;
			default:
				break;
			}
		}
	}
}

void PSA::reinflat_particles(int param_index)
{
	if (ModelParam::solverType == ModelParam::Adap_SS_NoWall && param_index > 8 || ModelParam::solverType == ModelParam::Adap_SS_Wall && param_index > 12)
	{
		std::cout << "PARAMETERS out of index" << endl;
		exit(EXIT_FAILURE);
	}
	
	this->particles.clear();
	this->particles.assign(this->step_len+1, particle_temp[0]);

	for (int i=0; i<=this->step_len; i++)
	{
		particle &p = this->particles.at(i);
		p.X[param_index] = particle_temp[2].X[param_index]+i*(particle_temp[1].X[param_index] - particle_temp[2].X[param_index])/this->step_len;
	}
}

void PSA::run_params_sensitivity_analysis(int param_index)
{
	for (int i=0; i<this->step_len; i++)
	{
		if(ModelParam::solverType==ModelParam::Adap_SS_Wall){
			// For wall adaptation
			AdapParam::PO2Ref[0]   = particles[i].X[0];
			AdapParam::kc[0]       = particles[i].X[1];
			AdapParam::kmd[0]      = particles[i].X[2];
			AdapParam::kmg[0]      = particles[i].X[3];
			AdapParam::ksd[0]      = particles[i].X[4];
			AdapParam::ksg[0]      = particles[i].X[5];
			AdapParam::kwtau[0]    = particles[i].X[6];
			AdapParam::kwsigma[0]  = particles[i].X[7];
			AdapParam::tauRef[0]   = particles[i].X[8];     
			AdapParam::sigmaRef[0] = particles[i].X[9];
			AdapParam::wRef[0]     = particles[i].X[10];
			AdapParam::J0[0]       = particles[i].X[11];
			AdapParam::LRef[0]     = particles[i].X[12];
		}
		else if(ModelParam::solverType==ModelParam::Adap_SS_NoWall){
			// For no wall adaptation
			AdapParam::kc2[0]      = particles[i].X[0];
			AdapParam::kp2[0]      = particles[i].X[1];
			AdapParam::km2[0]      = particles[i].X[2];
			AdapParam::ks2[0]      = particles[i].X[3];
			AdapParam::J02[0]      = particles[i].X[4];
			AdapParam::LRef2[0]    = particles[i].X[5];
			AdapParam::tauRef2[0]  = particles[i].X[6];
			AdapParam::QRef2[0]    = particles[i].X[7];
		}
		
		int errory_type = Adap_SS_Solver::AdapObjFunc();

		//if (AdapParam::errFlag == AdapParam::NO_ADAP_ERR)
		if (true)
		{
			// 将Good Particle及其对应的Fitness value写入文件记录
			AdapParam::adapLogFile << setw(15) << i << setw(15) << AdapParam::errFlag
				<< setprecision(6) << setw(15) << AdapParam::ErrorV
				<< setprecision(6) << setw(15) << AdapParam::ErrorD 
				<< setprecision(6) << setw(15) << AdapParam::ErrorQ;
			for(int j=0;j<PSO_Dim;j++)
				AdapParam::adapLogFile << setprecision(6) << setw(15) << particles[i].X[j];
			AdapParam::adapLogFile << endl;

			// 对第nIter次迭代的第i个粒子，记录仿真结果
			char str[50];
			sprintf(str,"PSAResult_ParamIndex%02d_Particle%02d.txt", param_index+1, i+1);
			AdapParam::adapHemoFile.open(str);
			for(int j=0;j<ModelParam::Ndoms;j++){
				AdapParam::adapHemoFile << setw(12) << j << setw(12) << sqrt(4*ModelParam::omega[j].A[0].h[0]/ModelParam::PI)*1e6 
					<< setw(12) << ModelParam::omega[j].WallTh*1e6 << setw(12) << Adap_SS_Solver::Debug_MeanP[j] 
				<< setw(12) << Adap_SS_Solver::Debug_MeanQ[j] << endl;
			}
			AdapParam::adapHemoFile.close();
		}

	}
}