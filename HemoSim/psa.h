#pragma once
#include "pso.h"
#include <vector>
using namespace std;

class PSA
{
public:
	int step_len;
	particle particle_temp[3];
	vector<particle> particles;
	PSA(void);
	~PSA(void);
	void set_step_len(int step_len);
	int get_step_len();
	void init_particle_temp();
	void reinflat_particles(int param_index);
	void run_params_sensitivity_analysis(int param_index);
}; 
