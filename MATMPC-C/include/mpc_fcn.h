//
// Created by pico on 03/03/22.
//

#ifndef MPC_FCN_H
#define MPC_FCN_H

//#define PI 3.14159
#define INF 1.e+10

void initSettings(TypeSettings *settings);
void initInput(TypeInput *input);
void initMem(TypeMem *mem);
//void setSimInput(TypeSimInput *sim_input, TypeInput *input, double* state_sim);
int qp_generation(TypeInput *input, TypeMem *mem, TypeSettings *settings);
int hpipm_sparse(TypeMem *mem, TypeSettings *settings);
int line_search(TypeInput *input, TypeMem *mem, TypeSettings *settings);
int solution_info(TypeInput *input, TypeMem *mem, TypeSettings *settings, double *eq_res_out,double *ineq_res_out,double* KKT_out,double* OBJ_out);
void simulate_system(TypeSimInput *sim_input, TypeMem *mem, TypeSettings *settings, TypeSimOutput *sim_output);

#endif //MPC_FCN_H
