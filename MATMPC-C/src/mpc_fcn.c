//
// Created by pico on 03/03/22.
//

#include <string.h>
#include "main_utils.h"
#include "mpc_fcn.h"

void initSettings(TypeSettings *settings){
    strncpy(settings->model,MODEL_NAME,100);
    settings->Ts_st = TS_ST;
    settings->nx = NX;
    settings->nu = NU;
    settings->nz = NZ;
    settings->ny = NY;
    settings->nyN = NYN;
    settings->np = NP;
    settings->nc = NC;
    settings->ncN = NCN;
    settings->nbx = NBX;
    settings->nbu = NBU;
    int arr_nbx[NBX] = NBX_IDX;
    setArrayIntValues(settings->nbx_idx,arr_nbx,NBX);
    int arr_nbu[NBU] = NBU_IDX;
    setArrayIntValues(settings->nbu_idx,arr_nbu,NBU);
    settings->N = NS;
    settings->N2 = NS2;
    settings->r = RS;
}

void initInput(TypeInput *input){
    /* USER DEFINED INITIALIZATION */
    // input.x0: to be read from a configuration file!
    double x0[NX] = {0, 3.14, 0, 0};
    setArrayDoubleValues(input->x0,x0,NX);
    // input.u0:  to be read from a configuration file!
    double u0[NU] = {0};
    setArrayDoubleValues(input->u0,u0,NU);
    // input.z0:  to be read from a configuration file!
    double z0[NZ] = {};
    setArrayDoubleValues(input->z0,z0,NZ);
    // input.para0: to be read from a configuration file!
    double para0[NP] = {};
    // cost function weights h and hN: to be read from a configuration file!
    double h[NY] = {1.0e+1, 1.0e+1, 1.0e-1, 1.0e-1, 1.0e-2};
    double hN[NYN] = {1.0e+1, 1.0e+1, 1.0e-1, 1.0e-1};
    // upper and lower bounds for states (=nbx)
    double lb_x[NBX] = { -2};
    double ub_x[NBX] = { +2};
    // upper and lower bounds for controls (=nbu)
    double lb_u[NBU] = { -20};
    double ub_u[NBU] = { +20};
    // upper and lower bounds for general constraints (=nc)
    // *b_g = [e_y + slack_lat]
    double lb_g[NC] = {};
    double ub_g[NC] = {};
    // *b_gN = [e_y]
    double lb_gN[NCN] = {};
    double ub_gN[NCN] = {};

    /* AUTOMATIC INITIALIZATION */
    // general constraints
    double lb[NC*(NS+1)];
    repeatDoubleArray(lb_g, lb, NC, NC*NS);
    repeatDoubleArray(lb_gN, &lb[NC*NS], NCN, NCN);
    setArrayDoubleValues(input->lb,lb,NC*(NS+1));
    double ub[NC*(NS+1)];
    repeatDoubleArray(ub_g, ub, NC, NC*NS);
    repeatDoubleArray(ub_gN, &ub[NC*NS], NCN, NCN);
    setArrayDoubleValues(input->ub,ub,NC*(NS+1));
    // input constraints
    double lbu_0[NU]  = {};
    double ubu_0[NU]  = {};
    double neg_inf = -INF;
    repeatDoubleArray(&neg_inf, lbu_0, 1, NU);
    double pos_inf = +INF;
    repeatDoubleArray(&pos_inf, ubu_0, 1, NU);
    int nbu_idx[NBU] = NBU_IDX;
    for (int i = 0; i < NBU; i++) {
        lbu_0[nbu_idx[i]-1] = lb_u[i];
        ubu_0[nbu_idx[i]-1] = ub_u[i];
    }
    double lbu[NU*NS];
    repeatDoubleArray(lbu_0, lbu, NU, NU*NS);
    setArrayDoubleValues(input->lbu,lbu,NU*NS);
    double ubu[NU*NS];
    repeatDoubleArray(ubu_0, ubu, NU, NU*NS);
    setArrayDoubleValues(input->ubu,ubu,NU*NS);
    // state constraints
    double lbx[NBX*NS];
    repeatDoubleArray(lb_x, lbx, NBX, NBX*NS);
    setArrayDoubleValues(input->lbx,lbx,NBX*NS);
    double ubx[NBX*NS];
    repeatDoubleArray(ub_x, ubx, NBX, NBX*NS);
    setArrayDoubleValues(input->ubx,ubx,NBX*NS);
    // states
    double x[NX*(NS+1)];
    repeatDoubleArray(x0, x, NX, NX*(NS+1));
    setArrayDoubleValues(input->x,x,NX*(NS+1));
    // inputs
    double u[NU*NS];
    repeatDoubleArray(u0, u, NU, NU*NS);
    setArrayDoubleValues(input->u,u,NU*NS);
    // algebraic states
    double z[NZ*NS];
    repeatDoubleArray(z0, z, NZ, NZ*NS);
    setArrayDoubleValues(input->z,z,NZ*NS);
    // params
    int params_num = NP*(NS+1);
    double params[params_num];
    repeatDoubleArray(para0, params, NP, params_num);
    setArrayDoubleValues(input->od,params,params_num);
    // cost function weights
    double Q[NY*NS];
    repeatDoubleArray(h, Q, NY, NY*NS);
    setArrayDoubleValues(input->W,Q,NY*NS);
    setArrayDoubleValues(input->WN,hN,NYN);
    // langrangian multipliers
    double lambda[NX*(NS+1)] = {0};
    setArrayDoubleValues(input->lambda,lambda,NX*(NS+1));
    double mu[NC*NS+NCN] = {};
    setArrayDoubleValues(input->mu,mu,NC*NS+NCN);
    double mu_u[NU*NS] = {0};
    setArrayDoubleValues(input->mu_u,mu_u,NU*NS);
    double mu_x[NBX*NS] = {0};
    setArrayDoubleValues(input->mu_x,mu_x,NBX*NS);
}

void initMem(TypeMem *mem){

    mem->warm_start=0;
    mem->hot_start=0; // see opt.hot_start

    // settings2 for opt.partial_condensing NOT IMPLEMENTED

    // solver initialization (see opt.qpsolver) IMPLEMENTED FOR HPIPM_SPARSE
    mem->mu0 = 1.0e4;
    mem->max_qp_it = 1e3;
    mem->pred_corr = 1;
    mem->cond_pred_corr = 1;
    mem->solver_mode = 3;

    // hessian setting IMPLEMENTED FOR GAUSS_NEWTON (see opt.hessian)
    mem->hessian = 0;

    // integrator settings IMPLEMENTED FOR ERK4 (see opt.integrator)
    mem->sim_method = 1;
//    double A_tab[16] = {  0,   0,   0,   0,
//                        0.5,   0,   0,   0,
//                          0, 0.5,  0,   0,
//                         0,  0,   1,  0};
    double A_tab[16] = {  0,   0.5, 0,   0,
                          0,   0,   0.5, 0,
                          0,   0,   0,  1,
                          0,  0,   0,  0};
    setArrayDoubleValues(mem->A_tab,A_tab,16);
    printf("Check row/column major for A_tab! Now it is column-major (as it seems to be Blasfeo). \n");
    double B_tab[4] = {1.0/6, 1.0/3, 1.0/3, 1.0/6};
    setArrayDoubleValues(mem->B_tab,B_tab,4);
    mem->num_steps = 2;
    mem->num_stages = 4;
    mem->h = TS_ST/mem->num_steps;
    mem->nx = NX;
    mem->nu = NU;
    mem->nz = NZ;
    double Sx[NX*NX] = {0};
    identityMatrixVector(Sx,NX);
    setArrayDoubleValues(mem->Sx ,Sx,NX*NX);
    double Su[NX*NU] = {0};
    setArrayDoubleValues(mem->Su ,Su,NX*NU);

    // RTI settings HARD CODED (see opt.RTI)
    mem->sqp_maxit = 1;
    mem->kkt_lim = 1.0e-2;

    // globalization
    mem->mu_merit = 0,
    mem->eta = 1.0e-8;
    mem->tau = 0.5;
    mem->mu_safty = 1.1;
    mem->rho = 0.5;
    mem->alpha = 1;
    mem->obj = 0;

    // initialize memory
    double A[NX*NX*NS] = {0};
    setArrayDoubleValues(mem->A ,A,NX*NX*NS);
//    double A[NX][NX*NS] = {0};
//    setMatrixDoubleValues(NX, NX*NS,mem->A ,A);
    double B[NX*NU*NS] = {0};
    setArrayDoubleValues(mem->B ,B,NX*NU*NS);
    double Cx[NBX*NX] = {0};
    setArrayDoubleValues(mem->Cx ,Cx,NBX*NX);
    double Cgx[NC*NX*NS] = {};
    setArrayDoubleValues(mem->Cgx ,Cgx,NC*NX*NS);
    double Cgu[NC*NU*NS] = {};
    setArrayDoubleValues(mem->Cgu ,Cgu,NC*NU*NS);
    double CgN[NCN*NX] = {};
    setArrayDoubleValues(mem->CgN ,CgN,NCN*NX);
    double gx[NX*(NS+1)] = {0};
    setArrayDoubleValues(mem->gx ,gx,NX*(NS+1));
    double gu[NU*NS] = {0};
    setArrayDoubleValues(mem->gu ,gu,NU*NS);
    double a[NX*NS] = {0};
    setArrayDoubleValues(mem->a ,a,NX*NS);
    double ds0[NX] = {0};
    setArrayDoubleValues(mem->ds0 ,ds0,NX);
    double lc[NS*NC+NCN] = {};
    setArrayDoubleValues(mem->lc ,lc,NS*NC+NCN);
    double uc[NS*NC+NCN] = {};
    setArrayDoubleValues(mem->uc ,uc,NS*NC+NCN);
    double lb_du[NS*NU] = {0};
    setArrayDoubleValues(mem->lb_du ,lb_du,NS*NU);
    double ub_du[NS*NU] = {0};
    setArrayDoubleValues(mem->ub_du ,ub_du,NS*NU);
    double lb_dx[NS*NBX] = {0};
    setArrayDoubleValues(mem->lb_dx ,lb_dx,NS*NBX);
    double ub_dx[NS*NBX] = {0};
    setArrayDoubleValues(mem->ub_dx ,ub_dx,NS*NBX);
    double Hc[NS*NU*NS*NU] = {0};
    setArrayDoubleValues(mem->Hc ,Hc,NS*NU*NS*NU);
    double Ccx[NS*NBX*NS*NU] = {0};
    setArrayDoubleValues(mem->Ccx ,Ccx,NS*NBX*NS*NU);
    double Ccg[(NS*NC+NCN)*NS*NU] = {};
    setArrayDoubleValues(mem->Ccg ,Ccg,(NS*NC+NCN)*NS*NU);
    double gc[NS*NU] = {0};
    setArrayDoubleValues(mem->gc ,gc,NS*NU);
    double lcc[NS*NC+NCN] = {};
    setArrayDoubleValues(mem->lcc ,lcc,NS*NC+NCN);
    double ucc[NS*NC+NCN] = {};
    setArrayDoubleValues(mem->ucc ,ucc,NS*NC+NCN);
    double lxc[NS*NBX] = {0};
    setArrayDoubleValues(mem->lxc ,lxc,NS*NBX);
    double uxc[NS*NBX] = {0};
    setArrayDoubleValues(mem->uxc ,uxc,NS*NBX);

    double dx[NX*(NS+1)] = {0};
    setArrayDoubleValues(mem->dx ,dx,NX*(NS+1));
    double du[NU*NS] = {0};
    setArrayDoubleValues(mem->du ,du,NU*NS);
    double lambda_new[NX*(NS+1)] = {0};
    setArrayDoubleValues(mem->lambda_new ,lambda_new,NX*(NS+1));
    double mu_new[NS*NC+NCN] = {};
    setArrayDoubleValues(mem->mu_new ,mu_new,NS*NC+NCN);
    double mu_x_new[NS*NX] = {0};
    setArrayDoubleValues(mem->mu_x_new ,mu_x_new,NS*NX);
    double mu_u_new[NS*NU] = {0};
    setArrayDoubleValues(mem->mu_u_new ,mu_u_new,NS*NU);
    double z_out[NS*NZ] = {};
    setArrayDoubleValues(mem->z_out ,z_out,NS*NZ);
    double q_dual[NX*(NS+1)] = {0};
    setArrayDoubleValues(mem->q_dual ,q_dual,NX*(NS+1));
    double dmu[NS*(NU+NBX+NC)+NCN] = {0};
    setArrayDoubleValues(mem->dmu ,dmu,NS*(NU+NBX+NC)+NCN);

    int nbx_idx[NBX] = NBX_IDX;
    for (int i=0;i<NBX;i++){
        mem->Cx[(nbx_idx[i]-1)*NBX+i]=1;
    }

    mem->reg = 1.0e-4;
    double Q[NX*NX*(NS+1)] = {0};
    setArrayDoubleValues(mem->Q ,Q,NX*NX*(NS+1));
    double S[NX*NU*NS] = {0};
    setArrayDoubleValues(mem->S ,S,NX*NU*NS);
    double R[NU*NU*NS] = {0};
    setArrayDoubleValues(mem->R ,R,NU*NU*NS);

    mem->iter=1;

    mem->gp_flag = 0;

    // Non-Uniform Grid NOT IMPLEMENTED YET (see opt.nonuniform_grid)
}



