#ifndef TYPEDEF_UTILS_H
#define TYPEDEF_UTILS_H

#include "mpc_set.h"


typedef struct
{
    double x0[NX];
    double u0[NU];
    double z0[NZ];
    double lb[NC*NS+NCN];
    double ub[NC*NS+NCN];
    double lbu[NU*NS];
    double ubu[NU*NS];
    double lbx[NBX*NS];
    double ubx[NBX*NS];
    double x[NX*(NS+1)];
    double u[NU*NS];
    double z[NZ*NS];
    double od[NP*(NS+1)];
    double W[NY*NS];
    double WN[NYN];
    double lambda[NX*(NS+1)];
    double mu[NS*NC+NCN];
    double mu_u[NS*NU];
    double mu_x[NS*NBX];
    double y[NY*NS];
    double yN[NYN];
} TypeInput;

typedef struct
{
    int warm_start;
    int hot_start;
    double mu0;
    int max_qp_it;
    int pred_corr;
    int cond_pred_corr;
    int solver_mode;
    int hessian;
    int sim_method;
    double A_tab[16];
    double B_tab[4];
    int num_steps;
    int num_stages;
    double h;
    double nx;
    double nu;
    double nz;
    double Sx[NX*NX];
    double Su[NX*NU];
    int sqp_maxit;
    double kkt_lim;
    double mu_merit;
    double eta;
    double tau;
    double mu_safty;
    double rho;
    double alpha;
    double obj;
    double A[NX*NX*NS];
//    double A[NX][NX*NS];
    double B[NX*NU*NS];
    double Cx[NBX*NX];
    double Cgx[NC*NX*NS];
    double Cgu[NC*NU*NS];
    double CgN[NCN*NX];
    double gx[NX*(NS+1)];
    double gu[NU*NS];
    double a[NX*NS];
    double ds0[NX];
    double lc[NS*NC+NCN];
    double uc[NS*NC+NCN];
    double lb_du[NS*NU];
    double ub_du[NS*NU];
    double lb_dx[NS*NBX];
    double ub_dx[NS*NBX];
    double Hc[NS*NU*NS*NU];
    double Ccx[NS*NBX*NS*NU];
    double Ccg[(NS*NC+NCN)*NS*NU];
    double gc[NS*NU];
    double lcc[NS*NC+NCN];
    double ucc[NS*NC+NCN];
    double lxc[NS*NBX];
    double uxc[NS*NBX];
    double dx[NX*(NS+1)];
    double du[NU*NS];
    double lambda_new[NX*(NS+1)];
    double mu_new[NS*NC+NCN];
    double mu_x_new[NS*NBX];
    double mu_u_new[NS*NU];
    double z_out[NZ*NS];
    double q_dual[NX*(NS+1)];
    double dmu[NS*(NU+NBX+NC)+NCN];
    double reg;
    double Q[NX*NX*(NS+1)];
    double S[NX*NU*NS];
    double R[NU*NU*NS];
    int iter;
    int sqp_it;
    int newton_iter;
    int gp_flag;
} TypeMem;

typedef struct
{
    char   model[100];
    double Ts_st;
    int nx;
    int nu;
    int nz;
    int ny;
    int nyN;
    int npODE;
    int T;
    int xGP;
    int yGP;
    int npGP;
    int np;
    int nc;
    int ncN;
    int nbx;
    int nbu;
    int nbx_idx[NBX];
    int nbu_idx[NBU];
    int N;
    int N2;
    int r;
    double ref_Ts;
} TypeSettings;

typedef struct
{
    double x[NX];
    double u[NU];
    double z[NZ];
    double p[NP*(NS+1)];
} TypeSimInput;

typedef struct
{
    double x[NX];
    double z[NZ];
} TypeSimOutput;

/*
typedef struct
{
    double H[6400];      //  H = mem.Hc;
    double g[80];        //  g = mem.gc;
    double A[6400+6480]; //  [mem.Ccx;mem.Ccg];
    double lb[80];       //  mem.lb_du;
    double ub[80];       //  mem.ub_du;
    double lbA[80];      //  [mem.lxc;mem.lcc]
    double ubA[80];      //  [mem.uxc;mem.ucc];
    double isSimplyBoundedQp;
    //TypeQpoases_opt qpoases_opt;
} TypeQPOASES_input;
 */

// typedef struct
// {
// 
// }TypeQPOASES_output;

#endif