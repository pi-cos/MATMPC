#include "string.h"

#include "typedef_utils.h"
#include "casadi_wrapper.h"
#include "sim.h"
#include "erk.h"
#include "mpc_common.h"
//#include "erk_gp.h"
//#include "irk_ode.h"
//#include "irk_dae.h"


#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

static sim_opts *opts = NULL;
static sim_in *in = NULL;
static sim_out *out = NULL;
static sim_erk_workspace *erk_workspace = NULL;
//static sim_erk_gp_workspace *erk_gp_workspace = NULL;
//static sim_irk_ode_workspace *irk_ode_workspace = NULL;
//static sim_irk_dae_workspace *irk_dae_workspace = NULL;
static bool mem_alloc = false;

static double *casadi_out[3];
static double *L = NULL;
static double *eq_res_vec = NULL;
static double *lu, *uu, *lx, *ux, *tmp;

void exitFcn_si(){   
    if (erk_workspace!=NULL)
        sim_erk_workspace_free(opts, erk_workspace);
//    if (erk_gp_workspace!=NULL)
//        sim_erk_gp_workspace_free(opts, erk_gp_workspace);
//    if (irk_ode_workspace!=NULL)
//        sim_irk_ode_workspace_free(opts, irk_ode_workspace);
//    if (irk_dae_workspace!=NULL)
//        sim_irk_dae_workspace_free(opts, irk_dae_workspace);
    if (opts!=NULL)
        sim_opts_free(opts);
    if (in!=NULL)
        sim_in_free(in);
    if (out!=NULL)
        sim_out_free(out);
    if (mem_alloc){
        free(casadi_out[1]);
        free(casadi_out[2]);
        free(L);
        free(eq_res_vec);
        free(lu);
        free(uu);
        free(lx);
        free(ux);
        free(tmp);
    }
}

int solution_info(TypeInput *input,TypeMem *mem, TypeSettings *settings, double *eq_res_out,double *ineq_res_out,double* KKT_out,double* OBJ_out)
{
    double *x = input[0].x;
    double *u = input[0].u;
    double *z = input[0].z;
    double *lambda = input[0].lambda;
    double *mu = input[0].mu;
    double *mu_u = input[0].mu_u;
    double *mu_x = input[0].mu_x;
    double *y = input[0].y;
    double *yN = input[0].yN;
    double *od = input[0].od;
    double *W = input[0].W;
    double *WN = input[0].WN;
    double *lb = input[0].lb;
    double *ub = input[0].ub;
    double *lbu = input[0].lbu;
    double *ubu = input[0].ubu;
    double *lbx = input[0].lbx;
    double *ubx = input[0].ubx;

    double *ds0 = mem[0].ds0;
    double *lc = mem[0].lc;
    double *uc = mem[0].uc;
    double *obj = &(mem[0].obj);
    double *z_out = mem[0].z_out;

    int nx =settings->nx;
    int nu =settings->nu;
    int nz =settings->nz;
    int np =settings->np; if(np==0) np++;
    int ny =settings->ny;
    int nyN =settings->nyN;
    int nc =settings->nc;
    int ncN =settings->ncN;
    int nbx =settings->nbx;
    int *nbx_idx = settings[0].nbx_idx;
    int N =settings->N;
    double Ts_st = settings->Ts_st;
    
    int sim_method = mem->sim_method;
          
    int i=0,j=0,l=0;
    int nv = nx+nu;
    int nw = N*nv+nx;
    int neq = (N+1)*nx;
    int nineq = N*nc+ncN;
    char *nTrans = "N", *Trans="T", *Norm="O";
    double one_d = 1.0, zero = 0.0;
    int one_i = 1;
    int idx;
    
    double *casadi_in[9];
        
    if (!mem_alloc){
        casadi_out[1] = (double *)malloc(nv * sizeof(double));
        //mexMakeMemoryPersistent(casadi_out[1]);
        casadi_out[2] = (double *)malloc(nv * sizeof(double));
        //mexMakeMemoryPersistent(casadi_out[2]);     
        L = (double *)malloc( nw * sizeof(double));
        //mexMakeMemoryPersistent(L);
        
        eq_res_vec = (double *)malloc( neq * sizeof(double));
        //mexMakeMemoryPersistent(eq_res_vec);
        
        lu = (double *)malloc( N*nu * sizeof(double));
        //mexMakeMemoryPersistent(lu);
        uu = (double *)malloc( N*nu * sizeof(double));
        //mexMakeMemoryPersistent(uu);
        lx = (double *)malloc( N*nbx * sizeof(double)); 
        //mexMakeMemoryPersistent(lx);
        ux = (double *)malloc( N*nbx * sizeof(double));
        //mexMakeMemoryPersistent(ux);

        tmp = (double *)malloc(nbx*sizeof(double));
        memset(tmp,0,nbx*sizeof(double));
        
        switch(sim_method){
            case 1:
                opts = sim_opts_create(mem);
                opts->forw_sens_flag = false;
                opts->adj_sens_flag = true;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                erk_workspace = sim_erk_workspace_create(opts);               
                sim_erk_workspace_init(opts, mem, erk_workspace);
                break;
            /*
            case 2:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens_flag = false;
                opts->adj_sens_flag = true;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                irk_ode_workspace = sim_irk_ode_workspace_create(opts);               
                sim_irk_ode_workspace_init(opts, prhs[2], irk_ode_workspace);
                break;
            case 3:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens_flag = false;
                opts->adj_sens_flag = true;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                irk_dae_workspace = sim_irk_dae_workspace_create(opts);               
                sim_irk_dae_workspace_init(opts, prhs[2], irk_dae_workspace);
                break;
            case 4:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens_flag = false;
                opts->adj_sens_flag = true;
                opts->gp_status_flag = false; // no gp info in output
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                erk_gp_workspace = sim_erk_gp_workspace_create(opts);               
                sim_erk_gp_workspace_init(opts, prhs[2], erk_gp_workspace); 
                break;
            */
            default:
                printf("Please choose a supported integrator");
                return -1;
//                break;
         
        }   
                
        mem_alloc = true;     
//        mexAtExit(exitFcn);
    }
    
    memcpy(eq_res_vec, ds0, nx*sizeof(double));
    
    double *work;
    double KKT=0, eq_res=0, ineq_res=0, OBJ=0;
    
    for (i=0;i<N;i++){
        casadi_in[0] = x+i*nx;
        casadi_in[1] = u+i*nu;
        casadi_in[2] = od+i*np;
        casadi_in[3] = y+i*ny;
        casadi_in[4] = W+i*ny;
        casadi_in[5] = lambda+(i+1)*nx;
        if (i==0){
            casadi_in[6] = tmp;
        }
        else{
            casadi_in[6] = mu_x+(i-1)*nbx;
        }
        casadi_in[7] = mu_u+i*nu;
        casadi_in[8] = mu+i*nc;

        // objective value      
        casadi_out[0] = obj;
        obji_Fun(casadi_in, casadi_out);
        OBJ += obj[0];
              
        // call integrator to compute xend and adjoint sensitivities
        switch(sim_method){
            case 1:      
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->lambda = lambda+(i+1)*nx;
                out->xn = eq_res_vec+(i+1)*nx;
                out->adj_sens = casadi_out[2];
                sim_erk(in, out, opts, erk_workspace);
                break;
            /*
            case 2:                         
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->z = z+i*nz;
                in->lambda = lambda+(i+1)*nx;
                out->xn = eq_res_vec+(i+1)*nx;
                out->adj_sens = casadi_out[2];
                sim_irk_ode(in, out, opts, irk_ode_workspace);
                break;
            case 3:                         
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->z = z+i*nz;
                in->lambda = lambda+(i+1)*nx;
                out->xn = eq_res_vec+(i+1)*nx;
                out->zn = z_out+i*nz;
                out->adj_sens = casadi_out[2];
                sim_irk_dae(in, out, opts, irk_dae_workspace);
                break;
            case 4:      
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->lambda = lambda+(i+1)*nx;
                out->xn = eq_res_vec+(i+1)*nx;
                out->adj_sens = casadi_out[2];
                opts->gp_status_flag = false; // no gp info in output
                sim_erk_gp(in, out, opts, erk_gp_workspace);
                break;
            */
            default :
                printf("Please choose a supported integrator");
                return -1;
//                break;
        }

        // compute gradient of Lagrangian
        if (i>0){
            for (j=0;j<nx;j++)
                casadi_out[2][j] -= lambda[i*nx+j];
        }else{
            for (j=0;j<nx;j++)
                casadi_out[2][j] += lambda[j];
        }
        casadi_out[0] = L+i*nv;
        adj_Fun(casadi_in, casadi_out);  // dojb                    
        daxpy(&nv, &one_d, casadi_out[1], &one_i, casadi_out[0], &one_i); // dojb+dB'*mu -> disequality constraints on x, u, general
        daxpy(&nv, &one_d, casadi_out[2], &one_i, casadi_out[0], &one_i); // dobj+dB'*mu+dG'*lambda -> equality constraints (dynamics)
        
        // equality constraint residual
        for (j=0;j<nx;j++)
            eq_res_vec[(i+1)*nx+j] -= x[(i+1)*nx+j];
        
        // control bound residual
        for (j=0;j<nu;j++){
            lu[i*nu+j] = lbu[i*nu+j] - u[i*nu+j];
            uu[i*nu+j] = ubu[i*nu+j] - u[i*nu+j];
        }
        
        // state bound residual
        for (j=0;j<nbx;j++){
            idx = (int)nbx_idx[j]-1;
            lx[i*nbx+j] = lbx[i*nbx+j] - x[(i+1)*nx+idx];
            ux[i*nbx+j] = ubx[i*nbx+j] - x[(i+1)*nx+idx];
        }
        
        // general constraint residual
        if (nc>0){
            casadi_out[0] = lc + i*nc;
            path_con_Fun(casadi_in, casadi_out);
            for (j=0;j<nc;j++){
                uc[i*nc+j] = ub[i*nc+j] - casadi_out[0][j];
                casadi_out[0][j] = lb[i*nc+j] - casadi_out[0][j];            
            }
        }
    }
    casadi_in[0] = x+N*nx;
    casadi_in[1] = od+N*np;
    casadi_in[2] = yN;
    casadi_in[3] = WN;
    casadi_in[4] = mu_x+(N-1)*nbx;
    casadi_in[5] = mu+N*nc;
    
    casadi_out[0] = L+N*nv;
    adjN_Fun(casadi_in, casadi_out);
    for (j=0;j<nx;j++)
        casadi_out[0][j] -= lambda[N*nx+j];
    
    daxpy(&nx, &one_d, casadi_out[1], &one_i, casadi_out[0], &one_i);
               
    if (ncN>0){
        casadi_out[0] = lc + N*nc;
        path_con_N_Fun(casadi_in, casadi_out); 
        for (j=0;j<ncN;j++){
            uc[N*nc+j] = ub[N*nc+j] - casadi_out[0][j];
            casadi_out[0][j] = lb[N*nc+j] - casadi_out[0][j];            
        }
    }
         
    eq_res = dlange(Norm, &neq, &one_i, eq_res_vec, &neq, work);
    KKT = dnrm2(&nw, L, &one_i);
    
    for (i=0;i<N*nu;i++)
        ineq_res += MAX(-1*uu[i],0) + MAX(lu[i],0);
    for (i=0;i<N*nbx;i++)
        ineq_res += MAX(-1*ux[i],0) + MAX(lx[i],0);
    for (i=0;i<nineq;i++)
        ineq_res += MAX(-1*uc[i],0) + MAX(lc[i],0);
    
    casadi_out[0] = obj;
    objN_Fun(casadi_in, casadi_out);
    OBJ += obj[0];
    
    OBJ *= Ts_st;

    memcpy(z, z_out, N*nz*sizeof(double));

    *eq_res_out = eq_res;
    *ineq_res_out = ineq_res;
    *KKT_out = KKT;
    *OBJ_out = OBJ;

//    exitFcn_si();
    return 1;

}