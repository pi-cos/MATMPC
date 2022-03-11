

#include <string.h>

#include "typedef_utils.h"
#include "sim.h"
#include "erk.h"
#include "main_utils.h"
//#include "irk_ode.h"
//#include "irk_dae.h"
#include "casadi_wrapper.h"
#include "mpc_common.h"

static sim_opts *opts = NULL;
static sim_in *in = NULL;
static sim_out *out = NULL;
static sim_erk_workspace *erk_workspace = NULL;
//static sim_irk_ode_workspace *irk_ode_workspace = NULL;
//static sim_irk_dae_workspace *irk_dae_workspace = NULL;
static bool mem_alloc = false;
static double *Hes[1];
static double *HesN[1];
static double *Jac[2];
static double *JacN[1];
static double *temp[3];

void exitFcn_qp(){
    if (erk_workspace!=NULL)
        sim_erk_workspace_free(opts, erk_workspace);
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
    if (Hes[0]!=NULL)
        free(Hes[0]);
    if (HesN[0]!=NULL)
        free(HesN[0]);
    if (Jac[0]!=NULL)
        free(Jac[0]);
    if (Jac[1]!=NULL)
        free(Jac[1]);
    if (JacN[0]!=NULL)
        free(JacN[0]);
    if (temp[0]!=NULL)
        free(temp[0]);
    if (temp[1]!=NULL)
        free(temp[1]);
    if (temp[2]!=NULL)
        free(temp[2]);
}

int qp_generation(TypeInput *input,TypeMem *mem, TypeSettings *settings)
{

    double *x = input[0].x;
    double *u = input[0].u;
    double *z = input[0].z;
    double *y = input[0].y;
    double *yN = input[0].yN;
    double *od = input[0].od;
    double *W = input[0].W;
    double *WN = input[0].WN;
    double *lb = input[0].lb;
    double *ub = input[0].ub;
    double *x0 = input[0].x0;
    double *lbu = input[0].lbu;
    double *ubu = input[0].ubu;
    double *lbx = input[0].lbx;
    double *ubx = input[0].ubx;

    int nx = settings->nx;
    int nu = settings->nu;
    int nz = settings->nz;
    int np = settings->np; if(np==0) np++;
    int ny = settings->ny;
    int nyN = settings->nyN;
    int nc = settings->nc;
    int ncN = settings->ncN;
    int nbx = settings->nbx;
    int *nbx_idx = settings[0].nbx_idx;
//    size_t N = settings->N;
    int N = settings->N;
    int sim_method = mem->sim_method;
    
    int i=0,j=0;
    char *nTrans = "N", *Trans="T", *UPLO="L";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    int one_i = 1;
    int idx;

    double *Q = mem[0].Q;
    double *S = mem[0].S;
    double *R = mem[0].R;
    double *A = mem[0].A;
    double *B = mem[0].B;
    double *Cx = mem[0].Cx;
    double *Cgx = mem[0].Cgx;
    double *Cgu = mem[0].Cgu;
    double *CgN = mem[0].CgN;
    double *gx = mem[0].gx;
    double *gu = mem[0].gu;
    double *a = mem[0].a;
    double *ds0 = mem[0].ds0;
    double *lc = mem[0].lc;
    double *uc = mem[0].uc;
    double *lb_du = mem[0].lb_du;
    double *ub_du = mem[0].ub_du;
    double *lb_dx = mem[0].lb_dx;
    double *ub_dx = mem[0].ub_dx;
    double *z_out = mem[0].z_out;
    
    double reg = mem->reg;
    int hessian_type = mem->hessian;
    
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - x[i];
    
    // allocate memory    
    double *Cons[2];
      
    double *casadi_in[5];
    double *casadi_out[2];    
          
    if (!mem_alloc){
        switch(sim_method){
            case 1:
                opts = sim_opts_create(mem);
                opts->forw_sens_flag = true;
                opts->adj_sens_flag = false;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                erk_workspace = sim_erk_workspace_create(opts);               
                sim_erk_workspace_init(opts, mem, erk_workspace);                
                break;
            /*
            case 2:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens_flag = true;
                opts->adj_sens_flag = false;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                irk_ode_workspace = sim_irk_ode_workspace_create(opts);               
                sim_irk_ode_workspace_init(opts, prhs[2], irk_ode_workspace);
                break;
            case 3:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens_flag = true;
                opts->adj_sens_flag = false;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                irk_dae_workspace = sim_irk_dae_workspace_create(opts);               
                sim_irk_dae_workspace_init(opts, prhs[2], irk_dae_workspace);
                break;
            */
            default:
                printf("Please choose a supported integrator");
                return -1;
                break;
        }  
        
        Hes[0] = (double *) malloc(ny*ny * sizeof(double));	
        //mexMakeMemoryPersistent(Hes[0]); 	
        HesN[0] = (double *) malloc(nyN*nyN * sizeof(double));	
        //mexMakeMemoryPersistent(HesN[0]); 	
                       
        Jac[0] = (double *) malloc(ny*nx * sizeof(double));	
        //mexMakeMemoryPersistent(Jac[0]); 	
        Jac[1] = (double *) malloc(ny*nu * sizeof(double));	
        //mexMakeMemoryPersistent(Jac[1]); 	
        JacN[0] = (double *) malloc(nyN*nx * sizeof(double));	
        //mexMakeMemoryPersistent(JacN[0]);       
        
        temp[0] = (double *) malloc(ny*nx * sizeof(double));	
        //mexMakeMemoryPersistent(temp[0]); 	
        temp[1] = (double *) malloc(ny*nu * sizeof(double));	
        //mexMakeMemoryPersistent(temp[1]);        
        temp[2] = (double *) malloc(nyN*nx * sizeof(double));	
        //mexMakeMemoryPersistent(temp[2]); 
        
        mem_alloc=true;
//        mexAtExit(exitFcn_qp);
    }
    
    // start loop
    for(i=0;i<N;i++){
        casadi_in[0] = x+i*nx;
        casadi_in[1] = u+i*nu;
        casadi_in[2] = od+i*np;
        casadi_in[3] = y+i*ny;
        casadi_in[4] = W+i*ny;
        
        // control bounds
        for (j=0;j<nu;j++){
            lb_du[i*nu+j] = lbu[i*nu+j]-u[i*nu+j];
            ub_du[i*nu+j] = ubu[i*nu+j]-u[i*nu+j];
        }
        
        // state bounds
        for (j=0;j<nbx;j++){
            idx = (int)nbx_idx[j]-1;
            lb_dx[i*nbx+j] = lbx[i*nbx+j]-x[(i+1)*nx+idx];
            ub_dx[i*nbx+j] = ubx[i*nbx+j]-x[(i+1)*nx+idx];
        }
        
        // integration                      
        switch(sim_method){
            case 1:
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                out->xn = a+i*nx;
                out->Sx = A + i*nx*nx;
                out->Su = B + i*nx*nu;
                sim_erk(in, out, opts, erk_workspace);
                break;
            /*
            case 2:
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->z = z+i*nz;
                out->xn = a+i*nx;
                out->Sx = A + i*nx*nx;
                out->Su = B + i*nx*nu;
                sim_irk_ode(in, out, opts, irk_ode_workspace);
                break;
            case 3:
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->z = z+i*nz;
                out->xn = a+i*nx;
                out->Sx = A + i*nx*nx;
                out->Su = B + i*nx*nu;
                out->zn = z_out + i*nz;
                sim_irk_dae(in, out, opts, irk_dae_workspace);
                break;
            */
            default:
                printf("Please choose a supported integrator");
                return -1;
                break;
        }
       
        // equality residual        
        for (j=0;j<nx;j++)
            a[i*nx+j] -= x[(i+1)*nx+j];
       
        // Hessian
        Ji_Fun(casadi_in, Jac);
        switch(hessian_type){                
            case 0:                   
                dgemm(Trans, nTrans, &nx, &nx, &ny, &one_d, Jac[0], &ny, Jac[0], &ny, &zero, Q+i*nx*nx, &nx);
                dgemm(Trans, nTrans, &nx, &nu, &ny, &one_d, Jac[0], &ny, Jac[1], &ny, &zero, S+i*nx*nu, &nx);
                dgemm(Trans, nTrans, &nu, &nu, &ny, &one_d, Jac[1], &ny, Jac[1], &ny, &zero, R+i*nu*nu, &nu);
                break;
            
            case 1:
                Hi_Fun(casadi_in, Hes);
                dgemm(Trans, nTrans, &nx, &ny, &ny, &one_d, Jac[0], &ny, Hes[0], &ny, &zero, temp[0], &nx);
                dgemm(nTrans, nTrans, &nx, &nx, &ny, &one_d, temp[0], &nx, Jac[0], &ny, &zero, Q+i*nx*nx, &nx);
                
                dgemm(nTrans, nTrans, &nx, &nu, &ny, &one_d, temp[0], &nx, Jac[1], &ny, &zero, S+i*nx*nu, &nx);
                
                dgemm(Trans, nTrans, &nu, &ny, &ny, &one_d, Jac[1], &ny, Hes[0], &ny, &zero, temp[1], &nu);
                dgemm(nTrans, nTrans, &nu, &nu, &ny, &one_d, temp[1], &nu, Jac[1], &ny, &zero, R+i*nu*nu, &nu);
                
                break;
            default:
                printf("Please choose a supported Hessian type");
                return -1;
                break;
                            
        }
        regularization(nx, Q+i*nx*nx, reg);
        regularization(nu, R+i*nu*nu, reg);
        
        // gradient
        casadi_out[0] = gx+i*nx;
        casadi_out[1] = gu+i*nu;
        gi_Fun(casadi_in, casadi_out);
                
        // constraint residual
        if (nc>0){  
            casadi_out[0] = lc + i*nc;
            path_con_Fun(casadi_in, casadi_out);
            for (j=0;j<nc;j++){
                uc[i*nc+j] = ub[i*nc+j] - casadi_out[0][j];
                casadi_out[0][j] = lb[i*nc+j] - casadi_out[0][j];            
            }
        
            // constraint Jacobian
            Cons[0] = Cgx+i*nc*nx;
            Cons[1] = Cgu+i*nc*nu;
            Ci_Fun(casadi_in, Cons);
        }
    }
    
    // terminal data
    casadi_in[0] = x+N*nx;
    casadi_in[1] = od+N*np;
    casadi_in[2] = yN;
    casadi_in[3] = WN;
    
    JN_Fun(casadi_in, JacN);
    switch(hessian_type){            
        case 0:
            dgemm(Trans, nTrans, &nx, &nx, &nyN, &one_d, JacN[0], &nyN, JacN[0], &nyN, &zero, Q+N*nx*nx, &nx);
            break;
        case 1:
            HN_Fun(casadi_in, HesN);
            dgemm(Trans, nTrans, &nx, &nyN, &nyN, &one_d, JacN[0], &nyN, HesN[0], &nyN, &zero, temp[2], &nx);
            dgemm(nTrans, nTrans, &nx, &nx, &nyN, &one_d, temp[2], &nx, JacN[0], &nyN, &zero, Q+N*nx*nx, &nx);
            break;
        default:
            printf("Please choose a supported Hessian type");
            return -1;
            break;
            
    }
    regularization(nx, Q+N*nx*nx, reg);
        
    casadi_out[0] = gx+N*nx;
    gN_Fun(casadi_in, casadi_out);

    if (ncN>0){
        casadi_out[0] = lc + N*nc;
        path_con_N_Fun(casadi_in, casadi_out);
        for (j=0;j<ncN;j++){
            uc[N*nc+j] = ub[N*nc+j] - casadi_out[0][j];
            casadi_out[0][j] = lb[N*nc+j] - casadi_out[0][j];            
        }

        CN_Fun(casadi_in, &CgN);
    }

//    exitFcn_qp();
    return 1;
    
}