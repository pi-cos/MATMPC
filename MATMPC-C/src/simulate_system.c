
//#include "mex.h"
#include "string.h"

#include "sim.h"
#include "erk.h"
#include "typedef_utils.h"
//#include "erk_gp.h"
//#include "irk_ode.h"
//#include "irk_dae.h"

static sim_opts *opts = NULL;
static sim_in *in = NULL;
static sim_out *out = NULL;
static sim_erk_workspace *erk_workspace = NULL;
//static sim_erk_gp_workspace *erk_gp_workspace = NULL;
//static sim_irk_ode_workspace *irk_ode_workspace = NULL;
//static sim_irk_dae_workspace *irk_dae_workspace = NULL;
static bool mem_alloc = false;

void exitFcn_ss(){
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
}

void simulate_system(TypeSimInput *sim_input, TypeMem *mem, TypeSettings *settings, TypeSimOutput *sim_output)
{
    double *x = sim_input[0].x;
    double *u = sim_input[0].u;
    double *z = sim_input[0].z;
    double *od = sim_input[0].p;

    int nx = settings->nx;
    int nz = settings->nz;
    int sim_method = mem->sim_method;

    double *x_out = sim_output[0].x;
    double *z_out = sim_output[0].z;
              
    if (!mem_alloc){
        opts = sim_opts_create(mem);
        opts->forw_sens_flag = false;
        opts->adj_sens_flag = false;
        in = sim_in_create(opts);              
        out = sim_out_create(opts);     
        switch(sim_method){
            case 1:                         
                erk_workspace = sim_erk_workspace_create(opts);               
                sim_erk_workspace_init(opts, mem, erk_workspace);
                break;
            /*
            case 2:          
                irk_ode_workspace = sim_irk_ode_workspace_create(opts);               
                sim_irk_ode_workspace_init(opts, prhs[4], irk_ode_workspace);
                break;
            case 3:
                irk_dae_workspace = sim_irk_dae_workspace_create(opts);               
                sim_irk_dae_workspace_init(opts, prhs[4], irk_dae_workspace);
                break;
            case 4:                         
                erk_gp_workspace = sim_erk_gp_workspace_create(opts);               
                sim_erk_gp_workspace_init(opts, prhs[4], erk_gp_workspace);                
                break;
            */
            default:
                printf("Please choose a supported integrator");
                break;
        }  
                
        mem_alloc=true;
//        mexAtExit(exitFcn);
    }
    
        
     // integration                      
    switch(sim_method){
        case 1:
            in->x = x;
            in->u = u;
            in->p = od;
            out->xn = x_out;
            sim_erk(in, out, opts, erk_workspace);
            break;
        /*
        case 2:
            in->x = x;
            in->u = u;
            in->p = od;
            in->z = z;
            out->xn = x_out;
            sim_irk_ode(in, out, opts, irk_ode_workspace);
            break;
        case 3:
            in->x = x;
            in->u = u;
            in->p = od;
            in->z = z;
            out->xn = x_out;
            out->zn = z_out;
            sim_irk_dae(in, out, opts, irk_dae_workspace);
            break;
        case 4:
            in->x = x;
            in->u = u;
            in->p = od;
            out->xn = x_out;
            opts->gp_status_flag = false; // no gp info in output
            sim_erk_gp(in, out, opts, erk_gp_workspace);
            break;
        */
        default:
            printf("Please choose a supported integrator");
            break;
    }

//    exitFcn_ss();
        
}