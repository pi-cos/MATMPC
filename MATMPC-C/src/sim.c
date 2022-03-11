#include <stdlib.h>

#include "typedef_utils.h"
#include "sim.h"
#include "stdio.h"
#include <stdbool.h>


sim_opts* sim_opts_create(const TypeMem *mem)
{
    sim_opts* opts = (sim_opts*)malloc(sizeof(sim_opts));
    //mexMakeMemoryPersistent(opts);
       
    opts->nx = mem[0].nx;
    opts->nu =mem[0].nu;
    opts->nz = mem[0].nz;
    opts->num_stages = mem[0].num_stages;
    opts->num_steps = mem[0].num_steps;
    opts->h = mem[0].h;
    
	return opts;
}

sim_in* sim_in_create(sim_opts *opts)
{
    sim_in* in = (sim_in*)malloc(sizeof(sim_in));
    //mexMakeMemoryPersistent(in);
    
    return in;
}

sim_out* sim_out_create(sim_opts *opts)
{
    sim_out* out = (sim_out*)malloc(sizeof(sim_out));
    //mexMakeMemoryPersistent(out);
    
    return out;
}

void sim_opts_free(sim_opts *opts)
{
    free(opts);
}

void sim_in_free(sim_in *in)
{
    free(in);
}

void sim_out_free(sim_out *out)
{
    free(out);
}