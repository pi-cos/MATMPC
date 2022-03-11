#include <time.h>

#include "typedef_utils.h"
#include "main_utils.h"
#include "mpc_fcn.h"

#define MAX(a,b) (((a)>(b))?(a):(b))


int nmpc_solver(TypeInput *input, TypeMem *mem, TypeSettings *settings);
int simulation_loop(TypeInput *input, TypeMem *mem, TypeSettings *settings, int sim_length, double reference[],
                    double state_sim[sim_length+1][settings->nx], double input_sim[sim_length][settings->nu]);
void write_results(int nx, int steps, double state_sim[steps+1][nx], int nu, double input_sim[steps][nu]);

/* main */
int main() {

    // strctures initialization
    TypeInput input;
    TypeMem mem;
    TypeSettings settings;

    initSettings(&settings);
    initInput(&input);
    initMem(&mem);

    // reference definition
    double ref[NY] = {0};

    // simulation
    int sim_length = 160;
    double state_sim[sim_length+1][settings.nx];
    setArrayDoubleValues(state_sim[0],input.x0,settings.nx);
    double input_sim[sim_length][settings.nu];

    mem.iter = 1;
    int result = simulation_loop(&input,&mem,&settings,sim_length,ref,state_sim,input_sim);

    write_results(settings.nx, sim_length, state_sim, settings.nu, input_sim);

    return result;
}

int simulation_loop(TypeInput *input, TypeMem *mem, TypeSettings *settings, int sim_length, double reference[],
                    double state_sim[sim_length+1][settings->nx], double input_sim[sim_length][settings->nu]){

    // results init
    int result_sl = 0;
    int result_nmpc;

    // sim struct definition
    TypeSimInput sim_input;
    TypeSimOutput sim_output;

    // simulation state
    printf("Starting simulation, system state is: ");
    printDoubleArrayLine(state_sim[0],settings->nx);
    printf("\n");

    // set time-invarying reference, params and constraints
    setArrayDoubleValues(input->y,reference,settings->ny);
    setArrayDoubleValues(input->yN,reference,settings->nyN);

    for (int i = 0; i<sim_length; i++){

        // set input.x0
        setArrayDoubleValues(input->x0,state_sim[i],settings->nx);

        // set time-varying reference, params and constraints (if any)


        // nmpc call
        clock_t start = clock(), diff;
        result_nmpc = nmpc_solver(input,mem,settings);
        diff = clock() - start;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Iteration no. %d,  CPT = %d ms, SQP_IT = %d", mem->iter, msec,mem->sqp_it);

        // set the simulation input
        setArrayDoubleValues(sim_input.x, state_sim[i], settings->nx);
        setArrayDoubleValues(sim_input.u, input->u, settings->nu);
        setArrayDoubleValues(sim_input.p, input->od, settings->np);
        // remember to set z if present!!
        // simulate the system
        simulate_system(&sim_input, mem, settings, &sim_output);
        // print in console
        printf(", u = ");
        printDoubleArrayLine(sim_input.u, settings->nu);
        printf(" x = ");
        printDoubleArrayLine(sim_output.x, settings->nx);
        printf("\n");

        // evolve state and save everything
        mem->iter++;
        setArrayDoubleValues(state_sim[i+1],sim_output.x,settings->nx);
        setArrayDoubleValues(input_sim[i],sim_input.u,settings->nu);

        if (result_nmpc != 1) {
            result_sl = -1;
            break;
        }
        else
            result_sl = +1;
    }

    return result_sl;

}

int nmpc_solver(TypeInput *input, TypeMem *mem, TypeSettings *settings){

    double eq_res_out = 0;
    double ineq_res_out = 0;
    double KKT_out= 0;
    double OBJ_out= 0;

    int result_qp = 0;
    int result_hpipm = 0;
    int result_ls = 0;
    int result_si = 0;

    mem->sqp_it=0;
    mem->alpha=1;
    double StopCrit = 2*mem->kkt_lim;

    while (mem->sqp_it < mem->sqp_maxit && StopCrit > mem->kkt_lim && mem->alpha>1.0e-8){

        result_qp = qp_generation(input,mem,settings);
        result_hpipm = hpipm_sparse(mem,settings);
        result_ls = line_search(input,mem,settings);
        result_si = solution_info(input,mem,settings,&eq_res_out,&ineq_res_out,&KKT_out,&OBJ_out);

        StopCrit = MAX(eq_res_out, ineq_res_out);
        StopCrit = MAX(StopCrit, KKT_out);

        mem->sqp_it=mem->sqp_it+1;

    }

    if (result_qp != 1 || result_hpipm != 1 || result_ls != 1 || result_si != 1)
        return -1;
    else
        return 1;


}

void write_results(int nx, int steps, double state_sim[steps+1][nx], int nu, double input_sim[steps][nu]){
    FILE *fpt;
    fpt = fopen("../csv/sim_results.csv", "w+");
    fprintf(fpt,"p theta v omega u\n");

    for (int i=0;i<steps+1;i++){
        for (int j=0; j<nx;j++){
            fprintf(fpt,"%f ", state_sim[i][j]);
        }
        if (i<steps) {
            for (int j = 0; j<nu-1; j++) {
                fprintf(fpt, "%f ", input_sim[i][j]);
            }
            fprintf(fpt,"%f\n", input_sim[i][nu-1]);
        }
        else
            fprintf(fpt,"\n");
    }
    fclose(fpt);
};
