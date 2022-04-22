#ifndef MPC_SET_H_
#define MPC_SET_H_

#define MODEL_NAME "InvertedPendulum" // Name of the compiled model from Matlab

#define NS 80 // No. of shooting points
#define NS2 10 // No. of horizon length after partial condensing (if used)
#define RS 10 // No. of input blocks for Non-Uniform Grid (if used)

#define NX 4 // No. of differential states
#define NU 1 // No. of controls
#define NZ 0 // No. of algebraic states
#define NY 5 // No. of outputs
#define NYN 4 // No. of outputs at the terminal point
#define NP 0 // No. of model parameters (NOTE: it must be NPODE+NPGP)
#define NC 0 // No. of general constraints
#define NCN 0 // No. of general constraints at the terminal point
#define NBX 1 // No. of bounds on states
#define NBU 1 // No. of bounds on controls
#define NBX_IDX  {1}// indexes of state bounds (NOTE: from 1 to NX)
#define NBU_IDX {1} // indexes of control bounds (NOTE: from 1 to NU)

#define TS_ST 0.025 // NMPC shooting interval


#endif