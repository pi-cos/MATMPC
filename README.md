# MATMPC-C

## C implementation of MATMPC

This tool aims at providing an easy-to-use NMPC implementation, based on MATMPC (https://github.com/chenyutao36/MATMPC) and hpipm (https://github.com/giaf/hpipm). The peculiar feature is that it is based on MATLAB-based MATMPC, in order to ease the prototyping and debugging, and allows to easily deploy the code from the MATLAB generated files.

## To use MATMPC-C, follow the steps below.

Once the NMPC controller has been developed and tested in MATMPC, it is possible to test the c implementation by:
Within MATMPC-C folder:

1. set all the MACROS in include/mpc_set.h

2. initialize the system in src/mpc_fcn.c -> initInput() function

3. set your simulation parameters in src/main

4. copy into src the file 'casadi_src.c' from MATMPC-MATLAB/mex_core


## To build and run

To build, run the following commands in the terminal:

1. cd build
2. cmake ..
3. cmake --build .

Run the executable by typing ./matmpc_c

To verify the results copy the file 'sim_results.csv' from csv folder into MATMPC-MATLAB/data and compare the obtained results launching MATMPC-MATLAB/draw_from_csv.m in MATLAB.


## MATMPC: MATLAB based nonlinear MPC tool

This tool aims at providing an easy-to-use nonlinear MPC implementation. The optimal control problem (OCP) that should be solved is transcribed by multiple shooting and the resulting nonlinear program (NLP) is solved by Sequential Quadratic Programming (SQP) method.

The tool supports fixed step (explicit/implicit) Runge-Kutta (RK) integrator for multiple shooting. The derivatives that are needed to perform optimization are obtained by CasADi (https://github.com/casadi/casadi/wiki), the state-of-the-art automatic/algorithmic differentiation toolbox. The Quadratic Programming (QP) problems can be solved by both dense and sparse solvers. By now, MATMPC supports interfaces with the following external solvers: qpOASES (https://projects.coin-or.org/qpOASES/wiki/QpoasesInstallation), Ipopt (https://projects.coin-or.org/Ipopt), hpipm (https://github.com/giaf/hpipm), osqp (https://osqp.org/) and qpalm (https://github.com/Benny44/QPALM).

**The most unique feature of MATMPC is that it does not require to install any external libraries. The users need not to understand how to make, compile and link any library. Except for external QP solvers, all algorithmic routines are written directly using MATLAB C API and can be compiled into independent MEX functions using compilers that belong to GCC class (e.g. GCC, MinGW and Clang). MATMPC employs MATLAB built-in linear algebra library provided by Intel MKL. Therefore, MATMPC is able to provide compatible runtime performance as other libraries written directly in C/C++.**

The mathematics and parameter tuning of MATMPC is explained at 
[matmpc-docs.readthedocs.io/en/latest/](https://matmpc-docs.readthedocs.io/en/latest/ "https://matmpc-docs.readthedocs.io/en/latest/")

### Windows:

Please install Matlab supported MinGW compiler at https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler

### Linux:

Install gcc by running "sudo apt-get install gcc"

**Note:** MATLAB is incompatible with the latest releases of gcc, It is recommended to download the gcc version that is compatible with your MATLAB version.  Go check [Compilers - MATLAB & Simulink (mathworks.com)](https://www.mathworks.com/support/requirements/supported-compilers.html).

### MacOS:

Install Xcode from app store

### To use MATMPC, follow the steps below.

1. Download and install CasADi-3.3.0. MATMPC currently is compatible with CasADi 3.3.0.

2. Write your own model using following the styles given by examples, e.g. Inverted Pendulum.

3. In the model file you created, set your own sampling time and multiple shooting time. 

4. Run Model_Generation.m

5. In InitData.m, set your own initialization data, e.g. data path, initial states, references and etc.

6. In Simulation.m, choose your integrator, set the prediction horizon and the solver options. You may also need to modify the reference type for constant or time-varying reference tracking problems.

7. In Draw.m, write your own plot functions to display your results.

8. Run Simulation.m and see the results!

### To use HPIPM as the QP solver
please read /doc/HPIPM-tutorial for the detailed installation process
