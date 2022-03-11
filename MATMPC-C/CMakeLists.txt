cmake_minimum_required(VERSION 3.16)
project(matmpc_c C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS -pthread)

add_executable(matmpc_c src/main.c src/main_utils.c src/mpc_fcn.c src/qp_generation.c
        src/sim.c src/erk.c src/casadi_wrapper.c src/mpc_common.c src/casadi_src.c src/hpipm_sparse.c src/line_search.c src/solution_info.c src/simulate_system.c)

include_directories(include /opt/blasfeo/include /opt/hpipm/include)

link_directories(${CMAKE_SOURCE_DIR}/../lib)
target_link_libraries(matmpc_c m libopenblas.a /opt/hpipm/lib/libhpipm.a /opt/blasfeo/lib/libblasfeo.a)
