#ifndef MPC_COMMON_H_
#define MPC_COMMON_H_

#include "stdlib.h"

extern void daxpy_(int*,double *,double*,int*,double*,int*);//(&nx, &a, K_lambda+(s+1)*nw, &one_i, lambda_t, &one_i);
extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern double dlange_(char*, int*, int*, double*, int*, double*);//dlange(Norm, &neq, &one_i, eq_res_vec, &neq, work)
extern void dgemv_(char*, int*, int*, double*,double*,int*,double*,int*,double*,double*,int*);//(nTrans,&nx,&nx,&one_d,Q+i*nx*nx,&nx,dx+i*nx,&one_i,&zero,tmp,&one_i);
extern double ddot_(int*,double*,int*,double*,int*);//(&nx, dx+i*nx, &one_i, tmp, &one_i);
extern double dnrm2_(int*,double*,int*);//(&nw, L, &one_i);
#define daxpy daxpy_
#define dgemm dgemm_
#define dlange dlange_
#define dgemv dgemv_
#define ddot ddot_
#define dnrm2 dnrm2_

void Block_Fill(size_t m, size_t n, double *Gi, double *G,
     size_t idm, size_t idn, size_t ldG);

void Block_Fill_Trans(size_t m, size_t n, double *Gi, double *G,
     size_t idm, size_t idn, size_t ldG);

void Block_Get(size_t m, size_t n, double *Gi, double *G, size_t idm, size_t idn, size_t ldG);

void set_zeros(size_t dim, double *A);

void print_matrix(double *A, size_t m, size_t n);

void print_vector(double *x, size_t m);

void regularization(size_t n, double *A, double reg);

// void product_mb(double *P,double *M, size_t a, double *index, size_t r,size_t nu);    // P = M*T,   M is (a x (N*nu))
// 
// void product_mb_Trans(double *P,double *M, size_t b, double *index, size_t r,size_t nu);     // P = T'* M, M is (N*nu x b) 

#endif