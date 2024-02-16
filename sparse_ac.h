#ifndef SPARSE_AC_H
#define SPARSE_AC_H

#include "gsl.h"
#include "../CXSparse-master/Include/cs.h"
#include <complex.h>

extern double complex* sparse_ac_b_vector;
extern gsl_vector_complex *sparse_gsl_ac_b_vector;

extern cs_ci *sparse_ac_A;    // Sparse array A in triplet form
extern cs_ci *sparse_ac_cc_A;    // Sparse array A in compressed column form

void ac_sweep_sparse();
void solve_ac_sweep_system_sparse(cs_ci *sparse_ac_cc_A);

cs_ci* form_ac_sparse_A(double omega);

void solve_sparse_lu_complex(cs_ci* A, double complex* b, double complex* cur_x);
void solve_sparse_bicg_complex(cs_ci* A, double complex* b, double complex* cur_x);

void sparse_ac_fill_with_r(cs_ci *sparse_ac_A, component *current, int *nonzero_counter);
void sparse_ac_fill_with_v(cs_ci *sparse_ac_A, component *current, int *nonzero_counter, int* sparse_ac_m2_count);
void sparse_ac_fill_with_c(cs_ci *sparse_ac_A, component *current, int *nonzero_counter, double omega);
void sparse_ac_fill_with_l(cs_ci *sparse_ac_A, component* current, int *nonzero_counter, int* sparse_ac_m2_count, double omega);
void sparse_ac_fill_with_i(cs_ci *sparse_ac_A, component *current);

int cs_gaxpy_complex_with_gsl_x (const cs_ci *A, gsl_vector_complex* cur_gsl_x, gsl_vector_complex *gsl_y);
int cs_gaxpy_complex_with_gsl_x_conjtrans (const cs_ci *A, gsl_vector_complex* cur_gsl_x, gsl_vector_complex *gsl_y);
double complex* find_diag_a_complex(const cs_ci *sparse_C);
double complex* find_diag_a_conj_complex(const cs_ci *sparse_C); 
gsl_vector_complex* preconditioner_solve_sparse_complex(gsl_vector_complex *gsl_r, gsl_vector_complex *gsl_z, const double complex *diag_a);

void add_to_sparse_ac_A(cs_ci *sparse_ac_A, int k, int i, int j, double complex x);

#endif