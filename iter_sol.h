#ifndef ITER_SOL_H
#define ITER_SOL_H

#define EPS 1e-14

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

void solve_cg(gsl_vector* cur_gsl_b, gsl_vector *cur_gsl_x);
void solve_bicg(gsl_vector* cur_gsl_b, gsl_vector* cur_gsl_x);
gsl_vector* preconditioner_solve(gsl_vector *gsl_r, gsl_vector *gsl_z);
double vector_inner_product(int m_row, double *first, double *second);
void daxpy(const double alpha, const double *x, double *y, const int size);
double euclidean_norm(const double *vector, const int vector_size);
void solve_sparse_cg();
void solve_sparse_bicg();

#endif
