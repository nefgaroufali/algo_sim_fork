#ifndef DIRECT_SOL_H
#define DIRECT_SOL_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>


#define LU_SOL 0
#define CHOL_SOL 1

extern gsl_matrix* gsl_A;
extern gsl_matrix* gsl_LU;
extern gsl_matrix* gsl_chol;
extern gsl_vector* gsl_b;
extern gsl_vector* gsl_x;
extern gsl_permutation *gsl_p;


void form_gsl_system();
void print_gsl_matrix(gsl_matrix *matrix, int dim);
void print_gsl_vector(gsl_vector *vector, int dim);
void free_gsl();
void form_LU();
void form_chol();
void solve_dc_system(int solver_type);

void gslErrorHandler(const char *reason, const char *file, int line, int gsl_errno);

#endif