#ifndef DIRECT_SOL_H
#define DIRECT_SOL_H

#include <gsl/gsl_matrix.h>

extern gsl_matrix* gsl_A;
extern gsl_matrix* gsl_LU;
extern gsl_matrix* gsl_chol;


void copy_to_A();
void print_gsl_matrix(gsl_matrix *matrix, int dim);
void free_gsl();
void form_LU();
void form_chol();

void gslErrorHandler(const char *reason, const char *file, int line, int gsl_errno);

#endif