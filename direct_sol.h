#ifndef DIRECT_SOL_H
#define DIRECT_SOL_H

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

extern gsl_matrix* gsl_LU;
extern gsl_matrix* gsl_chol;
extern gsl_permutation *gsl_p;

void form_LU();
void form_chol();

#endif