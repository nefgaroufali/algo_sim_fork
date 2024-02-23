#ifndef AC_H
#define AC_H

//#define M_PI 3.14159
#include <complex.h>
#include "structs.h"
#include "gsl.h"

extern double complex* ac_b_vector;
extern gsl_vector_complex *gsl_ac_b_vector;

void ac_sweep();
void fill_ac_A_array_temp(double complex* temp_ac_A_array, double omega);

void fill_ac_A_with_r(double complex* temp_ac_A_array, component* current);
void fill_ac_A_with_i(double complex* temp_ac_A_array, component* current);
void fill_ac_A_with_c(double complex* temp_ac_A_array, component* current, double omega);
void fill_ac_A_with_v(double complex* temp_ac_A_array, component* current, int *m2_count_ac);
void fill_ac_A_with_l(double complex* temp_ac_A_array, component* current, int *m2_count_ac, double omega);

void print_ac_A_array(double complex* temp_ac_A_array);

void solve_ac_sweep_system(gsl_matrix_complex *temp_gsl_ac_A_array, double freq)
;

#endif