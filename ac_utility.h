#ifndef AC_UTILITY_H
#define AC_UTILITY_H

#include "gsl.h"
#include <complex.h>
#include "structs.h"
#include "parse.h"
#include <math.h>

void fill_ac_b_vector(double complex* ac_b_vector);
double complex* alloc_ac_b_vector();
void print_ac_vector(double complex* ac_vector);

void fill_ac_bvector_with_r(double complex* ac_b_vector, component* current);
void fill_ac_bvector_with_i(double complex* ac_b_vector, component* current);
void fill_ac_bvector_with_c(double complex* ac_b_vector, component* current);
void fill_ac_bvector_with_v(double complex* ac_b_vector, component* current, int *m2_count_ac);
void fill_ac_bvector_with_l(double complex* ac_b_vector, component* current, int *m2_count_ac);

double complex convert_to_complex(double mag, double phase);
double get_magnitude(double complex complex_num, int ac_sweep_method);
double get_phase(double complex complex_num);

void create_ac_gnuplot(int i, int ac_sweep_method);



#endif