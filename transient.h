#ifndef MATHFUNCS_H
#define MATHFUNCS_H

#include "structs.h"

double get_pwl_val(pwl_spec *data, double t);
double get_exp_val(exp_spec *data, double t);
double get_sin_val(sin_spec *data, double t);
double get_pulse_val(pulse_spec *data, double t);

double get_comp_transient_val(component *curr, double t);

#endif