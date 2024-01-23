#include <stdio.h>
#include <math.h>
#include "transient.h"
#include "gsl.h"

// get the value of the exp transient function at time t
double get_exp_val(exp_spec *data, double t) {
	double res;

    // negative time -> return zero
	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}
    // t < td1
	else if (t < data->td1) {
		res = data->i1;
	}
    // td1 < t < td2
	else if (t < data->td2) {
		res = data->i1 + (data->i2 - data->i1)*(1.0 - exp(-1*(t - data->td1)/data->tc1));
	}
    // t > td2
	else {
		res = data->i1 + (data->i2 - data->i1)*(exp(-1*(t- data->td2)/data->tc2) - exp(-1*(t - data->td1)/data->tc1));
	}
	return res;
}


// get the value of the sin tansient function at time t
double get_sin_val(sin_spec *data, double t) {
	double res;

    // negative time -> return zero
	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}
    // t < td 
	else if (t < data->td) {
		res = data->i1 + data->ia*sin(2*M_PI*(data->ph)/360.0);
	}
    // t > td
	else {
		res = data->i1 + data->ia * sin(2*M_PI * data->fr*(t - data->td) + 2*M_PI * data->ph/360.0)*exp(-1*(t - data->td)*data->df);
	}
	return res;
}


// get the value of the pulse transient function at time t
double get_pulse_val(pulse_spec *data, double t) {
	double res;
	double y2, y1;
	double x1;


    // negative time -> return zero
	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}
    // t < td
	else if (t < data->td) {
		res = data->i1;
	}
	else {
		t = t - (int)((t-data->td)/data->per)*data->per;

		if (t < data->td + data->tr) {
			// i1->i2 linearly
			/*x2 = data->td + data->tr;*/
			x1 = data->td;
			y2 = data->i2;
			y1 = data->i1;
			/*res = ((y2-y1)/(x2-x1)) * (t - x1) + y1;*/
			res = ((y2 - y1)/data->tr)* (t - x1) + y1;
		}
		else if (t < data->td + data->tr + data->pw) {
			res = data->i2;
		}
		else if (t < data->td + data->tr + data->pw + data->tf) {
			// i2->i1 linearly
			/*x2 = data->td + data->tr + data->pw + data->tf;*/
			x1 = data->td + data->tr + data->pw;
			y2 = data->i1;
			y1 = data->i2;
			res = ((y2 - y1)/data->tf) * (t - x1) + y1;
		}
		else {
			res = data->i1;
		}
	}

	return res;
}



// get the value of pwl transient function at time t
double get_pwl_val(pwl_spec *data, double t) {
	double res;
	double y2, y1;
	double x2, x1;
	int idx;

	// negative time -> return zero
	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}

	// only one pair is equivalent to a line parallel to x axis (or a constant DC value)
	if (data->pairs == 1)
		return data->i[0];

	// t is smaller than the first pair time (return first value)
	if (t < data->t[0])
		return data->i[0];

	// t is bigger that the last pair time (return the last value)
	if (t > data->t[data->pairs-1])
		return data->i[data->pairs-1];


	// find the very FIRST pair with a time bigger than t
	idx = 0;
	while (1) {
		idx++;
		if (data->t[idx] > t)
			break;
	}

	x2 = data->t[idx];
	x1 = data->t[idx-1];
	y2 = data->i[idx];
	y1 = data->i[idx-1];
	res = ((y2-y1)/(x2-x1))*(t - x1) + y1;

	return res;
}

double get_comp_transient_val(component *curr, double t) {

	if (curr->spec_type == EXP_SPEC) return get_exp_val(&curr->spec->exp, t);
	else if (curr->spec_type == PWL_SPEC) return get_pwl_val(&curr->spec->pwl, t);
	else if (curr->spec_type == PULSE_SPEC) return get_pulse_val(&curr->spec->pulse, t);
	else if (curr->spec_type == SIN_SPEC) return get_sin_val(&curr->spec->sin, t);

	return 0;
}