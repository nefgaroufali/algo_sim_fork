#ifndef EQUATIONS_H
#define EQUATIONS_H

extern double** A_array;
extern double* b_array;
extern int count;

double** alloc_A_array();
double* alloc_b_array();
int create_arrays();
void print_arrays();
void create_equations();

#endif
