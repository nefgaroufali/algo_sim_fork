#ifndef EQUATIONS_H
#define EQUATIONS_H

extern double** A_array; // The square matrix A we want to fill
extern double* b_array; // The array b we want to fill
extern int count; // Variable to count how many elements V and L we have inserted

double** alloc_A_array();
double* alloc_b_array();
int create_arrays();
void print_arrays();
void create_equations();
void free_A_array();
void free_b_array();

#endif
