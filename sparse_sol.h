#ifndef SPARSE_SOL_H
#define SPARSE_SOL_H

#include "structs.h"
#include "csparse.h"

extern cs *sparse_A;    // Sparse array A in triplet form
extern cs *sparse_C;    // Sparse array A in compressed column form
extern int sparse_m2_count; // m2 count
extern double *b_array_sparse;  // b vector for sparse methods
extern double *x_array_sparse;  // x vector for sparse methods

void form_sparse();
void sparse_fill_with_r(component *current, int *nonzero_counter);
void sparse_fill_with_v(component *current, int *nonzero_counter);
void sparse_fill_with_l(component *current, int *nonzero_counter);
void sparse_fill_with_i(component *current);
void create_sparse_vectors();
void add_to_sparse_A(int k, int i, int j, double x);
void print_sparse_arrays();
void print_sparse_x();

#endif
