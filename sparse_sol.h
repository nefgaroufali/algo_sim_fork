#ifndef SPARSE_SOL_H
#define SPARSE_SOL_H

#include "structs.h"
#include "csparse.h"

extern cs *sparse_Α;
extern cs *sparse_C;
extern int sparse_m2_count;

void form_sparse();
void sparse_fill_with_r(component *current, int *nonzero_counter);
void sparse_fill_with_v(component *current, int *nonzero_counter);
void sparse_fill_with_l(component *current, int *nonzero_counter);
void add_to_sparse_A(int k, int i, int j, double x);

#endif
