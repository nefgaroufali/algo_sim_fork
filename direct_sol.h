#ifndef DIRECT_SOL_H
#define DIRECT_SOL_H

#define LU_SOL 0
#define CHOL_SOL 1

void form_LU();
void form_chol();
void solve_dc_system(int solver_type);
void dc_sweep();

#endif