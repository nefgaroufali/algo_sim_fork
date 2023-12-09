#include "direct_sol.h"
#include "gsl.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "mna.h"
#include "structs.h"
#include "parse.h"



// This function forms the LU decomposition. 
// The array A is copied to a new array LU, which will have its elements changed to match the LU decomp.
// A gsl permutation is used

void form_LU() {

    gsl_p = gsl_permutation_alloc(A_dim);
    int signnum;

    gsl_LU = gsl_matrix_alloc(A_dim, A_dim);
    gsl_matrix_memcpy(gsl_LU, gsl_A);

    gsl_linalg_LU_decomp(gsl_LU, gsl_p, &signnum);

}

// This function forms the cholesky decomposition
// IMPORTANT: The array that is produced, contains A's elements in the upper triangular section.
// So, only the lower triangular part (including the diagonal) contains elements of the new decomp.

void form_chol() {

    gsl_set_error_handler(&gslErrorHandler);    // Handle error for non-spd matrices

    int chol_success;

    gsl_chol = gsl_matrix_alloc(A_dim, A_dim);
    gsl_matrix_memcpy(gsl_chol, gsl_A);

    chol_success = gsl_linalg_cholesky_decomp1(gsl_chol);

    if (chol_success != GSL_EDOM) {
        spd = TRUE;
    }

}

// This function solves the  DC system, with LU or Cholesky, depending on whether the user has specified SPD circuit
// Then writes the dc operating point in a .op file

void solve_dc_system(int solver_type) {

    // Allocate the x vector
    gsl_x = gsl_vector_alloc(A_dim);

    // Solve the system depending on the solver type
    if (solver_type == LU_SOL) {
        gsl_linalg_LU_solve(gsl_LU, gsl_p, gsl_b, gsl_x);
    }
    else {
        gsl_linalg_cholesky_solve(gsl_chol, gsl_b, gsl_x);
    }

    printf("Vector x:\n");
    print_gsl_vector(gsl_x, A_dim);
    printf("\n");

    // Open the file in which the operating point will be printed
    FILE *op_file;
    op_file = fopen("dc_solution.op", "w");

    // Print the solution to the file
    // First n-1 elements are the voltages of the nodes
    // Last m2 elements are the currents of the V and L components

    for (int i = 0, j=0; i < A_dim; i++) {

        if(i < nodes_n-1){
            fprintf(op_file, "v(%s)   %.5lf\n", node_array[i+1], gsl_vector_get(gsl_x, i));    
        }
        else
        {
            fprintf(op_file, "i(%s)   %.5lf\n", m2_array[j], gsl_vector_get(gsl_x, i));
            j++;
        }

    }

    fclose(op_file);


}

// This function performs DC sweep: It solves the Ax =b system for all different values of a specific voltage or current source
// The name of this source is parsed earlier and its pointer is sweep_component
void dc_sweep() {

    double cur_value;
    double low = DC_arguments[0];
    double high = DC_arguments[1];
    double step = DC_arguments[2];

    int pos_sweep_node_i, neg_sweep_node_i, sweep_node_i;

    // Temp_gsl_b: Same as gsl_b, with only the sweep value changing

    gsl_vector *temp_gsl_b = gsl_vector_alloc(A_dim);
    gsl_vector_memcpy(temp_gsl_b, gsl_b);

    component *current = sweep_component;

    // If the component is a current source, find the index of both nodes in the hash table

    // Change the current value, adding the step in each iteration, then change the bvector
    // Then solve the linear system in each iteration, in function solve_dc_sweep_system

    // Only change the b vector if the node is not ground

    if (current->comp_type == 'i') {
        pos_sweep_node_i = find_hash_node(&node_hash_table, current->positive_node) - 1; // -1 because 0 is ground
        neg_sweep_node_i = find_hash_node(&node_hash_table, current->negative_node) - 1;

        cur_value = low;
        do {
            if (pos_sweep_node_i != -1) {
                gsl_vector_set(temp_gsl_b, pos_sweep_node_i, -cur_value);
            }
            if (neg_sweep_node_i != -1) {
                gsl_vector_set(temp_gsl_b, neg_sweep_node_i, cur_value);
            }

            solve_dc_sweep_system(temp_gsl_b, cur_value, 'i'); 
            cur_value = cur_value + step;
        } while (cur_value <= high);

    }

    // If the component is a voltage source, take its index in the m2 part of the b vector

    // Change the current value, adding the step in each iteration, then change the bvector
    // Then solve the linear system in each iteration, in function solve_dc_sweep_system

    if (current->comp_type == 'v') {
        sweep_node_i = nodes_n-1 +current->m2_i;

        cur_value = low;
        do {
            gsl_vector_set(temp_gsl_b, sweep_node_i, cur_value); 

            solve_dc_sweep_system(temp_gsl_b, cur_value, 'v'); 
            cur_value = cur_value + step;
        } while (cur_value <= high);
    }

    gsl_vector_free(temp_gsl_b);

}