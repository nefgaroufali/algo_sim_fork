#include "direct_sol.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "mna.h"
#include "structs.h"
#include "parse.h"

gsl_matrix* gsl_A = NULL;
gsl_matrix* gsl_LU = NULL;
gsl_matrix* gsl_chol = NULL;
gsl_vector *gsl_b = NULL;
gsl_vector *gsl_x = NULL;
gsl_permutation *gsl_p = NULL;

int spd = FALSE;

// This function copies the double array A made in mna.c, to a GSL array, to be used in the linear algebra routines
// It also copies the elements of the b vector

void form_gsl_system() {

    int i,j;
    i = j = A_dim;  // A_dim = n

    gsl_A = gsl_matrix_alloc(A_dim, A_dim);
    gsl_b = gsl_vector_alloc(A_dim);

    // Copying each element of A to the respective position of the GSL array

    for (i=0; i<A_dim; i++) {
        for (j=0; j<A_dim; j++) {
            gsl_matrix_set(gsl_A, i, j, A_array[i*A_dim + j]);
        }
    }

    // Copying each element of b to the respective position of the GSL Vector

    for (i=0; i<A_dim; i++) {
        gsl_vector_set(gsl_b, i, b_array[i]);
    }

    // print_gsl_matrix(gsl_A, A_dim);
    // printf("\n");
    // print_gsl_vector(gsl_b, A_dim);

    

}

// This function prints the specified gsl matrix. MUST be square
void print_gsl_matrix(gsl_matrix *matrix, int dim){

    int i,j;

    printf("\n");
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            printf ("%.2lf\t", gsl_matrix_get (matrix, i, j));
        }
        printf("\n");
    }

}

// This function prints the specified gsl vector

void print_gsl_vector(gsl_vector *vector, int dim){

    int i;

    printf("\n");
    for (i = 0; i < dim; i++) {
        printf ("%.2lf ", gsl_vector_get(vector, i));
    }
    printf("\n");

}

// This function frees all memory occupied by gsl arrays

void free_gsl() {

     gsl_matrix_free(gsl_A);
     gsl_vector_free(gsl_b);
     gsl_matrix_free(gsl_LU);
     gsl_matrix_free(gsl_chol);
}

// This function forms the LU decomposition. 
// The array A is copied to a new array LU, which will have its elements changed to match the LU decomp.
// A gsl permutation is used

void form_LU() {

    gsl_p = gsl_permutation_alloc(A_dim);
    int signnum;

    gsl_LU = gsl_matrix_alloc(A_dim, A_dim);
    gsl_matrix_memcpy(gsl_LU, gsl_A);

    gsl_linalg_LU_decomp(gsl_LU, gsl_p, &signnum);


    // print_gsl_matrix(gsl_LU, A_dim);


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
    // else {
    //     fprintf(stderr, "Cholesky decomposition failed: Matrix is not positive definite\n");
    // }


}

// Custom error handler function
void gslErrorHandler(const char *reason, const char *file, int line, int gsl_errno) {

    // Print an error message
    fprintf(stderr, "GSL Error: %s \n", reason);
    
    // Handle the error as needed
    // For example, you could set a flag or take corrective action
    
    // In this example, we exit the program, but you might choose to do something else
}

// This function solves the  DC system, with LU or Cholesky, depending on whether the user has specified SPD circuit
// Then writes the dc operating point in a .op file

void solve_dc_system(int solver_type) {

    gsl_x = gsl_vector_alloc(A_dim);

    if (solver_type == LU_SOL) {
        gsl_linalg_LU_solve(gsl_LU, gsl_p, gsl_b, gsl_x);
    }
    else {
        gsl_linalg_cholesky_solve(gsl_chol, gsl_b, gsl_x);
    }

    printf("Vector x:\n");
    print_gsl_vector(gsl_x, A_dim);
    printf("\n");

    FILE *op_file;

    op_file = fopen("dc_solution.op", "w");

    // ONLY for the DC system!
    for (int i = 0; i < A_dim; i++) {
        fprintf(op_file, "%.5lf\n",gsl_vector_get(gsl_x, i));
    }

    fclose(op_file);


}

void dc_sweep() {

    int i;
    int low = DC_arguments[0];
    int high = DC_arguments[1];
    int step = DC_arguments[2];

    for (i=low; i<=high; i=i+step) {

    }

}