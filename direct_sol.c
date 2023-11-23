#include "direct_sol.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "mna.h"

gsl_matrix* gsl_A = NULL;
gsl_matrix* gsl_LU = NULL;
gsl_matrix* gsl_chol = NULL;

// This function copies the double array A made in mna.c, to a GSL array, to be used in the linear algebra routines

void copy_to_A() {

    int i,j;
    i = j = A_dim;  // A_dim = n

    gsl_A = gsl_matrix_alloc(A_dim, A_dim);

    // Copying each element of A to the respective position of the GSL array

    for (i=0; i<A_dim; i++) {
        for (j=0; j<A_dim; j++) {
            gsl_matrix_set(gsl_A, i, j, A_array[i*A_dim + j]);
        }
    }

    print_gsl_matrix(gsl_A, A_dim);

    

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

// This function frees all memory occupied by gsl arrays

void free_gsl() {

     gsl_matrix_free(gsl_A);
     gsl_matrix_free(gsl_LU);
}

// This function forms the LU decomposition. 
// The array A is copied to a new array LU, which will have its elements changed to match the LU decomp.
// A gsl permutation is used

void form_LU() {

    gsl_permutation *permutation;
    permutation = gsl_permutation_alloc(A_dim);
    int signnum;

    gsl_LU = gsl_matrix_alloc(A_dim, A_dim);
    gsl_matrix_memcpy(gsl_LU, gsl_A);

    gsl_linalg_LU_decomp(gsl_LU, permutation, &signnum);


    print_gsl_matrix(gsl_LU, A_dim);

    for (int i=0; i<A_dim; i++) {
        printf("%ld\n", gsl_permutation_get(permutation, i));
    }

    gsl_permutation_free(permutation);



}

void form_chol() {

    gsl_set_error_handler(&gslErrorHandler);

    int chol_success;

    gsl_chol = gsl_matrix_alloc(A_dim, A_dim);
    gsl_matrix_memcpy(gsl_chol, gsl_A);

    chol_success = gsl_linalg_cholesky_decomp1(gsl_chol);

    if (chol_success == GSL_EDOM) {
        //fprintf(stderr, "Cholesky decomposition failed: Matrix is not positive definite\n");
    }
    else {
        print_gsl_matrix(gsl_chol, A_dim);
    }


}

// Custom error handler function
void gslErrorHandler(const char *reason, const char *file, int line, int gsl_errno) {

    // Print an error message
    fprintf(stderr, "GSL Error: %s \n", reason);
    
    // Handle the error as needed
    // For example, you could set a flag or take corrective action
    
    // In this example, we exit the program, but you might choose to do something else
}