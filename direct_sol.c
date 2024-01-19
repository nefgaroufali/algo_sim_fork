#include "direct_sol.h"
#include "gsl.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "mna.h"
#include "structs.h"
#include "parse.h"
#include "sparse_sol.h"

#include <time.h>


gsl_matrix* gsl_LU = NULL;
gsl_matrix* gsl_chol = NULL;
gsl_permutation *gsl_p = NULL;
css *css_S = NULL;
csn *csn_N = NULL;


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

// This function solves the LU system using sparse methods
void solve_sparse_lu(gsl_vector* cur_gsl_b, gsl_vector* cur_gsl_x){

    double* temp_b = malloc(A_dim*sizeof(double *));
    double* temp_x = malloc(A_dim*sizeof(double *));
    double *x = (double *) malloc(sizeof(double)* A_dim);

    // The function inputs are gsl vectors, so the b vector needs to be converted to double
    // so that it can be used by the cs functions
    gsl_to_double(cur_gsl_b, temp_b);

    // LU Solution

    css_S = cs_sqr(sparse_cc_A, 2, 0);
    csn_N = cs_lu(sparse_cc_A, css_S, 1);

    cs_ipvec(A_dim, csn_N->Pinv, temp_b, x);
    cs_lsolve(csn_N->L, x);
    cs_usolve(csn_N->U, x);
    cs_ipvec(A_dim, css_S->Q, x, temp_x);

    // The vector needs to be in gsl form, so it must be converted
    double_to_gsl(cur_gsl_x, temp_x);

    cs_sfree(css_S);
    cs_nfree(csn_N);
    free(x);
    free(temp_b);
    free(temp_x);
}

// This function solves the system with the Cholesky algorithm using sparse methods
void solve_sparse_chol(gsl_vector* cur_gsl_b, gsl_vector* cur_gsl_x){

    double* temp_b = malloc(A_dim*sizeof(double *));
    double* temp_x = malloc(A_dim*sizeof(double *));
    double * x= (double *) malloc(sizeof(double)*A_dim);

    // The function inputs are gsl vectors, so the b vector needs to be converted to double
    // so that it can be used by the cs functions
    gsl_to_double(cur_gsl_b, temp_b);

    css_S = cs_schol(sparse_cc_A, 1);
    csn_N = cs_chol(sparse_cc_A, css_S);

    cs_ipvec(A_dim, css_S->Pinv, temp_b, x);
    cs_lsolve(csn_N->L, x);
    cs_ltsolve(csn_N->L, x);
    cs_pvec(A_dim, css_S->Pinv, x, temp_x);

    // The vector needs to be in gsl form, so it must be converted
    double_to_gsl(cur_gsl_x, temp_x);

    cs_sfree(css_S);
    cs_nfree(csn_N);
    free(x);
    free(temp_b);
    free(temp_x);
}