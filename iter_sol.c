#include "direct_sol.h"
#include "gsl.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include "mna.h"
#include "structs.h"
#include "parse.h"
#include "iter_sol.h"

// Algorithm for the CG algorithm, which is used to solve SPD systems iteratively

void solve_cg(gsl_vector* cur_gsl_b, gsl_vector* cur_gsl_x) {

    int max_iter = A_dim;   // A_dim = n

    gsl_vector *gsl_r = gsl_vector_alloc(A_dim);    // residual vector
    gsl_vector *gsl_Ax = gsl_vector_alloc(A_dim);   // vector Ax
    gsl_vector *gsl_z = gsl_vector_alloc(A_dim);    // vector z
    gsl_vector *gsl_p = gsl_vector_alloc(A_dim);    // vector p
    gsl_vector *gsl_q = gsl_vector_alloc(A_dim);    // vector q

    gsl_blas_dgemv(CblasNoTrans, 1.0, gsl_A, cur_gsl_x, 0.0, gsl_Ax);   // Ax = A*x

    gsl_vector_memcpy(gsl_r, cur_gsl_b);
    gsl_vector_sub(gsl_r, gsl_Ax);      // r = b - Ax

    double b_norm = gsl_blas_dnrm2(cur_gsl_b);  // norm of vector b
    double r_norm ;                         // norm of vector r
    double rho;
    double rho1;
    double beta;
    double pq_dot;
    double alpha;

    if (b_norm == 0) b_norm = 1;            // If the norm of b is 0, make it 1

    for (int i = 0; i < max_iter; i++) {

        r_norm = gsl_blas_dnrm2(gsl_r);
        if (r_norm/b_norm <= itol) break;

        gsl_z = preconditioner_solve(gsl_r, gsl_z);    // Solve Mz = r

        gsl_blas_ddot(gsl_r, gsl_z, &rho);   // dot(r,z)

        if (i == 0) {
            gsl_vector_memcpy(gsl_p, gsl_z);    // p = z
        }
        else {
            beta = rho/rho1;
            gsl_vector_axpby(1, gsl_z, beta, gsl_p);    // p = z + beta*p

        }

        rho1 = rho;
        gsl_blas_dgemv(CblasNoTrans, 1.0, gsl_A, gsl_p, 0.0, gsl_q);   // q = A*p

        gsl_blas_ddot(gsl_p, gsl_q, &pq_dot);   // dot(p,q)

        alpha = rho/pq_dot; // alpha = rho/(dot(p,q))

        gsl_blas_daxpy(alpha, gsl_p, cur_gsl_x);        // x = x + alpha*p
        gsl_blas_daxpy(-alpha, gsl_q, gsl_r);           // r = r - alpha*q

    }

    gsl_vector_free(gsl_r);
    gsl_vector_free(gsl_Ax);
    gsl_vector_free(gsl_z);
    gsl_vector_free(gsl_p);
    gsl_vector_free(gsl_q);
    

}

void solve_bicg(gsl_vector* cur_gsl_b, gsl_vector* cur_gsl_x) {

    int max_iter = A_dim;   // A_dim = n

    gsl_vector *gsl_r = gsl_vector_alloc(A_dim);        // residual vector
    gsl_vector *gsl_r_bicg = gsl_vector_alloc(A_dim);   // bicg residual vector

    gsl_vector *gsl_Ax = gsl_vector_alloc(A_dim);       // vector Ax
    gsl_vector *gsl_z = gsl_vector_alloc(A_dim);        // vector z
    gsl_vector *gsl_z_bicg = gsl_vector_alloc(A_dim);   // vector z'
    gsl_vector *gsl_p = gsl_vector_alloc(A_dim);        // vector p
    gsl_vector *gsl_p_bicg = gsl_vector_alloc(A_dim);   // vector p'
    gsl_vector *gsl_q = gsl_vector_alloc(A_dim);        // vector q
    gsl_vector *gsl_q_bicg = gsl_vector_alloc(A_dim);   // vector q'

    gsl_blas_dgemv(CblasNoTrans, 1.0, gsl_A, cur_gsl_x, 0.0, gsl_Ax);   // Ax = A*x

    gsl_vector_memcpy(gsl_r, cur_gsl_b);
    gsl_vector_sub(gsl_r, gsl_Ax);              // r = b - Ax

    gsl_vector_memcpy(gsl_r_bicg, cur_gsl_b);
    gsl_vector_sub(gsl_r_bicg, gsl_Ax);         // r' = b - Ax

    double b_norm = gsl_blas_dnrm2(cur_gsl_b);  // norm of vector b
    double r_norm;                         // norm of vector r
    double rho;
    double rho1;
    double beta;
    double omega;
    double alpha;

    if (b_norm == 0) b_norm = 1;            // If the norm of b is 0, make it 1

    for (int i = 0; i < 2*max_iter; i++) {

        r_norm = gsl_blas_dnrm2(gsl_r);
        if (r_norm/b_norm <= itol) break;

        gsl_z = preconditioner_solve(gsl_r, gsl_z);    // Solve Mz = r
        gsl_z_bicg = preconditioner_solve(gsl_r_bicg, gsl_z_bicg);    // Solve Mz' = r'

        gsl_blas_ddot(gsl_r_bicg, gsl_z, &rho);   // dot(r',z)

        if (fabs(rho) < EPS) break;

        if (i == 0) {
            gsl_vector_memcpy(gsl_p, gsl_z);    // p = z
            gsl_vector_memcpy(gsl_p_bicg, gsl_z_bicg);    // p' = z'
        }
        else {
            beta = rho/rho1;
            gsl_vector_axpby(1, gsl_z, beta, gsl_p);    // p = z + beta*p
            gsl_vector_axpby(1, gsl_z_bicg, beta, gsl_p_bicg);    // p' = z' + beta*p'

        }

        rho1 = rho;
        gsl_blas_dgemv(CblasNoTrans, 1.0, gsl_A, gsl_p, 0.0, gsl_q);   // q = A*p
        gsl_blas_dgemv(CblasTrans, 1.0, gsl_A, gsl_p_bicg, 0.0, gsl_q_bicg);   // q' = A(T)*p'

        gsl_blas_ddot(gsl_p_bicg, gsl_q, &omega);   // omega = dot(p',q)

        if (fabs(omega) < EPS) break;

        alpha = rho/omega;
        
        gsl_blas_daxpy(alpha, gsl_p, cur_gsl_x);        // x = x + alpha*p
        gsl_blas_daxpy(-alpha, gsl_q, gsl_r);           // r = r - alpha*q
        gsl_blas_daxpy(-alpha, gsl_q_bicg, gsl_r_bicg); // r' = r' - alpha*q'

    }

    gsl_vector_free(gsl_r);
    gsl_vector_free(gsl_r_bicg);
    gsl_vector_free(gsl_Ax);
    gsl_vector_free(gsl_z);
    gsl_vector_free(gsl_z_bicg);
    gsl_vector_free(gsl_p);
    gsl_vector_free(gsl_p_bicg);
    gsl_vector_free(gsl_q);
    gsl_vector_free(gsl_q_bicg);

}

// This function executes the preconditioning process, filling vector gsl_z
// Element i of gsl_z[i] = r[i] * 1/A[i][i]
gsl_vector* preconditioner_solve(gsl_vector *gsl_r, gsl_vector *gsl_z) {

    //gsl_vector *gsl_z = gsl_vector_alloc(A_dim);
    double cur_A_element;
    double cur_z_element;

    for (int i=0; i<A_dim; i++) {

        cur_A_element = gsl_matrix_get(gsl_A, i, i);
        if (cur_A_element == 0) cur_A_element = 1;

        cur_z_element = gsl_vector_get(gsl_r, i) / cur_A_element;

        gsl_vector_set(gsl_z, i, cur_z_element);
    }

    return gsl_z;
}