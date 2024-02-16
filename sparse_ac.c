#include "ac.h"
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parse.h"
#include "direct_sol.h"
#include "iter_sol.h"
#include "sparse_ac.h"
#include "ac_utility.h"

double complex* sparse_ac_b_vector = NULL;
gsl_vector_complex *sparse_gsl_ac_b_vector;

void ac_sweep_sparse() {

    double omega;
    
    cs_ci *sparse_ac_cc_A = NULL;

    sparse_ac_b_vector = alloc_ac_b_vector();
    fill_ac_b_vector(sparse_ac_b_vector);

    double iter_num = 10;
    for (double freq=0; freq<iter_num; freq++) {

        omega=2*M_PI*freq;

        // fill AC array in double complex form
        sparse_ac_cc_A = form_ac_sparse_A(omega);

        // then copy in GSL form

        // for (int i=0; i<A_dim; i++) {
        //     for (int j=0; j<A_dim; j++) {
        //         gsl_matrix_complex_set(temp_gsl_ac_A_array, i, j, temp_ac_A_array[i*A_dim + j]);
        //     }
        // }
        
        solve_ac_sweep_system_sparse(sparse_ac_cc_A);

        //cs_ci_print(sparse_ac_cc_A, 0);

        cs_ci_spfree(sparse_ac_cc_A);
        sparse_ac_cc_A = NULL;
        // print_gsl_matrix_complex(temp_gsl_ac_A_array, A_dim);
        // memset(temp_ac_A_array, 0, sizeof(double complex) *A_dim*A_dim);




    }
}

cs_ci* form_ac_sparse_A(double omega) { 

    component* current = head;
    int nonzero_counter = 0;
    int sparse_m2_count_ac = 0;
    int nonzeros;

    // Since the c and l components are muitiplied with omega, we must account for than in the nonzero counting
    if (omega == 0) {
        nonzeros = nonzeros_A;
    }
    else {
        nonzeros = nonzeros_A + nonzeros_C;
    }

    // Allocate the sparse complex A array
    cs_ci* sparse_ac_A = cs_ci_spalloc(A_dim, A_dim, nonzeros,1,1);

    // For each component of the circuit fill the sparse complex array A

    while (current != NULL) {

        if(current->comp_type == 'r'){
            sparse_ac_fill_with_r(sparse_ac_A, current, &nonzero_counter);
        }
        else if(current->comp_type == 'v'){
            sparse_ac_fill_with_v(sparse_ac_A, current, &nonzero_counter, &sparse_m2_count_ac);
        }
        else if(current->comp_type == 'l'){
            sparse_ac_fill_with_l(sparse_ac_A, current, &nonzero_counter, &sparse_m2_count_ac, omega);
        }
        else if(current->comp_type == 'c'){
            sparse_ac_fill_with_c(sparse_ac_A, current, &nonzero_counter, omega);
        }
        else if(current->comp_type == 'i'){
            sparse_ac_fill_with_i(sparse_ac_A, current);
        }

        current = current->next;
    }

    sparse_ac_A->nz = nonzero_counter;

    // Form the compressed-column arrays cc_A and cc_C
    cs_ci* sparse_ac_cc_A = cs_ci_compress(sparse_ac_A);

    cs_ci_spfree(sparse_ac_A);

    // Remove the duplicates of the compressed-column arrays
    cs_ci_dupl(sparse_ac_cc_A);

    return sparse_ac_cc_A;


}

// Component: R ---> Fill A1*G*A1t array
void sparse_ac_fill_with_r(cs_ci* sparse_ac_A, component* current, int *nonzero_counter){

    int positive_node_i;
    int negative_node_i;
    int k = *nonzero_counter;

    double current_g;

    if (current->value == 0) current_g = 0;
    else current_g =  1/current->value; // g = 1/R

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: r 0 x value
    if (positive_node_i == 0) {

        add_to_sparse_ac_A(sparse_ac_A, k, negative_node_i-1, negative_node_i-1, current_g);
        *nonzero_counter += 1;
    }

    // Case: r x 0 value
    else if (negative_node_i == 0) {

        add_to_sparse_ac_A(sparse_ac_A, k, positive_node_i-1, positive_node_i-1, current_g);
        *nonzero_counter += 1;
    }

    // Case: x y value
    else {
        add_to_sparse_ac_A(sparse_ac_A, k, negative_node_i-1, negative_node_i-1, current_g);
        add_to_sparse_ac_A(sparse_ac_A, k+1, positive_node_i-1, positive_node_i-1, current_g);
        add_to_sparse_ac_A(sparse_ac_A, k+2, negative_node_i-1, positive_node_i-1, -current_g);
        add_to_sparse_ac_A(sparse_ac_A, k+3, positive_node_i-1, negative_node_i-1, -current_g);
        *nonzero_counter += 4;
    }

    // return;
}

// Component: C ---> Fill C array (Transient)
void sparse_ac_fill_with_c(cs_ci* sparse_ac_A, component* current, int *nonzero_counter, double omega){

    int positive_node_i;
    int negative_node_i;
    int k = *nonzero_counter;

    double current_c = omega*current->value;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));


    if (omega != 0) {
        // Case: r 0 x value
        if (positive_node_i == 0) {

            add_to_sparse_ac_A(sparse_ac_A, k, negative_node_i-1, negative_node_i-1, I*current_c);
            *nonzero_counter += 1;
        }

        // Case: r x 0 value
        else if (negative_node_i == 0) {

            add_to_sparse_ac_A(sparse_ac_A, k, positive_node_i-1, positive_node_i-1, I*current_c);
            *nonzero_counter += 1;
        }

        // Case: x y value
        else {
            add_to_sparse_ac_A(sparse_ac_A, k, negative_node_i-1, negative_node_i-1, I*current_c);
            add_to_sparse_ac_A(sparse_ac_A, k+1, positive_node_i-1, positive_node_i-1, I*current_c);
            add_to_sparse_ac_A(sparse_ac_A, k+2, negative_node_i-1, positive_node_i-1, -I*current_c);
            add_to_sparse_ac_A(sparse_ac_A, k+3, positive_node_i-1, negative_node_i-1, -I*current_c);
            *nonzero_counter += 4;
        }
    }

    // return;
}

// Component: I ---> Fill b array
void sparse_ac_fill_with_i(cs_ci* sparse_ac_A, component* current){

}

// Component: V ---> Fill A2t, A2 array and b array
void sparse_ac_fill_with_v(cs_ci* sparse_ac_A, component* current, int *nonzero_counter, int* sparse_ac_m2_count) {

    int positive_node_i;
    int negative_node_i;
    int sparse_m2_i = nodes_n - 1 + (*sparse_ac_m2_count);
    int k = *nonzero_counter;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: v 0 x value
    if(positive_node_i == 0){
        add_to_sparse_ac_A(sparse_ac_A, k, negative_node_i-1, sparse_m2_i, -1);
        add_to_sparse_ac_A(sparse_ac_A, k+1, sparse_m2_i, negative_node_i-1, -1);
        *nonzero_counter += 2;
    }
    // Case: v x 0 value
    else if(negative_node_i == 0){
        add_to_sparse_ac_A(sparse_ac_A, k, positive_node_i-1, sparse_m2_i, 1);
        add_to_sparse_ac_A(sparse_ac_A, k+1, sparse_m2_i, positive_node_i-1, 1);
        *nonzero_counter += 2;
    }
    // Case: v x y value
    else{
        add_to_sparse_ac_A(sparse_ac_A, k, negative_node_i-1, sparse_m2_i, -1);
        add_to_sparse_ac_A(sparse_ac_A, k+1, sparse_m2_i, negative_node_i-1, -1);
        add_to_sparse_ac_A(sparse_ac_A, k+2, positive_node_i-1, sparse_m2_i, 1);
        add_to_sparse_ac_A(sparse_ac_A, k+3, sparse_m2_i, positive_node_i-1, 1);
        *nonzero_counter += 4;
    }

    (*sparse_ac_m2_count)++;
}

// Component: L ---> Fill A2t, A2 array and b array, as well as C(transient) array
void sparse_ac_fill_with_l(cs_ci* sparse_ac_A, component* current, int *nonzero_counter, int* sparse_ac_m2_count, double omega) {

    int positive_node_i;
    int negative_node_i;
    int sparse_m2_i = nodes_n - 1 + (*sparse_ac_m2_count);
    int k = *nonzero_counter;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: v 0 x value
    if(positive_node_i == 0){
        add_to_sparse_ac_A(sparse_ac_A, k, negative_node_i-1, sparse_m2_i, -1);
        add_to_sparse_ac_A(sparse_ac_A, k+1, sparse_m2_i, negative_node_i-1, -1);
        *nonzero_counter += 2;
    }
    // Case: v x 0 value
    else if(negative_node_i == 0){
        add_to_sparse_ac_A(sparse_ac_A, k, positive_node_i-1, sparse_m2_i, 1);
        add_to_sparse_ac_A(sparse_ac_A, k+1, sparse_m2_i, positive_node_i-1, 1);
        *nonzero_counter += 2;
    }
    // Case: v x y value
    else{
        add_to_sparse_ac_A(sparse_ac_A, k, negative_node_i-1, sparse_m2_i, -1);
        add_to_sparse_ac_A(sparse_ac_A, k+1, sparse_m2_i, negative_node_i-1, -1);
        add_to_sparse_ac_A(sparse_ac_A, k+2, positive_node_i-1, sparse_m2_i, 1);
        add_to_sparse_ac_A(sparse_ac_A, k+3, sparse_m2_i, positive_node_i-1, 1);
        *nonzero_counter += 4;
    }

    // Must add the imaginary part!
      
    if (omega != 0) {
        add_to_sparse_ac_A(sparse_ac_A, (*nonzero_counter), sparse_m2_i, sparse_m2_i, -I*omega*current->value);
        *nonzero_counter += 1;
    }  

    (*sparse_ac_m2_count)++;
}

// This function adds an element to the sparse array A
void add_to_sparse_ac_A(cs_ci *sparse_ac_A, int k, int i, int j, double complex x) {
        sparse_ac_A->i[k] = i;
        sparse_ac_A->p[k] = j;
        sparse_ac_A->x[k] = x;
}

// This function solves the LU system using sparse methods
void solve_ac_sweep_system_sparse(cs_ci *sparse_ac_cc_A){

    double complex *temp_sparse_ac_x = (double complex*) calloc((A_dim), sizeof(double complex));

    //print_ac_vector(sparse_ac_b_vector);
    // Solve the system based on the solver flag, generated during parsing
    if (solver_type == SPARSE_LU_SOL) {
        solve_sparse_lu_complex(sparse_ac_cc_A, sparse_ac_b_vector, temp_sparse_ac_x);
        printf("x is\n");
        print_ac_vector(temp_sparse_ac_x);
    }


    else if (solver_type == SPARSE_BICG_SOL) {
        solve_sparse_bicg_complex(sparse_ac_cc_A, sparse_ac_b_vector, temp_sparse_ac_x);
        printf("x is\n");
        print_ac_vector(temp_sparse_ac_x);
    }
    // else if (solver_type == SPARSE_LU_SOL) {
    //     solve_sparse_lu_complex(temp_gsl_ac_A_array, gsl_ac_b_vector, temp_gsl_x);
    // }
    // else if (solver_type == SPARSE_BICG_SOL) {
    //     solve_sparse_bicg(sparse_cc_A, temp_gsl_b, temp_gsl_x);
    // }

    // int i;
    // int plot_node_i;
    // double b_vector_value, x_vector_value;

    // // Plot_node_i: The index of the node(s) that are to be plotted
    // // Sweep_node_i: The index of the node whose value in the b vector changes

    // // Get the values that will be printed to the files, then call add_to_plot_file
    // for (i = 0; i < plot_node_count; i++) {
    //     plot_node_i = plot_node_indexes[i];
    //     b_vector_value = cur_value;
    //     x_vector_value = gsl_vector_get(temp_gsl_x, plot_node_i);
    //     add_to_plot_file(b_vector_value, x_vector_value, i);

    // }

    // //gsl_vector_memcpy(gsl_x, temp_gsl_x); // copy the value to gsl_x, so that it is the initial value for the next cg/bicg
    // gsl_vector_free(temp_gsl_x);

}

// This function solves the LU system using sparse methods
void solve_sparse_lu_complex(cs_ci* A, double complex* b, double complex* cur_x){

    double complex* temp_b = (double complex*) malloc(A_dim*sizeof(double complex));
    double complex* temp_x = (double complex*) malloc(A_dim*sizeof(double complex));
    double complex* x = (double complex*) malloc(A_dim*sizeof(double complex));

    memcpy(temp_b, b, A_dim*sizeof(double complex));

    // LU Solution

    cs_cis *css_S = cs_ci_sqr(2, A, 0);
    cs_cin *csn_N = cs_ci_lu(A, css_S, 1); 

    cs_ci_ipvec(csn_N->pinv, temp_b, x, A_dim);
    cs_ci_lsolve(csn_N->L, x);
    cs_ci_usolve(csn_N->U, x);
    cs_ci_ipvec(css_S->q, x, temp_x, A_dim);

    memcpy(cur_x, temp_x, A_dim*sizeof(double complex));

    // printf("temp b and tempx are\n");
    // print_ac_vector(temp_b);
    // print_ac_vector(temp_x);

    cs_ci_sfree(css_S);
    cs_ci_nfree(csn_N);

    free(x);
    free(temp_b);
    free(temp_x);

}

// Sparse implementation of the BiCG algorithm, which is used to solve non-SPD systems iteratively
void solve_sparse_bicg_complex(cs_ci* A, double complex* b, double complex* cur_x) {

    int max_iter = A_dim;

    gsl_vector_complex *cur_gsl_x = gsl_vector_complex_alloc(A_dim);
    gsl_vector_complex *gsl_b = gsl_vector_complex_alloc(A_dim);

    // convert to gsl
    for (int i=0; i<A_dim; i++) {
        gsl_vector_complex_set(cur_gsl_x, i, cur_x[i]);
        gsl_vector_complex_set(gsl_b, i, b[i]);
    }

    gsl_vector_complex *gsl_r = gsl_vector_complex_alloc(A_dim);        // residual vector
    gsl_vector_complex *gsl_r_bicg = gsl_vector_complex_alloc(A_dim);   // bicg residual vector

    gsl_vector_complex *gsl_Ax = gsl_vector_complex_calloc(A_dim);       // vector Ax
    gsl_vector_complex *gsl_z = gsl_vector_complex_alloc(A_dim);        // vector z
    gsl_vector_complex *gsl_z_bicg = gsl_vector_complex_alloc(A_dim);   // vector z'
    gsl_vector_complex *gsl_p = gsl_vector_complex_alloc(A_dim);        // vector p
    gsl_vector_complex *gsl_p_bicg = gsl_vector_complex_alloc(A_dim);   // vector p'
    gsl_vector_complex *gsl_q = gsl_vector_complex_alloc(A_dim);        // vector q
    gsl_vector_complex *gsl_q_bicg = gsl_vector_complex_alloc(A_dim);   // vector q'

    gsl_complex complex_0 = gsl_complex_rect(0.0, 0.0);
    gsl_complex complex_1 = gsl_complex_rect(1.0, 0.0);

    cs_gaxpy_complex_with_gsl_x(A, cur_gsl_x, gsl_Ax);   // Ax = A*x

    gsl_vector_complex_memcpy(gsl_r, gsl_b);
    gsl_vector_complex_sub(gsl_r, gsl_Ax);              // r = b - Ax

    gsl_vector_complex_memcpy(gsl_r_bicg, gsl_b);
    gsl_vector_complex_sub(gsl_r_bicg, gsl_Ax);         // r' = b - Ax

    double b_norm = gsl_blas_dznrm2(gsl_b);  // norm of vector b
    double r_norm = gsl_blas_dznrm2(gsl_r);      // norm of vector r
    gsl_complex rho;
    gsl_complex rho1 = complex_1;
    gsl_complex beta;
    gsl_complex beta_conj;
    gsl_complex omega;
    gsl_complex alpha;
    gsl_complex alpha_conj;
    gsl_complex neg_alpha;
    gsl_complex neg_alpha_conj;
    int i = 0;

    double complex *diag_a = find_diag_a_complex(A);
    double complex *conj_diag_a = find_diag_a_conj_complex(A); //used for preconditioner! Also transforms into the conjugate

    if (b_norm == 0) b_norm = 1;            // If the norm of b is 0, make it 1

    while(r_norm/b_norm > itol && i<2*max_iter) {

        i++;

        gsl_z = preconditioner_solve_sparse_complex(gsl_r, gsl_z, diag_a);    // Solve Mz = r
        gsl_z_bicg = preconditioner_solve_sparse_complex(gsl_r_bicg, gsl_z_bicg, conj_diag_a);    // Solve M(H)z' = r'

        gsl_blas_zdotc(gsl_r_bicg, gsl_z, &rho);   // dot(r',z)

        if (gsl_complex_abs(rho) < EPS) break;

        if (i == 1) {
            gsl_vector_complex_memcpy(gsl_p, gsl_z);    // p = z
            gsl_vector_complex_memcpy(gsl_p_bicg, gsl_z_bicg);    // p' = z'
        }
        else {
            beta = gsl_complex_div(rho, rho1);
            beta_conj = gsl_complex_conjugate(beta);
            gsl_vector_complex_axpby(complex_1, gsl_z, beta, gsl_p);    // p = z + beta*p
            gsl_vector_complex_axpby(complex_1, gsl_z_bicg, beta_conj, gsl_p_bicg);    // p' = z' + beta*p'

        }

        rho1 = rho;

        gsl_vector_complex_set_zero(gsl_q); // initialize gsl_q to 0, necessary for the following axpy function

        cs_gaxpy_complex_with_gsl_x(A, gsl_p, gsl_q);  // q = A*p 
        cs_gaxpy_complex_with_gsl_x_conjtrans(A, gsl_p_bicg, gsl_q_bicg);  // q' = A(T)*p'

        gsl_blas_zdotc(gsl_p_bicg, gsl_q, &omega);   // omega = dot(p',q)

        if (gsl_complex_abs(omega) < EPS) break;

        alpha = gsl_complex_div(rho,omega);
        alpha_conj = gsl_complex_conjugate(alpha);
        neg_alpha = gsl_complex_negative(alpha);
        neg_alpha_conj = gsl_complex_negative(alpha_conj);
        
        gsl_blas_zaxpy(alpha, gsl_p, cur_gsl_x);        // x = x + alpha*p
        gsl_blas_zaxpy(neg_alpha, gsl_q, gsl_r);           // r = r - alpha*q
        gsl_blas_zaxpy(neg_alpha_conj, gsl_q_bicg, gsl_r_bicg); // r' = r' - alpha*q'

        r_norm = gsl_blas_dznrm2(gsl_r);

    }

    // convert back
    for (int i=0; i<A_dim; i++) {
        cur_x[i] = gsl_vector_complex_get(cur_gsl_x, i);
    }

    gsl_vector_complex_free(gsl_r);
    gsl_vector_complex_free(gsl_r_bicg);
    gsl_vector_complex_free(gsl_Ax);
    gsl_vector_complex_free(gsl_z);
    gsl_vector_complex_free(gsl_z_bicg);
    gsl_vector_complex_free(gsl_p);
    gsl_vector_complex_free(gsl_p_bicg);
    gsl_vector_complex_free(gsl_q);
    gsl_vector_complex_free(gsl_q_bicg);
    gsl_vector_complex_free(cur_gsl_x);
    gsl_vector_complex_free(gsl_b);

    free(diag_a);
    free(conj_diag_a);

}

// This function performs the axpy operation y = A*x
int cs_gaxpy_complex_with_gsl_x(const cs_ci *A, gsl_vector_complex* cur_gsl_x, gsl_vector_complex *gsl_y) {
    int p, j, n, *Ap, *Ai;
    double complex *Ax;

    n = A->n;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;

    for (j = 0; j < n; j++) {
        int start = Ap[j];
        int end = Ap[j + 1];
        double complex xj = gsl_vector_complex_get(cur_gsl_x, j);

        for (p = start; p < end; p++) {
            gsl_vector_complex_set(gsl_y, Ai[p], gsl_vector_complex_get(gsl_y, Ai[p]) + Ax[p] * xj);
        }
    }

    return 1;
}

// This function performs the axpy operation y = Atranspose*x
int cs_gaxpy_complex_with_gsl_x_conjtrans(const cs_ci *A, gsl_vector_complex* cur_gsl_x, gsl_vector_complex *gsl_y) {
    int p, j, n, *Ap, *Ai;
    double complex *Ax;

    n = A->n;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;

    for (j = 0; j < n; j++) {
        int start = Ap[j];
        int end = Ap[j + 1];
        double complex yj = 0.0;

        for (p = start; p < end; p++) {
            yj += gsl_complex_conjugate(Ax[p]) * gsl_vector_complex_get(cur_gsl_x, Ai[p]);
        }

        gsl_vector_complex_set(gsl_y, j, yj);
    }

    return 1;
}


// This function calculates the diagonal of A
double complex* find_diag_a_complex(const cs_ci *sparse_C) {

    int N = sparse_C->n;
    double complex *diag_a = (double complex*)malloc(N * sizeof(double complex));

    for (int j = 0; j < N; j++) {
        diag_a[j] = gsl_complex_rect(0.0, 0.0);  // Default value for missing diagonal elements

        for (int p = sparse_C->p[j]; p < sparse_C->p[j + 1]; p++) {
            if (sparse_C->i[p] == j) {
                diag_a[j] = sparse_C->x[p];
                break;  // Found the diagonal element, no need to continue
            }
        }
    }

    return diag_a;
}

// This function calculates the diagonal of A and also the conjugate of it
double complex* find_diag_a_conj_complex(const cs_ci *sparse_C) {

    int N = sparse_C->n;
    double complex *diag_a = (double complex*)malloc(N * sizeof(double complex));

    for (int j = 0; j < N; j++) {
        diag_a[j] = gsl_complex_rect(0.0, 0.0);  // Default value for missing diagonal elements

        for (int p = sparse_C->p[j]; p < sparse_C->p[j + 1]; p++) {
            if (sparse_C->i[p] == j) {
                diag_a[j] = conj(sparse_C->x[p]);
                break;  // Found the diagonal element, no need to continue
            }
        }
    }

    return diag_a;
}

gsl_vector_complex* preconditioner_solve_sparse_complex(gsl_vector_complex *gsl_r, gsl_vector_complex *gsl_z, const double complex *diag_a) {

    //gsl_vector *gsl_z = gsl_vector_alloc(A_dim);
    gsl_complex cur_A_element;
    gsl_complex cur_z_element;

    gsl_complex complex_1 = gsl_complex_rect(1.0, 0.0);

    for (int i=0; i<A_dim; i++) {

        cur_A_element = gsl_complex_rect(creal(diag_a[i]), cimag(diag_a[i]));

        if (GSL_REAL(cur_A_element) == 0 && GSL_IMAG(cur_A_element) == 0) cur_A_element = complex_1;

        cur_z_element = gsl_complex_div(gsl_vector_complex_get(gsl_r, i), cur_A_element);

        gsl_vector_complex_set(gsl_z, i, cur_z_element);
    }

    return gsl_z;
}