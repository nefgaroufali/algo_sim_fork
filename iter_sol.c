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
#include "sparse_sol.h"
#include "csparse.h"
#include "math.h"
#include "string.h"


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
    //printf("b_norm = %lf\n", b_norm);
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

        //printf("r_norm = %lf\n", r_norm);
        //printf("b_norm = %lf\n", b_norm);

        gsl_z = preconditioner_solve(gsl_r, gsl_z);    // Solve Mz = r

        gsl_blas_ddot(gsl_r, gsl_z, &rho);   // dot(r,z)
        //printf("rho = %lf\n", rho);

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

//////////////////////////////////////////////////////////////////////
/////////////////////////// Sparse Solvers ///////////////////////////
//////////////////////////////////////////////////////////////////////


double sparse_norm(int m_row, int m_col, cs *m, double *v){
    
    double* between = (double *) calloc(m_row, sizeof(double));
    cs_gaxpy(m,v,between);
    double between2 = vector_inner_product(m_row, v, between);
   
    double absolute = fabs(between2);
    double res = sqrt(absolute);
    free(between);
    return(res);
}


double *solve_diagonal(int m_row, double *left, double *right){    
    
    int i; 
    double *res = (double *) calloc(m_row, sizeof(double));
    
    for(i=0; i<m_row; i++){
        res[i] = right[i]*left[i];
    }
    
    return(res);
    
}

double vector_inner_product(int m_row, double *first, double *second){

  int i; 
  double res;
  
  for(i=0; i<m_row; i++){
     res += first[i] * second[i];
  }
  
  return(res);
}

double *vector_sub(int m_row, double *first, double *second){

  int i; 
  double *res = (double *) calloc(m_row, sizeof(double));
  
  for(i=0; i<m_row; i++){
     res[i] = first[i] - second[i];
  }
  
  return(res);
}

double *vector_add(int m_row, double *first, double *second){

  int i; 
  double *res = (double *) calloc(m_row, sizeof(double));
  
  for(i=0; i<m_row; i++){
     res[i] = first[i] + second[i];
  }
  
  return(res);
}

double *vector_copy(int m_row, double *initial){

  int i; 
  double *res = (double *) calloc(m_row, sizeof(double));
  
  for(i=0; i<m_row; i++){
     res[i] = initial[i];
  }
  
  return(res);
}

double *constant_x_vector(int m_row, double *vector, double constant){
  int i; 
  double *res = (double *) calloc(m_row, sizeof(double));

  for(i=0; i<m_row; i++){
      res[i] = constant * vector[i];
  }
  
  return(res);

}


void solve_sparse_cg() {

    double* r = calloc(N, sizeof(double));    // residual vector
    double* Ax = calloc(N, sizeof(double));   // vector Ax
    double* z = calloc(N, sizeof(double));    // vector z
    double* p = calloc(N, sizeof(double));    // vector p
    double* q = calloc(N, sizeof(double));    // vector q

    double b_norm;
    double r_norm ;                         // norm of vector r
    double rho;
    double rho1;
    double beta;
    double alpha;
    int iter = 0;

    double* M = (double *) calloc(N, sizeof(double));

    // Calculate M   
    for(int j = 0 ; j < N ; j++){
        for (int k = sparse_C->p[j] ; k < sparse_C->p[j+1] ; k++){
            if (sparse_C->i[k] == j){
                M[j] += sparse_C->x[k];
            }
        }
    }
    
    for(int k=0;k<N;k++){
       if(M[k]==0){
         M[k] = 1;  
       }
       else{
         M[k] = 1/M[k];    
       }
    }

    Ax = (double *) calloc(N, sizeof(double));

    r = vector_sub(N, b_array_sparse, Ax);    // r = b - Ax
    
    r_norm = b_norm = sparse_norm(N, N, sparse_C, b_array_sparse);  // norm of vector b
    
    if(isnan(fabs(b_norm)) || b_norm == 0){
      r_norm = b_norm = 1;  
    }
    
    while(((r_norm/b_norm)>itol) && (iter<N)){
        iter++;
        
        z = solve_diagonal(N, M, r);        // Solve Mz = r
        rho = vector_inner_product(N, r, z);    // dot(r,z)
        
        if(iter == 1){
           p = vector_copy(N, z);       // p = z
        }
        else{
           beta = rho/rho1;
           double *temp = constant_x_vector(N, p, beta);    // beta*p
           free(p);
           p = vector_add(N, z, temp);      // p = z + beta*p
           free(temp);
        }
        
        rho1 = rho;
        for(int i=0;i<N;i++){
            q[i]=0;
        }
        cs_gaxpy(sparse_C,p,q);     // q = A*p
        
        alpha = rho/vector_inner_product(N, p, q);      // alpha = rho/(dot(p,q))
        double *temp1 = constant_x_vector(N, p, alpha);     // alpha*p
        double *temp2 = constant_x_vector(N, q, alpha);     // alpha*q
        double *temp3 = x_array_sparse;
        double *temp4 = r;
        x_array_sparse = vector_add(N, temp3, temp1);       // x = x + alpha*p
        r = vector_sub(N, temp4, temp2);        // r = r - alpha*q
        free(temp1);
        free(temp2);
        free(temp3);
        free(temp4);
        r_norm = sparse_norm(N,N,sparse_C,r);   // norm of vector r
        free(z);
    }

    free(Ax);
    free(r);
    //free(z);
    free(M);
    free(p);
    free(q);

}

void solve_sparse_bicg() {

    double *Ax, *r, *z, *M, *p, *q, *r_bicg, *z_bicg, *q_bicg, *p_bicg;
    double b_norm, r_norm, rho, rho1 = 1, beta, alpha, omega;
    int iter=0;
    
    q = (double *) calloc(N, sizeof(double));
    q_bicg = (double *) calloc(N, sizeof(double));
    M = (double *) calloc(N, sizeof(double));
    
    for(int j = 0 ; j < N ; j++){
        for (int k = sparse_C->p[j] ; k < sparse_C->p[j+1] ; k++){
            if (sparse_C->i[k] == j){
                M[j] += sparse_C->x[k];
            }
        }
    }
    
    for(int k=0;k<N;k++){
       if(M[k]==0){
         M[k] = 1;  
       }
       else{
         M[k] = 1/M[k];    
       }
    }
    
    Ax = (double *) calloc(N, sizeof(double));
    r = vector_sub(N, b_array_sparse, Ax);      // r = b - Ax
    r_bicg = vector_sub(N, b_array_sparse, Ax);     // r' = b - Ax
    
    r_norm = b_norm = sparse_norm(N, N, sparse_C, b_array_sparse);  // norm of vector b
    
    if(isnan(fabs(b_norm))|| b_norm == 0){
      r_norm = b_norm = 1;  
    }
    
    
    while(((r_norm/b_norm)>itol) && (iter<N)){
        iter++;
        
        z = solve_diagonal(N, M, r);    // Solve Mz = r
        z_bicg = solve_diagonal(N, M, r_bicg);  // Solve Mz' = r'
        rho = vector_inner_product(N, r_bicg, z);   // dot(r',z)
        
        if(isnan(fabs(rho))){
           break; 
        }
        if(iter == 1){
           p = vector_copy(N, z);   // p = z
           p_bicg = vector_copy(N, z_bicg);   // p' = z'
        }
        else{
           beta = rho/rho1;
           double *temp1 = constant_x_vector(N, p, beta);   // beta*p
           double *temp2 = constant_x_vector(N, p_bicg, beta);  // beta*p'
           free(p);
           free(p_bicg);
           p = vector_add(N, z, temp1);  // p = z + beta*p
           p_bicg = vector_add(N, z_bicg, temp2);   // p' = z' + beta*p'
           free(temp1);
           free(temp2);
        }
        
        rho1 = rho;
        for(int i=0;i<N;i++){
            q[i]=0;
        }
        cs_gaxpy(sparse_C,p,q);     // q = A*p

        // q' = A*p'
        for (int i = 0 ; i < N ; i++){
            q_bicg[i] = 0.0;
            for (int j = sparse_C->p[i] ; j < sparse_C->p[i+1] ; j++){
                q_bicg[i] = q_bicg[i] + sparse_C->x[j] * p_bicg[sparse_C->i[j]];
            }
        }
        
        // omega = dot(p',q)
        omega = vector_inner_product(N, p_bicg, q);
        if(isnan(fabs(omega))){
           break; 
        }
        alpha = rho/omega;
        double *temp1 = constant_x_vector(N, p, alpha);     // alpha*p
        double *temp2 = constant_x_vector(N, q, alpha);     // alpha*q
        double *temp3 = constant_x_vector(N, q_bicg, alpha);    // alpha*q'
        double *temp4 = x_array_sparse;
        x_array_sparse = vector_add(N, temp4, temp1);       // x = x + alpha*p
        double *temp5 = r;
        double *temp6 = r_bicg;
        r = vector_sub(N, temp5, temp2);    // r = r - alpha*q
        r_bicg = vector_sub(N, temp6, temp3);   // r' = r' - alpha*q'
        r_norm = sparse_norm(N,N,sparse_C,r);   // norm of vector r
        free(temp1);
        free(temp2);
        free(temp3);
        free(temp4);
        free(temp5);
        free(temp6);
        free(z);
        free(z_bicg);
    }
    
    free(p);
    free(p_bicg);
    free(r);
    free(r_bicg);
    free(q);
    free(q_bicg);
    free(M);
    free(Ax);
    
}
