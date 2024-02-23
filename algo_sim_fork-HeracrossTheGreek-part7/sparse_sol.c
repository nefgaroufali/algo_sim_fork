#include <stdio.h>
#include <stdlib.h>
#include "sparse_sol.h"
#include "mna.h"
#include "parse.h"
#include "structs.h"
#include "gsl.h"

cs *sparse_A;
cs *sparse_cc_A;
cs *sparse_C;
cs *sparse_cc_C;

double *b_array_sparse = NULL;
double *x_array_sparse = NULL;
int sparse_m2_count=0;

// This function fills the sparse array A, which is in triplet form(see Timothy Davies book)
void form_sparse() { 

    component* current = head;
    int nonzero_counter_A = 0;
    int nonzero_counter_C = 0;

    // Allocate the sparse A array
    sparse_A = cs_spalloc(A_dim, A_dim, nonzeros_A,1,1);
    sparse_C = cs_spalloc(A_dim, A_dim, nonzeros_C,1,1);

    // create vectors for x and b
    create_sparse_vectors();

    // For each component of the circuit fill the sparse array

    while (current != NULL) {

        if(current->comp_type == 'r'){
            sparse_fill_with_r(current, &nonzero_counter_A);
        }
        else if(current->comp_type == 'v'){
            sparse_fill_with_v(current, &nonzero_counter_A);
        }
        else if(current->comp_type == 'l'){
            sparse_fill_with_l(current, &nonzero_counter_A, nonzero_counter_C);
            nonzero_counter_C++; // always added 1 element for l in C
        }
        else if(current->comp_type == 'c'){
            sparse_fill_with_c(current, &nonzero_counter_C);
        }
        else if(current->comp_type == 'i'){
            sparse_fill_with_i(current);
        }

        current = current->next;
    }

    sparse_A->nz = nonzeros_A;
    sparse_C->nz = nonzeros_C;

    // Form the compressed-column arrays cc_A and cc_C
    sparse_cc_A = cs_compress(sparse_A);
    sparse_cc_C = cs_compress(sparse_C);

    cs_spfree(sparse_A);
    cs_spfree(sparse_C);

    // Remove the duplicates of the compressed-column arrays
    cs_dupl(sparse_cc_A);
    cs_dupl(sparse_cc_C);

    // We also create gsl b for iter solver
    gsl_b = gsl_vector_alloc(A_dim);

    // Copy the array b to a gsl vector, which will be used for sparse methods
    double_to_gsl(gsl_b, b_array_sparse);


}

// Component: R ---> Fill A1*G*A1t array
void sparse_fill_with_r(component* current, int *nonzero_counter){

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

        add_to_sparse_A(k, negative_node_i-1, negative_node_i-1, current_g);
        *nonzero_counter += 1;
    }

    // Case: r x 0 value
    else if (negative_node_i == 0) {

        add_to_sparse_A(k, positive_node_i-1, positive_node_i-1, current_g);
        *nonzero_counter += 1;
    }

    // Case: x y value
    else {
        add_to_sparse_A(k, negative_node_i-1, negative_node_i-1, current_g);
        add_to_sparse_A(k+1, positive_node_i-1, positive_node_i-1, current_g);
        add_to_sparse_A(k+2, negative_node_i-1, positive_node_i-1, -current_g);
        add_to_sparse_A(k+3, positive_node_i-1, negative_node_i-1, -current_g);
        *nonzero_counter += 4;
    }

    // return;
}

// Component: C ---> Fill C array (Transient)
void sparse_fill_with_c(component* current, int *nonzero_counter){

    int positive_node_i;
    int negative_node_i;
    int k = *nonzero_counter;


    double current_c = current->value;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: r 0 x value
    if (positive_node_i == 0) {

        add_to_sparse_C(k, negative_node_i-1, negative_node_i-1, current_c);
        *nonzero_counter += 1;
    }

    // Case: r x 0 value
    else if (negative_node_i == 0) {

        add_to_sparse_C(k, positive_node_i-1, positive_node_i-1, current_c);
        *nonzero_counter += 1;
    }

    // Case: x y value
    else {
        add_to_sparse_C(k, negative_node_i-1, negative_node_i-1, current_c);
        add_to_sparse_C(k+1, positive_node_i-1, positive_node_i-1, current_c);
        add_to_sparse_C(k+2, negative_node_i-1, positive_node_i-1, -current_c);
        add_to_sparse_C(k+3, positive_node_i-1, negative_node_i-1, -current_c);
        *nonzero_counter += 4;
    }

    // return;
}

// Component: I ---> Fill b array
void sparse_fill_with_i(component* current){

    int positive_node_i;
    int negative_node_i;

    double current_i = current->value;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: i 0 x value
    if(positive_node_i == 0){
        b_array_sparse[negative_node_i-1] += current_i;
    }
    // Case: i x 0 value
    else if(negative_node_i == 0) {
        b_array_sparse[positive_node_i-1] -= current_i;
    }
    // Case: i x y value
    else{
        b_array_sparse[negative_node_i-1] += current_i;
        b_array_sparse[positive_node_i-1] -= current_i;
    }

    return;
}

// Component: V ---> Fill A2t, A2 array and b array
void sparse_fill_with_v(component* current, int *nonzero_counter) {

    int positive_node_i;
    int negative_node_i;
    int sparse_m2_i = nodes_n - 1 + sparse_m2_count;
    int k = *nonzero_counter;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: v 0 x value
    if(positive_node_i == 0){
        add_to_sparse_A(k, negative_node_i-1, sparse_m2_i, -1);
        add_to_sparse_A(k+1, sparse_m2_i, negative_node_i-1, -1);
        *nonzero_counter += 2;
    }
    // Case: v x 0 value
    else if(negative_node_i == 0){
        add_to_sparse_A(k, positive_node_i-1, sparse_m2_i, 1);
        add_to_sparse_A(k+1, sparse_m2_i, positive_node_i-1, 1);
        *nonzero_counter += 2;
    }
    // Case: v x y value
    else{
        add_to_sparse_A(k, negative_node_i-1, sparse_m2_i, -1);
        add_to_sparse_A(k+1, sparse_m2_i, negative_node_i-1, -1);
        add_to_sparse_A(k+2, positive_node_i-1, sparse_m2_i, 1);
        add_to_sparse_A(k+3, sparse_m2_i, positive_node_i-1, 1);
        *nonzero_counter += 4;
    }

    b_array_sparse[sparse_m2_i] = current->value;

    sparse_m2_count++;
}

// Component: L ---> Fill A2t, A2 array and b array, as well as C(transient) array
void sparse_fill_with_l(component* current, int *nonzero_counter_A, int nonzero_counter_C) {

    int positive_node_i;
    int negative_node_i;
    int sparse_m2_i = nodes_n - 1 + sparse_m2_count;
    int k = *nonzero_counter_A;
    int l = nonzero_counter_C;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: v 0 x value
    if(positive_node_i == 0){
        add_to_sparse_A(k, negative_node_i-1, sparse_m2_i, -1);
        add_to_sparse_A(k+1, sparse_m2_i, negative_node_i-1, -1);
        *nonzero_counter_A += 2;
    }
    // Case: v x 0 value
    else if(negative_node_i == 0){
        add_to_sparse_A(k, positive_node_i-1, sparse_m2_i, 1);
        add_to_sparse_A(k+1, sparse_m2_i, positive_node_i-1, 1);
        *nonzero_counter_A += 2;
    }
    // Case: v x y value
    else{
        add_to_sparse_A(k, negative_node_i-1, sparse_m2_i, -1);
        add_to_sparse_A(k+1, sparse_m2_i, negative_node_i-1, -1);
        add_to_sparse_A(k+2, positive_node_i-1, sparse_m2_i, 1);
        add_to_sparse_A(k+3, sparse_m2_i, positive_node_i-1, 1);
        *nonzero_counter_A += 4;
    }

    // Always add to C (Transient)
    add_to_sparse_C(l, sparse_m2_i, sparse_m2_i, -current->value);   

    b_array_sparse[sparse_m2_i] = 0;

    sparse_m2_count++;
}

// This function adds an element to the sparse array A
void add_to_sparse_A(int k, int i, int j, double x) {
        sparse_A->i[k] = i;
        sparse_A->p[k] = j;
        sparse_A->x[k] = x;
}

void add_to_sparse_C(int k, int i, int j, double x) {
        sparse_C->i[k] = i;
        sparse_C->p[k] = j;
        sparse_C->x[k] = x;
}

// This function fills the vectors b and x, used in sparse operations
void create_sparse_vectors() {

    // Allocate memory for the array of pointers to rows
    b_array_sparse = (double*)calloc(A_dim, sizeof(double));
    x_array_sparse = (double*)calloc(A_dim, sizeof(double));

    if(x_array_sparse == NULL){
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Check if memory allocation was successful
    if (b_array_sparse == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
}

// This function prints arrays A, b and x, used in sparse operations
void print_sparse_arrays() {

    printf("\n *** Printing A array *** \n");

    cs_print(sparse_cc_A, 0);

    printf("\n *** Printing b array *** \n");

    for (int k = 0; k < A_dim; k++) {
        printf("%.2lf\t", b_array_sparse[k]);
    }

    printf("\n *** Printing C array *** \n");

    cs_print(sparse_cc_C, 0);

    printf("\n *** Printing x array *** \n");

    for (int l = 0; l < A_dim; l++) {
        printf("%.2lf\t", x_array_sparse[l]);
    }

    printf("\n");
}

// This function only prints vector x
void print_sparse_x() {

    for (int l = 0; l < A_dim; l++) {
        printf("%.2lf\t", x_array_sparse[l]);
    }

    printf("\n");
}