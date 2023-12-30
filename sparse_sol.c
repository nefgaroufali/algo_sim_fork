#include <stdio.h>
#include "sparse_sol.h"
#include "mna.h"
#include "parse.h"
#include "structs.h"

cs *sparse_A;
cs *sparse_C;

int sparse_m2_count=0;

void form_sparse() { 

    component* current = head;
    int nonzero_counter = 0;

    sparse_A = cs_spalloc(A_dim, A_dim, nonzeros,1,1);

    // For each component of the circuit fill the sparse array
    // Also increment the nonzero counter

    while (current != NULL) {

        if(current->comp_type == 'r'){
            sparse_fill_with_r(current, &nonzero_counter);
        }
        else if(current->comp_type == 'v'){
            sparse_fill_with_v(current, &nonzero_counter);
        }
        else if(current->comp_type == 'l'){
            sparse_fill_with_v(current, &nonzero_counter);
        }

        current = current->next;
    }

    sparse_A->nz = nonzeros;

    printf("Before compression:\n");
    cs_print(sparse_A, 0);

    sparse_C = cs_compress(sparse_A);

    // cs_spfree(sparse_A);
    // cs_dupl(sparse_C);

    printf("After compression:\n");
    cs_print(sparse_C, 0);

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

// Component: V or L ---> Fill A2t, A2 array and b array
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

    sparse_m2_count++;
}

void add_to_sparse_A(int k, int i, int j, double x) {
        sparse_A->i[k] = i;
        sparse_A->p[k] = j;
        sparse_A->x[k] = x;
}
