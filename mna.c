#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mna.h"
#include "structs.h"
#include "parse.h"


double* A_array = NULL;
double* b_array = NULL;
double *C_array = NULL;
int m2_count=0;

double* alloc_A_array() {

    // Allocate memory for the array
    double* new_array = (double*)calloc((A_dim*A_dim), sizeof(double));

    // Check if memory allocation was successful
    if (new_array == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Return a pointer to the created array
    A_array = new_array;

    return new_array;
}

double* alloc_b_array() {

    // Allocate memory for the array of pointers to rows
    double* new_array = (double*)calloc((A_dim), sizeof(double));

    // Check if memory allocation was successful
    if (new_array == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Return a pointer to the created array
    b_array = new_array;

    return new_array;
}

double* alloc_C_array() {

    // Allocate memory for the array
    double* new_array = (double*)calloc((A_dim*A_dim), sizeof(double));

    // Check if memory allocation was successful
    if (new_array == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Return a pointer to the created array
    C_array = new_array;

    return new_array;
}


void fill_arrays(){

    component* current = head;

    // For each component of the circuit fill the arrays
    while (current != NULL) {
        
        if(current->comp_type == 'r'){
            fill_with_r(current);
        }
        else if(current->comp_type == 'i'){
            fill_with_i(current);
        }
        else if(current->comp_type == 'v'){
            fill_with_v(current);
        }
        else if(current->comp_type == 'c'){
            fill_with_c(current);
        }
        else if(current->comp_type == 'l'){
            fill_with_l(current);
        }

        current = current->next;
    }
}

// Component: R ---> Fill A1*G*A1t array
void fill_with_r(component* current){

    int positive_node_i;
    int negative_node_i;

    double current_g = 1/current->value; // g = 1/R

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: r 0 x value
    if (positive_node_i == 0) {
        A_array[(negative_node_i-1) * A_dim + (negative_node_i-1)] += current_g;
    }

    // Case: r x 0 value
    else if (negative_node_i == 0) {
        A_array[(positive_node_i-1) * A_dim + (positive_node_i-1)] += current_g;
    }

    // Case: x y value
    else {
        A_array[(negative_node_i-1) * A_dim + (negative_node_i-1)] += current_g;
        A_array[(positive_node_i-1) * A_dim + (positive_node_i-1)] += current_g;
        A_array[(negative_node_i-1) * A_dim + (positive_node_i-1)] -= current_g;
        A_array[(positive_node_i-1) * A_dim + (negative_node_i-1)] -= current_g;
    }


    return;
}


// Component: I ---> Fill b array
void fill_with_i(component* current){

    int positive_node_i;
    int negative_node_i;

    double current_i = current->value;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: i 0 x value
    if(positive_node_i == 0){
        b_array[negative_node_i-1] += current_i;
    }
    // Case: i x 0 value
    else if(negative_node_i == 0) {
        b_array[positive_node_i-1] -= current_i;
    }
    // Case: i x y value
    else{
        b_array[negative_node_i-1] += current_i;
        b_array[positive_node_i-1] -= current_i;
    }

    return;
}

// Component: V ---> Fill A2t, A2 array and b array
void fill_with_v(component* current){

    int positive_node_i;
    int negative_node_i;
    int m2_i = nodes_n - 1 + m2_count;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: v 0 x value
    if(positive_node_i == 0){
        A_array[(negative_node_i-1)*A_dim +m2_i] = -1;
        A_array[m2_i*A_dim +(negative_node_i-1)] = -1;
        b_array[m2_i] = current->value;
    }
    // Case: v x 0 value
    else if(negative_node_i == 0){
        A_array[(positive_node_i-1)*A_dim+m2_i] = 1;
        A_array[m2_i*A_dim +(positive_node_i-1)] = 1;
        b_array[m2_i] = current->value;
    }
    // Case: v x y value
    else{
        A_array[(negative_node_i-1)*A_dim +m2_i] = -1;
        A_array[m2_i*A_dim +(negative_node_i-1)] = -1;
        A_array[(positive_node_i-1)*A_dim+m2_i] = 1;
        A_array[m2_i*A_dim +(positive_node_i-1)] = 1;
        b_array[m2_i] = current->value;
    }

    m2_count++;

    return;
}

// ONLY fill the C (transient) array 
void fill_with_c(component* current) {

    int positive_node_i;
    int negative_node_i;

    double current_c = current->value; // g = 1/R

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: r 0 x value
    if (positive_node_i == 0) {
        C_array[(negative_node_i-1) * A_dim + (negative_node_i-1)] += current_c;
    }

    // Case: r x 0 value
    else if (negative_node_i == 0) {
        C_array[(positive_node_i-1) * A_dim + (positive_node_i-1)] += current_c;
    }

    // Case: x y value
    else {
        C_array[(negative_node_i-1) * A_dim + (negative_node_i-1)] += current_c;
        C_array[(positive_node_i-1) * A_dim + (positive_node_i-1)] += current_c;
        C_array[(negative_node_i-1) * A_dim + (positive_node_i-1)] -= current_c;
        C_array[(positive_node_i-1) * A_dim + (negative_node_i-1)] -= current_c;
    }

}

// Component: L ---> Fill A2t, A2 array and b array (V node+ node- 0)
void fill_with_l(component* current){

    int positive_node_i;
    int negative_node_i;
    int m2_i = nodes_n - 1 + m2_count;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: l 0 x value
    if(positive_node_i == 0){
        A_array[(negative_node_i-1)*A_dim +m2_i] = -1;
        A_array[m2_i*A_dim +(negative_node_i-1)] = -1;
        b_array[m2_i] = 0;
    }
    // Case: l x 0 value
    else if(negative_node_i == 0){
        A_array[(positive_node_i-1)*A_dim+m2_i] = 1;
        A_array[m2_i*A_dim +(positive_node_i-1)] = 1;
        b_array[m2_i] = 0;
    }
    // Case: l x y value
    else{
        A_array[(negative_node_i-1)*A_dim + m2_i] = -1;
        A_array[m2_i*A_dim +(negative_node_i-1)] = -1;
        A_array[(positive_node_i-1)*A_dim + m2_i] = 1;
        A_array[m2_i*A_dim +(positive_node_i-1)] = 1;
        b_array[m2_i] = 0;
    }

    // Always fill the C array (for Transient analysis)
    C_array[m2_i*A_dim + m2_i] = - current->value;
 
    m2_count++;

    return;
}

// print the arrays A and b
void print_arrays(){

    
    printf("\n *** Printing A array *** \n");

    for (int i = 0; i < (A_dim); i++) {
        for (int j = 0; j < (A_dim); j++) {
            printf("%.2lf\t", A_array[i*A_dim + j]);
        }
        printf("\n");
    }

    printf("\n *** Printing b array *** \n");

    for (int k = 0; k < A_dim; k++) {
        printf("%.2lf\t", b_array[k]);
    }

    printf("\n *** Printing C array *** \n");

    for (int i = 0; i < (A_dim); i++) {
        for (int j = 0; j < (A_dim); j++) {
            printf("%.2lf\t", C_array[i*A_dim + j]);
        }
        printf("\n");
    }
    printf("\n");

    return;
}

// Function to free the memory allocated for a (n+m) x (n+m) array
void free_A_array() {

    free(A_array);
}

// Function to free the memory allocated for a 1D array
void free_b_array() {

    free(b_array);
}

void free_C_array() {

    free(C_array);
}

// This function creates and fills the necessary arrays to solve the MNA system
void create_equations(){

    alloc_A_array();
    alloc_b_array();
    alloc_C_array();
    fill_arrays();
    //print_arrays();

    return;
}