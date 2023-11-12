#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "equations.h"
#include "structs.h"
#include "parse.h"


double** A_array;
double* b_array;
int m2_count = 0;

// Function prototypes for filling the arrays
void fill_with_i(component* current);
void fill_with_r(component* current);
void fill_with_v(component* current);
void fill_with_c(component* current);
void fill_with_l(component* current);

double** alloc_A_array() {

    // Allocate memory for the array of pointers to rows
    double** new_array = (double**)calloc((nodes_n-1 + m2), sizeof(double*));

    // Check if memory allocation was successful
    if (new_array == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each row
    for (int i = 0; i < (nodes_n-1 + m2); i++) {
        new_array[i] = (double*)calloc((nodes_n-1 + m2), sizeof(double));

        // Check if memory allocation was successful
        if (new_array[i] == NULL) {
            fprintf(stderr, "Memory allocation failed\n");

            // If allocation fails, free previously allocated memory and exit
            for (int j = 0; j < i; j++) {
                free(new_array[j]);
            }
            free(new_array);

            exit(EXIT_FAILURE);
        }
    }

    // Return a pointer to the created array
    A_array = new_array;
    return new_array;
}

double* alloc_b_array() {

    // Allocate memory for the array of pointers to rows
    double* new_array = (double*)calloc((nodes_n-1 + m2), sizeof(double));

    // Check if memory allocation was successful
    if (new_array == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Return a pointer to the created array
    b_array = new_array;
    return new_array;
}


void create_arrays(){

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
        // else if(current->comp_type == 'c'){
        //     current = current->next;
        //     continue;
        // }
        else if(current->comp_type == 'l'){
            fill_with_l(current);
        }

        current = current->next;
    }
    return;
}

// Component: R ---> Fill A1*G*A1t array
void fill_with_r(component* current){

    int positive_node;
    int negative_node;

    // Take the node_index of each node
    positive_node = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: r 0 x value
    if(positive_node == 0){
        A_array[negative_node-1][negative_node-1] += 1/current->value;
    }
    // Case: r x 0 value
    else if(negative_node == 0){
        A_array[positive_node-1][positive_node-1] += 1/current->value;
    }
    // Case: x y value
    else{
        A_array[negative_node-1][negative_node-1] += 1/current->value;
        A_array[positive_node-1][positive_node-1] += 1/current->value;
        A_array[negative_node-1][positive_node-1] += -1/current->value;
        A_array[positive_node-1][negative_node-1] += -1/current->value;
    }

    return;
}

// Component: I ---> Fill b array
void fill_with_i(component* current){

    int positive_node;
    int negative_node;

    // Take the node_index of each node
    positive_node = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: i 0 x value
    if(positive_node == 0){
        b_array[negative_node-1] += current->value;
    }
    // Case: i x 0 value
    else if(negative_node == 0){
        b_array[positive_node-1] += -current->value;
    }
    // Case: i x y value
    else{
        b_array[negative_node-1] += current->value;
        b_array[positive_node-1] += -current->value;
    }

    return;
}

// Component: V ---> Fill A2t, A2 array and b array
void fill_with_v(component* current){

    int positive_node;
    int negative_node;

    // Take the node_index of each node
    positive_node = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: v 0 x value
    if(positive_node == 0){
        A_array[negative_node-1][nodes_n-1+m2_count] = -1;
        A_array[nodes_n-1+m2_count][negative_node-1] = -1;
        b_array[nodes_n-1+m2_count] = current->value;
    }
    // Case: v x 0 value
    else if(negative_node == 0){
        A_array[positive_node-1][nodes_n-1+m2_count] = 1;
        A_array[nodes_n-1+m2_count][positive_node-1] = 1;
        b_array[nodes_n-1+m2_count] = current->value;
    }
    // Case: v x y value
    else{
        A_array[negative_node-1][nodes_n-1+m2_count] = -1;
        A_array[nodes_n-1+m2_count][negative_node-1] = -1;
        A_array[positive_node-1][nodes_n-1+m2_count] = 1;
        A_array[nodes_n-1+m2_count][positive_node-1] = 1;
        b_array[nodes_n-1+m2_count] = current->value;
    }

    m2_count++;

    return;
}

// Component: L ---> Fill A2t, A2 array and b array (V node+ node- 0)
void fill_with_l(component* current){

    int positive_node;
    int negative_node;

    // Take the node_index of each node
    positive_node = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // Case: l 0 x value
    if(positive_node == 0){
        A_array[negative_node-1][nodes_n-1+m2_count] = -1;
        A_array[nodes_n-1+m2_count][negative_node-1] = -1;
        b_array[nodes_n-1+m2_count] = 0;
    }
    // Case: l x 0 value
    else if(negative_node == 0){
        A_array[positive_node-1][nodes_n-1+m2_count] = 1;
        A_array[nodes_n-1+m2_count][positive_node-1] = 1;
        b_array[nodes_n-1+m2_count] = 0;
    }
    // Case: l x y value
    else{
        A_array[negative_node-1][nodes_n-1+m2_count] = -1;
        A_array[nodes_n-1+m2_count][negative_node-1] = -1;
        A_array[positive_node-1][nodes_n-1+m2_count] = 1;
        A_array[nodes_n-1+m2_count][positive_node-1] = 1;
        b_array[nodes_n-1+m2_count] = 0;
    }

    m2_count++;

    return;
}

// print the arrays A and b
void print_arrays(){

    
    printf("\n *** Printing A array *** \n");

    for (int i = 0; i < (nodes_n-1 + m2); i++) {
        for (int j = 0; j < (nodes_n-1 + m2); j++) {
            printf("%.2lf\t", A_array[i][j]);
        }
        printf("\n");
    }

    printf("\n *** Printing b array *** \n");

    for (int k = 0; k < nodes_n-1 + m2; k++) {
        printf("%.2lf\t", b_array[k]);
    }
    printf("\n");

    return;
}

// Function to free the memory allocated for a (n+m) x (n+m) array
void free_A_array() {
    for (int i = 0; i < nodes_n-1 + m2; i++) {
        free(A_array[i]);
    }
    free(A_array);
}

// Function to free the memory allocated for a 1D array
void free_b_array() {
    free(b_array);
}

void create_equations(){

    alloc_A_array();
    alloc_b_array();
    create_arrays();
    print_arrays();

    return;
}