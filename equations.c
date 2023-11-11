#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "equations.h"
#include "structs.h"
#include "parse.h"


double** A_array;
double* b_array;
int count = 0;

int fill_with_i(component* current);
int fill_with_r(component* current);
int fill_with_v(component* current);
int fill_with_c(component* current);
int fill_with_l(component* current);

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

    // Return a pointer to the created matrix
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

    // Return a pointer to the created matrix
    b_array = new_array;
    return new_array;
}


int create_arrays(){

    component* current = head;

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
            current = current->next;
            continue;
        }
        else if(current->comp_type == 'l'){
            fill_with_l(current);
        }

        current = current->next;
    }
    return(1);
}

int fill_with_r(component* current){

    int positive_node;
    int negative_node;

    positive_node = atoi(current->positive_node);
    negative_node = atoi(current->negative_node);

    if(positive_node == 0){
        A_array[negative_node-1][negative_node-1] += 1/current->value;
    }
    else if(negative_node == 0){
        A_array[positive_node-1][positive_node-1] += 1/current->value;
    }
    else{
        A_array[negative_node-1][negative_node-1] += 1/current->value;
        A_array[positive_node-1][positive_node-1] += 1/current->value;
        A_array[negative_node-1][positive_node-1] += -1/current->value;
        A_array[positive_node-1][negative_node-1] += -1/current->value;
    }

    return(1);
}

int fill_with_i(component* current){

    int positive_node;
    int negative_node;

    positive_node = atoi(current->positive_node);
    negative_node = atoi(current->negative_node);

    if(positive_node == 0){
        b_array[negative_node-1] += current->value;
    }
    else if(negative_node == 0){
        b_array[positive_node-1] += -current->value;
    }
    else{
        b_array[negative_node-1] += current->value;
        b_array[positive_node-1] += -current->value;
    }

    return(1);
}

int fill_with_v(component* current){

    int positive_node;
    int negative_node;

    positive_node = atoi(current->positive_node);
    negative_node = atoi(current->negative_node);

    if(positive_node == 0){
        A_array[negative_node-1][nodes_n-1+count] = -1;
        A_array[nodes_n-1+count][negative_node-1] = -1;
        b_array[nodes_n-1+count] = current->value;
    }
    else if(negative_node == 0){
        A_array[positive_node-1][nodes_n-1+count] = 1;
        A_array[nodes_n-1+count][positive_node-1] = 1;
        b_array[nodes_n-1+count] = current->value;
    }
    else{
        A_array[negative_node-1][nodes_n-1+count] = -1;
        A_array[nodes_n-1+count][negative_node-1] = -1;
        A_array[positive_node-1][nodes_n-1+count] = 1;
        A_array[nodes_n-1+count][positive_node-1] = 1;
        b_array[nodes_n-1+count] = current->value;
    }

    count++;

    return(1);
}

int fill_with_l(component* current){

    int positive_node;
    int negative_node;

    positive_node = atoi(current->positive_node);
    negative_node = atoi(current->negative_node);

    if(positive_node == 0){
        A_array[negative_node-1][nodes_n-1+count] = -1;
        A_array[nodes_n-1+count][negative_node-1] = -1;
        b_array[nodes_n-1+count] = 0;
    }
    else if(negative_node == 0){
        A_array[positive_node-1][nodes_n-1+count] = 1;
        A_array[nodes_n-1+count][positive_node-1] = 1;
        b_array[nodes_n-1+count] = 0;
    }
    else{
        A_array[negative_node-1][nodes_n-1+count] = -1;
        A_array[nodes_n-1+count][negative_node-1] = -1;
        A_array[positive_node-1][nodes_n-1+count] = 1;
        A_array[nodes_n-1+count][positive_node-1] = 1;
        b_array[nodes_n-1+count] = 0;
    }

    count++;

    return(1);
}

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

void create_equations(){

    alloc_A_array();
    alloc_b_array();
    create_arrays();
    print_arrays();

    return;
}