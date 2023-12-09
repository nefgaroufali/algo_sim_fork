#include "direct_sol.h"
#include "gsl.h"
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "mna.h"
#include "structs.h"
#include "parse.h"

gsl_matrix* gsl_A = NULL;
gsl_matrix* gsl_LU = NULL;
gsl_matrix* gsl_chol = NULL;
gsl_vector *gsl_b = NULL;
gsl_vector *gsl_x = NULL;
gsl_permutation *gsl_p = NULL;


FILE *filePointers[5] = {NULL};

int spd = FALSE;
int *plot_node_indexes = NULL;
int plot_node_count = 0;


// This function copies the double array A made in mna.c, to a GSL array, to be used in the linear algebra routines
// It also copies the elements of the b vector

void form_gsl_system() {

    int i,j;
    i = j = A_dim;  // A_dim = n

    gsl_A = gsl_matrix_alloc(A_dim, A_dim);
    gsl_b = gsl_vector_alloc(A_dim);

    // Copying each element of A to the respective position of the GSL array

    for (i=0; i<A_dim; i++) {
        for (j=0; j<A_dim; j++) {
            gsl_matrix_set(gsl_A, i, j, A_array[i*A_dim + j]);
        }
    }

    // Copying each element of b to the respective position of the GSL Vector

    for (i=0; i<A_dim; i++) {
        gsl_vector_set(gsl_b, i, b_array[i]);
    }
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

// This function prints the specified gsl vector
void print_gsl_vector(gsl_vector *vector, int dim){

    int i;

    printf("\n");
    for (i = 0; i < dim; i++) {
        printf ("%.5lf\n", gsl_vector_get(vector, i));
    }
    printf("\n");

}

// This function frees all memory occupied by gsl arrays

void free_gsl() {

    gsl_vector_free(gsl_x);
    gsl_vector_free(gsl_b);
    gsl_matrix_free(gsl_A);
    gsl_matrix_free(gsl_LU);
    gsl_matrix_free(gsl_chol);
    gsl_permutation_free(gsl_p);
}

// Custom error handler function
void gslErrorHandler(const char *reason, const char *file, int line, int gsl_errno) {

    // Print an error message
    fprintf(stderr, "GSL Error: %s \n", reason);

}

// This function solves the linear system with the altered b vector
// Then prints the value of each plot node in a separate file, that corresponds to each different value of the sweep component

void solve_dc_sweep_system(gsl_vector *temp_gsl_b, double cur_value, char type) {

    gsl_vector *temp_gsl_x = gsl_vector_alloc(A_dim);

    if (solver_type == NON_SPD_SOL) {
        gsl_linalg_LU_solve(gsl_LU, gsl_p, temp_gsl_b, temp_gsl_x);
    }
    else {
        gsl_linalg_cholesky_solve(gsl_chol, temp_gsl_b, temp_gsl_x);
    }

    int i;
    int plot_node_i;
    double b_vector_value, x_vector_value;

    // Plot_node_i: The index of the node(s) that are to be plotted
    // Sweep_node_i: The index of the node whose value in the b vector changes

    // Get the values that will be printed to the files, then call add_to_plot_file
    for (i = 0; i < plot_node_count; i++) {
        plot_node_i = plot_node_indexes[i];
        b_vector_value = cur_value;
        x_vector_value = gsl_vector_get(temp_gsl_x, plot_node_i);

        add_to_plot_file(b_vector_value, x_vector_value, i);

    }

    gsl_vector_free(temp_gsl_x);

}

// Adds the node to be plotted to a list of all the nodes to be plotted
// IMPORTANT: These indexes start at 0.
void add_plot_node(int node_i) {

    if (plot_node_indexes== NULL) {
        plot_node_indexes= (int *)malloc(sizeof(int) * 1);
    } else {
        plot_node_indexes= (int *)realloc(plot_node_indexes, sizeof(int) * (plot_node_count + 1));
    }

    // Check if memory allocation failed
    if (plot_node_indexes== NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    plot_node_indexes[plot_node_count] = node_i;
    plot_node_count++;
}

// This function keeps as many file pointers as the plot commands
// Then prints the corresponding values to each one

void add_to_plot_file(double b_vector_value, double x_vector_value, int i) {

    // Check if the file corresponding to i is already open
    if (filePointers[i] == NULL) {
        char filename[20];
        snprintf(filename, sizeof(filename), "file_%d.txt", i);

        // Open the file in write mode
        filePointers[i] = fopen(filename, "w");

        // Check if the file is opened successfully
        if (filePointers[i] == NULL) {
            printf("Error opening file %d\n", i);
            return;
        }
    }

    // Write the values to the file
    fprintf(filePointers[i], "%lf %lf\n", b_vector_value, x_vector_value);

    // Flush the file buffer to ensure data is written immediately
    fflush(filePointers[i]);

    // Code to generate the GNU Plot command to create the plot
    char plot_command[100];
    snprintf(plot_command, sizeof(plot_command), "set terminal png; set output 'file_%d.png'; plot 'file_%d.txt' with lines", i, i);

    // Generate a temporary script file to run the GNU Plot commands
    FILE *gnuplotScript = fopen("plot_script.gnu", "w");
    if (gnuplotScript == NULL) {
        printf("Error creating GNU Plot script file\n");
        return;
    }

    // Write the GNU Plot command to the script file
    fprintf(gnuplotScript, "%s", plot_command);
    fclose(gnuplotScript);
    

    // // Execute GNU Plot using the script file
    // system("gnuplot plot_script.gnu");

    // // Clean up: remove the temporary script file
    // remove("plot_script.gnu");
}

// This function frees the memory occupied by the plot node indexes
void free_plot_node() {

    free(plot_node_indexes);
}