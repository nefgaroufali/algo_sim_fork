#include "direct_sol.h"
#include "iter_sol.h"
#include "gsl.h"
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "mna.h"
#include "structs.h"
#include "parse.h"
#include "sparse_sol.h"
#include "csparse.h"
#include "time.h"
#include <math.h>

gsl_matrix* gsl_A = NULL;
gsl_vector *gsl_b = NULL;
gsl_matrix *gsl_C = NULL;
gsl_vector *gsl_x = NULL;

FILE *filePointers[5] = {NULL};

int spd = FALSE;
int *plot_node_indexes = NULL;
int plot_node_count = 0;

// This function copies the double array A made in mna.c, to a GSL array, to be used in the linear algebra routines
// It also copies the elements of the b vector

// This function cretes the gsl arrays A and b and C, and fills them
void form_gsl_system() {

    int i,j;
    i = j = A_dim;  // A_dim = n

    gsl_A = gsl_matrix_alloc(A_dim, A_dim);
    gsl_b = gsl_vector_alloc(A_dim);
    gsl_C = gsl_matrix_alloc(A_dim, A_dim);

    // Copying each element of A and C to the respective position of the GSL array

    for (i=0; i<A_dim; i++) {
        for (j=0; j<A_dim; j++) {
            gsl_matrix_set(gsl_A, i, j, A_array[i*A_dim + j]);
            gsl_matrix_set(gsl_C, i, j, C_array[i*A_dim + j]);
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
    gsl_matrix_free(gsl_C);
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

void solve_dc_sweep_system(gsl_vector *temp_gsl_b, double cur_value) {

    gsl_vector *temp_gsl_x;


    temp_gsl_x = gsl_vector_calloc(A_dim);


    //gsl_vector_memcpy(temp_gsl_x, gsl_x);   // for the iterative methods, temp_gsl_x receives the value of the previous solution



    // Solve the system based on the solver flag, generated during parsing

    if (solver_type == LU_SOL) {
        gsl_linalg_LU_solve(gsl_LU, gsl_p, temp_gsl_b, temp_gsl_x);
    }
    else if (solver_type == CHOL_SOL) {
        gsl_linalg_cholesky_solve(gsl_chol, temp_gsl_b, temp_gsl_x);
    }
    else if (solver_type == CG_SOL) {
        solve_cg(temp_gsl_b, temp_gsl_x);
    }
    else if (solver_type == BICG_SOL) {
        solve_bicg(temp_gsl_b, temp_gsl_x);
    }
    else if (solver_type == SPARSE_LU_SOL) {
        solve_sparse_lu(temp_gsl_b, temp_gsl_x);
    }
    else if (solver_type == SPARSE_CHOL_SOL) {
        solve_sparse_chol(temp_gsl_b, temp_gsl_x);
    }
    if (solver_type == SPARSE_CG_SOL) {
        solve_sparse_cg(temp_gsl_b, temp_gsl_x);
    }
    else if (solver_type == SPARSE_BICG_SOL) {
        solve_sparse_bicg(temp_gsl_b, temp_gsl_x);
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

    //gsl_vector_memcpy(gsl_x, temp_gsl_x); // copy the value to gsl_x, so that it is the initial value for the next cg/bicg
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
        char filename[60];
        snprintf(filename, sizeof(filename), "output/%s_%d.txt", circuit_name, i);

        // Open the file in write mode inside the "output" subdirectory
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
    char plot_command[150];
    snprintf(plot_command, sizeof(plot_command), "set terminal png; set output 'output/%s_%d.png'; plot 'output/%s_%d.txt' with lines", circuit_name, i, circuit_name, i);

    // Generate a temporary script file to run the GNU Plot commands
    FILE *gnuplotScript = fopen("plot_script.gnu", "w");
    if (gnuplotScript == NULL) {
        printf("Error creating GNU Plot script file\n");
        return;
    }

    // Write the GNU Plot command to the script file
    fprintf(gnuplotScript, "%s", plot_command);
    fclose(gnuplotScript);
    

    // // // Execute GNU Plot using the script file
    // system("gnuplot plot_script.gnu");

    // // // Clean up: remove the temporary script file
    // remove("plot_script.gnu");
}

// This function frees the memory occupied by the plot node indexes
void free_plot_node() {

    free(plot_node_indexes);
}

// This function performs DC sweep: It solves the Ax =b system for all different values of a specific voltage or current source
// The name of this source is parsed earlier and its pointer is sweep_component
void dc_sweep() {

    double cur_value;
    double low = DC_arguments[0];
    double high = DC_arguments[1];
    double step = DC_arguments[2];

    int pos_sweep_node_i, neg_sweep_node_i, sweep_node_i;

    // Temp_gsl_b: Same as gsl_b, with only the sweep value changing

    gsl_vector *temp_gsl_b = gsl_vector_alloc(A_dim);
    gsl_vector_memcpy(temp_gsl_b, gsl_b);

    component *current = sweep_component;

    // If the component is a current source, find the index of both nodes in the hash table

    // Change the current value, adding the step in each iteration, then change the bvector
    // Then solve the linear system in each iteration, in function solve_dc_sweep_system

    // Only change the b vector if the node is not ground

    if (current->comp_type == 'i') {
        pos_sweep_node_i = find_hash_node(&node_hash_table, current->positive_node) - 1; // -1 because 0 is ground
        neg_sweep_node_i = find_hash_node(&node_hash_table, current->negative_node) - 1;

        cur_value = low;
        do {
            if (pos_sweep_node_i != -1) {
                gsl_vector_set(temp_gsl_b, pos_sweep_node_i, -cur_value);
            }
            if (neg_sweep_node_i != -1) {
                gsl_vector_set(temp_gsl_b, neg_sweep_node_i, cur_value);
            }

            solve_dc_sweep_system(temp_gsl_b, cur_value); 

            cur_value = cur_value + step;
        } while (cur_value <= high);

    }

    // If the component is a voltage source, take its index in the m2 part of the b vector

    // Change the current value, adding the step in each iteration, then change the bvector
    // Then solve the linear system in each iteration, in function solve_dc_sweep_system

    if (current->comp_type == 'v') {
        sweep_node_i = nodes_n-1 +current->m2_i;

        cur_value = low;
        do {
            gsl_vector_set(temp_gsl_b, sweep_node_i, cur_value); 
            solve_dc_sweep_system(temp_gsl_b, cur_value); 
            cur_value = cur_value + step;
        } while (cur_value <= high);
    }

    gsl_vector_free(temp_gsl_b);

}


// This function solves the  DC system,
// Then writes the dc operating point in a .op file

void solve_dc_system(int solver_type) {

    clock_t t1, t2;
    double time;

    gsl_x = gsl_vector_calloc(A_dim);   // fill it with zeros, useful for the iterative methods

    t1 = clock();

    // Solve the system depending on the solver type
    if (solver_type == LU_SOL) {
        gsl_linalg_LU_solve(gsl_LU, gsl_p, gsl_b, gsl_x);
    }
    else if (solver_type == CHOL_SOL) {
        gsl_linalg_cholesky_solve(gsl_chol, gsl_b, gsl_x);
    }
    else if (solver_type == CG_SOL){
        solve_cg(gsl_b, gsl_x);
    }
    else if (solver_type == BICG_SOL){
        solve_bicg(gsl_b, gsl_x);
    }
    else if (solver_type == SPARSE_LU_SOL){
        solve_sparse_lu(gsl_b, gsl_x);
    }
    else if (solver_type == SPARSE_CHOL_SOL){
        solve_sparse_chol(gsl_b, gsl_x);
    }
    else if (solver_type == SPARSE_CG_SOL){
        solve_sparse_cg(gsl_b, gsl_x);
    }
    else if (solver_type == SPARSE_BICG_SOL){
        solve_sparse_bicg(gsl_b, gsl_x);
    }

    t2 = clock();
    time = (double) (t2 - t1) / CLOCKS_PER_SEC;
    printf("Completed in %f seconds\n", time);

    // Open the file in which the operating point will be printed
    FILE *op_file;
    char op_file_name[40];

    snprintf(op_file_name, sizeof(op_file_name), "%s%s%s", "output/", circuit_name, ".op");
    op_file = fopen(op_file_name, "w");

    // Print the solution to the file
    // First n-1 elements are the voltages of the nodes
    // Last m2 elements are the currents of the V and L components

    for (int i = 0, j=0; i < A_dim; i++) {

        if(i < nodes_n-1){
            fprintf(op_file, "%s\t\t%.5e\n", node_array[i+1], gsl_vector_get(gsl_x, i));    
        }
        else
        {
            fprintf(op_file, "%s\t\t%.5e\n", m2_array[j], gsl_vector_get(gsl_x, i));
            j++;
        }
    }

    fclose(op_file);

}

// This function converts a double to a gsl vector
void double_to_gsl(gsl_vector *gsl_v, double* v) {

    int i;

    for (i=0; i<A_dim; i++) {
        gsl_vector_set(gsl_v, i, v[i]);
    }

}

// This function converts a gsl vector to a double 
void gsl_to_double(gsl_vector *gsl_v, double* v) {

    int i;

    for (i=0; i<A_dim; i++) {
        v[i] = gsl_vector_get(gsl_v, i);

    }

}


// get the value of the exp transient function at time t
double get_exp_val(exp_spec *data, double t) {
	double res;

    // negative time -> return zero
	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}
    // t < td1
	else if (t < data->td1) {
		res = data->i1;
	}
    // td1 < t < td2
	else if (t < data->td2) {
		res = data->i1 + (data->i2 - data->i1)*(1.0 - exp(-1*(t - data->td1)/data->tc1));
	}
    // t > td2
	else {
		res = data->i1 + (data->i2 - data->i1)*(exp(-1*(t- data->td2)/data->tc2) - exp(-1*(t - data->td1)/data->tc1));
	}
	return res;
}


// get the value of the sin tansient function at time t
double get_sin_val(sin_spec *data, double t) {
	double res;

    // negative time -> return zero
	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}
    // t < td 
	else if (t < data->td) {
		res = data->i1 + data->ia*sin(2*M_PI/360.0);
	}
    // t > td
	else {
		res = data->i1 + data->ia * sin(2*M_PI * data->fr*(t - data->td) + 2*M_PI * data->ph/360.0)*exp(-1*(t - data->td)*data->df);
	}
	return res;
}


// get the value of the pulse transient function at time t
double get_pulse_val(pulse_spec *data, double t) {
	double res;
	double y2, y1;
	double x1;


    // negative time -> return zero
	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}
    // t < td
	else if (t < data->td) {
		res = data->i1;
	}
	else {
		t = t - (int)((t-data->td)/data->per)*data->per;

		if (t < data->td + data->tr) {
			// i1->i2 linearly
			/*x2 = data->td + data->tr;*/
			x1 = data->td;
			y2 = data->i2;
			y1 = data->i1;
			/*res = ((y2-y1)/(x2-x1)) * (t - x1) + y1;*/
			res = ((y2 - y1)/data->tr)* (t - x1) + y1;
		}
		else if (t < data->td + data->tr + data->pw) {
			res = data->i2;
		}
		else if (t < data->td + data->tr + data->pw + data->tf) {
			// i2->i1 linearly
			/*x2 = data->td + data->tr + data->pw + data->tf;*/
			x1 = data->td + data->tr + data->pw;
			y2 = data->i1;
			y1 = data->i2;
			res = ((y2 - y1)/data->tf) * (t - x1) + y1;
		}
		else {
			res = data->i1;
		}
	}

	return res;
}



// get the value of pwl transient function at time t
double get_pwl_val(pwl_spec *data, double t) {
	double res;
	double y2, y1;
	double x2, x1;
	int idx;

	// negative time -> return zero
	if (t < 0) {
		printf("Transient Spec Function: Got negative time. Returning zero.\n");
		res = 0.0;
	}

	// only one pair is equivalent to a line parallel to x axis (or a constant DC value)
	if (data->pairs == 1)
		return data->i[0];

	// t is smaller than the first pair time (return first value)
	if (t < data->t[0])
		return data->i[0];

	// t is bigger that the last pair time (return the last value)
	if (t > data->t[data->pairs-1])
		return data->i[data->pairs-1];


	// find the very FIRST pair with a time bigger than t
	idx = 0;
	while (1) {
		idx++;
		if (data->t[idx] > t)
			break;
	}

	x2 = data->t[idx];
	x1 = data->t[idx-1];
	y2 = data->i[idx];
	y1 = data->i[idx-1];
	res = ((y2-y1)/(x2-x1))*(t - x1) + y1;

	return res;
}

