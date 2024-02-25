#include "ac_utility.h"


// The b-vector in the AC analysis is always the same, so it must only be created once

void fill_ac_b_vector(double complex* ac_b_vector) {

    component* current = head;
    int m2_count_ac = 0;

    // For each component of the circuit fill the b vector
    while (current != NULL) {
        
        if(current->comp_type == 'r'){
            fill_ac_bvector_with_r(ac_b_vector, current);
        }
        else if(current->comp_type == 'i'){
            fill_ac_bvector_with_i(ac_b_vector, current);
        }
        else if(current->comp_type == 'c'){
            fill_ac_bvector_with_c(ac_b_vector, current);
        }
        else if(current->comp_type == 'v'){
            fill_ac_bvector_with_v(ac_b_vector, current, &m2_count_ac);
        }
        else if(current->comp_type == 'l'){
            fill_ac_bvector_with_l(ac_b_vector, current, &m2_count_ac);
        }

        current = current->next;
    }


}

double complex* alloc_ac_b_vector(){

    // Allocate memory for the array of pointers to rows
    double complex* new_vector = (double complex*)calloc((A_dim), sizeof(double complex));

    // Check if memory allocation was successful
    if (new_vector == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(0);
    }

    // Return a pointer to the created array
    return new_vector;

}

void fill_ac_bvector_with_r(double complex* ac_b_vector, component *current) {

}

void fill_ac_bvector_with_i(double complex* ac_b_vector, component *current) {

    int positive_node_i;
    int negative_node_i;

    // Take the node_index of each node
    positive_node_i = find_hash_node(&node_hash_table, str_tolower(current->positive_node));
    negative_node_i = find_hash_node(&node_hash_table, str_tolower(current->negative_node));

    // If the I source is not AC, the value is 0
    if (current->spec_type != AC_SPEC) {
        return;
    }

    double complex current_i = convert_to_complex(current->spec->ac.mag, current->spec->ac.phase);

    // Case: i 0 x value
    if(positive_node_i == 0){
        ac_b_vector[negative_node_i-1] += current_i;
    }
    // Case: i x 0 value
    else if(negative_node_i == 0) {
        ac_b_vector[positive_node_i-1] -= current_i;
    }
    // Case: i x y value
    else{
        ac_b_vector[negative_node_i-1] += current_i;
        ac_b_vector[positive_node_i-1] -= current_i;
    }

    return;
}

void fill_ac_bvector_with_c(double complex* ac_b_vector, component *current) {
    
}

void fill_ac_bvector_with_v(double complex* ac_b_vector, component *current, int *m2_count_ac) {

    int m2_i = nodes_n - 1 + (*m2_count_ac);

    // If the I source is not AC, the value is 0
    if (current->spec_type != AC_SPEC) {
        return;
    }

    double complex current_v = convert_to_complex(current->spec->ac.mag, current->spec->ac.phase);

    ac_b_vector[m2_i] = current_v;

    (*m2_count_ac)++;
}

void fill_ac_bvector_with_l(double complex* ac_b_vector, component *current, int *m2_count_ac) {
    (*m2_count_ac)++;
}

// This function converts a complex number represented with the magnitude and phase, to its real and imaginary part
double complex convert_to_complex(double mag, double phase) {

    double phase_rad = phase*M_PI/180.0;

    double real = mag*cos(phase_rad);
    double imaginary = mag*sin(phase_rad);

    return real + I*imaginary;

}

// get magnitude (in 20log) and phase
double get_magnitude(double complex complex_num, int ac_sweep_method) {

    double magnitude = sqrt(creal(complex_num)*creal(complex_num) + cimag(complex_num)*cimag(complex_num));
    if (ac_sweep_method == AC_METHOD_LIN) {
        return magnitude;
    }
    else {
        return 20*log10(magnitude);
    }
}

double get_phase(double complex complex_num) {
    
    double phase_rad = atan2(cimag(complex_num), creal(complex_num));

    return phase_rad * 180 / M_PI;

    //return phase_rad;
}

void print_ac_vector(double complex* ac_vector) {

    int i;

    printf("\n");
    for (i=0; i<A_dim; i++) {
        printf("%.2lf + i*%.2lf\n", creal(ac_vector[i]), cimag(ac_vector[i]));
    }
    printf("\n");


}

void create_ac_gnuplot(int i, int ac_sweep_method) {


    // Generate a temporary script file to run the GNU Plot commands
    FILE *gnuplotScript = fopen("plot_script_ac.gnu", "w");
    if (gnuplotScript == NULL) {
        printf("Error creating GNU Plot script file\n");
        return;
    }

    // Write the GNU Plot command to the script file
    fprintf(gnuplotScript, "set terminal png; set output 'output/%s_%d_db.png'\n", circuit_name, i);

    if (ac_sweep_method == AC_METHOD_LOG) fprintf(gnuplotScript, "set logscale x\n");

    fprintf(gnuplotScript, "set xlabel \"Frequency (Hz)\"\n");
    fprintf(gnuplotScript, "set ylabel \"Magnitude (dB))\"\n");
    fprintf(gnuplotScript, "set title \"Bode diagram\"\n\n");
    fprintf(gnuplotScript, "plot 'output/%s_%d.txt' using 1:2 with lines title \"Magnitude\"\n", circuit_name, i);

    fprintf(gnuplotScript, "set terminal png; set output 'output/%s_%d_phase.png'\n", circuit_name, i);

    if (ac_sweep_method == AC_METHOD_LOG) fprintf(gnuplotScript, "set logscale x\n");

    fprintf(gnuplotScript, "set ylabel \"Phase (degrees))\"\n");
    fprintf(gnuplotScript, "set title \"Phase diagram\"\n\n");
    fprintf(gnuplotScript, "plot 'output/%s_%d.txt' using 1:3 with lines title \"Phase\"\n", circuit_name, i);
    
    fclose(gnuplotScript);


}