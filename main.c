#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "parse.h"
#include "structs.h"
#include "mna.h"
#include "direct_sol.h"
#include "gsl.h"
#include "iter_sol.h"
#include "sparse_sol.h"

int main(int argc, char* argv[]) {

    char file_name[40];

    // If there is no argument (or multiple arguments) exit
    if (argc!=2) {
        printf("Error! No file specified. Use ./spice_sim <file_name>\n");
        return 0;
    }

    strcpy(file_name, argv[1]);
    parse(file_name);

    // printf("DC Parameters: %f %f %f\n", DC_arguments[0], DC_arguments[1], DC_arguments[2]);

    // print_hash_table(&node_hash_table);
    // print_comp_list(head);

    switch (solver_type) {
        case LU_SOL:
            create_equations();
            form_gsl_system();
            form_LU();
            solve_dc_system(LU_SOL);
            break;

        case CHOL_SOL:
            create_equations();
            form_gsl_system();
            form_chol();
            if (spd == TRUE) {
                solve_dc_system(CHOL_SOL);
            }
            break;

        case CG_SOL:
            create_equations();
            form_gsl_system();
            solve_dc_system(CG_SOL);
            break;

        case BICG_SOL:
            create_equations();
            form_gsl_system();
            solve_dc_system(BICG_SOL);
            break;

        case SPARSE_LU_SOL:
            form_sparse();
            solve_dc_system(SPARSE_LU_SOL);
            // Print the solver type for the extra cases
            printf("Solver Type: %d\n", solver_type);
            break;
        case SPARSE_CHOL_SOL:
            form_sparse();
            solve_dc_system(SPARSE_CHOL_SOL);
            // Print the solver type for the extra cases
            printf("Solver Type: %d\n", solver_type);
            break;
        case SPARSE_CG_SOL:
            form_sparse();
            solve_dc_system(SPARSE_CG_SOL);
            // Print the solver type for the extra cases
            printf("Solver Type: %d\n", solver_type);
            break;
        case SPARSE_BICG_SOL:
            form_sparse();
            solve_dc_system(SPARSE_BICG_SOL);
            // Print the solver type for the extra cases
            printf("Solver Type: %d\n", solver_type);
            break;

        default:
            // Handle any other cases if needed
            break;
    }

    print_arrays();
    print_sparse_arrays();
    if (sweep_flag == TRUE) {
        dc_sweep();
    }


    // Close all open files before exiting the program
    for (int i = 0; i < 5; i++) {
        if (filePointers[i] != NULL) {
            fclose(filePointers[i]);
            filePointers[i] = NULL;
        }
    }

    free_gsl();
    free_hash_table(&node_hash_table);
    free_comp_list(head);
    free_A_array();
    free_b_array();
    free_node_array();
    free_m2_array();
    free_plot_node();   
    
    return 0;

}