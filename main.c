#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "parse.h"
#include "structs.h"
#include "mna.h"
#include "direct_sol.h"

int main(int argc, char* argv[]) {

    char file_name[40];

    // If there is no argument (or multiple arguments) exit
    if (argc!=2) {
        printf("Error! No file specified. Use ./spice_sim <file_name>\n");
        return 0;
    }

    strcpy(file_name, argv[1]);
    parse(file_name);

    print_hash_table(&node_hash_table);
    // print_comp_list(head);

    // printf("NUmber of nodes (INCLUDING ground) is %d and number of V or L branches is %d\n", nodes_n, m2);

    create_equations();

    copy_to_A();
    form_LU();
    //form_chol();

        free_gsl();

    free_hash_table(&node_hash_table);
    free_comp_list(head);
    free_A_array();
    free_b_array();
    
    return 0;

}