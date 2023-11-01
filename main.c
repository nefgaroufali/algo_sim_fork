#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "parse.h"
#include "structs.h"

int main(int argc, char* argv[]) {

    char file_name[20];

    // If there is no argument (or multiple arguments) exit
    if (argc!=2) {
        printf("Error! No file specified. Use ./spice_sim <file_name>\n");
        return 0;
    }

    strcpy(file_name, argv[1]);
    parse(file_name);

    print_hash_table(&node_hash_table);

    free_hash_table(&node_hash_table);
    free_comp_list(head);
    
    return 0;

}