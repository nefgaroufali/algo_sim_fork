#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "parse.h"
#include "structs.h"

char valid_comp_list[16] = {'V', 'I', 'R', 'C', 'L', 'D', 'M', 'Q', 
                           'v', 'i', 'r', 'c', 'l', 'd', 'm', 'q'};

int nodes_n = 0;
int m2 = 0;

void parse(char* file_name) {

    FILE *input_file;
    char line[LINE_MAX];
    int parse_line_result;

    input_file = fopen(file_name, "r");

    // Check if the file is valid
    if (input_file == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    // Parse every line until end of file
    while (fgets(line, LINE_MAX, input_file) != NULL) {

        if (is_blank_line(line) == IS_BLANK) {
            continue;
        }
        parse_line_result = parse_line(line);

        // If the parsing fails, free the memory and exit.
        if (parse_line_result == PARSING_ERROR) {
            printf("Invalid Syntax!!\n");
            free_comp_list(head);
            free_hash_table(&node_hash_table);
            fclose(input_file);
            exit(1);
        }
    }

    fclose(input_file);
}

int parse_line(char* line) {

    char *token = NULL;
    char comp_type;
    char* comp_name = NULL;
    char* positive_node = NULL;
    char *negative_node = NULL;
    double value;
    int node_found;


    // String 1: Type of circuit element & name
    token = strtok(line, " \t");

    // If the first character is an asterisk, the line is a comment
    if (strcmp(token, "*") == 0) {
        return 1; 
    }

    // Spice command: For now this will be treated as a comment
    if (strcmp(token, ".") == 0) {
        return 1;
    }

    if ((token == NULL) || (valid_comp_type(token[0]) == FALSE)) {
        return PARSING_ERROR;
    }

    // Store the component type and name in lowercase

    comp_type = tolower(token[0]);

    if ((check_for_V_or_L(comp_type)) == TRUE) {
        m2++;
    }

    comp_name = str_tolower(token);

    // String 2: Positive node name
    token = strtok(NULL, " \t");
    if (token == NULL) {
        return PARSING_ERROR;
    }

    positive_node = token;

    // Search  for the node in the hash table. If it doesnt exist, add it in lowercase

    node_found = find_hash_node(&node_hash_table, str_tolower(positive_node));

    if (node_found == NOT_FOUND) {
        insert_node(&node_hash_table, str_tolower(positive_node));
        nodes_n++;
    }  

    // String 3: Negative node name
    token = strtok(NULL, " \t");
    if (token == NULL) {
        return PARSING_ERROR;
    }

    negative_node = token;

    // Search  for the node in the hash table. If it doesnt exist, add it in lowercase

    node_found = find_hash_node(&node_hash_table, str_tolower(negative_node));

    if (node_found == NOT_FOUND) {
        insert_node(&node_hash_table, str_tolower(negative_node));
        nodes_n++;
    } 
 
    // String 4: Numeric value
    token = strtok(NULL, " \t");
    if (token == NULL) {
        return PARSING_ERROR;
    }

    // strtod converts the string to a float
    value = strtod(token, NULL);

    // Add the component to the linked list
    append_component(&head, &tail, comp_type, comp_name, positive_node, negative_node, value);

    return PARSING_SUCCESSFUL;

}

// This function checks if the component type is valid, i.e. if it starts with one of the 8 key letters
int valid_comp_type(char comp_type) {

    int i;
    for (i = 0; i < 16; i++) {
        if (comp_type == valid_comp_list[i]) {
            return TRUE;
        }
    }

    return FALSE;
}

// This function checks if the line is blank
int is_blank_line(const char *line) {

    for (int i = 0; line[i] != '\0'; i++) {
        if (!isspace(line[i])) {
            return ISNOT_BLANK;
        }
    }
    return IS_BLANK;
}

// This function turns a string to lowercase
char* str_tolower(char *str) {

    for (int i = 0; str[i] != '\0'; i++) {
        str[i] = tolower(str[i]);
    }

    return str;
}

// This function checks if the given node is ground or not
int isnot_ground(char *node) {

    if ((node[0] == '0') && (node[1] == '\0')) {
        return 0;
    }

    return 1;
}

// This function checks if the component type is a voltage source or an inductor
int check_for_V_or_L(char comp_type) {

    if (comp_type == 'v' || comp_type == 'l') {
        return 1;
    }

    return 0;

}




