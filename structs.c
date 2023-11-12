#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"

component* head = NULL;
component* tail = NULL;
int index_counter = 1;
hash_table node_hash_table;


// ------- Component Structure ------ //


// This function creates a new component
component* create_component(char comp_type, const char* comp_name, const char* positive_node, const char* negative_node, double value) {

    component* new_node = (component*) malloc(sizeof(component));
    if (new_node == NULL) {
        printf("Error! Not enough memory.\n");
        exit(1);
    }

    new_node->comp_type = comp_type;

    new_node->comp_name = (char*) malloc(strlen(comp_name) + 1);
    new_node->positive_node = (char*) malloc(strlen(positive_node) + 1);
    new_node->negative_node = (char*) malloc(strlen(negative_node) + 1);

    if (new_node->comp_name == NULL || new_node->positive_node == NULL || new_node->negative_node == NULL) {
        printf("Error! Not enough memory.\n");
        exit(1);
    }

    strcpy(new_node->comp_name, comp_name);
    strcpy(new_node->positive_node, positive_node);
    strcpy(new_node->negative_node, negative_node);
    new_node->value = value;
    new_node->next = NULL;

    return new_node;
}

// This functions appends a component to the end of the list
void append_component(component** head, component **tail, char comp_type, const char* comp_name, const char* positive_node, const char* negative_node, double value) {

    component* new_node = create_component(comp_type, comp_name, positive_node, negative_node, value);

    // If the head is NULL, the lsit is empty, so we add it in both ends. Else, we only update the tail
	if (*head == NULL) {
        *head = new_node;
        *tail = new_node;
    } else {
        (*tail)->next = new_node;
        *tail = new_node;
    }
}

// This function prints all the components in the list
void print_comp_list(component* head) {

    component* current = head;

    while (current != NULL) {
        printf("Type: %c, Name: %s, Pos Node: %s, Neg Node: %s, Value: %lf\n",
               current->comp_type, current->comp_name, current->positive_node, current->negative_node, current->value);
        current = current->next;
    }
}

// This function frees the memory used in the list
void free_comp_list(component* head) {

    component* current = head;

    while (current != NULL) {
        component* temp = current;
        current = current->next;
        free(temp->comp_name);
        free(temp->positive_node);
        free(temp->negative_node);
        free(temp);
    }
}


// ------- Hash Table Structure ------ //

// djb2 hash function
unsigned long hash(char *str) {
    int hash = 5381;
    int c;
    while ((c = *str++)) {
        hash = ((hash << 5) + hash) + c;
    }
    return (hash);
}

// This function creates a new hash node.
hash_node *create_hash_node(char *node_str) {
    hash_node *new_node = (hash_node *) malloc(sizeof(hash_node));
    if (new_node == NULL) {
        printf("Eror! Not enough memory.\n");
        exit(EXIT_FAILURE);
    }

    new_node->node_str = malloc(strlen(node_str) + 1);
    strcpy(new_node->node_str, node_str);

    // Node index: The order in which each node is added. Irrelevant to the hash key.
    // Reserve the node index value of "0" for the ground node
    if(strcmp(new_node->node_str,"0")==0){
        new_node->node_index = 0;
        new_node->next = NULL;

        return new_node;    
    }
    new_node->node_index = index_counter;

    new_node->next = NULL;
    index_counter++;

    return new_node;
}

// This function inserts a string into the hash table
void insert_node(hash_table *ht, char *node_str) {
    
    int index = hash((char *)node_str)%HASH_TABLE_SIZE;

    hash_node *new_node = create_hash_node(node_str);

    // Insert the node at the beginning of the linked list (separate chaining)
    new_node->next = ht->table[index];
    ht->table[index] = new_node;
}

// This function checks if a string is in the hash table
int find_hash_node(hash_table *ht, char *node_str) {

    int index = hash((char *)node_str)%HASH_TABLE_SIZE;

    // Search the linked list at the calculated index
    hash_node *current = ht->table[index];

    while (current != NULL) {
        if (strcmp(node_str, current->node_str) == 0) {
            return current->node_index; // Found
        }
        current = current->next;
    }

    return NOT_FOUND; // Not found
}

// This function prints the hash table elements
void print_hash_table(hash_table *ht) {

    for (int i = 0; i < HASH_TABLE_SIZE; i++) {

        hash_node *current = ht->table[i];

        while (current != NULL) {
            printf("Hash index %d: %s with node index %d\n", i, current->node_str, current->node_index);
            current = current->next;
        }
    }
}

// This function frees the nodes in the hash table
void free_hash_table(hash_table *ht) {

    for (int i = 0; i < HASH_TABLE_SIZE; i++) {

        hash_node *current = ht->table[i];
        while (current != NULL) {

            hash_node *temp = current;
            current = current->next;
            free(temp->node_str);
            free(temp);
        }
        ht->table[i] = NULL; 
    }
}
