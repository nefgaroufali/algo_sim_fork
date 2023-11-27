#ifndef STRUCTS_H
#define STRUCTS_H

#define HASH_TABLE_SIZE 500
#define NOT_FOUND -1
#define FOUND 1
#define NOT_V_OR_I 0

// The structure for a component
typedef struct component {
    char comp_type;
    char* comp_name;
    char* positive_node;
    char* negative_node;
    double value;
    int m2_i;
    struct component* next;
} component;

// The structure for a circuit node
typedef struct hash_node {
    char *node_str;
    int node_index;
    struct hash_node *next;
} hash_node;

// The structure for the hash table
typedef struct {
    hash_node **table;
    int size;
} hash_table;

extern component* head;
extern component* tail;
extern component *sweep_component;
extern int index_counter;
extern hash_table node_hash_table;

extern char** node_array;
extern int node_array_index_counter;
extern char** m2_array;

// Function declarations
component* create_component(char comp_type, const char* comp_name, const char* positive_node, const char* negative_node, double value);
void append_component(component** head, component **tail, char comp_type, const char* comp_name, const char* positive_node, const char* negative_node, double value);
void print_comp_list(component* head);
void free_comp_list(component* head);
int find_component(const char* comp_name);

hash_node *create_hash_node(char *node_str);
void create_hash_table(int size);
void insert_node(hash_table *ht, char *node_str);
int find_hash_node(hash_table *ht, char *node_str);
void add_node_array(char* node_name);
void print_node_array();
void add_m2_array(char* comp_name);
void print_m2_array();

void print_hash_table(hash_table *ht);
void free_hash_table(hash_table *ht);

#endif
