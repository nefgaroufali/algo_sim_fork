

char** node_array;
int index_counter = 0;

void add_node_array(char* node_name)
{
    // if it is the first time of adding a node allocate the array
    if (node_array == NULL)
    {
        node_array = (char**) malloc(sizeof(char*));
        node_array[0] = (char*) malloc(strlen(node_name) + 1);
        strcpy(node_array[0], node_name);
        return;
    }
    else
    {
        // realloc the array and add the new node
        node_array = (char**) realloc(node_array, sizeof(char*) * (index_counter + 1));
        node_array[index_counter] = (char*) malloc(strlen(node_name) + 1);
        strcpy(node_array[index_counter], node_name);
        index_counter++;
    }
}

print_node_array()
{
    for (int i = 0; i < index_counter; i++)
    {
        printf("%s\n", node_array[i]);
    }
}

main (int argc, char* argv[]){
}