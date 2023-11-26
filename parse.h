#ifndef PARSE_H
#define PARSE_H

#define LINE_MAX 500
#define TRUE 1
#define FALSE 0

#define NOT_FOUND -1
#define PARSING_ERROR -1
#define PARSING_SUCCESSFUL 1
#define ISNOT_BLANK 0
#define IS_BLANK 1

extern char valid_comp_list[16];
extern int nodes_n; // circuit nodes INCLUDING ground
extern int m2;      // circuit branches that correspond to voltage sources or inductors
extern int m2_i;    // index for V or L components
extern int lines;
extern int solver_type;

extern double DC_arguments[3];

void number_of_lines(char *file_name);
void parse(char *file_name);
int parse_line(char *line);
int valid_comp_type(char comp_type);
int is_blank_line(const char *line);
char* str_tolower(char *str);
int isnot_ground(char *node);
int check_for_V_or_L(char comp_type);
int parse_plot_arg(char *token);

#endif