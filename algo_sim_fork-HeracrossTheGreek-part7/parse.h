#ifndef PARSE_H
#define PARSE_H

#include "structs.h"

#define LINE_MAX 500
#define TRUE 1
#define FALSE 0

#define NOT_FOUND -1
#define PARSING_ERROR -1
#define PARSING_SUCCESSFUL 1
#define ISNOT_BLANK 0
#define IS_BLANK 1

#define LU_SOL 0
#define CHOL_SOL 1
#define BICG_SOL 2
#define CG_SOL 3
#define SPARSE_LU_SOL 4
#define SPARSE_CHOL_SOL 5
#define SPARSE_BICG_SOL 6
#define SPARSE_CG_SOL 7

#define TR 0
#define BE 1

#define DELIMITERS " /t/n/r"

extern char valid_comp_list[16];
extern int nodes_n; // circuit nodes INCLUDING ground
extern int m2;      // circuit branches that correspond to voltage sources or inductors
extern int m2_i;    // Index for circuit branches of type V or L
extern int lines;   // The number of lines in the input file
extern int solver_type; // The type of solver (LU or Chol)
extern double DC_arguments[3];  // The three arguments to .dc: Low, High, Step
extern int dc_sweep_flag; // When this flag is raised, do dc_sweep
extern int tran_sweep_flag; // When this flag is raised, do tran_sweep
extern float itol; // tolerance as exit condition for iterative solver
extern int nonzeros_A; // Number of nonzeros, used by sparse structures
extern int nonzeros_C; // Number of nonzeros, used by sparse structures, for array C(transient)
extern char circuit_name[30];   // Name of the circuit, based on the file name
extern int spd_flag;    // A   flag that enables/disables SPD methods
extern int A_dim;       // The dimension of the square A array, equals n-1+m2
extern double tran_time_step; // time_step for .tran
extern double tran_fin_time;  // fin_time for .tran
extern int tran_method; // TR (trapezoidal) or BE (backward Euler)
extern double start_freq; // start frequency for .ac
extern double end_freq;   // end frequency for .ac
extern int ac_sweep_flag; // When this flag is raised, do ac_sweep
extern int points; // Number of points for .ac
extern int point_num; // Current point number
extern double step; // Step for .ac
extern int ac_sweep_method; // Method for ac: 0 for linear, 1 for logarithmic
extern double* sweep_points; // Array of points for .ac



void number_of_lines(char *file_name);
void parse(char *file_name);
int parse_line(char *line);
int valid_comp_type(char comp_type);
int is_blank_line(const char *line);
char* str_tolower(char *str);
int is_ground(char *node);
int check_for_V_or_L(char comp_type);
int parse_plot_arg(char *token);
int parse_spice_command(char* token);
int option_command(char* token);
int dc_command(char* token);
int plot_command(char* token) ;
int tran_command(char *token);
int ac_command(char *token);
int increment_nonzeros_A(char comp_type, char* positive_node, char* negative_node);
int increment_nonzeros_C(char comp_type, char* positive_node, char* negative_node);
transient_spec_type parse_spec_type(char *token);

#endif