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

void parse(char *file_name);
int parse_line(char *line);
int valid_comp_type(char comp_type);
int is_blank_line(const char *line);
char* str_tolower(char *str);

#endif