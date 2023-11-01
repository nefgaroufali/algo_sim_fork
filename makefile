# Compiler
CC = gcc -ggdb

# Compiler flags
CFLAGS = -Wall  -std=c11

# Source files
SRC = main.c parse.c structs.c

# Object files
OBJ = $(SRC:.c=.o)

# Executable name
EXEC = spice_sim

# Include directories
INCLUDES = -I.

# Libraries
LIBS =

# Main target
all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(EXEC) $(OBJ) $(LIBS)

# Compile source files into object files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# Clean up object files and executable
clean:
	rm -f $(OBJ) $(EXEC)
