# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -g -std=c11

# Include directories
INCLUDES = -I/home/makaragiannis/gsl/include

# Libraries
LIBS = -L/home/makaragiannis/gsl/lib -lgsl -lgslcblas -lm

# Source files
SRC = main.c parse.c structs.c mna.c direct_sol.c

# Object files

# Executable name
EXEC = spice_sim

# Main target
all: $(EXEC)

# Compile source files into object files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# Link object files to create the executable
$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(EXEC) $(OBJ) $(LIBS)

# Clean up object files and executables
clean:
	rm -f $(OBJ) $(EXEC)
