# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -std=c99 -O2

# Linker flags
LDFLAGS = -I . -lm

# Executable name
TARGET = minlexrot

# Source files
SRC = generalised_min_lex_rot.c minlexrot.c

# Default target
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o $@

# Clean up compiled file
clean:
	rm -f $(TARGET)
