# Compiler settings - Can be customized.
CXX = g++
CXXFLAGS = -Wall -std=c++17

# Define directories for headers and source files
LIB_DIR = ./lib
MAT_DIR = ./materials
SRC_DIR = .

# Include directories for header files
INCLUDES = -I$(LIB_DIR) -I$(SRC_DIR) -I$(MAT_DIR)

# Source files
SOURCES = $(SRC_DIR)/main.cpp $(SRC_DIR)/loader.cpp $(SRC_DIR)/scene.cpp 

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Define the output binary
BIN = main

# Default target
all: $(BIN)

$(BIN): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^

# To obtain object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	rm -f $(SRC_DIR)/*.o $(BIN)

.PHONY: all clean
