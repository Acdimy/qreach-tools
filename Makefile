# Project Name (executable)
PROJECT = libqreach.so
# Test Project Name (executable)
TEST = test_qreach
# Compiler
CC = g++

# Run Options       
COMMANDLINE_OPTIONS = #/dev/ttyS0

# Compiler options during compilation
COMPILE_OPTIONS = -g -O3 -std=c++2a -w -Wall -Wextra -DHAVE_CONFIG_H -Werror -Wunused-but-set-variable -fPIC
# -ansi -pedantic -Wall 

#Header include directories
HEADERS = -I. -I $(BOOST_PATH) -I.cflobdd/CFLOBDD -I.cflobdd/CFLOBDD/Solver/uwr/bit_vector/ -I.cflobdd/CFLOBDD/Solver/uwr/assert/ -I.cflobdd/CFLOBDD/Solver/uwr/matrix/ -I.cflobdd/CFLOBDD/Solver/uwr/parsing/

# Dependency options
DEPENDENCY_OPTIONS = -MM


# Subdirs to search for additional source files
SOURCE_FILES := $(shell ls *.cpp)
SOURCE_FILES += $(shell ls cflobdd/CFLOBDD/Solver/uwr/bit_vector/*.cpp)
SOURCE_FILES += $(shell ls cflobdd/CFLOBDD/Solver/uwr/parsing/*.cpp)
SOURCE_FILES += $(shell find cflobdd/CFLOBDD -maxdepth 1 -mindepth 1 -name \*.cpp -a -not -name main.cpp)
SOURCE_FILES := $(filter-out quantum_circuit.cpp, $(SOURCE_FILES))
# [SOURCE_FILES] -= $(shell ls cflobdd/CFLOBDD/main.cpp)
# $(info $(SOURCE_FILES))

# Create an object file of every cpp file
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCE_FILES))

# Dependencies
DEPENDENCIES = $(patsubst %.cpp, %.d, $(SOURCE_FILES))


# Make $(PROJECT) the default target
all: $(PROJECT)
#$(DEPENDENCIES) -shared 
$(PROJECT): $(OBJECTS)
	$(CC) -shared -o $(PROJECT) $(OBJECTS) -fPIC

# Include dependencies (if there are any)
# ifneq "$(strip $(DEPENDENCIES))" ""
#   include $(DEPENDENCIES)
# endif

$(TEST): $(OBJECTS)
	$(CC) -o $(TEST) $(OBJECTS)

test: $(TEST)
	@echo "Target $(TEST) is built."

# Compile every cpp file to an object
# %.cpp 
%.o: %.cpp
	$(CC) -c $(COMPILE_OPTIONS) -o $@ $^ $(HEADERS)

# Build & Run Project
run: $(PROJECT)
	./$(PROJECT) $(COMMANDLINE_OPTIONS)

# Clean & Debug
.PHONY: makefile-debug
makefile-debug:

.PHONY: clean
clean:
	rm -f $(PROJECT) $(OBJECTS)

.PHONY: depclean
depclean:
	rm -f $(DEPENDENCIES)

clean-all: clean depclean

# -include $(DEPENDENCIES)
