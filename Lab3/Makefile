# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

# Target executable
TARGET = Legalizer

# Source files
SRCS = Lab3.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default rule
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to build object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run the program with testcase
run: $(TARGET)
	@echo "Running testcase..."
	@for lgfile in testcase/*.lg; do \
		optfile="$${lgfile%.lg}.opt"; \
		outfile="$${lgfile%.lg}_post.lg"; \
		if [ -f "$$optfile" ]; then \
			echo "Processing $$lgfile with $$optfile -> $$outfile"; \
			./$(TARGET) "$$lgfile" "$$optfile" "$$outfile"; \
		else \
			echo "Warning: No corresponding .opt file for $$lgfile"; \
		fi \
	done

# Clean up build files and generated output
clean:
	rm -f $(TARGET) $(OBJS) testcase/*_post.lg

# Phony targets
.PHONY: all clean run
