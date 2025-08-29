CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -pthread
TARGET = D2DGRter

SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)

TESTCASES := $(notdir $(wildcard publicCase/*))

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "Compiling $(TARGET)..."
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@echo "Cleaning up..."
	rm -f $(OBJS) $(TARGET)
	find publicCase -name "*.lg" -delete

.PHONY: $(TESTCASES)
$(TESTCASES): $(TARGET)
	@echo "=== Running $@ ==="
	@start_time=$$(date +%s%N); \
	./$(TARGET) \
	publicCase/$@/$@.gmp \
	publicCase/$@/$@.gcl \
	publicCase/$@/$@.cst \
	publicCase/$@/$@.lg; \
	end_time=$$(date +%s%N); \
	duration=$$(echo "scale=9; ($$end_time - $$start_time) / 1000000000" | bc); \
	printf "Execution time: %.9f seconds\n" $$duration

run: $(TARGET) $(TESTCASES)

.PHONY: all clean run