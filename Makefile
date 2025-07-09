# Makefile for FDTD Laser Simulation

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -march=native
DEBUG_FLAGS = -g -O0 -DDEBUG
LDFLAGS = -lm

# Target executable
TARGET = laser_simulation

# Source files
SOURCES = Laser.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Build the main executable
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS)
	@echo "Build complete: $(TARGET)"

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Debug build
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(TARGET)
	@echo "Debug build complete: $(TARGET)"

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *.dat *.txt *.png *.gif
	@echo "Clean complete"

# Clean only simulation output files
clean-output:
	rm -f Ez_field_*.dat geometry_params.txt *.png *.gif
	@echo "Output files cleaned"

# Run the simulation
run: $(TARGET)
	./$(TARGET)

# Run with timing
time: $(TARGET)
	time ./$(TARGET)

# Fast build (optimized for speed)
fast: CXXFLAGS = -std=c++11 -Ofast -march=native -funroll-loops -ffast-math
fast: $(TARGET)
	@echo "Fast build complete: $(TARGET)"

# Profile build
profile: CXXFLAGS += -pg -O2
profile: LDFLAGS += -pg
profile: $(TARGET)
	@echo "Profile build complete: $(TARGET)"

# Install dependencies (placeholder - adjust for your system)
install-deps:
	@echo "Installing dependencies..."
	@echo "Make sure you have g++ compiler installed"
	@echo "On Ubuntu/Debian: sudo apt-get install g++ build-essential"
	@echo "On CentOS/RHEL: sudo yum install gcc-c++ make"
	@echo "On Windows: Install MinGW-w64 or use Visual Studio"

# Help target
help:
	@echo "Available targets:"
	@echo "  all          - Build the simulation (default)"
	@echo "  debug        - Build with debug symbols"
	@echo "  fast         - Build with maximum optimization"
	@echo "  profile      - Build with profiling enabled"
	@echo "  run          - Build and run the simulation"
	@echo "  time         - Build and run with timing"
	@echo "  clean        - Remove all build artifacts and output files"
	@echo "  clean-output - Remove only simulation output files"
	@echo "  install-deps - Show dependency installation instructions"
	@echo "  help         - Show this help message"

# Declare phony targets
.PHONY: all debug fast profile run time clean clean-output install-deps help

# Advanced targets for analysis
check-syntax:
	$(CXX) $(CXXFLAGS) -fsyntax-only $(SOURCES)

# Static analysis (if cppcheck is available)
analyze:
	@if command -v cppcheck >/dev/null 2>&1; then \
		cppcheck --enable=all --std=c++11 $(SOURCES); \
	else \
		echo "cppcheck not found. Install it for static analysis."; \
	fi

# Format code (if clang-format is available)
format:
	@if command -v clang-format >/dev/null 2>&1; then \
		clang-format -i $(SOURCES); \
		echo "Code formatted"; \
	else \
		echo "clang-format not found. Install it for code formatting."; \
	fi
