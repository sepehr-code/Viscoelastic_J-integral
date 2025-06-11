# Makefile for Viscoelastic J-Integral Calculator
# ==============================================

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g
LIBS = -lm

# Project name
TARGET = viscoelastic_j_integral

# Source files
SOURCES = main.c \
          memory_management.c \
          contour_operations.c \
          field_initialization.c \
          viscoelastic_model.c \
          j_integral_calculation.c \
          time_integration.c \
          utilities.c

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Header files
HEADERS = viscoelastic_j_integral.h

# Default target
all: $(TARGET)

# Build the main executable
$(TARGET): $(OBJECTS)
	@echo "Linking $(TARGET)..."
	$(CC) $(OBJECTS) -o $(TARGET) $(LIBS)
	@echo "Build successful! Run with: ./$(TARGET)"

# Compile object files
%.o: %.c $(HEADERS)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(OBJECTS) $(TARGET)
	rm -f *.dat *.txt  # Remove generated output files

# Clean everything including output files
distclean: clean
	@echo "Cleaning all generated files..."
	rm -f demo*.dat demo*.txt *.report

# Install (copy to /usr/local/bin)
install: $(TARGET)
	@echo "Installing $(TARGET) to /usr/local/bin..."
	sudo cp $(TARGET) /usr/local/bin/

# Uninstall
uninstall:
	@echo "Removing $(TARGET) from /usr/local/bin..."
	sudo rm -f /usr/local/bin/$(TARGET)

# Run demonstrations
demo: $(TARGET)
	@echo "Running all demonstrations..."
	./$(TARGET) 0

# Run specific demo (usage: make demo1, demo2, etc.)
demo1: $(TARGET)
	./$(TARGET) 1

demo2: $(TARGET)
	./$(TARGET) 2

demo3: $(TARGET)
	./$(TARGET) 3

demo4: $(TARGET)
	./$(TARGET) 4

demo5: $(TARGET)
	./$(TARGET) 5

# Debug build
debug: CFLAGS += -DDEBUG -g3 -O0
debug: $(TARGET)

# Release build
release: CFLAGS += -DNDEBUG -O3
release: clean $(TARGET)

# Check for memory leaks with valgrind
memcheck: $(TARGET)
	@echo "Running memory check with valgrind..."
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./$(TARGET) 1

# Static analysis with cppcheck
check:
	@echo "Running static analysis..."
	cppcheck --enable=all --std=c99 --suppress=missingIncludeSystem $(SOURCES)

# Count lines of code
count:
	@echo "Lines of code:"
	wc -l $(SOURCES) $(HEADERS)

# Show help
help:
	@echo "Available targets:"
	@echo "  all       - Build the main executable (default)"
	@echo "  clean     - Remove build artifacts"
	@echo "  distclean - Remove all generated files"
	@echo "  install   - Install to /usr/local/bin"
	@echo "  uninstall - Remove from /usr/local/bin"
	@echo "  demo      - Run all demonstrations"
	@echo "  demo[1-5] - Run specific demonstration"
	@echo "  debug     - Build debug version"
	@echo "  release   - Build optimized version"
	@echo "  memcheck  - Run with valgrind memory check"
	@echo "  check     - Run static analysis"
	@echo "  count     - Count lines of code"
	@echo "  help      - Show this help message"

# Mark phony targets
.PHONY: all clean distclean install uninstall demo demo1 demo2 demo3 demo4 demo5 debug release memcheck check count help

# Dependencies (automatically generated)
-include $(OBJECTS:.o=.d)

# Generate dependency files
%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$ 