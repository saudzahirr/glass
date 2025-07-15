# Makefile for Gauss-Seidel solvers (Fortran + C++)

# Compiler settings
FC      = f77
CXX     = g++
FFLAGS  = -O2 -c
CXXFLAGS = -O2 -std=c++11
LDFLAGS = -lgfortran

# Sources and object files
FORTRAN_SOURCES = sgssv.f dgssv.f cgssv.f zgssv.f
F77_OBJS = $(FORTRAN_SOURCES:.f=.o)

CPP_SRC  = test.cpp
CPP_OBJ  = test.o

TARGET = solver.out

# Default rule
all: $(TARGET)

# Link final executable
$(TARGET): $(F77_OBJS) $(CPP_OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Compile Fortran sources
%.o: %.f
	$(FC) $(FFLAGS) $< -o $@

# Compile C++ test file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -f *.o $(TARGET)

# Useful phony targets
.PHONY: all clean
