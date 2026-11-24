# Makefile for PV-LAND model
# Fortran 90 compilation

# Compiler
#FC = gfortran
FC = ifx

# Compiler flags
#FFLAGS = -O2 -fdefault-real-8 -fdefault-double-8 -ffree-form -ffree-line-length-none -fno-align-commons
#FFLAGS = -g -fdefault-real-8 -fdefault-double-8 -ffree-form -ffree-line-length-none -fno-align-commons
FFLAGS = -g -traceback -check bounds

# Linker flags
# ifx : LAPACK/BLAS fournis par Intel MKL (gfortran : -llapack -lblas)
MKLROOT ?= /opt/intel/oneapi/mkl/latest
INTEL_LIB ?= /opt/intel/oneapi/compiler/latest/lib
LDFLAGS = -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(INTEL_LIB) -qmkl=sequential

# Executable name
EXEC = pv_land.exe

# Source files (order matters for dependencies!)
SOURCES = \
	view_factors.f90 \
	solar_radiation.f90 \
	radiosity_SW.f90 \
	radiosity_LW.f90 \
	profil.f90 \
	pv_thermal.f90 \
	cdh_stability.f90 \
	dragcoef_sfcflx.f90 \
	soilebm3d.f90 \
	pv_convection.f90 \
	compute_all_forcings.f90 \
	pv_land.f90 \
	pv_land_main.f90 \
	thermolib3d.f90 \
	tridiag2.f90

# Object files
OBJECTS = $(SOURCES:.f90=.o)

# Include files needed
INCLUDES = \
	timelvls.inc \
	globcst.inc \
	phycst.inc \
	sfcphycst.inc \
	grid.inc \
	soilcst.inc

# Default target
all: check_includes $(EXEC)

# Check if include files exist
check_includes:
	@echo "Checking for required include files..."
	@for inc in $(INCLUDES); do \
		if [ ! -f $$inc ]; then \
			echo "WARNING: Missing include file: $$inc"; \
		fi; \
	done
	@echo ""

# Link executable
$(EXEC): $(OBJECTS)
	@echo "Linking executable..."
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)
	@echo "Build successful! Executable: $(EXEC)"
	@echo ""

# Compile Fortran files
%.o: %.f90
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -c $<

# Clean build files
clean:
	@echo "Cleaning build files..."
	rm -f *.o *.mod $(EXEC)
	@echo "Clean complete."

# Clean and rebuild
rebuild: clean all

# Run the program
run: $(EXEC)
	@echo "Running PV-LAND model..."
	./$(EXEC)

# Help
help:
	@echo "PV-LAND Makefile"
	@echo "----------------"
	@echo "Targets:"
	@echo "  all (default) - Build the executable"
	@echo "  clean         - Remove object files and executable"
	@echo "  rebuild       - Clean and rebuild"
	@echo "  run           - Build and run the model"
	@echo "  help          - Show this help message"
	@echo ""
	@echo "Required files:"
	@echo "  Source files: $(SOURCES)"
	@echo "  Include files: $(INCLUDES)"
	@echo ""
	@echo "Note: You need the .inc files from ARPS!"

# Dependencies (important for correct compilation order)
view_factors.o: 
solar_radiation.o: view_factors.o
radiosity_SW.o: 
radiosity_LW.o: 
profil.o: 
pv_thermal.o: 
cdh_stability.o:
dragcoef_sfcflx.o: 
soilebm3d.o: thermolib3d.o tridiag2.o
pv_convection.o: dragcoef_sfcflx.o
compute_all_forcings.o: view_factors.o solar_radiation.o radiosity_SW.o radiosity_LW.o \
                        soilebm3d.o dragcoef_sfcflx.o pv_thermal.o thermolib3d.o
pv_land.o: compute_all_forcings.o profil.o pv_thermal.o pv_convection.o dragcoef_sfcflx.o \
						cdh_stability.o thermolib3d.o tridiag2.o
pv_land_main.o: pv_land.o

.PHONY: all clean rebuild run help check_includes
