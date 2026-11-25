# Makefile for the PV-LAND offline model
#
# Source layout (see README.md):
#   core/      PV-LAND physics, exact copies of ARPS src/arps/pvland*.f90 (do not edit)
#   arps_lib/  ARPS routines and include files used by PV-LAND (do not edit)
#   offline/   offline main program (CSV forcing, outputs) and park map
#   scripts/   sync_from_arps.py (copy/generate from ARPS, or check)
#
# Targets:
#   make              build pv_land.exe
#   make run          run with the full forcing
#   make run-short    run with the 5-day forcing (~1 s)
#   make check-arps   check that core/ and arps_lib/ are identical to ARPS
#                     ($(ARPS_DIR))
#   make sync-arps    copy them again from ARPS
#   make clean / rebuild
#
# Reproducible build (to compare two versions bit for bit):
#   make repro                    build pv_land_source.exe (objects in build_source/)
#   make run FPMODEL=source       same, and run it (any target accepts FPMODEL)
# FPMODEL=source compiles with -fp-model source, as ARPS: floating-point
# expressions are evaluated as written, so a pure code reorganisation gives
# bit-identical outputs (about 20 % slower). The default build lets the
# compiler reorder floating-point operations (faster, last-digit changes).

# Compiler
#FC = gfortran
FC = ifx

# Compiler flags
#FFLAGS = -O2 -fdefault-real-8 -fdefault-double-8 -ffree-form -ffree-line-length-none -fno-align-commons
#FFLAGS = -g -fdefault-real-8 -fdefault-double-8 -ffree-form -ffree-line-length-none -fno-align-commons
# ifx : debug (lent)
#FFLAGS = -g -traceback -check bounds
# ifx : optimisé (meme precision que la version debug : REAL simple precision)
FFLAGS = -O2 -traceback

# Floating-point model: empty = compiler default (fast), source = reproducible
FPMODEL ?=

# Linker flags
# ifx : LAPACK/BLAS fournis par Intel MKL (gfortran : -llapack -lblas)
MKLROOT ?= /opt/intel/oneapi/mkl/latest
INTEL_LIB ?= /opt/intel/oneapi/compiler/latest/lib
LDFLAGS = -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(INTEL_LIB) -qmkl=sequential

# ARPS PV-LAND directory (reference for core/ and arps_lib/)
ARPS_DIR ?= ../arps5.3.4_PVLAND_v2
PYTHON   ?= python3

# Executable name and build directory (objects and .mod files); a
# separate executable and build directory for each floating-point model
ifeq ($(FPMODEL),)
EXEC     = pv_land.exe
BUILDDIR = build
else
FFLAGS  += -fp-model $(FPMODEL)
EXEC     = pv_land_$(FPMODEL).exe
BUILDDIR = build_$(FPMODEL)
endif

# Include files (ARPS common blocks)
INCFLAGS = -Iarps_lib/include

# Module file location (ifx: -module, gfortran: -J)
MODFLAG  = -module $(BUILDDIR)
#MODFLAG = -J$(BUILDDIR)

# Source files (module_pvland first: other files USE its modules)
CORE_SRC = module_pvland.f90 \
           pvland_profil.f90 \
           pvland_view_factors.f90 \
           pvland_solar_radiation.f90 \
           pvland_radiosity_SW.f90 \
           pvland_radiosity_LW.f90 \
           pvland.f90 \
           pvland_driver.f90 \
           pvland_obukhov.f90 \
           pvland_soil.f90 \
           pvland_dragcoef.f90
LIB_SRC  = sfcphy3d_extract.f90 \
           tmix3d_extract.f90 \
           thermolib3d.f90 \
           soilebm3d.f90
OFF_SRC  = pvland_park_map.f90 \
           pvland_offline_main.f90

SOURCES  = $(CORE_SRC) $(LIB_SRC) $(OFF_SRC)
OBJECTS  = $(addprefix $(BUILDDIR)/,$(SOURCES:.f90=.o))

vpath %.f90 core arps_lib offline

# Default target
all: $(EXEC)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Link executable
$(EXEC): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)
	@echo "Build successful! Executable: $(EXEC)"

# Compile Fortran files
$(BUILDDIR)/%.o: %.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) $(INCFLAGS) $(MODFLAG) -c $< -o $@

# Every PV-LAND file USEs the modules of module_pvland.f90
$(filter-out $(BUILDDIR)/module_pvland.o,$(OBJECTS)): $(BUILDDIR)/module_pvland.o

# Include files: rebuild everything when they change
$(OBJECTS): $(wildcard arps_lib/include/*.inc)

# Reproducible build (-fp-model source)
repro:
	$(MAKE) FPMODEL=source

# Clean build files (all floating-point models)
clean:
	rm -rf build build_* pv_land.exe pv_land_*.exe

rebuild: clean all

# Run the model
run: $(EXEC)
	./$(EXEC)

run-short: $(EXEC)
	./$(EXEC) pvland.nml recorded_data_2025_short_5days.csv

# Consistency with ARPS
check-arps:
	$(PYTHON) -I scripts/sync_from_arps.py --check --arps $(ARPS_DIR)

sync-arps:
	$(PYTHON) -I scripts/sync_from_arps.py --arps $(ARPS_DIR)

help:
	@sed -n '1,/^$$/p' Makefile | sed -n '/^#/p'

.PHONY: all repro clean rebuild run run-short check-arps sync-arps help
