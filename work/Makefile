# Some versions of make default FC to f77
ifeq ($(FC),f77)
    FC = gfortran
endif
FC ?= gfortran
FCFLAG=
FPPFLAG=

# for gfortran
ifeq (${FC},gfortran)
ifeq (${MPI},yes)
override FC := mpifort
FPPFLAG += -DMPI
endif
MODFLAG= -J
FCFLAG:= $(CFLAG) -fopenmp -ffree-line-length-512 #-fmax-stack-var-size=32768
DBFLAG= -g -ffpe-trap=invalid,zero,overflow -fbacktrace -O0 -Wall -Wextra -fcheck=all
PRFLAG= -pg -save-temps
endif

# for Intel fortran compiler
ifeq (${FC},ifx)
ifeq (${MPI},yes)
override FC := mpiifx
FPPFLAG += -DMPI
endif
MODFLAG= -module # trailing space is important
FCFLAG:= $(CFLAG) -qopenmp -fp-model precise #-heap-arrays
DBFLAG= -traceback -fpe0 -CB -warn all -O0 -std #-pg #-check all
PRFLAG= -pg
endif

# debug mode
ifeq (${DEBUG},yes)
FCFLAG+=$(DBFLAG)
endif

# profile mode
ifeq (${PROFILE},yes)
FCFLAG+=$(PRFLAG)
endif

# PETSc settings
ifeq (${USE_PETSC},yes)
# Useful as an escape hatch if libpetsc.so was not compiled properly
# and has missing symbols from external libraries.
ifeq ($(PETSC_PC_STATIC), yes)
STATIC = "--static"
else
STATIC =
endif

PETSC_PC = $(shell pkg-config --path petsc || pkg-config --path PETSc || if [ -n "${PETSC_DIR}" ] && [ -d "${PETSC_DIR}" ]; then find "${PETSC_DIR}" -type f -iname petsc.pc 2>/dev/null | head -n 1; fi)
PETSC_INCLUDE ?= $(shell pkg-config --cflags $(PETSC_PC) 2>/dev/null)
PETSC_LDFLAGS ?= $(shell pkg-config $(STATIC) --libs $(PETSC_PC) 2>/dev/null)

FCFLAG += $(PETSC_INCLUDE)
LDFLAGS += $(PETSC_LDFLAGS)

FPPFLAG += -DUSE_PETSC
endif

# directories
OBJ_DIR = obj
SRC_DIR = ../src/main
INI_DIR = ../src/init
POS_DIR = ../src/post
BIN_DIR = $(HORMONE_DIR)/bin

TARGET = hormone

TARGETA = ascii_output

TARGETL = lightcurve

# compile modules that are being depended on first
MOD = mpi_utils.F90 modules.f90 profiler.f90 mpi_domain.F90 io.F90 utils.f90 conserve.f90 ionization.f90 eos.f90 miccg.f90 fluxlimiter.f90 hlldflux.f90 dirichlet.f90 particles.f90 opacity.f90 radiation_utils.f90 radiation.f90 cooling.f90 sinks.f90 externalforce.f90 composition.f90 star.f90 shockfind.f90 output.f90 smear.f90 input.f90 setup.f90 readbin.f90  restart.f90 fluxboundary.f90 rungekutta.f90 timestep.f90 tests.f90 gridset.f90 gravbound.f90 gravity_hyperbolic.f90 petsc_solver.F90 matrix_coeffs.f90 matrix_solver.F90 gravity_miccg.f90 interpolation.f90 numflux.f90 source.f90
MODF= $(patsubst %,$(SRC_DIR)/%,$(MOD))

# Initial condition routines should be compiled before initialcondition.f90
INI = $(sort $(wildcard $(INI_DIR)/*.f90))

# Post-processing routines
POS = utils_analysis.f90
POSF = $(patsubst %,$(POS_DIR)/%,$(POS))

TARGETF = $(patsubst %,$(SRC_DIR)/%.f90,$(TARGET))
TARGETAF= $(patsubst %,$(POS_DIR)/%.f90,$(TARGETA))
TARGETLF= $(patsubst %,$(POS_DIR)/%.f90,$(TARGETL))
SRC = $(MODF) $(INI) \
      $(filter-out $(MODF) $(TARGETF), \
      $(sort $(wildcard $(SRC_DIR)/*.f90 $(SRC_DIR)/*.F90)))
SRCPOS = $(POSF) $(filter-out $(POSF), $(sort $(wildcard $(POS_DIR)/*.f90)))
OBJ  = $(patsubst $(SRC_DIR)/%.f90,%.o,$(SRC))
OBJ := $(patsubst $(INI_DIR)/%.f90,%.o,$(OBJ))
OBJ := $(patsubst $(SRC_DIR)/%.F90,%.o,$(OBJ))
OBJ_BASE = $(patsubst %.o,$(OBJ_DIR)/%.o,$(OBJ))
OBJ_POS  = $(filter-out $(TARGETAF) $(TARGETLF),$(SRCPOS))
OBJ_POS := $(patsubst $(POS_DIR)/%.f90,$(OBJ_DIR)/%.o,$(OBJ_POS))

OBJ_F  = $(OBJ_BASE) $(OBJ_DIR)/$(TARGET).o
OBJA_F = $(OBJ_BASE) $(OBJ_DIR)/$(TARGETA).o
OBJL_F = $(OBJ_BASE) $(OBJ_POS) $(OBJ_DIR)/$(TARGETL).o

all: $(TARGET)

debug: FCFLAG+=$(DBFLAG)
debug: $(TARGET)

profile: FCFLAG+=$(PRFLAG)
profile: $(TARGET)

ascii: $(TARGETA)

lc: $(TARGETL)

$(TARGET): $(OBJ_F)
	$(FC) $(FCFLAG) $(LDFLAGS) -o $@ $(OBJ_F)

$(TARGETA): $(OBJA_F)
	$(FC) $(FCFLAG) -o $@ $(OBJA_F)

$(TARGETL): $(OBJL_F)
	$(FC) $(FCFLAG) -o $@ $(OBJL_F)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90
	$(FC) $(FPPFLAG) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

$(OBJ_DIR)/%.o: $(INI_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

$(OBJ_DIR)/%.o: $(POS_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

.PHONY: clean
clean:
	rm -f *~ $(SRC_DIR)/*~ $(INI_DIR)/*~ $(POS_DIR)/*~ $(TARGET) $(TARGETA) $(TARGETL) $(OBJ_DIR)/*
