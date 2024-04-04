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
DBFLAG= -ffpe-trap=invalid,zero,overflow -fbacktrace -O0
PRFLAG= -pg -save-temps
endif

# for ifort
ifeq (${FC},ifort)
ifeq (${MPI},yes)
$(error "MPI not supported for ifort")
endif
MODFLAG= -module # trailing space is important
FCFLAG:= $(CFLAG) -qopenmp -fp-model precise #-heap-arrays
DBFLAG= -traceback -fpe0 -CB -warn all -O0 -std #-pg #-check all
PRFLAG= -pg
endif

# directories
OBJ_DIR = obj
SRC_DIR = ../src
BIN_DIR = $(HORMONE_DIR)/bin

TARGET = hormone

# compile modules that are being depended on first
MOD = mpi_utils.F90 modules.f90 profiler.f90 utils.f90 conserve.f90 ionization.f90 eos.f90 miccg.f90 fluxlimiter.f90 hlldflux.f90 dirichlet.f90 particles.f90 radiation.f90 cooling.f90 sinks.f90 composition.f90 star.f90 shockfind.f90 output.f90 smear.f90 input.f90 setup.f90 readbin.f90 tests.f90 restart.f90 fluxboundary.f90 rungekutta.f90 timestep.f90 gridset.f90

# Initial condition routines should be compiled before initialcondition.f90
INI = eostest.f90 shocktube.f90 orszagtang.f90 sedov.f90 KHtests.f90 star_init.f90 rsg.f90 agndisk.f90 windtunnel.f90 rad_box.f90 stellarcollision.f90 modify.f90

MODF= $(patsubst %,$(SRC_DIR)/%,$(MOD))
INIF= $(patsubst %,$(SRC_DIR)/%,$(INI))
TARGETF= $(patsubst %,$(SRC_DIR)/%.f90,$(TARGET))
SRC = $(MODF) $(INIF) \
      $(filter-out $(MODF) $(INIF) $(TARGETF), $(sort $(wildcard $(SRC_DIR)/*.f90 $(SRC_DIR)/*.F90))) \
      $(TARGETF)
_OBJ = $(patsubst $(SRC_DIR)/%.f90,%.o,$(SRC))
OBJ = $(patsubst $(SRC_DIR)/%.F90,%.o,$(_OBJ))
OBJ_F = $(patsubst %.o,$(OBJ_DIR)/%.o,$(OBJ))

all: $(TARGET)

debug: FCFLAG+=$(DBFLAG)
debug: $(TARGET)

profile: FCFLAG+=$(PRFLAG)
profile: $(TARGET)

$(TARGET): $(OBJ_F)
	$(FC) $(FCFLAG) -o $@ $(OBJ_F)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90
	$(FC) $(FPPFLAG) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

.PHONY: clean
clean:
	rm -f *~ $(SRC_DIR)/*~ hormone $(OBJ_DIR)/*
