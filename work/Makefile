# Some versions of make default FC to f77
ifeq ($(FC),f77)
    FC = gfortran
endif
FC ?= gfortran
FCFLAG=

# for gfortran
ifeq (${FC},gfortran)
MODFLAG= -J
FCFLAG+= -fopenmp -ffree-line-length-512 #-fmax-stack-var-size=32768
DBFLAG= -ffpe-trap=invalid,zero,overflow -fbacktrace
PRFLAG= -pg
endif

# for ifort
ifeq (${FC},ifort)
MODFLAG= -module 
FCFLAG+= -qopenmp #-heap-arrays
DBFLAG= -traceback -fpe0 -CB #-pg #-check all
PRFLAG= -pg
endif

# directories
OBJ_DIR = obj
SRC_DIR = ../src
BIN_DIR = $(HORMONE_DIR)/bin

TARGET = hormone

# compile modules that are being depended on first
MOD = modules.f90 utils.f90 setup.f90 ionization.f90 eos.f90 fluxlimiter.f90 hlldflux.f90 dirichlet.f90 particles.f90 cooling.f90 composition.f90 star.f90 output.f90 smear.f90 readbin.f90 input.f90
MODF= $(patsubst %,$(SRC_DIR)/%,$(MOD))
TARGETF= $(patsubst %,$(SRC_DIR)/%.f90,$(TARGET))
SRC = $(MODF) \
      $(filter-out $(MODF) $(TARGETF), $(sort $(wildcard $(SRC_DIR)/*.f90))) \
      $(TARGETF)
OBJ = $(patsubst $(SRC_DIR)/%.f90,%.o,$(SRC))
OBJ_F = $(patsubst %.o,$(OBJ_DIR)/%.o,$(OBJ))

all: $(TARGET)

debug: FCFLAG+=$(DBFLAG)
debug: $(TARGET)

profile: FCFLAG+=$(PRFLAG)
profile: $(TARGET)

$(TARGET): $(OBJ_F)
	$(FC) $(FCFLAG) -o $@ $(OBJ_F)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

.PHONY: clean
clean:
	rm -f *~ $(SRC_DIR)/*~ hormone $(OBJ_DIR)/*
