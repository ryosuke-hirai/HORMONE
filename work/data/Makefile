# Some versions of make default FC to f77
ifeq ($(FC),f77)
    FC = gfortran
endif
FC ?= gfortran
FCFLAG=

# for gfortran
ifeq (${FC},gfortran)
MODFLAG= -J
FCFLAG+=  -ffree-line-length-512 -fopenmp #-fmax-stack-var-size=32768
DBFLAG= -ffpe-trap=invalid,zero,overflow -fbacktrace 
endif

# for ifort
ifeq (${FC},ifort)
MODFLAG= -module 
FCFLAG+= -qopenmp #-heap-arrays
DBFLAG= -traceback -fpe0 -CB #-pg -check all
endif

# directories
OBJ_DIR = obj
SRC_DIR = ../../src

LCTARGET = lightcurve

# compile modules that are being depended on first
MOD = modules.f90 recombination.f90 eos.f90 composition.f90 checksetup.f90 dirichlet.f90 allocations.f90 readbin.f90

MODF= $(patsubst %,$(SRC_DIR)/%,$(MOD)) 
TARGETF= $(patsubst %,%.f90,$(LCTARGET))
SRC = $(MODF)
#SRC = $(MODF) \
#      $(filter-out $(MODF) $(TARGETF), $(sort $(wildcard $(SRC_DIR)/*.f90))) 
#      $(TARGETF)
OBJ = $(patsubst $(SRC_DIR)/%.f90,%.o,$(SRC))
OBJLC=$(patsubst %,$(OBJ_DIR)/%.o,$(LCTARGET))
OBJ_F = $(patsubst %.o,$(OBJ_DIR)/%.o,$(OBJ)) $(OBJLC)

all: $(LCTARGET)

debug: FCFLAG+=$(DBFLAG)
debug: $(LCTARGET)

#$(TARGET): $(OBJ_F)
#	$(FC) $(FCFLAG) -o $@ $(OBJ_F)

$(LCTARGET): $(OBJ_F)
	$(FC) $(FCFLAG) -o $@ $(OBJ_F)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

$(OBJLC): $(TARGETF)
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

.PHONY: clean
clean:
	rm -f *~ $(SRC_DIR)/*~ lightcurve $(OBJ_DIR)/*
